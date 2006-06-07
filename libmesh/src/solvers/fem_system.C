#include "dof_map.h"
#include "elem.h"
#include "fe_base.h"
#include "equation_systems.h"
#include "fem_system.h"
#include "libmesh_logging.h"
#include "mesh.h"
#include "numeric_vector.h"
#include "quadrature.h"
#include "sparse_matrix.h"
#include "time_solver.h"



FEMSystem::FEMSystem (EquationSystems& es,
                      const std::string& name,
                      const unsigned int number)
  : Parent(es, name, number)
{
}


FEMSystem::~FEMSystem ()
{
  this->clear();
}



void FEMSystem::clear()
{
  // We don't want to store AutoPtrs in STL containers, but we don't
  // want to leak memory either
  for (std::map<FEType, FEBase *>::iterator i = element_fe.begin();
       i != element_fe.end(); ++i)
    delete i->second;
  for (std::map<FEType, FEBase *>::iterator i = side_fe.begin();
       i != side_fe.end(); ++i)
    delete i->second;
  delete element_qrule;
  delete side_qrule;
  for (unsigned int i=0; i != elem_subsolutions.size(); ++i)
    {
      delete elem_subsolutions[i];
      delete elem_subresiduals[i];

      for (unsigned int j=0; j != elem_subjacobians[i].size(); ++j)
        delete elem_subjacobians[i][j];
    }
}



void FEMSystem::init_data ()
{
  // First initialize LinearImplicitSystem data
  Parent::init_data();

  const Mesh &mesh = this->get_mesh();

  unsigned int dim = mesh.mesh_dimension();

  // We need to know which of our variables has the hardest
  // shape functions to numerically integrate.

  unsigned int n_vars = this->n_vars();

  assert (n_vars);
  FEType hardest_fe_type = this->variable_type(0);

  for (unsigned int i=0; i != n_vars; ++i)
    {
      FEType fe_type = this->variable_type(i);

      // FIXME - we don't yet handle mixed finite elements from
      // different families 
      assert (fe_type.family == hardest_fe_type.family);
      if (fe_type.order > hardest_fe_type.order)
        hardest_fe_type = fe_type;
    }

  // Create an adequate quadrature rule
  element_qrule =
    hardest_fe_type.default_quadrature_rule(dim).release();
  side_qrule =
    hardest_fe_type.default_quadrature_rule(dim-1).release();

  // Next, create finite element objects
  for (unsigned int i=0; i != n_vars; ++i)
    {
      FEType fe_type = this->variable_type(i);
      if (element_fe[fe_type] == NULL)
        {
          element_fe[fe_type] = FEBase::build(dim, fe_type).release();
          element_fe[fe_type]->attach_quadrature_rule(element_qrule);
          side_fe[fe_type] = FEBase::build(dim, fe_type).release();
          side_fe[fe_type]->attach_quadrature_rule(side_qrule);
        }
    }
}



void FEMSystem::assembly (bool get_residual, bool get_jacobian)
{
  START_LOG("assembly", "FEMSystem");
  assert(get_residual || get_jacobian);

  const Mesh& mesh = this->get_mesh();

  const DofMap& dof_map = this->get_dof_map();

  // Global nonlinear solution
  NumericVector<Number> &nonlinear_solution =
    this->get_vector("_nonlinear_solution");

  // Is this definitely necessary? [RHS]
  if (get_jacobian)
    matrix->zero();
  if (get_residual)
    rhs->zero();

  // In time-dependent problems, the nonlinear function we're trying
  // to solve at each timestep may depend on the particular solver
  // we're using
  assert (time_solver.get() != NULL);

  MeshBase::const_element_iterator el =
    mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =
    mesh.active_local_elements_end();

  for ( ; el != end_el; ++el)
    {
      elem = *el;

      dof_map.dof_indices (elem, dof_indices);
      unsigned int n_dofs = dof_indices.size();

      elem_solution.resize(n_dofs);
      for (unsigned int i=0; i != n_dofs; ++i)
        elem_solution(i) = nonlinear_solution(dof_indices[i]);

      elem_residual.resize(n_dofs);
      elem_jacobian.resize(n_dofs, n_dofs);

      unsigned int sub_dofs = 0;

      for (unsigned int i=0; i != this->n_vars(); ++i)
        {
          dof_map.dof_indices (elem, dof_indices_var[i], i);

          elem_subsolutions[i]->reposition
            (sub_dofs, dof_indices_var[i].size());

          elem_subresiduals[i]->reposition
            (sub_dofs, dof_indices_var[i].size());

          for (unsigned int j=0; j != i; ++j)
            {
              elem_subjacobians[i][j]->reposition
                (sub_dofs, elem_subresiduals[j]->i_off(),
                 dof_indices_var[i].size(),
                 dof_indices_var[j].size());
              elem_subjacobians[j][i]->reposition
                (elem_subresiduals[j]->i_off(), sub_dofs,
                 dof_indices_var[j].size(),
                 dof_indices_var[i].size());
            }
          elem_subjacobians[i][i]->reposition
            (sub_dofs, sub_dofs,
             dof_indices_var[i].size(),
             dof_indices_var[i].size());
          sub_dofs += dof_indices_var[i].size();
        }
      assert(sub_dofs == n_dofs);

      std::map<FEType, FEBase *>::iterator fe_end = element_fe.end();
      for (std::map<FEType, FEBase *>::iterator i = element_fe.begin();
           i != fe_end; ++i)
        {
          i->second->reinit(elem);
        }
      
      bool jacobian_computed = time_solver->element_residual(get_jacobian);

      for (side = 0; side != elem->n_sides(); ++side)
        {
          // Don't compute on non-boundary sides unless requested
          if (!compute_internal_sides && elem->neighbor(side) != NULL)
            continue;

          fe_end = side_fe.end();
          for (std::map<FEType, FEBase *>::iterator i = side_fe.begin();
               i != fe_end; ++i)
            {
              i->second->reinit(elem, side);
            }
          jacobian_computed = time_solver->side_residual(get_jacobian)
                              && jacobian_computed;
        }

      if (get_jacobian && !jacobian_computed)
        {
          PAUSE_LOG("assembly", "FEMSystem");
          this->compute_numerical_jacobian();
          RESTART_LOG("assembly", "FEMSystem");
        }

      if (get_residual && get_jacobian)
        this->get_dof_map().constrain_element_matrix_and_vector
          (elem_jacobian, elem_residual, dof_indices);
      else if (get_residual)
        this->get_dof_map().constrain_element_vector
          (elem_residual, dof_indices);
      else if (get_jacobian)
        this->get_dof_map().constrain_element_matrix
          (elem_jacobian, dof_indices);

      if (get_jacobian)
        this->matrix->add_matrix (elem_jacobian, dof_indices);
      if (get_residual)
        this->rhs->add_vector (elem_residual, dof_indices);
    }


  if (get_residual && print_residual_norms)
    {
      this->rhs->close();
      std::cerr << "|F| = " << this->rhs->l1_norm() << std::endl;
    }
  if (get_residual && print_residuals)
    {
      this->rhs->close();
      unsigned int old_precision = std::cerr.precision();
      std::cerr.precision(15);
      std::cerr << "F = [" << *(this->rhs) << "];" << std::endl;
      std::cerr.precision(old_precision);
    }
  if (get_jacobian && print_jacobian_norms)
    {
      this->matrix->close();
      std::cerr << "|J| = " << this->matrix->l1_norm() << std::endl;
    }
  if (get_jacobian && print_jacobians)
    {
      this->matrix->close();
      unsigned int old_precision = std::cerr.precision();
      std::cerr.precision(15);
      std::cerr << "J = [" << *(this->matrix) << "];" << std::endl;
      std::cerr.precision(old_precision);
    }
  STOP_LOG("assembly", "FEMSystem");
}



void FEMSystem::compute_numerical_jacobian ()
{
  START_LOG("compute_numerical_jacobian", "FEMSystem");
  Real h = 1e-6;  // FIXME: How do we scale this?

  DenseVector<Number> original_residual(elem_residual);
  DenseVector<Number> backwards_residual(elem_residual);
  DenseMatrix<Number> numerical_jacobian(elem_jacobian);

  for (unsigned int j = 0; j != dof_indices.size(); ++j)
    {
      Number original_solution = elem_solution(j);
      elem_solution(j) -= h;
      elem_residual.zero();
      time_solver->element_residual(false);

      for (side = 0; side != elem->n_sides(); ++side)
        {
          // Don't compute on non-boundary sides unless requested
          if (!compute_internal_sides && elem->neighbor(side) != NULL)
            continue;

          std::map<FEType, FEBase *>::iterator fe_end = side_fe.end();
          for (std::map<FEType, FEBase *>::iterator i = side_fe.begin();
               i != fe_end; ++i)
            {
              i->second->reinit(elem, side);
            }
          time_solver->side_residual(false);
        }
      backwards_residual = elem_residual;

      elem_solution(j) = original_solution + h;
      elem_residual.zero();
      time_solver->element_residual(false);

      for (side = 0; side != elem->n_sides(); ++side)
        {
          // Don't compute on non-boundary sides unless requested
          if (!compute_internal_sides && elem->neighbor(side) != NULL)
            continue;

          std::map<FEType, FEBase *>::iterator fe_end = side_fe.end();
          for (std::map<FEType, FEBase *>::iterator i = side_fe.begin();
               i != fe_end; ++i)
            {
              i->second->reinit(elem, side);
            }
          time_solver->side_residual(false);
        }

      for (unsigned int i = 0; i != dof_indices.size(); ++i)
        {
          numerical_jacobian(i,j) =
            (elem_residual(i) - backwards_residual(i)) / 2 / h;
        }
      elem_solution(j) = original_solution;
    }

  elem_residual = original_residual;
  elem_jacobian = numerical_jacobian;
  STOP_LOG("compute_numerical_jacobian", "FEMSystem");
}



void FEMSystem::time_evolving (unsigned int var)
{
  // Call the parent function
  Parent::time_evolving(var);

  // Then make sure we're prepared to do mass integration
  element_fe[this->variable_type(var)]->get_JxW();
  element_fe[this->variable_type(var)]->get_phi();
}



bool FEMSystem::mass_residual (bool request_jacobian)
{
  unsigned int n_qpoints = element_qrule->n_points();

  for (unsigned int var = 0; var != this->n_vars(); ++var)
    {
      if (!_time_evolving[var])
        continue;

      const std::vector<Real> &JxW = 
        element_fe[this->variable_type(var)]->get_JxW();

      const std::vector<std::vector<Real> > &phi =
        element_fe[this->variable_type(var)]->get_phi();

      const unsigned int n_dofs = dof_indices_var[var].size();

      DenseSubVector<Number> &Fu = *elem_subresiduals[var];
      DenseSubMatrix<Number> &Kuu = *elem_subjacobians[var][var];

      for (unsigned int qp = 0; qp != n_qpoints; ++qp)
        {
          Number u = 0.;
          for (unsigned int l = 0; l != n_dofs; ++l)
            {
              u += phi[l][qp] * (*elem_subsolutions[var])(l);
            }
          for (unsigned int i = 0; i != n_dofs; ++i)
            {
              Fu(i) += JxW[qp] * u * phi[i][qp];
              if (request_jacobian)
                for (unsigned int j = 0; j != n_dofs; ++j)
                  {
                    Kuu(i,j) += JxW[qp] * phi[i][qp] * phi[j][qp];
                  }
            }
        }
    }

  return request_jacobian;
}
