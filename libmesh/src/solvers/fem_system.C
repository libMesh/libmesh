#include "dof_map.h"
#include "elem.h"
#include "fe_base.h"
#include "equation_systems.h"
#include "fem_system.h"
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

  // Finally, create subvector and submatrix objects
  elem_subsolutions.clear();
  elem_subsolutions.reserve(n_vars);
  elem_subresiduals.clear();
  elem_subresiduals.reserve(n_vars);
  elem_subjacobians.clear();
  elem_subjacobians.resize(n_vars);
  for (unsigned int i=0; i != n_vars; ++i)
    {
      elem_subsolutions.push_back(new DenseSubVector<Number>(elem_solution));
      elem_subresiduals.push_back(new DenseSubVector<Number>(elem_residual));
      elem_subjacobians[i].clear();
      elem_subjacobians[i].reserve(n_vars);

      for (unsigned int j=0; j != n_vars; ++j)
        {
          elem_subjacobians[j].push_back
            (new DenseSubMatrix<Number>(elem_jacobian));
        }
    }
}



void FEMSystem::assembly (bool get_residual, bool get_jacobian)
{
  assert(get_residual || get_jacobian);

  const Mesh& mesh = this->get_mesh();

  const DofMap& dof_map = this->get_dof_map();

  // Is this definitely necessary? [RHS]
  if (get_jacobian)
    matrix->zero();
  if (get_residual)
    rhs->zero();

  // Global nonlinear solution
  NumericVector<Number> &nonlinear_solution =
    this->get_vector("_nonlinear_solution");

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

      elem_residual.resize(dof_indices.size());
      elem_jacobian.resize(dof_indices.size(), dof_indices.size());

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
          for (std::map<FEType, FEBase *>::iterator i = element_fe.begin();
               i != fe_end; ++i)
            {
              i->second->reinit(elem, side);
            }
          jacobian_computed = time_solver->side_residual(get_jacobian)
                              && jacobian_computed;
        }

      if (get_jacobian && !jacobian_computed)
        this->compute_numerical_jacobian();

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
}



void FEMSystem::compute_numerical_jacobian ()
{
  Real h = 1e-6;  // FIXME: How do we scale this?

  DenseVector<Number> original_residual(elem_residual);

  for (unsigned int j = 0; j != dof_indices.size(); ++j)
    {
      Number original_solution = elem_solution(j);
      elem_solution(j) += h;
      time_solver->element_residual(false);
      for (unsigned int i = 0; i != dof_indices.size(); ++i)
        {
          elem_jacobian(i,j) =
            (elem_residual(i) - original_residual(i)) / h;
        }
      elem_solution(j) = original_solution;
    }

  elem_residual = original_residual;
}
