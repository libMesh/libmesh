#include "dof_map.h"
#include "elem.h"
#include "fe_base.h"
#include "fe_interface.h"
#include "equation_systems.h"
#include "fem_system.h"
#include "libmesh_logging.h"
#include "mesh.h"
#include "numeric_vector.h"
#include "quadrature.h"
#include "sparse_matrix.h"
#include "time_solver.h"
#include "unsteady_solver.h" // For euler_residual



FEMSystem::FEMSystem (EquationSystems& es,
                      const std::string& name,
                      const unsigned int number)
  : Parent(es, name, number),
    fe_reinit_during_postprocess(true),
    extra_quadrature_order(0),
    numerical_jacobian_h(1.e-6),
    verify_analytic_jacobians(0.0),
    element_qrule(NULL), side_qrule(NULL), elem(NULL),
    _mesh_sys(libMesh::invalid_uint),
    _mesh_x_var(libMesh::invalid_uint),
    _mesh_y_var(libMesh::invalid_uint),
    _mesh_z_var(libMesh::invalid_uint)
{
}


FEMSystem::~FEMSystem ()
{
  this->clear();
}



Number FEMSystem::interior_value(unsigned int var, unsigned int qp)
{
  // Get local-to-global dof index lookup
  libmesh_assert (dof_indices.size() > var);
  const unsigned int n_dofs = dof_indices_var[var].size();

  // Get current local coefficients
  libmesh_assert (elem_subsolutions.size() > var);
  libmesh_assert (elem_subsolutions[var] != NULL);
  DenseSubVector<Number> &coef = *elem_subsolutions[var];

  // Get shape function values at quadrature point
  const std::vector<std::vector<Real> > &phi =
    element_fe_var[var]->get_phi();

  // Accumulate solution value
  Number u = 0.;

  for (unsigned int l=0; l != n_dofs; l++)
    u += phi[l][qp] * coef(l);

  return u;
}



Gradient FEMSystem::interior_gradient(unsigned int var, unsigned int qp)
{
  // Get local-to-global dof index lookup
  libmesh_assert (dof_indices.size() > var);
  const unsigned int n_dofs = dof_indices_var[var].size();

  // Get current local coefficients
  libmesh_assert (elem_subsolutions.size() > var);
  libmesh_assert (elem_subsolutions[var] != NULL);
  DenseSubVector<Number> &coef = *elem_subsolutions[var];

  // Get shape function values at quadrature point
  const std::vector<std::vector<RealGradient> > &dphi =
    element_fe_var[var]->get_dphi();

  // Accumulate solution derivatives
  Gradient du;

  for (unsigned int l=0; l != n_dofs; l++)
    du.add_scaled(dphi[l][qp], coef(l));

  return du;
}



#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
Tensor FEMSystem::interior_hessian(unsigned int var, unsigned int qp)
{
  // Get local-to-global dof index lookup
  libmesh_assert (dof_indices.size() > var);
  const unsigned int n_dofs = dof_indices_var[var].size();

  // Get current local coefficients
  libmesh_assert (elem_subsolutions.size() > var);
  libmesh_assert (elem_subsolutions[var] != NULL);
  DenseSubVector<Number> &coef = *elem_subsolutions[var];

  // Get shape function values at quadrature point
  const std::vector<std::vector<RealTensor> > &d2phi =
    element_fe_var[var]->get_d2phi();

  // Accumulate solution second derivatives
  Tensor d2u;

  for (unsigned int l=0; l != n_dofs; l++)
    d2u.add_scaled(d2phi[l][qp], coef(l));

  return d2u;
}
#endif // ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES



Number FEMSystem::side_value(unsigned int var, unsigned int qp)
{
  // Get local-to-global dof index lookup
  libmesh_assert (dof_indices.size() > var);
  const unsigned int n_dofs = dof_indices_var[var].size();

  // Get current local coefficients
  libmesh_assert (elem_subsolutions.size() > var);
  libmesh_assert (elem_subsolutions[var] != NULL);
  DenseSubVector<Number> &coef = *elem_subsolutions[var];

  // Get shape function values at quadrature point
  const std::vector<std::vector<Real> > &phi =
    side_fe_var[var]->get_phi();

  // Accumulate solution value
  Number u = 0.;

  for (unsigned int l=0; l != n_dofs; l++)
    u += phi[l][qp] * coef(l);

  return u;
}



Gradient FEMSystem::side_gradient(unsigned int var, unsigned int qp)
{
  // Get local-to-global dof index lookup
  libmesh_assert (dof_indices.size() > var);
  const unsigned int n_dofs = dof_indices_var[var].size();

  // Get current local coefficients
  libmesh_assert (elem_subsolutions.size() > var);
  libmesh_assert (elem_subsolutions[var] != NULL);
  DenseSubVector<Number> &coef = *elem_subsolutions[var];

  // Get shape function values at quadrature point
  const std::vector<std::vector<RealGradient> > &dphi =
    side_fe_var[var]->get_dphi();

  // Accumulate solution derivatives
  Gradient du;

  for (unsigned int l=0; l != n_dofs; l++)
    du.add_scaled(dphi[l][qp], coef(l));

  return du;
}



#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
Tensor FEMSystem::side_hessian(unsigned int var, unsigned int qp)
{
  // Get local-to-global dof index lookup
  libmesh_assert (dof_indices.size() > var);
  const unsigned int n_dofs = dof_indices_var[var].size();

  // Get current local coefficients
  libmesh_assert (elem_subsolutions.size() > var);
  libmesh_assert (elem_subsolutions[var] != NULL);
  DenseSubVector<Number> &coef = *elem_subsolutions[var];

  // Get shape function values at quadrature point
  const std::vector<std::vector<RealTensor> > &d2phi =
    side_fe_var[var]->get_d2phi();

  // Accumulate solution second derivatives
  Tensor d2u;

  for (unsigned int l=0; l != n_dofs; l++)
    d2u.add_scaled(d2phi[l][qp], coef(l));

  return d2u;
}
#endif // ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES



Number FEMSystem::point_value(unsigned int var, Point &p)
{
  // Get local-to-global dof index lookup
  libmesh_assert (dof_indices.size() > var);
  const unsigned int n_dofs = dof_indices_var[var].size();

  
  // Get current local coefficients
  libmesh_assert (elem_subsolutions.size() > var);
  libmesh_assert (elem_subsolutions[var] != NULL);
  DenseSubVector<Number> &coef = *elem_subsolutions[var];

  Number u = 0.;

  unsigned int dim = get_mesh().mesh_dimension();
  FEType fe_type = element_fe_var[var]->get_fe_type();
  Point p_master = FEInterface::inverse_map(dim, fe_type, elem, p);

  for (unsigned int l=0; l != n_dofs; l++)
    u += FEInterface::shape(dim, fe_type, elem, l, p_master)
         * coef(l);

  return u;
}



Number FEMSystem::fixed_interior_value(unsigned int var, unsigned int qp)
{
  // Get local-to-global dof index lookup
  libmesh_assert (dof_indices.size() > var);
  const unsigned int n_dofs = dof_indices_var[var].size();

  // Get current local coefficients
  libmesh_assert (elem_fixed_subsolutions.size() > var);
  libmesh_assert (elem_fixed_subsolutions[var] != NULL);
  DenseSubVector<Number> &coef = *elem_fixed_subsolutions[var];

  // Get shape function values at quadrature point
  const std::vector<std::vector<Real> > &phi =
    element_fe_var[var]->get_phi();

  // Accumulate solution value
  Number u = 0.;

  for (unsigned int l=0; l != n_dofs; l++)
    u += phi[l][qp] * coef(l);

  return u;
}



Gradient FEMSystem::fixed_interior_gradient(unsigned int var, unsigned int qp)
{
  // Get local-to-global dof index lookup
  libmesh_assert (dof_indices.size() > var);
  const unsigned int n_dofs = dof_indices_var[var].size();

  // Get current local coefficients
  libmesh_assert (elem_fixed_subsolutions.size() > var);
  libmesh_assert (elem_fixed_subsolutions[var] != NULL);
  DenseSubVector<Number> &coef = *elem_fixed_subsolutions[var];

  // Get shape function values at quadrature point
  const std::vector<std::vector<RealGradient> > &dphi =
    element_fe_var[var]->get_dphi();

  // Accumulate solution derivatives
  Gradient du;

  for (unsigned int l=0; l != n_dofs; l++)
    du.add_scaled(dphi[l][qp], coef(l));

  return du;
}



#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
Tensor FEMSystem::fixed_interior_hessian(unsigned int var, unsigned int qp)
{
  // Get local-to-global dof index lookup
  libmesh_assert (dof_indices.size() > var);
  const unsigned int n_dofs = dof_indices_var[var].size();

  // Get current local coefficients
  libmesh_assert (elem_fixed_subsolutions.size() > var);
  libmesh_assert (elem_fixed_subsolutions[var] != NULL);
  DenseSubVector<Number> &coef = *elem_fixed_subsolutions[var];

  // Get shape function values at quadrature point
  const std::vector<std::vector<RealTensor> > &d2phi =
    element_fe_var[var]->get_d2phi();

  // Accumulate solution second derivatives
  Tensor d2u;

  for (unsigned int l=0; l != n_dofs; l++)
    d2u.add_scaled(d2phi[l][qp], coef(l));

  return d2u;
}
#endif // ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES



Number FEMSystem::fixed_side_value(unsigned int var, unsigned int qp)
{
  // Get local-to-global dof index lookup
  libmesh_assert (dof_indices.size() > var);
  const unsigned int n_dofs = dof_indices_var[var].size();

  // Get current local coefficients
  libmesh_assert (elem_fixed_subsolutions.size() > var);
  libmesh_assert (elem_fixed_subsolutions[var] != NULL);
  DenseSubVector<Number> &coef = *elem_fixed_subsolutions[var];

  // Get shape function values at quadrature point
  const std::vector<std::vector<Real> > &phi =
    side_fe_var[var]->get_phi();

  // Accumulate solution value
  Number u = 0.;

  for (unsigned int l=0; l != n_dofs; l++)
    u += phi[l][qp] * coef(l);

  return u;
}



Gradient FEMSystem::fixed_side_gradient(unsigned int var, unsigned int qp)
{
  // Get local-to-global dof index lookup
  libmesh_assert (dof_indices.size() > var);
  const unsigned int n_dofs = dof_indices_var[var].size();

  // Get current local coefficients
  libmesh_assert (elem_fixed_subsolutions.size() > var);
  libmesh_assert (elem_fixed_subsolutions[var] != NULL);
  DenseSubVector<Number> &coef = *elem_fixed_subsolutions[var];

  // Get shape function values at quadrature point
  const std::vector<std::vector<RealGradient> > &dphi =
    side_fe_var[var]->get_dphi();

  // Accumulate solution derivatives
  Gradient du;

  for (unsigned int l=0; l != n_dofs; l++)
    du.add_scaled(dphi[l][qp], coef(l));

  return du;
}



#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
Tensor FEMSystem::fixed_side_hessian(unsigned int var, unsigned int qp)
{
  // Get local-to-global dof index lookup
  libmesh_assert (dof_indices.size() > var);
  const unsigned int n_dofs = dof_indices_var[var].size();

  // Get current local coefficients
  libmesh_assert (elem_fixed_subsolutions.size() > var);
  libmesh_assert (elem_fixed_subsolutions[var] != NULL);
  DenseSubVector<Number> &coef = *elem_fixed_subsolutions[var];

  // Get shape function values at quadrature point
  const std::vector<std::vector<RealTensor> > &d2phi =
    side_fe_var[var]->get_d2phi();

  // Accumulate solution second derivatives
  Tensor d2u;

  for (unsigned int l=0; l != n_dofs; l++)
    d2u.add_scaled(d2phi[l][qp], coef(l));

  return d2u;
}
#endif // ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES



Number FEMSystem::fixed_point_value(unsigned int var, Point &p)
{
  // Get local-to-global dof index lookup
  libmesh_assert (dof_indices.size() > var);
  const unsigned int n_dofs = dof_indices_var[var].size();

  
  // Get current local coefficients
  libmesh_assert (elem_fixed_subsolutions.size() > var);
  libmesh_assert (elem_fixed_subsolutions[var] != NULL);
  DenseSubVector<Number> &coef = *elem_fixed_subsolutions[var];

  Number u = 0.;

  unsigned int dim = get_mesh().mesh_dimension();
  FEType fe_type = element_fe_var[var]->get_fe_type();
  Point p_master = FEInterface::inverse_map(dim, fe_type, elem, p);

  for (unsigned int l=0; l != n_dofs; l++)
    u += FEInterface::shape(dim, fe_type, elem, l, p_master)
         * coef(l);

  return u;
}



void FEMSystem::clear()
{
  this->clear_fem_ptrs();

  Parent::clear();
}



void FEMSystem::clear_fem_ptrs()
{
  // We don't want to store AutoPtrs in STL containers, but we don't
  // want to leak memory either
  for (std::map<FEType, FEBase *>::iterator i = element_fe.begin();
       i != element_fe.end(); ++i)
    delete i->second;
  element_fe.clear();

  for (std::map<FEType, FEBase *>::iterator i = side_fe.begin();
       i != side_fe.end(); ++i)
    delete i->second;
  side_fe.clear();

  delete element_qrule;
  element_qrule = NULL;

  delete side_qrule;
  side_qrule = NULL;
}



void FEMSystem::init_data ()
{
  // We may have already been initialized once
  this->clear_fem_ptrs();

  // First initialize LinearImplicitSystem data
  Parent::init_data();

  const MeshBase &mesh = this->get_mesh();

  unsigned int dim = mesh.mesh_dimension();

  // We need to know which of our variables has the hardest
  // shape functions to numerically integrate.

  unsigned int n_vars = this->n_vars();

  libmesh_assert (n_vars);
  FEType hardest_fe_type = this->variable_type(0);

  for (unsigned int i=0; i != n_vars; ++i)
    {
      FEType fe_type = this->variable_type(i);

      // FIXME - we don't yet handle mixed finite elements from
      // different families which require different quadrature rules
      // libmesh_assert (fe_type.family == hardest_fe_type.family);

      if (fe_type.order > hardest_fe_type.order)
        hardest_fe_type = fe_type;
    }

  // Create an adequate quadrature rule
  element_qrule = hardest_fe_type.default_quadrature_rule
    (dim, extra_quadrature_order).release();
  side_qrule = hardest_fe_type.default_quadrature_rule
    (dim-1, extra_quadrature_order).release();

  // Next, create finite element objects
  element_fe_var.resize(n_vars);
  side_fe_var.resize(n_vars);
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
      element_fe_var[i] = element_fe[fe_type];
      side_fe_var[i] = side_fe[fe_type];
    }
}


void FEMSystem::elem_reinit(Real theta)
{
  // Handle a moving element if necessary
  if (_mesh_sys != libMesh::invalid_uint)
    {
      elem_position_set(theta);
      elem_fe_reinit();
    }
}


void FEMSystem::elem_side_reinit(Real theta)
{
  // Handle a moving element if necessary
  if (_mesh_sys != libMesh::invalid_uint)
    {
      elem_position_set(theta);
      elem_side_fe_reinit();
    }
}


void FEMSystem::elem_fe_reinit ()
{
  // Initialize all the interior FE objects on elem.
  // Logging of FE::reinit is done in the FE functions
  std::map<FEType, FEBase *>::iterator fe_end = element_fe.end();
  for (std::map<FEType, FEBase *>::iterator i = element_fe.begin();
       i != fe_end; ++i)
    {
      i->second->reinit(elem);
    }
}


void FEMSystem::elem_side_fe_reinit ()
{
  // Initialize all the interior FE objects on elem/side.
  // Logging of FE::reinit is done in the FE functions
  std::map<FEType, FEBase *>::iterator fe_end = side_fe.end();
  for (std::map<FEType, FEBase *>::iterator i = side_fe.begin();
       i != fe_end; ++i)
    {
      i->second->reinit(elem, side);
    }
}


void FEMSystem::assembly (bool get_residual, bool get_jacobian)
{
  libmesh_assert(get_residual || get_jacobian);
  std::string log_name;
  if (get_residual && get_jacobian)
    log_name = "assembly()";
  else if (get_residual)
    log_name = "assembly(get_residual)";
  else
    log_name = "assembly(get_jacobian)";

  START_LOG(log_name, "FEMSystem");

  const MeshBase& mesh = this->get_mesh();

  const DofMap& dof_map = this->get_dof_map();

//  this->get_vector("_nonlinear_solution").localize
//    (*current_local_nonlinear_solution,
//     dof_map.get_send_list());
  this->update();

  if (print_solution_norms)
    {
//      this->get_vector("_nonlinear_solution").close();
      this->solution->close();
      std::cout << "|U| = "
//                << this->get_vector("_nonlinear_solution").l1_norm()
                << this->solution->l1_norm()
                << std::endl;
    }
  if (print_solutions)
    {
      unsigned int old_precision = std::cout.precision();
      std::cout.precision(16);
//      std::cout << "U = [" << this->get_vector("_nonlinear_solution")
      std::cout << "U = [" << *(this->solution)
                << "];" << std::endl;
      std::cout.precision(old_precision);
    }

  // Is this definitely necessary? [RHS]
  if (get_jacobian)
    matrix->zero();
  if (get_residual)
    rhs->zero();

  // Stupid C++ lets you set *Real* verify_analytic_jacobians = true!
  if (verify_analytic_jacobians > 0.5)
    {
      std::cerr << "WARNING!  verify_analytic_jacobians was set "
                << "to absurdly large value of "
                << verify_analytic_jacobians << std::endl;
      std::cerr << "Resetting to 1e-6!" << std::endl;
      verify_analytic_jacobians = 1e-6;
    }

  // In time-dependent problems, the nonlinear function we're trying
  // to solve at each timestep may depend on the particular solver
  // we're using
  libmesh_assert (time_solver.get() != NULL);

  // Build the residual and jacobian contributions on every active
  // mesh element on this processor
  MeshBase::const_element_iterator el =
    mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =
    mesh.active_local_elements_end();

  for ( ; el != end_el; ++el)
    {
      elem = *el;

      // Initialize the per-element data for elem.
      dof_map.dof_indices (elem, dof_indices);
      unsigned int n_dofs = dof_indices.size();

      elem_solution.resize(n_dofs);
      if (use_fixed_solution)
        elem_fixed_solution.resize(n_dofs);

      for (unsigned int i=0; i != n_dofs; ++i)
//        elem_solution(i) = current_nonlinear_solution(dof_indices[i]);
        elem_solution(i) = current_solution(dof_indices[i]);

      // These resize calls also zero out the residual and jacobian
      elem_residual.resize(n_dofs);
      elem_jacobian.resize(n_dofs, n_dofs);

      // Initialize the per-variable data for elem.
      unsigned int sub_dofs = 0;
      for (unsigned int i=0; i != this->n_vars(); ++i)
        {
          dof_map.dof_indices (elem, dof_indices_var[i], i);

          elem_subsolutions[i]->reposition
            (sub_dofs, dof_indices_var[i].size());

          if (use_fixed_solution)
            elem_fixed_subsolutions[i]->reposition
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
      libmesh_assert(sub_dofs == n_dofs);

      // Moving the mesh may be necessary
      if (_mesh_sys != libMesh::invalid_uint)
        elem_position_set(1.);
      // Reinitializing the FE objects is definitely necessary
      this->elem_fe_reinit();

      bool jacobian_computed = time_solver->element_residual(get_jacobian);

      // Compute a numeric jacobian if we have to
      if (get_jacobian && !jacobian_computed)
        {
          // Make sure we didn't compute a jacobian and lie about it
          libmesh_assert(elem_jacobian.l1_norm() == 0.0);
          // Logging of numerical jacobians is done separately
          this->numerical_elem_jacobian();
        }

      // Compute a numeric jacobian if we're asked to verify the
      // analytic jacobian we got
      if (get_jacobian && jacobian_computed &&
          verify_analytic_jacobians != 0.0)
        {
          DenseMatrix<Number> analytic_jacobian(elem_jacobian);

          elem_jacobian.zero();
          // Logging of numerical jacobians is done separately
          this->numerical_elem_jacobian();

          Real analytic_norm = analytic_jacobian.l1_norm();
          Real numerical_norm = elem_jacobian.l1_norm();

          // If we can continue, we'll probably prefer the analytic jacobian
          analytic_jacobian.swap(elem_jacobian);

          // The matrix "analytic_jacobian" will now hold the error matrix
          analytic_jacobian.add(-1.0, elem_jacobian);
          Real error_norm = analytic_jacobian.l1_norm();

          Real relative_error = error_norm /
                                std::max(analytic_norm, numerical_norm);

          if (relative_error > verify_analytic_jacobians)
            {
              std::cerr << "Relative error " << relative_error
                        << " detected in analytic jacobian on element "
                        << elem->id() << '!' << std::endl;

              unsigned int old_precision = std::cout.precision();
              std::cout.precision(16);
	      std::cout << "J_analytic " << elem->id() << " = "
                        << elem_jacobian << std::endl;
              analytic_jacobian.add(1.0, elem_jacobian);
	      std::cout << "J_numeric " << elem->id() << " = "
                        << analytic_jacobian << std::endl;

              std::cout.precision(old_precision);

              libmesh_error();
            }
        }

      for (side = 0; side != elem->n_sides(); ++side)
        {
          // Don't compute on non-boundary sides unless requested
          if (!compute_internal_sides && elem->neighbor(side) != NULL)
            continue;

          // Any mesh movement has already been done (and restored,
          // if the TimeSolver isn't broken), but
          // reinitializing the side FE objects is still necessary
          this->elem_side_fe_reinit();

          DenseMatrix<Number> old_jacobian;
          // If we're in DEBUG mode, we should always verify that the
	  // user's side_residual function doesn't alter our existing
          // jacobian and then lie about it
#ifndef DEBUG
	  // Even if we're not in DEBUG mode, when we're verifying
	  // analytic jacobians we'll want to verify each side's
          // jacobian contribution separately
          if (verify_analytic_jacobians != 0.0 && get_jacobian)
#endif // ifndef DEBUG
            {
              old_jacobian = elem_jacobian;
              elem_jacobian.zero();
            }
          jacobian_computed = time_solver->side_residual(get_jacobian);

          // Compute a numeric jacobian if we have to
          if (get_jacobian && !jacobian_computed)
            {
	      // In DEBUG mode, we've already set elem_jacobian == 0,
	      // so we can make sure side_residual didn't compute a
              // jacobian and lie about it
#ifdef DEBUG
              libmesh_assert(elem_jacobian.l1_norm() == 0.0);
#endif
              // Logging of numerical jacobians is done separately
              this->numerical_side_jacobian();

              // If we're in DEBUG mode or if
	      // verify_analytic_jacobians is on, we've moved
              // elem_jacobian's accumulated values into old_jacobian.
              // Now let's add them back.
#ifndef DEBUG
              if (verify_analytic_jacobians != 0.0)
#endif // ifndef DEBUG
                elem_jacobian += old_jacobian;
            }

          // Compute a numeric jacobian if we're asked to verify the
          // analytic jacobian we got
	  if (get_jacobian && jacobian_computed &&
              verify_analytic_jacobians != 0.0)
            {
              DenseMatrix<Number> analytic_jacobian(elem_jacobian);

              elem_jacobian.zero();
              // Logging of numerical jacobians is done separately
              this->numerical_side_jacobian();

              Real analytic_norm = analytic_jacobian.l1_norm();
              Real numerical_norm = elem_jacobian.l1_norm();

              // If we can continue, we'll probably prefer the analytic jacobian
              analytic_jacobian.swap(elem_jacobian);

              // The matrix "analytic_jacobian" will now hold the error matrix
              analytic_jacobian.add(-1.0, elem_jacobian);
              Real error_norm = analytic_jacobian.l1_norm();

              Real relative_error = error_norm /
                                    std::max(analytic_norm, numerical_norm);

              if (relative_error > verify_analytic_jacobians)
                {
                  std::cerr << "Relative error " << relative_error
                            << " detected in analytic jacobian on element "
                            << elem->id()
			    << ", side " << side << '!' << std::endl;

                  unsigned int old_precision = std::cout.precision();
                  std::cout.precision(16);
	          std::cout << "J_analytic " << elem->id() << " = "
                            << elem_jacobian << std::endl;
                  analytic_jacobian.add(1.0, elem_jacobian);
	          std::cout << "J_numeric " << elem->id() << " = "
                            << analytic_jacobian << std::endl;
                  std::cout.precision(old_precision);

                  libmesh_error();
                }
              // Once we've verified a side, we'll want to add back the
              // rest of the accumulated jacobian
              elem_jacobian += old_jacobian;
            }
	  // In DEBUG mode, we've set elem_jacobian == 0, and we
          // may still need to add the old jacobian back
#ifdef DEBUG
	  if (get_jacobian && jacobian_computed &&
              verify_analytic_jacobians == 0.0)
            {
              elem_jacobian += old_jacobian;
            }
#endif // ifdef DEBUG
        }

#ifdef LIBMESH_ENABLE_AMR
      // We turn off the asymmetric constraint application;
      // enforce_constraints_exactly() should be called in the solver
      if (get_residual && get_jacobian)
        this->get_dof_map().constrain_element_matrix_and_vector
          (elem_jacobian, elem_residual, dof_indices, false);
      else if (get_residual)
        this->get_dof_map().constrain_element_vector
          (elem_residual, dof_indices, false);
      else if (get_jacobian)
        this->get_dof_map().constrain_element_matrix
          (elem_jacobian, dof_indices, false);
#endif // #ifdef LIBMESH_ENABLE_AMR

      if (get_jacobian && print_element_jacobians)
        {
          unsigned int old_precision = std::cout.precision();
          std::cout.precision(16);
	  std::cout << "J_elem " << elem->id() << " = "
                    << elem_jacobian << std::endl;
          std::cout.precision(old_precision);
        }

      if (get_jacobian)
        this->matrix->add_matrix (elem_jacobian, dof_indices);
      if (get_residual)
        this->rhs->add_vector (elem_residual, dof_indices);
    }


  if (get_residual && print_residual_norms)
    {
      this->rhs->close();
      std::cout << "|F| = " << this->rhs->l1_norm() << std::endl;
    }
  if (get_residual && print_residuals)
    {
      this->rhs->close();
      unsigned int old_precision = std::cout.precision();
      std::cout.precision(16);
      std::cout << "F = [" << *(this->rhs) << "];" << std::endl;
      std::cout.precision(old_precision);
    }
  if (get_jacobian && print_jacobian_norms)
    {
      this->matrix->close();
      std::cout << "|J| = " << this->matrix->l1_norm() << std::endl;
    }
  if (get_jacobian && print_jacobians)
    {
      this->matrix->close();
      unsigned int old_precision = std::cout.precision();
      std::cout.precision(16);
      std::cout << "J = [" << *(this->matrix) << "];" << std::endl;
      std::cout.precision(old_precision);
    }
  STOP_LOG(log_name, "FEMSystem");
}



void FEMSystem::postprocess ()
{
  START_LOG("postprocess()", "FEMSystem");

  const MeshBase& mesh = this->get_mesh();

  const DofMap& dof_map = this->get_dof_map();

//  this->get_vector("_nonlinear_solution").localize
//    (*current_local_nonlinear_solution,
//     dof_map.get_send_list());
  this->update();

  // Loop over every active mesh element on this processor
  MeshBase::const_element_iterator el =
    mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =
    mesh.active_local_elements_end();

  for ( ; el != end_el; ++el)
    {
      elem = *el;

      // Initialize the per-element data for elem.
      dof_map.dof_indices (elem, dof_indices);
      unsigned int n_dofs = dof_indices.size();

      elem_solution.resize(n_dofs);
      // This resize call also zeros out the residual
      elem_residual.resize(n_dofs);

      for (unsigned int i=0; i != n_dofs; ++i)
//        elem_solution(i) = current_nonlinear_solution(dof_indices[i]);
        elem_solution(i) = current_solution(dof_indices[i]);

      // Initialize the per-variable data for elem.
      unsigned int sub_dofs = 0;
      for (unsigned int i=0; i != this->n_vars(); ++i)
        {
          dof_map.dof_indices (elem, dof_indices_var[i], i);

          elem_subsolutions[i]->reposition
            (sub_dofs, dof_indices_var[i].size());

          if (use_fixed_solution)
            elem_fixed_subsolutions[i]->reposition
              (sub_dofs, dof_indices_var[i].size());

          elem_subresiduals[i]->reposition
            (sub_dofs, dof_indices_var[i].size());

          sub_dofs += dof_indices_var[i].size();
        }
      libmesh_assert(sub_dofs == n_dofs);

      // Optionally initialize all the interior FE objects on elem.
      if (fe_reinit_during_postprocess)
        this->elem_fe_reinit();
      
      this->element_postprocess();

      for (side = 0; side != elem->n_sides(); ++side)
        {
          // Don't compute on non-boundary sides unless requested
          if (!postprocess_sides ||
              (!compute_internal_sides &&
               elem->neighbor(side) != NULL))
            continue;

          // Optionally initialize all the interior FE objects on elem/side.
          // Logging of FE::reinit is done in the FE functions
          if (fe_reinit_during_postprocess)
            {
              std::map<FEType, FEBase *>::iterator fe_end = element_fe.end();
              fe_end = side_fe.end();
              for (std::map<FEType, FEBase *>::iterator i = side_fe.begin();
                   i != fe_end; ++i)
                {
                  i->second->reinit(elem, side);
                }
            }

          this->side_postprocess();
        }
    }

  STOP_LOG("postprocess()", "FEMSystem");
}



void FEMSystem::numerical_elem_jacobian ()
{
  START_LOG("numerical_elem_jacobian()", "FEMSystem");

  DenseVector<Number> original_residual(elem_residual);
  DenseVector<Number> backwards_residual(elem_residual);
  DenseMatrix<Number> numerical_jacobian(elem_jacobian);

  for (unsigned int j = 0; j != dof_indices.size(); ++j)
    {
      // Take the "minus" side of a central differenced first derivative
      Number original_solution = elem_solution(j);
      elem_solution(j) -= numerical_jacobian_h;

      // Make sure to catch any moving mesh terms
      // FIXME - this could be less ugly
      Real *coord = NULL;
      if (_mesh_sys == this->number())
        {
          if (_mesh_x_var != libMesh::invalid_uint)
            for (unsigned int k = 0;
                 k != dof_indices_var[_mesh_x_var].size(); ++k)
              if (dof_indices_var[_mesh_x_var][k] == dof_indices[j])
                coord = &(elem->point(k)(0));
          if (_mesh_y_var != libMesh::invalid_uint)
            for (unsigned int k = 0;
                 k != dof_indices_var[_mesh_y_var].size(); ++k)
              if (dof_indices_var[_mesh_y_var][k] == dof_indices[j])
                coord = &(elem->point(k)(1));
          if (_mesh_z_var != libMesh::invalid_uint)
            for (unsigned int k = 0;
                 k != dof_indices_var[_mesh_z_var].size(); ++k)
              if (dof_indices_var[_mesh_z_var][k] == dof_indices[j])
                coord = &(elem->point(k)(2));
        }
      if (coord)
        *coord = libmesh_real(elem_solution(j));

      elem_residual.zero();
      time_solver->element_residual(false);
      backwards_residual = elem_residual;

      // Take the "plus" side of a central differenced first derivative
      elem_solution(j) = original_solution + numerical_jacobian_h;
      if (coord)
        *coord = libmesh_real(elem_solution(j));
      elem_residual.zero();
      time_solver->element_residual(false);

      for (unsigned int i = 0; i != dof_indices.size(); ++i)
        {
          numerical_jacobian(i,j) =
            (elem_residual(i) - backwards_residual(i)) /
            2. / numerical_jacobian_h;
        }
      elem_solution(j) = original_solution;
      if (coord)
        *coord = libmesh_real(elem_solution(j));
    }

  elem_residual = original_residual;
  elem_jacobian = numerical_jacobian;
  STOP_LOG("numerical_elem_jacobian()", "FEMSystem");
}



void FEMSystem::numerical_side_jacobian ()
{
  START_LOG("numerical_side_jacobian()", "FEMSystem");

  DenseVector<Number> original_residual(elem_residual);
  DenseVector<Number> backwards_residual(elem_residual);
  DenseMatrix<Number> numerical_jacobian(elem_jacobian);

  for (unsigned int j = 0; j != dof_indices.size(); ++j)
    {
      // Take the "minus" side of a central differenced first derivative
      Number original_solution = elem_solution(j);
      elem_solution(j) -= numerical_jacobian_h;

      // Make sure to catch any moving mesh terms
      // FIXME - this could be less ugly
      Real *coord = NULL;
      if (_mesh_sys == this->number())
        {
          if (_mesh_x_var != libMesh::invalid_uint)
            for (unsigned int k = 0;
                 k != dof_indices_var[_mesh_x_var].size(); ++k)
              if (dof_indices_var[_mesh_x_var][k] == dof_indices[j])
                coord = &(elem->point(k)(0));
          if (_mesh_y_var != libMesh::invalid_uint)
            for (unsigned int k = 0;
                 k != dof_indices_var[_mesh_y_var].size(); ++k)
              if (dof_indices_var[_mesh_y_var][k] == dof_indices[j])
                coord = &(elem->point(k)(1));
          if (_mesh_z_var != libMesh::invalid_uint)
            for (unsigned int k = 0;
                 k != dof_indices_var[_mesh_z_var].size(); ++k)
              if (dof_indices_var[_mesh_z_var][k] == dof_indices[j])
                coord = &(elem->point(k)(2));
        }
      if (coord)
        *coord = libmesh_real(elem_solution(j));

      elem_residual.zero();
      time_solver->side_residual(false);
      backwards_residual = elem_residual;

      // Take the "plus" side of a central differenced first derivative
      elem_solution(j) = original_solution + numerical_jacobian_h;
      if (coord)
        *coord = libmesh_real(elem_solution(j));

      elem_residual.zero();
      time_solver->side_residual(false);

      for (unsigned int i = 0; i != dof_indices.size(); ++i)
        {
          numerical_jacobian(i,j) +=
            (elem_residual(i) - backwards_residual(i))
            / 2. / numerical_jacobian_h;
        }
      elem_solution(j) = original_solution;
      if (coord)
        *coord = libmesh_real(elem_solution(j));
    }

  elem_residual = original_residual;
  elem_jacobian = numerical_jacobian;
  STOP_LOG("numerical_side_jacobian()", "FEMSystem");
}



void FEMSystem::time_evolving (unsigned int var)
{
  // Call the parent function
  Parent::time_evolving(var);

  // Then make sure we're prepared to do mass integration
  element_fe_var[var]->get_JxW();
  element_fe_var[var]->get_phi();
}



void FEMSystem::mesh_x_position (unsigned int sysnum, unsigned int var)
{
  if (_mesh_sys != libMesh::invalid_uint && _mesh_sys != sysnum)
    libmesh_error();
  if (sysnum != this->number())
    libmesh_not_implemented();
  _mesh_sys = sysnum;
  _mesh_x_var = var;
}



void FEMSystem::mesh_y_position (unsigned int sysnum, unsigned int var)
{
  if (_mesh_sys != libMesh::invalid_uint && _mesh_sys != sysnum)
    libmesh_error();
  if (sysnum != this->number())
    libmesh_not_implemented();
  _mesh_sys = sysnum;
  _mesh_y_var = var;
}



void FEMSystem::mesh_z_position (unsigned int sysnum, unsigned int var)
{
  if (_mesh_sys != libMesh::invalid_uint && _mesh_sys != sysnum)
    libmesh_error();
  if (sysnum != this->number())
    libmesh_not_implemented();
  _mesh_sys = sysnum;
  _mesh_z_var = var;
}



void FEMSystem::initialize_mesh_variables ()
{
  // This function makes no sense unless we've already picked out some
  // variable(s) to reflect mesh position coordinates
  if (_mesh_sys == libMesh::invalid_uint)
    libmesh_error();

  // We currently assume mesh variables are in our own system
  if (_mesh_sys != this->number())
    libmesh_not_implemented();

  // Loop over every active mesh element on this processor
  const MeshBase& mesh = this->get_mesh();

  const DofMap& dof_map = this->get_dof_map();

  MeshBase::const_element_iterator el =
    mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =
    mesh.active_local_elements_end();

  for ( ; el != end_el; ++el)
    {
      elem = *el;

      // Initialize the per-variable data for elem.
      for (unsigned int i=0; i != this->n_vars(); ++i)
        {
          dof_map.dof_indices (elem, dof_indices_var[i], i);
        }

      this->elem_position_get();
    }
}



void FEMSystem::elem_position_get()
{
  // This is too expensive to call unless we've been asked to move the mesh
  libmesh_assert (_mesh_sys != libMesh::invalid_uint);

  // If the coordinate data is in our own system, it's already
  // been set up for us
  if (_mesh_sys == this->number())
    {
      unsigned int n_nodes = elem->n_nodes();
      // For simplicity we demand that mesh coordinates be stored
      // in a format that allows a direct copy
      libmesh_assert(_mesh_x_var == libMesh::invalid_uint ||
                     (element_fe_var[_mesh_x_var]->get_fe_type().family
                      == LAGRANGE &&
                      elem_subsolutions[_mesh_x_var]->size() == n_nodes));
      libmesh_assert(_mesh_y_var == libMesh::invalid_uint ||
                     (element_fe_var[_mesh_y_var]->get_fe_type().family
                      == LAGRANGE &&
                      elem_subsolutions[_mesh_y_var]->size() == n_nodes));
      libmesh_assert(_mesh_z_var == libMesh::invalid_uint ||
                     (element_fe_var[_mesh_z_var]->get_fe_type().family
                      == LAGRANGE &&
                      elem_subsolutions[_mesh_z_var]->size() == n_nodes));

      // Get degree of freedom coefficients from point coordinates
      if (_mesh_x_var != libMesh::invalid_uint)
        for (unsigned int i=0; i != n_nodes; ++i)
          this->solution->set(dof_indices_var[_mesh_x_var][i],
                              elem->point(i)(0));

      if (_mesh_y_var != libMesh::invalid_uint)
        for (unsigned int i=0; i != n_nodes; ++i)
          this->solution->set(dof_indices_var[_mesh_y_var][i],
                             elem->point(i)(1));

      if (_mesh_z_var != libMesh::invalid_uint)
        for (unsigned int i=0; i != n_nodes; ++i)
          this->solution->set(dof_indices_var[_mesh_z_var][i],
                              elem->point(i)(2));
    }
  // FIXME - If the coordinate data is not in our own system, someone
  // had better get around to implementing that... - RHS
  else
    {
      libmesh_not_implemented();
    }
}



void FEMSystem::elem_position_set(Real)
{
  // This is too expensive to call unless we've been asked to move the mesh
  libmesh_assert (_mesh_sys != libMesh::invalid_uint);

  // If the coordinate data is in our own system, it's already
  // been set up for us, and we can ignore our input parameter theta
  if (_mesh_sys == this->number())
    {
      unsigned int n_nodes = elem->n_nodes();
      // For simplicity we demand that mesh coordinates be stored
      // in a format that allows a direct copy
      libmesh_assert(_mesh_x_var == libMesh::invalid_uint ||
                     (element_fe_var[_mesh_x_var]->get_fe_type().family
                      == LAGRANGE &&
                      elem_subsolutions[_mesh_x_var]->size() == n_nodes));
      libmesh_assert(_mesh_y_var == libMesh::invalid_uint ||
                     (element_fe_var[_mesh_y_var]->get_fe_type().family
                      == LAGRANGE &&
                      elem_subsolutions[_mesh_y_var]->size() == n_nodes));
      libmesh_assert(_mesh_z_var == libMesh::invalid_uint ||
                     (element_fe_var[_mesh_z_var]->get_fe_type().family
                      == LAGRANGE &&
                      elem_subsolutions[_mesh_z_var]->size() == n_nodes));

      // Set the new point coordinates
      if (_mesh_x_var != libMesh::invalid_uint)
        for (unsigned int i=0; i != n_nodes; ++i)
          elem->point(i)(0) =
            libmesh_real((*elem_subsolutions[_mesh_x_var])(i));

      if (_mesh_y_var != libMesh::invalid_uint)
        for (unsigned int i=0; i != n_nodes; ++i)
          elem->point(i)(1) =
            libmesh_real((*elem_subsolutions[_mesh_y_var])(i));

      if (_mesh_z_var != libMesh::invalid_uint)
        for (unsigned int i=0; i != n_nodes; ++i)
          elem->point(i)(2) =
            libmesh_real((*elem_subsolutions[_mesh_z_var])(i));
    }
  // FIXME - If the coordinate data is not in our own system, someone
  // had better get around to implementing that... - RHS
  else
    {
      libmesh_not_implemented();
    }
}



bool FEMSystem::eulerian_residual (bool request_jacobian)
{
  // Only calculate a mesh movement residual if it's necessary
  if (!_mesh_sys)
    return request_jacobian;

  // This function only supports fully coupled mesh motion for now
  libmesh_assert(_mesh_sys == this->number());

  unsigned int n_qpoints = element_qrule->n_points();

  const unsigned int n_x_dofs = (_mesh_x_var == libMesh::invalid_uint) ?
                                0 : dof_indices_var[_mesh_x_var].size();
  const unsigned int n_y_dofs = (_mesh_y_var == libMesh::invalid_uint) ?
                                0 : dof_indices_var[_mesh_y_var].size();
  const unsigned int n_z_dofs = (_mesh_z_var == libMesh::invalid_uint) ?
                                0 : dof_indices_var[_mesh_z_var].size();

  const unsigned int mesh_xyz_var = n_x_dofs ? _mesh_x_var :
                                   (n_y_dofs ? _mesh_y_var :
                                   (n_z_dofs ? _mesh_z_var :
                                    libMesh::invalid_uint));

  // If we're our own _mesh_sys, we'd better be in charge of
  // at least one coordinate, and we'd better have the same
  // FE type for all coordinates we are in charge of
  libmesh_assert(mesh_xyz_var != libMesh::invalid_uint);
  libmesh_assert(!n_x_dofs || element_fe_var[_mesh_x_var] ==
                              element_fe_var[mesh_xyz_var]);
  libmesh_assert(!n_y_dofs || element_fe_var[_mesh_y_var] ==
                              element_fe_var[mesh_xyz_var]);
  libmesh_assert(!n_z_dofs || element_fe_var[_mesh_z_var] ==
                              element_fe_var[mesh_xyz_var]);

  const std::vector<std::vector<Real> >     &psi =
    element_fe_var[mesh_xyz_var]->get_phi();

  for (unsigned int var = 0; var != this->n_vars(); ++var)
    {
      // Mesh motion only affects time-evolving variables
      if (!_time_evolving[var])
        continue;

      // The mesh coordinate variables themselves are Lagrangian,
      // not Eulerian, and no convective term is desired.
      if (_mesh_sys == this->number() &&
          (var == _mesh_x_var ||
           var == _mesh_y_var ||
           var == _mesh_z_var))
        continue;

      // Some of this code currently relies on the assumption that
      // we can pull mesh coordinate data from our own system
      if (_mesh_sys != this->number())
        libmesh_not_implemented();

      // This residual should only be called by unsteady solvers:
      // if the mesh is steady, there's no mesh convection term!
      UnsteadySolver *unsteady =
        dynamic_cast<UnsteadySolver *>(this->time_solver.get());
      if (!unsteady)
        return request_jacobian;

      const std::vector<Real> &JxW = 
        element_fe_var[var]->get_JxW();

      const std::vector<std::vector<Real> >     &phi =
        element_fe_var[var]->get_phi();

      const std::vector<std::vector<RealGradient> > &dphi =
        element_fe_var[var]->get_dphi();

      const unsigned int n_u_dofs = dof_indices_var[var].size();

      DenseSubVector<Number> &Fu = *elem_subresiduals[var];
      DenseSubMatrix<Number> &Kuu = *elem_subjacobians[var][var];

      DenseSubMatrix<Number> *Kux = n_x_dofs ?
        elem_subjacobians[var][_mesh_x_var] : NULL;
      DenseSubMatrix<Number> *Kuy = n_y_dofs ?
        elem_subjacobians[var][_mesh_y_var] : NULL;
      DenseSubMatrix<Number> *Kuz = n_z_dofs ?
        elem_subjacobians[var][_mesh_z_var] : NULL;

      std::vector<Real> delta_x(n_x_dofs, 0.);
      std::vector<Real> delta_y(n_y_dofs, 0.);
      std::vector<Real> delta_z(n_z_dofs, 0.);

      for (unsigned int i = 0; i != n_x_dofs; ++i)
        {
          unsigned int j = dof_indices_var[_mesh_x_var][i];
          delta_x[i] = libmesh_real(this->current_solution(j)) -
                       libmesh_real(unsteady->old_nonlinear_solution(j));
        }

      for (unsigned int i = 0; i != n_y_dofs; ++i)
        {
          unsigned int j = dof_indices_var[_mesh_y_var][i];
          delta_y[i] = libmesh_real(this->current_solution(j)) -
                       libmesh_real(unsteady->old_nonlinear_solution(j));
        }

      for (unsigned int i = 0; i != n_z_dofs; ++i)
        {
          unsigned int j = dof_indices_var[_mesh_z_var][i];
          delta_z[i] = libmesh_real(this->current_solution(j)) -
                       libmesh_real(unsteady->old_nonlinear_solution(j));
        }

      for (unsigned int qp = 0; qp != n_qpoints; ++qp)
        {
          Gradient grad_u = interior_gradient(var, qp);
          RealGradient convection(0.);

          libmesh_error();
          for (unsigned int i = 0; i != n_x_dofs; ++i)
	    convection(0) += delta_x[i] * phi[i][qp];
          for (unsigned int i = 0; i != n_y_dofs; ++i)
	    convection(1) += delta_y[i] * phi[i][qp];
          for (unsigned int i = 0; i != n_z_dofs; ++i)
	    convection(2) += delta_z[i] * phi[i][qp];

          for (unsigned int i = 0; i != n_u_dofs; ++i)
            {
              Number JxWxPhiI = JxW[qp] * phi[i][qp];
              Fu(i) += (convection * grad_u) * JxWxPhiI;
              if (request_jacobian)
                {
                  Number JxWxPhiI = JxW[qp] * phi[i][qp];
                  for (unsigned int j = 0; j != n_u_dofs; ++j)
                    Kuu(i,j) += JxWxPhiI * (convection * dphi[j][qp]);

                  Number JxWxPhiIoverDT = JxWxPhiI/this->deltat;

                  Number JxWxPhiIxDUDXoverDT = JxWxPhiIoverDT * grad_u(0);
                  for (unsigned int j = 0; j != n_x_dofs; ++j)
                    (*Kux)(i,j) += JxWxPhiIxDUDXoverDT * psi[j][qp];

                  Number JxWxPhiIxDUDYoverDT = JxWxPhiIoverDT * grad_u(1);
                  for (unsigned int j = 0; j != n_y_dofs; ++j)
                    (*Kuy)(i,j) += JxWxPhiIxDUDYoverDT * psi[j][qp];

                  Number JxWxPhiIxDUDZoverDT = JxWxPhiIoverDT * grad_u(2);
                  for (unsigned int j = 0; j != n_z_dofs; ++j)
                    (*Kuz)(i,j) += JxWxPhiIxDUDZoverDT * psi[j][qp];
                }
            }
        }
    }

  return request_jacobian;
}

      

bool FEMSystem::mass_residual (bool request_jacobian)
{
  unsigned int n_qpoints = element_qrule->n_points();

  for (unsigned int var = 0; var != this->n_vars(); ++var)
    {
      if (!_time_evolving[var])
        continue;

      const std::vector<Real> &JxW = 
        element_fe_var[var]->get_JxW();

      const std::vector<std::vector<Real> > &phi =
        element_fe_var[var]->get_phi();

      const unsigned int n_dofs = dof_indices_var[var].size();

      DenseSubVector<Number> &Fu = *elem_subresiduals[var];
      DenseSubMatrix<Number> &Kuu = *elem_subjacobians[var][var];

      for (unsigned int qp = 0; qp != n_qpoints; ++qp)
        {
          Number u = interior_value(var, qp);
          Number JxWxU = JxW[qp] * u;
          for (unsigned int i = 0; i != n_dofs; ++i)
            {
              Fu(i) += JxWxU * phi[i][qp];
              if (request_jacobian && elem_solution_derivative)
                {
                  libmesh_assert (elem_solution_derivative == 1.0);

                  Number JxWxPhiI = JxW[qp] * phi[i][qp];
                  Kuu(i,i) += JxWxPhiI * phi[i][qp];
                  for (unsigned int j = i+1; j < n_dofs; ++j)
                    {
                      Number Kij = JxWxPhiI * phi[j][qp];
                      Kuu(i,j) += Kij;
                      Kuu(j,i) += Kij;
                    }
                }
            }
        }
    }

  return request_jacobian;
}
