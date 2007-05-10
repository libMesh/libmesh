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



FEMSystem::FEMSystem (EquationSystems& es,
                      const std::string& name,
                      const unsigned int number)
  : Parent(es, name, number),
    fe_reinit_during_postprocess(true),
    extra_quadrature_order(0),
    numerical_jacobian_h(1.e-6),
    verify_analytic_jacobians(0.0),
    element_qrule(NULL), side_qrule(NULL), elem(NULL)
{
}


FEMSystem::~FEMSystem ()
{
  this->clear();
}



Number FEMSystem::interior_value(unsigned int var, unsigned int qp)
{
  // Get local-to-global dof index lookup
  assert (dof_indices.size() > var);
  const unsigned int n_dofs = dof_indices_var[var].size();

  // Get current local coefficients
  assert (elem_subsolutions.size() > var);
  assert (elem_subsolutions[var] != NULL);
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
  assert (dof_indices.size() > var);
  const unsigned int n_dofs = dof_indices_var[var].size();

  // Get current local coefficients
  assert (elem_subsolutions.size() > var);
  assert (elem_subsolutions[var] != NULL);
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



#ifdef ENABLE_SECOND_DERIVATIVES
Tensor FEMSystem::interior_hessian(unsigned int var, unsigned int qp)
{
  // Get local-to-global dof index lookup
  assert (dof_indices.size() > var);
  const unsigned int n_dofs = dof_indices_var[var].size();

  // Get current local coefficients
  assert (elem_subsolutions.size() > var);
  assert (elem_subsolutions[var] != NULL);
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
#endif // ifdef ENABLE_SECOND_DERIVATIVES



Number FEMSystem::side_value(unsigned int var, unsigned int qp)
{
  // Get local-to-global dof index lookup
  assert (dof_indices.size() > var);
  const unsigned int n_dofs = dof_indices_var[var].size();

  // Get current local coefficients
  assert (elem_subsolutions.size() > var);
  assert (elem_subsolutions[var] != NULL);
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
  assert (dof_indices.size() > var);
  const unsigned int n_dofs = dof_indices_var[var].size();

  // Get current local coefficients
  assert (elem_subsolutions.size() > var);
  assert (elem_subsolutions[var] != NULL);
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



#ifdef ENABLE_SECOND_DERIVATIVES
Tensor FEMSystem::side_hessian(unsigned int var, unsigned int qp)
{
  // Get local-to-global dof index lookup
  assert (dof_indices.size() > var);
  const unsigned int n_dofs = dof_indices_var[var].size();

  // Get current local coefficients
  assert (elem_subsolutions.size() > var);
  assert (elem_subsolutions[var] != NULL);
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
#endif // ifdef ENABLE_SECOND_DERIVATIVES



Number FEMSystem::point_value(unsigned int var, Point &p)
{
  // Get local-to-global dof index lookup
  assert (dof_indices.size() > var);
  const unsigned int n_dofs = dof_indices_var[var].size();

  
  // Get current local coefficients
  assert (elem_subsolutions.size() > var);
  assert (elem_subsolutions[var] != NULL);
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



void FEMSystem::clear()
{
  Parent::clear();

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
  this->clear();

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
      // different families which require different quadrature rules
      // assert (fe_type.family == hardest_fe_type.family);

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



void FEMSystem::assembly (bool get_residual, bool get_jacobian)
{
  assert(get_residual || get_jacobian);
  std::string log_name;
  if (get_residual && get_jacobian)
    log_name = "assembly()";
  else if (get_residual)
    log_name = "assembly(get_residual)";
  else
    log_name = "assembly(get_jacobian)";

  START_LOG(log_name, "FEMSystem");

  const Mesh& mesh = this->get_mesh();

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
  assert (time_solver.get() != NULL);

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

      // Initialize all the interior FE objects on elem.
      // Logging of FE::reinit is done in the FE functions
      PAUSE_LOG(log_name, "FEMSystem");
      std::map<FEType, FEBase *>::iterator fe_end = element_fe.end();
      for (std::map<FEType, FEBase *>::iterator i = element_fe.begin();
           i != fe_end; ++i)
        {
          i->second->reinit(elem);
        }
      RESTART_LOG(log_name, "FEMSystem");
      
      bool jacobian_computed = time_solver->element_residual(get_jacobian);

      // Compute a numeric jacobian if we have to
      if (get_jacobian && !jacobian_computed)
        {
          // Make sure we didn't compute a jacobian and lie about it
          assert(elem_jacobian.l1_norm() == 0.0);
          // Logging of numerical jacobians is done separately
          PAUSE_LOG(log_name, "FEMSystem");
          this->numerical_elem_jacobian();
          RESTART_LOG(log_name, "FEMSystem");
        }

      // Compute a numeric jacobian if we're asked to verify the
      // analytic jacobian we got
      if (get_jacobian && jacobian_computed &&
          verify_analytic_jacobians != 0.0)
        {
          DenseMatrix<Number> analytic_jacobian(elem_jacobian);

          elem_jacobian.zero();
          // Logging of numerical jacobians is done separately
          PAUSE_LOG(log_name, "FEMSystem");
          this->numerical_elem_jacobian();
          RESTART_LOG(log_name, "FEMSystem");

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

              error();
            }
        }

      for (side = 0; side != elem->n_sides(); ++side)
        {
          // Don't compute on non-boundary sides unless requested
          if (!compute_internal_sides && elem->neighbor(side) != NULL)
            continue;

          // Initialize all the interior FE objects on elem/side.
          // Logging of FE::reinit is done in the FE functions
          PAUSE_LOG(log_name, "FEMSystem");
          fe_end = side_fe.end();
          for (std::map<FEType, FEBase *>::iterator i = side_fe.begin();
               i != fe_end; ++i)
            {
              i->second->reinit(elem, side);
            }
          RESTART_LOG(log_name, "FEMSystem");

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
              assert(elem_jacobian.l1_norm() == 0.0);
#endif
              // Logging of numerical jacobians is done separately
              PAUSE_LOG(log_name, "FEMSystem");
              this->numerical_side_jacobian();
              RESTART_LOG(log_name, "FEMSystem");

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
              PAUSE_LOG(log_name, "FEMSystem");
              this->numerical_side_jacobian();
              RESTART_LOG(log_name, "FEMSystem");

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

                  error();
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

  const Mesh& mesh = this->get_mesh();

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

          sub_dofs += dof_indices_var[i].size();
        }
      assert(sub_dofs == n_dofs);

      // Optionally initialize all the interior FE objects on elem.
      // Logging of FE::reinit is done in the FE functions
      if (fe_reinit_during_postprocess)
        {
          PAUSE_LOG("postprocess()", "FEMSystem");
          std::map<FEType, FEBase *>::iterator fe_end = element_fe.end();
          for (std::map<FEType, FEBase *>::iterator i = element_fe.begin();
               i != fe_end; ++i)
            {
              i->second->reinit(elem);
            }
          RESTART_LOG("postprocess()", "FEMSystem");
        }
      
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
              PAUSE_LOG("postprocess()", "FEMSystem");
              std::map<FEType, FEBase *>::iterator fe_end = element_fe.end();
              fe_end = side_fe.end();
              for (std::map<FEType, FEBase *>::iterator i = side_fe.begin();
                   i != fe_end; ++i)
                {
                  i->second->reinit(elem, side);
                }
              RESTART_LOG("postprocess()", "FEMSystem");
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
      Number original_solution = elem_solution(j);
      elem_solution(j) -= numerical_jacobian_h;
      elem_residual.zero();
      time_solver->element_residual(false);

      backwards_residual = elem_residual;
      elem_solution(j) = original_solution + numerical_jacobian_h;
      elem_residual.zero();
      time_solver->element_residual(false);

      for (unsigned int i = 0; i != dof_indices.size(); ++i)
        {
          numerical_jacobian(i,j) =
            (elem_residual(i) - backwards_residual(i)) /
            2. / numerical_jacobian_h;
        }
      elem_solution(j) = original_solution;
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
      Number original_solution = elem_solution(j);
      elem_solution(j) -= numerical_jacobian_h;
      elem_residual.zero();
      time_solver->side_residual(false);

      backwards_residual = elem_residual;
      elem_solution(j) = original_solution + numerical_jacobian_h;
      elem_residual.zero();
      time_solver->side_residual(false);

      for (unsigned int i = 0; i != dof_indices.size(); ++i)
        {
          numerical_jacobian(i,j) +=
            (elem_residual(i) - backwards_residual(i))
            / 2. / numerical_jacobian_h;
        }
      elem_solution(j) = original_solution;
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
                  assert (elem_solution_derivative == 1.0);

                  Number JxWxPhiI = JxW[qp] * phi[i][qp];
                  Kuu(i,i) += JxWxPhiI * phi[i][qp];
                  for (unsigned int j = i+1; j != n_dofs; ++j)
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
