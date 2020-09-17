// The libMesh Finite Element Library.
// Copyright (C) 2002-2020 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


#include "libmesh/diff_system.h"
#include "libmesh/euler2_solver.h"

#include "libmesh/error_vector.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/adjoint_refinement_estimator.h"

namespace libMesh
{



Euler2Solver::Euler2Solver (sys_type & s)
  : FirstOrderUnsteadySolver(s), theta(1.)
{
}



Euler2Solver::~Euler2Solver ()
{
}



Real Euler2Solver::error_order() const
{
  if (theta == 0.5)
    return 2.;
  return 1.;
}




bool Euler2Solver::element_residual (bool request_jacobian,
                                     DiffContext & context)
{
  bool compute_second_order_eqns = !this->_system.get_physics()->get_second_order_vars().empty();

  return this->_general_residual(request_jacobian,
                                 context,
                                 &DifferentiablePhysics::mass_residual,
                                 &DifferentiablePhysics::damping_residual,
                                 &DifferentiablePhysics::_eulerian_time_deriv,
                                 &DifferentiablePhysics::element_constraint,
                                 &DiffContext::elem_reinit,
                                 compute_second_order_eqns);
}



bool Euler2Solver::side_residual (bool request_jacobian,
                                  DiffContext & context)
{
  return this->_general_residual(request_jacobian,
                                 context,
                                 &DifferentiablePhysics::side_mass_residual,
                                 &DifferentiablePhysics::side_damping_residual,
                                 &DifferentiablePhysics::side_time_derivative,
                                 &DifferentiablePhysics::side_constraint,
                                 &DiffContext::elem_side_reinit,
                                 false);
}



bool Euler2Solver::nonlocal_residual (bool request_jacobian,
                                      DiffContext & context)
{
  bool compute_second_order_eqns = this->_system.have_second_order_scalar_vars();

  return this->_general_residual(request_jacobian,
                                 context,
                                 &DifferentiablePhysics::nonlocal_mass_residual,
                                 &DifferentiablePhysics::nonlocal_damping_residual,
                                 &DifferentiablePhysics::nonlocal_time_derivative,
                                 &DifferentiablePhysics::nonlocal_constraint,
                                 &DiffContext::nonlocal_reinit,
                                 compute_second_order_eqns);
}



bool Euler2Solver::_general_residual (bool request_jacobian,
                                      DiffContext & context,
                                      ResFuncType mass,
                                      ResFuncType damping,
                                      ResFuncType time_deriv,
                                      ResFuncType constraint,
                                      ReinitFuncType reinit_func,
                                      bool compute_second_order_eqns)
{
  unsigned int n_dofs = context.get_elem_solution().size();

  // Local nonlinear solution at old timestep
  DenseVector<Number> old_elem_solution(n_dofs);
  for (unsigned int i=0; i != n_dofs; ++i)
    old_elem_solution(i) =
      old_nonlinear_solution(context.get_dof_indices()[i]);

  // Local time derivative of solution
  context.get_elem_solution_rate() = context.get_elem_solution();
  context.get_elem_solution_rate() -= old_elem_solution;
  context.elem_solution_rate_derivative = 1 / _system.deltat;
  context.get_elem_solution_rate() *=
    context.elem_solution_rate_derivative;

  // Our first evaluations are at the final elem_solution
  context.elem_solution_derivative = 1.0;

  // If a fixed solution is requested, we'll use the elem_solution
  // at the new timestep
  // FIXME - should this be the theta solution instead?
  if (_system.use_fixed_solution)
    context.get_elem_fixed_solution() = context.get_elem_solution();

  context.fixed_solution_derivative = 1.0;

  // We need to save the old jacobian and old residual since we'll be
  // multiplying some of the new contributions by theta or 1-theta
  DenseMatrix<Number> old_elem_jacobian(n_dofs, n_dofs);
  DenseVector<Number> old_elem_residual(n_dofs);
  old_elem_residual.swap(context.get_elem_residual());
  if (request_jacobian)
    old_elem_jacobian.swap(context.get_elem_jacobian());

  // Local time derivative of solution
  context.get_elem_solution_rate() = context.get_elem_solution();
  context.get_elem_solution_rate() -= old_elem_solution;
  context.elem_solution_rate_derivative = 1 / _system.deltat;
  context.get_elem_solution_rate() *=
    context.elem_solution_rate_derivative;

  // If we are asked to compute residuals for second order variables,
  // we also populate the acceleration part so the user can use that.
  if (compute_second_order_eqns)
    this->prepare_accel(context);

  // Move the mesh into place first if necessary, set t = t_{n+1}
  (context.*reinit_func)(1.);

  // First, evaluate time derivative at the new timestep.
  // The element should already be in the proper place
  // even for a moving mesh problem.
  bool jacobian_computed =
    (_system.get_physics()->*time_deriv)(request_jacobian, context);

  // Next, evaluate the mass residual at the new timestep

  jacobian_computed = (_system.get_physics()->*mass)(jacobian_computed, context) &&
    jacobian_computed;

  // If we have second-order variables, we need to get damping terms
  // and the velocity equations
  if (compute_second_order_eqns)
    {
      jacobian_computed = (_system.get_physics()->*damping)(jacobian_computed, context) &&
        jacobian_computed;

      jacobian_computed = this->compute_second_order_eqns(jacobian_computed, context) &&
        jacobian_computed;
    }

  // Add the constraint term
  jacobian_computed = (_system.get_physics()->*constraint)(jacobian_computed, context) &&
    jacobian_computed;

  // The new solution's contribution is scaled by theta
  context.get_elem_residual() *= theta;
  context.get_elem_jacobian() *= theta;

  // Save the new solution's term
  DenseMatrix<Number> elem_jacobian_newterm(n_dofs, n_dofs);
  DenseVector<Number> elem_residual_newterm(n_dofs);
  elem_residual_newterm.swap(context.get_elem_residual());
  if (request_jacobian)
    elem_jacobian_newterm.swap(context.get_elem_jacobian());

  // Add the time-dependent term for the old solution

  // Make sure elem_solution is set up for elem_reinit to use
  // Move elem_->old_, old_->elem_
  context.get_elem_solution().swap(old_elem_solution);
  context.elem_solution_derivative = 0.0;

  // Move the mesh into place if necessary, set t = t_{n}
  (context.*reinit_func)(0.);

  jacobian_computed =
    (_system.get_physics()->*time_deriv)(jacobian_computed, context) &&
    jacobian_computed;

  // Add the mass residual term for the old solution

  // Evaluating the mass residual at both old and new timesteps will be
  // redundant in most problems but may be necessary for time accuracy
  // or stability in moving mesh problems or problems with user-overridden
  // mass_residual functions

  jacobian_computed =
    (_system.get_physics()->*mass)(jacobian_computed, context) &&
    jacobian_computed;

  // If we have second-order variables, we need to get damping terms
  // and the velocity equations
  if (compute_second_order_eqns)
    {
      jacobian_computed = (_system.get_physics()->*damping)(jacobian_computed, context) &&
        jacobian_computed;

      jacobian_computed = this->compute_second_order_eqns(jacobian_computed, context) &&
        jacobian_computed;
    }

  // The old solution's contribution is scaled by (1-theta)
  context.get_elem_residual() *= (1-theta);
  context.get_elem_jacobian() *= (1-theta);

  // Restore the elem_solution
  // Move elem_->elem_, old_->old_
  context.get_elem_solution().swap(old_elem_solution);
  context.elem_solution_derivative = 1;

  // Restore the elem position if necessary, set t = t_{n+1}
  (context.*reinit_func)(1.);

  // Add back (or restore) the old residual/jacobian
  context.get_elem_residual() += old_elem_residual;
  if (request_jacobian)
    {
      if (jacobian_computed)
        context.get_elem_jacobian() += old_elem_jacobian;
      else
        context.get_elem_jacobian().swap(old_elem_jacobian);
    }

  // Add the saved new-solution terms
  context.get_elem_residual() += elem_residual_newterm;
  if (jacobian_computed)
    context.get_elem_jacobian() += elem_jacobian_newterm;

  return jacobian_computed;
}

void Euler2Solver::integrate_qoi_timestep()
{
  // We are using the trapezoidal rule to integrate each timestep
  // (f(t_j) + f(t_j+1))/2 (t_j+1 - t_j)

  // Zero out the system.qoi vector
  for (auto j : make_range(_system.n_qois()))
  {
    (_system.qoi)[j] = 0.0;
  }

  // Left and right side contributions
  std::vector<Number> left_contribution(_system.qoi.size(), 0.0);
  Number time_left = 0.0;
  std::vector<Number> right_contribution(_system.qoi.size(), 0.0);
  Number time_right = 0.0;

  time_left = _system.time;

  // Base class assumes a direct steady evaluation
  this->_system.assemble_qoi();

  // Also get the spatially integrated errors for all the QoIs in the QoI set
  for (auto j : make_range(_system.n_qois()))
  {
    left_contribution[j] = (_system.qoi)[j];
  }

  // Advance to t_j+1
  _system.time = _system.time + _system.deltat;

  time_right = _system.time;

  // Load the solution at the next timestep
  retrieve_timestep();

  // Zero out the system.qoi vector
  for (auto j : make_range(_system.n_qois()))
  {
    (_system.qoi)[j] = 0.0;
  }

  // Base class assumes a direct steady evaluation
  this->_system.assemble_qoi();

  for(auto j : make_range(_system.n_qois()))
  {
    right_contribution[j] = (_system.qoi)[j];
  }

  // Combine the left and right side contributions as per the specified theta
  // theta = 0.5 (Crank-Nicholson) gives the trapezoidal rule.
  for (auto j : make_range(_system.n_qois()))
  {
    (_system.qoi)[j] = ( ( ((1.0 - theta)*left_contribution[j]) + (theta*right_contribution[j]) )/2.0 )*(time_right - time_left);
  }
}

void Euler2Solver::integrate_adjoint_refinement_error_estimate(AdjointRefinementEstimator & adjoint_refinement_error_estimator, ErrorVector & QoI_elementwise_error, std::vector<Real *> QoI_time_instant)
{
  // Make sure the system::qoi_error_estimates vector is of the same size as system::qoi
  if(_system.qoi_error_estimates.size() != _system.qoi.size())
    _system.qoi_error_estimates.resize(_system.qoi.size());

  // There are two possibilities regarding the integration rule we need to use for time integration.
  // If we have a instantaneous QoI, then we need to use a left sided Riemann sum, otherwise the trapezoidal rule for temporally smooth QoIs.

  // Create left and right error estimate vectors of the right size
  std::vector<Number> qoi_error_estimates_left(_system.qoi.size());
  std::vector<Number> qoi_error_estimates_right(_system.qoi.size());

  // Get t_j
  Real time_left = _system.time;

  // Get f(t_j)
  ErrorVector QoI_elementwise_error_left;

  // If we are at the very initial step, the error contribution is zero,
  // otherwise the old ajoint vector has been filled and we are the left end
  // of a subsequent timestep or sub-timestep
  if(old_adjoints[0] != nullptr)
  {
    // For evaluating the residual, we need to use the deltat that was used
    // to get us to this solution, so we save the current deltat as next_step_deltat
    // and set _system.deltat to the last completed deltat.
    next_step_deltat = _system.deltat;
    _system.deltat = last_step_deltat;

    // The adjoint error estimate expression for a backwards facing step
    // scheme needs the adjoint for the last time instant, so save the current adjoint for future use
    for (auto j : make_range(_system.n_qois()))
    {
      // Swap for residual weighting
      _system.get_adjoint_solution(j).swap(*old_adjoints[j]);
    }

    _system.update();

    // The residual has to be evaluated at the last time
    _system.time = _system.time - last_step_deltat;

    adjoint_refinement_error_estimator.estimate_error(_system, QoI_elementwise_error_left);

    // Shift the time back
    _system.time = _system.time + last_step_deltat;

    // Swap back the current and old adjoints
    for (auto j : make_range(_system.n_qois()))
    {
      _system.get_adjoint_solution(j).swap(*old_adjoints[j]);
    }

    // Set the system deltat back to what it should be to march to the next time
    _system.deltat = next_step_deltat;

  }
  else
  {
    for(unsigned int i = 0; i < QoI_elementwise_error.size(); i++)
      QoI_elementwise_error_left[i] = 0.0;
  }

  // Also get the left side contributions for the spatially integrated errors for all the QoIs in the QoI set
  for (auto j : make_range(_system.n_qois()))
  {
    // Skip this QoI if not in the QoI Set
    if (adjoint_refinement_error_estimator.qoi_set().has_index(j))
    {
      // If we are at the initial time, the error contribution is zero
      if(std::abs(_system.time) > TOLERANCE*sqrt(TOLERANCE))
      {
        qoi_error_estimates_left[j] = adjoint_refinement_error_estimator.get_global_QoI_error_estimate(j);
      }
      else
      {
        qoi_error_estimates_left[j] = 0.0;
      }
    }
  }

  // Advance to t_j+1
  _system.time = _system.time + _system.deltat;

  // Get t_j+1
  Real time_right = _system.time;

  // We will need to use the last step deltat for the weighted residual evaluation
  last_step_deltat = _system.deltat;

  // The adjoint error estimate expression for a backwards facing step
  // scheme needs the adjoint for the last time instant, so save the current adjoint for future use
  for (auto j : make_range(_system.n_qois()))
  {
    old_adjoints[j] = _system.get_adjoint_solution(j).clone();
  }

  // Retrieve the state and adjoint vectors for the next time instant
  retrieve_timestep();

  // Swap for residual weighting
  for (auto j : make_range(_system.n_qois()))
  {
   _system.get_adjoint_solution(j).swap(*old_adjoints[j]);
  }

  // Swap out the deltats as we did for the left side
  next_step_deltat = _system.deltat;
  _system.deltat = last_step_deltat;

  // Get f(t_j+1)
  ErrorVector QoI_elementwise_error_right;

  _system.update();

  // The residual has to be evaluated at the last time
  _system.time = _system.time - last_step_deltat;

  adjoint_refinement_error_estimator.estimate_error(_system, QoI_elementwise_error_right);

  // Shift the time back
  _system.time = _system.time + last_step_deltat;

  // Set the system deltat back to what it needs to be able to march to the next time
  _system.deltat = next_step_deltat;

  // Swap back now that the residual weighting is done
  for (auto j : make_range(_system.n_qois()))
  {
   _system.get_adjoint_solution(j).swap(*old_adjoints[j]);
  }

  // Also get the right side contributions for the spatially integrated errors for all the QoIs in the QoI set
  for (auto j : make_range(_system.n_qois()))
  {
    // Skip this QoI if not in the QoI Set
    if (adjoint_refinement_error_estimator.qoi_set().has_index(j))
    {
      qoi_error_estimates_right[j] = adjoint_refinement_error_estimator.get_global_QoI_error_estimate(j);
    }
  }

  // Error contribution from this timestep
  for(unsigned int i = 0; i < QoI_elementwise_error.size(); i++)
    QoI_elementwise_error[i] = ((QoI_elementwise_error_right[i] + QoI_elementwise_error_left[i])/2.)*(time_right - time_left);

  // QoI set spatially integrated errors contribution from this timestep
  for (auto j : make_range(_system.n_qois()))
  {
    // Skip this QoI if not in the QoI Set
    if (adjoint_refinement_error_estimator.qoi_set().has_index(j))
    {
      if(QoI_time_instant[j] == NULL)
      {
        (_system.qoi_error_estimates)[j] = ( (1.0 - theta)*qoi_error_estimates_left[j] + theta*qoi_error_estimates_right[j] )*last_step_deltat;
      }
      else if(time_right <= *(QoI_time_instant[j]) + TOLERANCE)
      {
        (_system.qoi_error_estimates)[j] = ( (1.0 - theta)*qoi_error_estimates_left[j] + theta*qoi_error_estimates_right[j] )*last_step_deltat;
      }
      else
      {
        (_system.qoi_error_estimates)[j] = 0.0;
      }
    }
  }

}

} // namespace libMesh
