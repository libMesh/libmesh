// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
  bool compute_second_order_eqns = !this->_system.get_second_order_vars().empty();

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
  if(compute_second_order_eqns)
    this->prepare_accel(context);

  // Move the mesh into place first if necessary, set t = t_{n+1}
  (context.*reinit_func)(1.);

  // First, evaluate time derivative at the new timestep.
  // The element should already be in the proper place
  // even for a moving mesh problem.
  bool jacobian_computed =
    (_system.*time_deriv)(request_jacobian, context);

  // Next, evaluate the mass residual at the new timestep

  jacobian_computed = (_system.*mass)(jacobian_computed, context) &&
    jacobian_computed;

  // If we have second-order variables, we need to get damping terms
  // and the velocity equations
  if(compute_second_order_eqns)
    {
      jacobian_computed = (_system.*damping)(jacobian_computed, context) &&
        jacobian_computed;

      jacobian_computed = this->compute_second_order_eqns(jacobian_computed, context) &&
        jacobian_computed;
    }

  // Add the constraint term
  jacobian_computed = (_system.*constraint)(jacobian_computed, context) &&
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
    (_system.*time_deriv)(jacobian_computed, context) &&
    jacobian_computed;

  // Add the mass residual term for the old solution

  // Evaluating the mass residual at both old and new timesteps will be
  // redundant in most problems but may be necessary for time accuracy
  // or stability in moving mesh problems or problems with user-overridden
  // mass_residual functions

  jacobian_computed =
    (_system.*mass)(jacobian_computed, context) &&
    jacobian_computed;

  // If we have second-order variables, we need to get damping terms
  // and the velocity equations
  if(compute_second_order_eqns)
    {
      jacobian_computed = (_system.*damping)(jacobian_computed, context) &&
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



} // namespace libMesh
