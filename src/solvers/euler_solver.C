// The libMesh Finite Element Library.
// Copyright (C) 2002-2014 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/euler_solver.h"

namespace libMesh
{



EulerSolver::EulerSolver (sys_type& s)
  : FirstOrderUnsteadySolver(s), theta(1.)
{
}



EulerSolver::~EulerSolver ()
{
}



Real EulerSolver::error_order() const
{
  if (theta == 0.5)
    return 2.;
  return 1.;
}




bool EulerSolver::element_residual (bool request_jacobian,
                                    DiffContext &context)
{
  return this->_general_residual(request_jacobian,
                                 context,
                                 &DifferentiablePhysics::mass_residual,
                                 &DifferentiablePhysics::_eulerian_time_deriv,
                                 &DifferentiablePhysics::element_constraint,
                                 &DiffContext::elem_reinit);
}



bool EulerSolver::side_residual (bool request_jacobian,
                                 DiffContext &context)
{
  return this->_general_residual(request_jacobian,
                                 context,
                                 &DifferentiablePhysics::side_mass_residual,
                                 &DifferentiablePhysics::side_time_derivative,
                                 &DifferentiablePhysics::side_constraint,
                                 &DiffContext::elem_side_reinit);
}



bool EulerSolver::nonlocal_residual (bool request_jacobian,
                                     DiffContext &context)
{
  return this->_general_residual(request_jacobian,
                                 context,
                                 &DifferentiablePhysics::nonlocal_mass_residual,
                                 &DifferentiablePhysics::nonlocal_time_derivative,
                                 &DifferentiablePhysics::nonlocal_constraint,
                                 &DiffContext::nonlocal_reinit);
}



bool EulerSolver::_general_residual (bool request_jacobian,
                                     DiffContext &context,
                                     ResFuncType mass,
                                     ResFuncType time_deriv,
                                     ResFuncType constraint,
                                     ReinitFuncType reinit_func)
{
  unsigned int n_dofs = context.get_elem_solution().size();

  // We might need to save the old jacobian in case one of our physics
  // terms later is unable to update it analytically.
  DenseMatrix<Number> old_elem_jacobian(n_dofs, n_dofs);
  if (request_jacobian)
    old_elem_jacobian.swap(context.get_elem_jacobian());

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

  // Local nonlinear solution at time t_theta
  DenseVector<Number> theta_solution(context.get_elem_solution());
  theta_solution *= theta;
  theta_solution.add(1. - theta, old_elem_solution);

  context.elem_solution_derivative = theta;
  context.fixed_solution_derivative = theta;

  // If a fixed solution is requested, we'll use theta_solution
  if (_system.use_fixed_solution)
    context.get_elem_fixed_solution() = theta_solution;

  // Move theta_->elem_, elem_->theta_
  context.get_elem_solution().swap(theta_solution);

  // Move the mesh into place first if necessary
  (context.*reinit_func)(theta);

  // Get the time derivative at t_theta
  bool jacobian_computed =
    (_system.*time_deriv)(request_jacobian, context);

  jacobian_computed = (_system.*mass)(jacobian_computed, context) &&
    jacobian_computed;

  // Restore the elem position if necessary
  (context.*reinit_func)(1);

  // Move elem_->elem_, theta_->theta_
  context.get_elem_solution().swap(theta_solution);
  context.elem_solution_derivative = 1;

  // Add the constraint term
  jacobian_computed = (_system.*constraint)(jacobian_computed, context) &&
    jacobian_computed;

  // Add back (or restore) the old jacobian
  if (request_jacobian)
    {
      if (jacobian_computed)
        context.get_elem_jacobian() += old_elem_jacobian;
      else
        context.get_elem_jacobian().swap(old_elem_jacobian);
    }

  return jacobian_computed;
}


} // namespace libMesh
