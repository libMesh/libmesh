// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


#include "diff_system.h"
#include "euler_solver.h"

namespace libMesh
{



EulerSolver::EulerSolver (sys_type& s)
 : UnsteadySolver(s), theta(1.)
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
  unsigned int n_dofs = context.elem_solution.size();

  // Local nonlinear solution at old timestep
  DenseVector<Number> old_elem_solution(n_dofs);
  for (unsigned int i=0; i != n_dofs; ++i)
    old_elem_solution(i) =
      old_nonlinear_solution(context.dof_indices[i]);

  // Local nonlinear solution at time t_theta
  DenseVector<Number> theta_solution(context.elem_solution);
  theta_solution *= theta;
  theta_solution.add(1. - theta, old_elem_solution);

  // Technically the elem_solution_derivative is either theta
  // or -1.0 in this implementation, but we scale the former part
  // ourselves
  context.elem_solution_derivative = 1.0;

// Technically the fixed_solution_derivative is always theta,
// but we're scaling a whole jacobian by theta after these first
// evaluations
  context.fixed_solution_derivative = 1.0;

  // If a fixed solution is requested, we'll use theta_solution
  if (_system.use_fixed_solution)
    context.elem_fixed_solution = theta_solution;

  // Move theta_->elem_, elem_->theta_
  context.elem_solution.swap(theta_solution);

  // Move the mesh into place first if necessary
  context.elem_reinit(theta);

  // We're going to compute just the change in elem_residual
  // (and possibly elem_jacobian), then add back the old values
  DenseVector<Number> old_elem_residual(context.elem_residual);
  DenseMatrix<Number> old_elem_jacobian;
  if (request_jacobian)
    {
      old_elem_jacobian = context.elem_jacobian;
      context.elem_jacobian.zero();
    }
  context.elem_residual.zero();

  // Get the time derivative at t_theta
  bool jacobian_computed =
    _system.element_time_derivative(request_jacobian, context);

  // For a moving mesh problem we may need the pseudoconvection term too
  jacobian_computed =
    _system.eulerian_residual(jacobian_computed, context) && jacobian_computed;

  // Scale the time-dependent residual and jacobian correctly
  context.elem_residual *= _system.deltat;
  if (jacobian_computed)
    context.elem_jacobian *= (theta * _system.deltat);

  // The fixed_solution_derivative is always theta,
  // and now we're done scaling jacobians
  context.fixed_solution_derivative = theta;

  // We evaluate mass_residual with the change in solution
  // to get the mass matrix, reusing old_elem_solution to hold that
  // delta_solution.  We're solving dt*F(u) - du = 0, so our
  // delta_solution is old_solution - new_solution.
  // We're still keeping elem_solution in theta_solution for now
  old_elem_solution -= theta_solution;

  // Move old_->elem_, theta_->old_
  context.elem_solution.swap(old_elem_solution);

  // We do a trick here to avoid using a non-1
  // elem_solution_derivative:
  context.elem_jacobian *= -1.0;
  jacobian_computed = _system.mass_residual(jacobian_computed, context) &&
    jacobian_computed;
  context.elem_jacobian *= -1.0;

  // Move elem_->elem_, old_->theta_
  context.elem_solution.swap(theta_solution);

  // Restore the elem position if necessary
  context.elem_reinit(1.);

  // Add the constraint term
  jacobian_computed = _system.element_constraint(jacobian_computed, context) &&
    jacobian_computed;

  // Add back the old residual and jacobian
  context.elem_residual += old_elem_residual;
  if (request_jacobian)
    {
      if (jacobian_computed)
        context.elem_jacobian += old_elem_jacobian;
      else
        context.elem_jacobian.swap(old_elem_jacobian);
    }

  return jacobian_computed;
}



bool EulerSolver::side_residual (bool request_jacobian,
                                 DiffContext &context)
{
  unsigned int n_dofs = context.elem_solution.size();

  // Local nonlinear solution at old timestep
  DenseVector<Number> old_elem_solution(n_dofs);
  for (unsigned int i=0; i != n_dofs; ++i)
    old_elem_solution(i) =
      old_nonlinear_solution(context.dof_indices[i]);

  // Local nonlinear solution at time t_theta
  DenseVector<Number> theta_solution(context.elem_solution);
  theta_solution *= theta;
  theta_solution.add(1. - theta, old_elem_solution);

  // Technically the elem_solution_derivative is either theta
  // or 1.0 in this implementation, but we scale the former part
  // ourselves
  context.elem_solution_derivative = 1.0;

// Technically the fixed_solution_derivative is always theta,
// but we're scaling a whole jacobian by theta after these first
// evaluations
  context.fixed_solution_derivative = 1.0;

  // If a fixed solution is requested, we'll use theta_solution
  if (_system.use_fixed_solution)
    context.elem_fixed_solution = theta_solution;

  // Move theta_->elem_, elem_->theta_
  context.elem_solution.swap(theta_solution);

  // Move the mesh into place first if necessary
  context.elem_side_reinit(theta);

  // We're going to compute just the change in elem_residual
  // (and possibly elem_jacobian), then add back the old values
  DenseVector<Number> old_elem_residual(context.elem_residual);
  DenseMatrix<Number> old_elem_jacobian;
  if (request_jacobian)
    {
      old_elem_jacobian = context.elem_jacobian;
      context.elem_jacobian.zero();
    }
  context.elem_residual.zero();

  // Get the time derivative at t_theta
  bool jacobian_computed =
    _system.side_time_derivative(request_jacobian, context);

  // Scale the time-dependent residual and jacobian correctly
  context.elem_residual *= _system.deltat;
  if (jacobian_computed)
    context.elem_jacobian *= (theta * _system.deltat);

  // The fixed_solution_derivative is always theta,
  // and now we're done scaling jacobians
  context.fixed_solution_derivative = theta;

  // We evaluate side_mass_residual with the change in solution
  // to get the mass matrix, reusing old_elem_solution to hold that
  // delta_solution.  We're solving dt*F(u) - du = 0, so our
  // delta_solution is old_solution - new_solution.
  // We're still keeping elem_solution in theta_solution for now
  old_elem_solution -= theta_solution;

  // Move old_->elem_, theta_->old_
  context.elem_solution.swap(old_elem_solution);

  // We do a trick here to avoid using a non-1
  // elem_solution_derivative:
  context.elem_jacobian *= -1.0;
  jacobian_computed = _system.side_mass_residual(jacobian_computed, context) &&
    jacobian_computed;
  context.elem_jacobian *= -1.0;

  // Move elem_->elem_, old_->theta_
  context.elem_solution.swap(theta_solution);

  // Restore the elem position if necessary
  context.elem_side_reinit(1.);

  // Add the constraint term
  jacobian_computed = _system.side_constraint(jacobian_computed, context) &&
    jacobian_computed;

  // Add back the old residual and jacobian
  context.elem_residual += old_elem_residual;
  if (request_jacobian)
    {
      if (jacobian_computed)
        context.elem_jacobian += old_elem_jacobian;
      else
        context.elem_jacobian.swap(old_elem_jacobian);
    }

  return jacobian_computed;
}

} // namespace libMesh
