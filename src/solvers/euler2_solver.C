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
#include "libmesh/euler2_solver.h"

namespace libMesh
{



Euler2Solver::Euler2Solver (sys_type& s)
  : UnsteadySolver(s), theta(1.)
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
                                     DiffContext &context)
{
  return this->_general_residual(request_jacobian,
                                 context,
                                 &DifferentiablePhysics::mass_residual,
                                 &DifferentiablePhysics::_eulerian_time_deriv,
                                 &DifferentiablePhysics::element_constraint,
                                 &DiffContext::elem_reinit);
}



bool Euler2Solver::side_residual (bool request_jacobian,
                                  DiffContext &context)
{
  return this->_general_residual(request_jacobian,
                                 context,
                                 &DifferentiablePhysics::side_mass_residual,
                                 &DifferentiablePhysics::side_time_derivative,
                                 &DifferentiablePhysics::side_constraint,
                                 &DiffContext::elem_side_reinit);
}



bool Euler2Solver::nonlocal_residual (bool request_jacobian,
                                      DiffContext &context)
{
  return this->_general_residual(request_jacobian,
                                 context,
                                 &DifferentiablePhysics::nonlocal_mass_residual,
                                 &DifferentiablePhysics::nonlocal_time_derivative,
                                 &DifferentiablePhysics::nonlocal_constraint,
                                 &DiffContext::nonlocal_reinit);
}



bool Euler2Solver::_general_residual (bool request_jacobian,
                                      DiffContext &context,
                                      ResFuncType mass,
                                      ResFuncType time_deriv,
                                      ResFuncType constraint,
                                      ReinitFuncType reinit_func)
{
  unsigned int n_dofs = context.get_elem_solution().size();

  // Local nonlinear solution at old timestep
  DenseVector<Number> old_elem_solution(n_dofs);
  for (unsigned int i=0; i != n_dofs; ++i)
    old_elem_solution(i) =
      old_nonlinear_solution(context.get_dof_indices()[i]);

  // We evaluate mass_residual with the change in solution
  // to get the mass matrix, reusing old_elem_solution to hold that
  // delta_solution.
  DenseVector<Number> delta_elem_solution(context.get_elem_solution());
  delta_elem_solution -= old_elem_solution;

  // Our first evaluations are at the true elem_solution
  context.elem_solution_derivative = 1.0;

  // If a fixed solution is requested, we'll use the elem_solution
  // at the new timestep
  if (_system.use_fixed_solution)
    context.get_elem_fixed_solution() = context.get_elem_solution();

  context.fixed_solution_derivative = 1.0;

  // We're going to compute just the change in elem_residual
  // (and possibly elem_jacobian), then add back the old values
  DenseVector<Number> total_elem_residual(context.get_elem_residual());
  DenseMatrix<Number> old_elem_jacobian, total_elem_jacobian;
  if (request_jacobian)
    {
      old_elem_jacobian = context.get_elem_jacobian();
      total_elem_jacobian = context.get_elem_jacobian();
      context.get_elem_jacobian().zero();
    }
  context.get_elem_residual().zero();

  // First, evaluate time derivative at the new timestep.
  // The element should already be in the proper place
  // even for a moving mesh problem.
  bool jacobian_computed =
    (_system.*time_deriv)(request_jacobian, context);

  // Scale the time-dependent residual and jacobian correctly
  context.get_elem_residual() *= (theta * _system.deltat);
  total_elem_residual += context.get_elem_residual();
  context.get_elem_residual().zero();

  if (jacobian_computed)
    {
      context.get_elem_jacobian() *= (theta * _system.deltat);
      total_elem_jacobian += context.get_elem_jacobian();
      context.get_elem_jacobian().zero();
    }

  // Next, evaluate the mass residual at the new timestep,
  // with the delta_solution.
  // Evaluating the mass residual at both old and new timesteps will be
  // redundant in most problems but may be necessary for time accuracy
  // or stability in moving mesh problems or problems with user-overridden
  // mass_residual functions

  // Move elem_->delta_, delta_->elem_
  context.get_elem_solution().swap(delta_elem_solution);

  jacobian_computed = (_system.*mass)(jacobian_computed, context) &&
    jacobian_computed;

  context.get_elem_residual() *= -theta;
  total_elem_residual += context.get_elem_residual();
  context.get_elem_residual().zero();

  if (jacobian_computed)
    {
      // The minus sign trick here is to avoid using a non-1
      // elem_solution_derivative:
      context.get_elem_jacobian() *= -theta;
      total_elem_jacobian += context.get_elem_jacobian();
      context.get_elem_jacobian().zero();
    }

  // Add the time-dependent term for the old solution

  // Make sure elem_solution is set up for elem_reinit to use
  // Move delta_->old_, old_->elem_
  context.get_elem_solution().swap(old_elem_solution);

  // Move the mesh into place first if necessary
  (context.*reinit_func)(0.);

  if (_system.use_fixed_solution)
    {
      context.elem_solution_derivative = 0.0;
      jacobian_computed =
        (_system.*time_deriv)(jacobian_computed, context) &&
        jacobian_computed;
      context.elem_solution_derivative = 1.0;
      context.get_elem_residual() *= ((1. - theta) * _system.deltat);
      total_elem_residual += context.get_elem_residual();
      if (jacobian_computed)
        {
          context.get_elem_jacobian() *= ((1. - theta) * _system.deltat);
          total_elem_jacobian += context.get_elem_jacobian();
          context.get_elem_jacobian().zero();
        }
    }
  else
    {
      // FIXME - we should detect if foo_time_derivative() edits
      // elem_jacobian and lies about it!
      (_system.*time_deriv)(false, context);
      context.get_elem_residual() *= ((1. - theta) * _system.deltat);
      total_elem_residual += context.get_elem_residual();
    }

  context.get_elem_residual().zero();

  // Add the mass residual term for the old solution

  // Move old_->old_, delta_->elem_
  context.get_elem_solution().swap(old_elem_solution);

  jacobian_computed =
    (_system.*mass)(jacobian_computed, context) &&
    jacobian_computed;

  context.get_elem_residual() *= -(1. - theta);
  total_elem_residual += context.get_elem_residual();
  context.get_elem_residual().zero();

  if (jacobian_computed)
    {
      // The minus sign trick here is to avoid using a non-1
      // *_solution_derivative:
      context.get_elem_jacobian() *= -(1. - theta);
      total_elem_jacobian += context.get_elem_jacobian();
      context.get_elem_jacobian().zero();
    }

  // Restore the elem_solution
  // Move elem_->elem_, delta_->delta_
  context.get_elem_solution().swap(delta_elem_solution);

  // Restore the elem position if necessary
  (context.*reinit_func)(1.);

  // Add the constraint term
  jacobian_computed = (_system.*constraint)(jacobian_computed, context) &&
    jacobian_computed;

  // Add back the previously accumulated residual and jacobian
  context.get_elem_residual() += total_elem_residual;
  if (request_jacobian)
    {
      if (jacobian_computed)
        context.get_elem_jacobian() += total_elem_jacobian;
      else
        context.get_elem_jacobian().swap(old_elem_jacobian);
    }

  return jacobian_computed;
}



} // namespace libMesh
