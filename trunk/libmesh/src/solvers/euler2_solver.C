
#include "diff_system.h"
#include "euler2_solver.h"
#include "numeric_vector.h"



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




bool Euler2Solver::element_residual (bool request_jacobian)
{
  unsigned int n_dofs = _system.elem_solution.size();

  // Local nonlinear solution at old timestep
  DenseVector<Number> old_elem_solution(n_dofs);
  for (unsigned int i=0; i != n_dofs; ++i)
    old_elem_solution(i) =
      old_nonlinear_solution(_system.dof_indices[i]);

  // If a fixed solution is requested, we'll use the elem_solution
  // at the new timestep
  if (_system.use_fixed_solution)
    {
      _system.elem_fixed_solution = _system.elem_solution;

      _system.fixed_solution_derivative = 1.0;
    }

  // We're going to store old values, since the new changes in
  // elem_residual and elem_jacobian need to be scaled.
  DenseVector<Number> old_elem_residual(_system.elem_residual);
  DenseMatrix<Number> old_elem_jacobian;
  if (request_jacobian)
    {
      old_elem_jacobian = _system.elem_jacobian;
      _system.elem_jacobian.zero();
    }
  _system.elem_residual.zero();

  // First, evaluate time derivative at the new timestep.
  // The element should already be in the proper place
  // even for a moving mesh problem.
  bool jacobian_computed =
    _system.element_time_derivative(request_jacobian);

  // Scale the new time-dependent residual and jacobian correctly
  _system.elem_residual *= (theta * _system.deltat);
  old_elem_residual += _system.elem_residual;
  _system.elem_residual.zero();
  if (jacobian_computed)
    {
      _system.elem_jacobian *= (theta * _system.deltat);
      old_elem_jacobian += _system.elem_jacobian;
      _system.elem_jacobian.zero();
    }

  // We'll add the time-dependent term for the old solution next
  _system.elem_solution.swap(old_elem_solution);

  // Move the mesh into place first if necessary
  _system.elem_reinit(0.);

  // Add the time-dependent term for the old solution
  if (_system.use_fixed_solution)
    {
      _system.elem_solution_derivative = 0.0;
      jacobian_computed = _system.element_time_derivative(jacobian_computed) &&
        jacobian_computed;
      _system.elem_solution_derivative = 1.0;
      _system.elem_residual *= ((1. - theta) * _system.deltat);
      old_elem_residual += _system.elem_residual;
      _system.elem_residual.zero();
      if (jacobian_computed)
        {
          _system.elem_jacobian *= ((1. - theta) * _system.deltat);
          old_elem_jacobian += _system.elem_jacobian;
          _system.elem_jacobian.zero();
        }
    }
  else
    {
      // FIXME - we should detect if element_time_derivative() edits
      // elem_jacobian and lies about it!
      _system.element_time_derivative(false);
      _system.elem_residual *= ((1. - theta) * _system.deltat);
      old_elem_residual += _system.elem_residual;
      _system.elem_residual.zero();
    }

  // Add the mass term for the old solution
  if (_system.use_fixed_solution)
    {
      _system.elem_solution_derivative = 0.0;
      jacobian_computed = _system.mass_residual(jacobian_computed) &&
        jacobian_computed;
      _system.elem_solution_derivative = 1.0;
    }
  else
    {
      // FIXME - we should detect if mass_residual() edits
      // elem_jacobian and lies about it!
      _system.mass_residual(false);
    }

  // Restore the elem_solution
  _system.elem_solution.swap(old_elem_solution);

  // Restore the elem position if necessary
  _system.elem_reinit(1.);

  // Subtract the mass term for the new solution
  if (jacobian_computed)
    _system.elem_jacobian *= -1.0;
  _system.elem_residual *= -1.0;
  jacobian_computed = _system.mass_residual(jacobian_computed) &&
    jacobian_computed;
  if (jacobian_computed)
    _system.elem_jacobian *= -1.0;
  _system.elem_residual *= -1.0;

  // Add the constraint term
  jacobian_computed = _system.element_constraint(jacobian_computed) &&
    jacobian_computed;

  // Add back the old residual and jacobian
  _system.elem_residual += old_elem_residual;
  if (request_jacobian)
    {
      if (jacobian_computed)
        _system.elem_jacobian += old_elem_jacobian;
      else
        _system.elem_jacobian.swap(old_elem_jacobian);
    }

  return jacobian_computed;
}



bool Euler2Solver::side_residual (bool request_jacobian)
{
  unsigned int n_dofs = _system.elem_solution.size();

  // Local nonlinear solution at old timestep
  DenseVector<Number> old_elem_solution(n_dofs);
  for (unsigned int i=0; i != n_dofs; ++i)
    old_elem_solution(i) =
      old_nonlinear_solution(_system.dof_indices[i]);

  // If a fixed solution is requested, we'll use the elem_solution
  // at the new timestep
  if (_system.use_fixed_solution)
    {
      _system.elem_fixed_solution = _system.elem_solution;

      _system.fixed_solution_derivative = 1.0;
    }

  // We're going to store old values, since the new changes in
  // elem_residual and elem_jacobian need to be scaled.
  DenseVector<Number> old_elem_residual(_system.elem_residual);
  DenseMatrix<Number> old_elem_jacobian;
  if (request_jacobian)
    {
      old_elem_jacobian = _system.elem_jacobian;
      _system.elem_jacobian.zero();
    }
  _system.elem_residual.zero();

  // First, evaluate time derivative at the new timestep.
  // The element should already be in the proper place
  // even for a moving mesh problem.
  bool jacobian_computed =
    _system.side_time_derivative(request_jacobian);

  // Scale the time-dependent residual and jacobian correctly
  _system.elem_residual *= (theta * _system.deltat);
  old_elem_residual += _system.elem_residual;
  _system.elem_residual.zero();
  if (jacobian_computed)
    {
      _system.elem_jacobian *= (theta * _system.deltat);
      old_elem_jacobian += _system.elem_jacobian;
      _system.elem_jacobian.zero();
    }

  // Add the time-dependent term for the old solution
  _system.elem_solution.swap(old_elem_solution);

  // Move the mesh into place first if necessary
  _system.elem_side_reinit(0.);

  if (_system.use_fixed_solution)
    {
      _system.elem_solution_derivative = 0.0;
      jacobian_computed = _system.side_time_derivative(jacobian_computed) &&
        jacobian_computed;
      _system.elem_solution_derivative = 1.0;
      _system.elem_residual *= ((1. - theta) * _system.deltat);
      old_elem_residual += _system.elem_residual;
      _system.elem_residual.zero();
      if (jacobian_computed)
        {
          _system.elem_jacobian *= ((1. - theta) * _system.deltat);
          old_elem_jacobian += _system.elem_jacobian;
          _system.elem_jacobian.zero();
        }
    }
  else
    {
      // FIXME - we should detect if side_time_derivative() edits
      // elem_jacobian and lies about it!
      _system.side_time_derivative(false);
      _system.elem_residual *= ((1. - theta) * _system.deltat);
      old_elem_residual += _system.elem_residual;
      _system.elem_residual.zero();
    }

  // Add the mass term for the old solution
  if (_system.use_fixed_solution)
    {
      _system.elem_solution_derivative = 0.0;
      jacobian_computed = _system.side_mass_residual(jacobian_computed) &&
        jacobian_computed;
      _system.elem_solution_derivative = 1.0;
    }
  else
    {
      // FIXME - we should detect if side_mass_residual() edits
      // elem_jacobian and lies about it!
      _system.side_mass_residual(false);
    }

  // Restore the elem_solution
  _system.elem_solution.swap(old_elem_solution);

  // Restore the elem position if necessary
  _system.elem_side_reinit(1.);

  // Subtract the mass term for the new solution
  if (jacobian_computed)
    _system.elem_jacobian *= -1.0;
  _system.elem_residual *= -1.0;
  jacobian_computed = _system.side_mass_residual(jacobian_computed) &&
    jacobian_computed;
  if (jacobian_computed)
    _system.elem_jacobian *= -1.0;
  _system.elem_residual *= -1.0;

  // Add the constraint term
  jacobian_computed = _system.side_constraint(jacobian_computed) &&
    jacobian_computed;

  // Add back the old residual and jacobian
  _system.elem_residual += old_elem_residual;
  if (request_jacobian)
    {
      if (jacobian_computed)
        _system.elem_jacobian += old_elem_jacobian;
      else
        _system.elem_jacobian.swap(old_elem_jacobian);
    }

  return jacobian_computed;
}
