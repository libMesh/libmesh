
#include "diff_system.h"
#include "euler_solver.h"
#include "numeric_vector.h"



EulerSolver::EulerSolver (sys_type& s)
 : TimeSolver(s), theta(1.)
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




bool EulerSolver::element_residual (bool request_jacobian)
{
  unsigned int n_dofs = _system.elem_solution.size();

  // Local nonlinear solution at old timestep
  DenseVector<Number> old_elem_solution(n_dofs);
  for (unsigned int i=0; i != n_dofs; ++i)
    old_elem_solution(i) =
      old_nonlinear_solution(_system.dof_indices[i]);

  // Local nonlinear solution at time t_theta
  DenseVector<Number> theta_solution(_system.elem_solution);
  theta_solution *= theta;
  theta_solution.add(1. - theta, old_elem_solution);

  // If a fixed solution is requested, we'll use theta_solution
  if (_system.use_fixed_solution)
    {
      _system.elem_fixed_solution = theta_solution;
      _system.fixed_solution_derivative = theta;
    }

  // Temporarily replace elem_solution with theta_solution
  _system.elem_solution.swap(theta_solution);

  // We're going to compute just the change in elem_residual
  // (and possibly elem_jacobian), then add back the old values
  DenseVector<Number> old_elem_residual(_system.elem_residual);
  DenseMatrix<Number> old_elem_jacobian;
  if (request_jacobian)
    {
      old_elem_jacobian = _system.elem_jacobian;
      _system.elem_jacobian.zero();
    }
  _system.elem_residual.zero();

  bool jacobian_computed =
    _system.element_time_derivative(request_jacobian);

  // Scale the time-dependent residual and jacobian correctly
  _system.elem_residual *= _system.deltat;
  if (jacobian_computed)
    _system.elem_jacobian *= (theta * _system.deltat);

  // Add the mass term for the old solution
  _system.elem_solution.swap(old_elem_solution);

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
  _system.elem_solution.swap(theta_solution);

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



bool EulerSolver::side_residual (bool request_jacobian)
{
  unsigned int n_dofs = _system.elem_solution.size();

  // Local nonlinear solution at old timestep
  DenseVector<Number> old_elem_solution(n_dofs);
  for (unsigned int i=0; i != n_dofs; ++i)
    old_elem_solution(i) =
      old_nonlinear_solution(_system.dof_indices[i]);

  // Local nonlinear solution at time t_theta
  DenseVector<Number> theta_solution(_system.elem_solution);
  theta_solution *= theta;
  theta_solution.add(1. - theta, old_elem_solution);

  // Temporarily replace elem_solution with theta_solution
  _system.elem_solution.swap(theta_solution);

  // If a fixed solution is requested, we'll use theta_solution
  if (_system.use_fixed_solution)
    {
      _system.elem_fixed_solution = theta_solution;
      _system.fixed_solution_derivative = theta;
    }

  // We're going to compute just the change in elem_residual,
  // then add back the old elem_residual.
  DenseVector<Number> old_elem_residual(_system.elem_residual);
  DenseMatrix<Number> old_elem_jacobian;
  if (request_jacobian)
    {
      old_elem_jacobian = _system.elem_jacobian;
      _system.elem_jacobian.zero();
    }
  _system.elem_residual.zero();

  bool jacobian_computed =
    _system.side_time_derivative(request_jacobian);

  // Scale the time-dependent residual and jacobian correctly
  _system.elem_residual *= _system.deltat;
  if (jacobian_computed)
    _system.elem_jacobian *= (theta * _system.deltat);

  // Restore the elem_solution
  _system.elem_solution.swap(theta_solution);

  // Add the constraint term (we shouldn't need a mass term on sides)
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
