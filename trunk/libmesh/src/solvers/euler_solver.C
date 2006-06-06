
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



bool EulerSolver::element_residual (bool request_jacobian)
{
  // Global nonlinear solution at old timestep
  NumericVector<Number> &old_nonlinear_solution =
    _system.get_vector("_old_nonlinear_solution");

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

  // We're going to compute just the change in elem_residual,
  // then add back the old elem_residual.
  DenseVector<Number> old_elem_residual(_system.elem_residual);
  DenseMatrix<Number> old_elem_jacobian;
  if (request_jacobian)
    old_elem_jacobian = _system.elem_jacobian;
  _system.elem_residual.zero();

  bool jacobian_computed =
    _system.element_time_derivative(request_jacobian);

  // Scale the time-dependent residual and jacobian correctly
  _system.elem_residual *= _system.deltat;
  if (jacobian_computed)
    _system.elem_jacobian *= (theta * _system.deltat);

  // Add the mass term for the old solution
  _system.elem_solution.swap(old_elem_solution);
  // FIXME - this will break if mass_residual() edits
  // elem_jacobian and lies about it!
  _system.mass_residual(false);

  // Restore the elem_solution
  _system.elem_solution.swap(theta_solution);

  // Add the mass term for the new solution
  jacobian_computed = _system.mass_residual(request_jacobian) &&
    jacobian_computed;

  // Add the constraint term
  jacobian_computed = _system.element_constraint(request_jacobian) &&
    jacobian_computed;

  // Add back the old residual and jacobian
  _system.elem_residual += old_elem_residual;
  if (request_jacobian)
    _system.elem_jacobian += old_elem_jacobian;

  return jacobian_computed;
}



bool EulerSolver::side_residual (bool request_jacobian)
{
  // Global nonlinear solution at old timestep
  NumericVector<Number> &old_nonlinear_solution =
    _system.get_vector("_old_nonlinear_solution");

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

  // We're going to compute just the change in elem_residual,
  // then add back the old elem_residual.
  DenseVector<Number> old_elem_residual(_system.elem_residual);
  DenseMatrix<Number> old_elem_jacobian;
  if (request_jacobian)
    old_elem_jacobian = _system.elem_jacobian;
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
  jacobian_computed = _system.side_constraint(request_jacobian) &&
    jacobian_computed;

  // Add back the old residual and jacobian
  _system.elem_residual += old_elem_residual;
  if (request_jacobian)
    _system.elem_jacobian += old_elem_jacobian;

  return jacobian_computed;
}
