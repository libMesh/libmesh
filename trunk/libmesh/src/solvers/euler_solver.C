
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

  DenseVector<Number> old_elem_solution(n_dofs);

  // Local nonlinear solution at old timestep
  for (unsigned int i=0; i != n_dofs; ++i)
    old_elem_solution(i) =
      old_nonlinear_solution(_system.dof_indices[i]);

  // Local nonlinear solution at time t_theta
  DenseVector<Number> theta_solution(_system.elem_solution);
  theta_solution *= theta;
  old_elem_solution *= (1. - theta);
  theta_solution += old_elem_solution;

  // Temporarily replace elem_solution with theta_solution
  old_elem_solution = _system.elem_solution;
  _system.elem_solution = theta_solution;

  bool jacobian_computed =
    _system.element_time_derivative(request_jacobian);

  // Scale the time-dependent residual and jacobian correctly
  _system.elem_residual *= _system.deltat;
  if (jacobian_computed)
    _system.elem_jacobian *= theta * _system.deltat;

  // Add the mass term and constraint term
  jacobian_computed = _system.mass_residual(request_jacobian) &&
    jacobian_computed;

  jacobian_computed = _system.element_constraint(request_jacobian) &&
    jacobian_computed;

  // Restore the elem_solution
  _system.elem_solution = old_elem_solution;

  return jacobian_computed;
}



bool EulerSolver::side_residual (bool request_jacobian)
{
  // Global nonlinear solution at old timestep
  NumericVector<Number> &old_nonlinear_solution =
    _system.get_vector("_old_nonlinear_solution");

  unsigned int n_dofs = _system.elem_solution.size();

  DenseVector<Number> old_elem_solution(n_dofs);

  // Local nonlinear solution at old timestep
  for (unsigned int i=0; i != n_dofs; ++i)
    old_elem_solution(i) =
      old_nonlinear_solution(_system.dof_indices[i]);

  // Local nonlinear solution at time t_theta
  DenseVector<Number> theta_solution(_system.elem_solution);
  theta_solution *= theta;
  old_elem_solution *= (1. - theta);
  theta_solution += old_elem_solution;

  // Temporarily replace elem_solution with theta_solution
  old_elem_solution = _system.elem_solution;
  _system.elem_solution = theta_solution;

  bool jacobian_computed =
    _system.side_time_derivative(request_jacobian);

  // Scale the time-dependent residual and jacobian correctly
  _system.elem_residual *= _system.deltat;
  if (jacobian_computed)
    _system.elem_jacobian *= theta * _system.deltat;

  // Add the constraint term (we shouldn't need a mass term on sides)
  jacobian_computed = _system.side_constraint(request_jacobian) &&
    jacobian_computed;

  // Restore the elem_solution
  _system.elem_solution = old_elem_solution;

  return jacobian_computed;
}
