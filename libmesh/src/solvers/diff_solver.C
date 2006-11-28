#include "diff_solver.h"
#include "newton_solver.h"



DiffSolver::DiffSolver (sys_type& s)
    : quiet(true),
      max_linear_iterations(1000),
      max_nonlinear_iterations(100),
      continue_after_max_iterations(true),
      absolute_residual_tolerance(0.),
      relative_residual_tolerance(0.),
      absolute_step_tolerance(0.),
      relative_step_tolerance(0.),
      initial_linear_tolerance(1e-12),
      max_solution_norm(0.),
      max_residual_norm(0.),
      _system (s)
{
}



AutoPtr<DiffSolver> DiffSolver::build (sys_type& s)
{
  return AutoPtr<DiffSolver>(new NewtonSolver(s));
}



void DiffSolver::reinit ()
{
  // Reset the max_step_size and max_residual_norm for a new mesh
  max_solution_norm = 0.;
  max_residual_norm = 0.;
}



void DiffSolver::init ()
{
  // Reset the max_step_size and max_residual_norm for a new problem
  max_solution_norm = 0.;
  max_residual_norm = 0.;
}
