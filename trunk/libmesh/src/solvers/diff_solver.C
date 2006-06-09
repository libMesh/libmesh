#include "diff_solver.h"
#include "newton_solver.h"



AutoPtr<DiffSolver> DiffSolver::build (sys_type& s)
{
  return AutoPtr<DiffSolver>(new NewtonSolver(s));
}



void DiffSolver::init ()
{
  // Reset the max_step_size and max_residual_norm for a new problem
  max_solution_norm = 0.;
  max_residual_norm = 0.;
}
