#include "diff_solver.h"
#include "newton_solver.h"



AutoPtr<DiffSolver> DiffSolver::build (sys_type& s)
{
  return AutoPtr<DiffSolver>(new NewtonSolver(s));
}
