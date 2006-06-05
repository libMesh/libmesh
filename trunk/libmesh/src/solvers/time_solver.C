
#include "diff_solver.h"
#include "time_solver.h"



TimeSolver::TimeSolver (sys_type& s) :
  _system (s),
  diff_solver(DiffSolver::build(s)) {}



TimeSolver::~TimeSolver ()
{
}
