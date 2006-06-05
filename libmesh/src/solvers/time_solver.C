
#include "diff_solver.h"
#include "time_solver.h"



TimeSolver::TimeSolver (sys_type& s) :
  time(0.),
  deltat(1.),
  _system (s),
  diff_solver(DiffSolver::build(s)) {}



TimeSolver::~TimeSolver ()
{
}



void TimeSolver::init ()
{
  diff_solver->init();
}



void TimeSolver::solve ()
{
  diff_solver->solve();
}
