
#include "diff_solver.h"
#include "diff_system.h"
#include "numeric_vector.h"
#include "time_solver.h"



TimeSolver::TimeSolver (sys_type& s)
  : diff_solver(DiffSolver::build(s)),
    _system(s) {}



TimeSolver::~TimeSolver ()
{
}



void TimeSolver::init ()
{
  diff_solver->init();

  _system.add_vector("_old_nonlinear_solution");
}



void TimeSolver::solve ()
{
  NumericVector<Number> &old_nonlinear_solution =
  _system.get_vector("_old_nonlinear_solution");
  NumericVector<Number> &nonlinear_solution =
    _system.get_vector("_nonlinear_solution");

  old_nonlinear_solution = nonlinear_solution;

  diff_solver->solve();

  _system.time += _system.deltat;
}
