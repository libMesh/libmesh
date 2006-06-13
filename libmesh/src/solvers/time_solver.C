
#include "diff_solver.h"
#include "diff_system.h"
#include "numeric_vector.h"
#include "time_solver.h"



TimeSolver::TimeSolver (sys_type& s)
  : diff_solver(NULL),
    _system(s),
    first_solve(true)
{
}



TimeSolver::~TimeSolver ()
{
}



void TimeSolver::reinit ()
{
  diff_solver->reinit();
}



void TimeSolver::init ()
{
  // If the user hasn't given us a solver to use,
  // just build a default solver
  if (diff_solver.get() == NULL)
    diff_solver = DiffSolver::build(_system);
  diff_solver->init();

  _system.add_vector("_old_nonlinear_solution");
}



void TimeSolver::solve ()
{
  if (first_solve)
    {
      advance_timestep();
      first_solve = false;
    }

  diff_solver->solve();
}



void TimeSolver::advance_timestep ()
{
  NumericVector<Number> &old_nonlinear_solution =
  _system.get_vector("_old_nonlinear_solution");
  NumericVector<Number> &nonlinear_solution =
    _system.get_vector("_nonlinear_solution");

  old_nonlinear_solution = nonlinear_solution;

  if (!first_solve)
    _system.time += _system.deltat;
}
