
#include "diff_solver.h"
#include "diff_system.h"
#include "dof_map.h"
#include "numeric_vector.h"
#include "time_solver.h"



TimeSolver::TimeSolver (sys_type& s)
  : quiet(true),
    _diff_solver                 (NULL),
    _system                      (s)
{
}



TimeSolver::~TimeSolver ()
{
}



void TimeSolver::reinit ()
{
  _diff_solver->reinit();
}



void TimeSolver::init ()
{
  // If the user hasn't given us a solver to use,
  // just build a default solver
  if (_diff_solver.get() == NULL)
    _diff_solver = DiffSolver::build(_system);
  _diff_solver->init();
}



void TimeSolver::solve ()
{
  _diff_solver->solve();
}



void TimeSolver::advance_timestep ()
{
}



AutoPtr<DiffSolver> & TimeSolver::diff_solver()
{
  return _diff_solver;
}
