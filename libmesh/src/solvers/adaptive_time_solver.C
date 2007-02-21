
#include "adaptive_time_solver.h"
#include "diff_system.h"
#include "euler_solver.h"
#include "numeric_vector.h"



AdaptiveTimeSolver::AdaptiveTimeSolver (sys_type& s)
 : TimeSolver(s),
   core_time_solver(new EulerSolver(s)),
   tolerance(1.e-3)
{
  // We start with a reasonable time solver: implicit Euler
}



AdaptiveTimeSolver::~AdaptiveTimeSolver ()
{
}



void AdaptiveTimeSolver::init()
{
  assert(core_time_solver.get());

  // We override this because our core_time_solver is the one that
  // needs to handle new vectors, diff_solver->init(), etc
  core_time_solver->init();
}



void AdaptiveTimeSolver::reinit()
{
  assert(core_time_solver.get());

  // We override this because our core_time_solver is the one that
  // needs to handle new vectors, diff_solver->reinit(), etc
  core_time_solver->reinit();
}



void AdaptiveTimeSolver::solve()
{
  assert(core_time_solver.get());

  // The core_time_solver will handle any first_solve actions
  first_solve = false;

  // Use the double-length timestep first (so the old_nonlinear_solution
  // won't have to change)
  core_time_solver->solve();

  // Save a copy of the double-length nonlinear solution
  // and the old nonlinear solution
  AutoPtr<NumericVector<Number> > double_solution =
    _system.solution->clone();
  AutoPtr<NumericVector<Number> > old_solution =
    _system.get_vector("_old_nonlinear_solution").clone();

  // Then reset the solution for our single-length calcs
  *(_system.solution) = _system.get_vector("_old_nonlinear_solution");

  // Call two single-length timesteps
  Real old_deltat = _system.deltat;
  _system.deltat /= 2;
  core_time_solver->solve();
  core_time_solver->advance_timestep();
  core_time_solver->solve();

  // But then back off just in case our advance_timestep() isn't
  // called - this probably doesn't work with multistep methods
  _system.get_vector("_old_nonlinear_solution") = *old_solution;
  _system.time -= _system.deltat;
  _system.deltat = old_deltat;

  // Find the relative error
  Real double_norm = double_solution->l2_norm();
  Real single_norm = _system.solution->l2_norm();
  *double_solution -= *(_system.solution);
  Real error_norm = double_solution->l2_norm();
  Real relative_error = error_norm / std::max(double_norm, single_norm);

  // If the relative error makes no sense, we're done
  if (!double_norm && !single_norm)
    return;
  
  // Otherwise, compare the relative error to the tolerance
  // and adjust deltat

  if (relative_error > tolerance)
    _system.deltat *= (tolerance / relative_error);
  // Use a little hysteresis when enlarging deltat
  else if (relative_error < tolerance)
    _system.deltat *= std::sqrt(tolerance / relative_error);
}



bool AdaptiveTimeSolver::element_residual (bool request_jacobian)
{
  assert(core_time_solver.get());

  return core_time_solver->element_residual(request_jacobian);
}



bool AdaptiveTimeSolver::side_residual (bool request_jacobian)
{
  assert(core_time_solver.get());

  return core_time_solver->element_residual(request_jacobian);
}



AutoPtr<DiffSolver> & AdaptiveTimeSolver::diff_solver()
{
  return core_time_solver->diff_solver();
}
