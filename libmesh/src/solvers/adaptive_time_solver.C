
#include "adaptive_time_solver.h"
#include "diff_system.h"
#include "euler_solver.h"
#include "numeric_vector.h"



AdaptiveTimeSolver::AdaptiveTimeSolver (sys_type& s)
 : TimeSolver(s),
   core_time_solver(new EulerSolver(s)),
   target_tolerance(1.e-3), upper_tolerance(1.e-2),
   max_deltat(0.),
   global_tolerance(true)
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
  _system.deltat *= 0.5;
  core_time_solver->solve();
  core_time_solver->advance_timestep();
  core_time_solver->solve();

  // But then back off just in case our advance_timestep() isn't
  // called - this probably doesn't work with multistep methods
  _system.get_vector("_old_nonlinear_solution") = *old_solution;
  _system.time -= _system.deltat;
  _system.deltat = old_deltat;

  // Find the relative error
  const Real double_norm = calculate_norm(_system, *double_solution);
  const Real single_norm = calculate_norm(_system, *_system.solution);
  *double_solution -= *(_system.solution);
  const Real error_norm  = calculate_norm(_system, *double_solution);
  Real relative_error = error_norm / _system.deltat /
    std::max(double_norm, single_norm);

  // If the relative error makes no sense, we're done
  if (!double_norm && !single_norm)
    return;

  if (!quiet)
    {
      std::cout << "\n === Computing adaptive timestep === " << std::endl;
      std::cout << "Double norm = " << double_norm << std::endl;
      std::cout << "Single norm = " << single_norm << std::endl;
      std::cout << "Error norm = " << error_norm << std::endl;
      std::cout << "Local relative error = "
		<< (error_norm / std::max(double_norm, single_norm)) << std::endl;
      std::cout << "Global relative error = "
		<< (error_norm / _system.deltat / std::max(double_norm, single_norm)) << std::endl;
      std::cout << "old delta t = " << _system.deltat << std::endl;
    }

  
  // Otherwise, compare the relative error to the tolerance
  // and adjust deltat
  last_deltat = _system.deltat;

  const Real global_shrink_or_growth_factor =
    std::pow(target_tolerance / relative_error,
	     1. / core_time_solver->error_order());

  const Real local_shrink_or_growth_factor =
    std::pow(target_tolerance / (error_norm/std::max(double_norm, single_norm)),
	     1. / (core_time_solver->error_order()+1.));

  if (!quiet)
    {
      std::cout << "The global growth/shrink factor is: "
		<< global_shrink_or_growth_factor << std::endl;
      std::cout << "The local growth/shrink factor is: "
		<< local_shrink_or_growth_factor << std::endl;
    }

  // The local s.o.g. factor is based on the expected **local**
  // truncation error for the timestepping method, the global
  // s.o.g. factor is based on the method's **global** truncation
  // error.  You can shrink/grow the timestep to attempt to satisfy
  // either a global or local time-discretization error tolerance.
  //bool use_global_shrink_or_growth_factor=false;
  
  if (relative_error > target_tolerance)
    {
      if (global_tolerance)
	_system.deltat *= global_shrink_or_growth_factor;
      else
	_system.deltat *= local_shrink_or_growth_factor;
    }
  
  // We might eventually want to use a little hysteresis when growing deltat
  // One way to achieve hysteresis is to use a different power when computing the
  // growth factor.
  else if (relative_error < target_tolerance)
    {
      if (global_tolerance)
	_system.deltat *= global_shrink_or_growth_factor;
      else
	_system.deltat *= local_shrink_or_growth_factor;
    }

  if ((this->max_deltat != 0.0) && (_system.deltat > this->max_deltat))
    {
      if (!quiet)
	{
	  std::cout << "delta t is constrained by maximum-allowable delta t." << std::endl;
	}
      _system.deltat = this->max_deltat;
    }
  
  if (!quiet)
    {
      std::cout << "new delta t = " << _system.deltat << std::endl;
    }

}



void AdaptiveTimeSolver::advance_timestep ()
{
  NumericVector<Number> &old_nonlinear_solution =
  _system.get_vector("_old_nonlinear_solution");
  NumericVector<Number> &nonlinear_solution =
    *(_system.solution);
//    _system.get_vector("_nonlinear_solution");

  old_nonlinear_solution = nonlinear_solution;

  if (!first_solve)
    _system.time += last_deltat;
}



Real AdaptiveTimeSolver::error_order () const
{
  assert(core_time_solver.get());

  return core_time_solver->error_order();
}



bool AdaptiveTimeSolver::element_residual (bool request_jacobian)
{
  assert(core_time_solver.get());

  return core_time_solver->element_residual(request_jacobian);
}



bool AdaptiveTimeSolver::side_residual (bool request_jacobian)
{
  assert(core_time_solver.get());

  return core_time_solver->side_residual(request_jacobian);
}



AutoPtr<DiffSolver> & AdaptiveTimeSolver::diff_solver()
{
  return core_time_solver->diff_solver();
}



Real AdaptiveTimeSolver::calculate_norm(System &s,
                                        NumericVector<Number> &v)
{
  return s.calculate_norm(v, component_norm, component_scale);
}
