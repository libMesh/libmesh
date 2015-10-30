
#include "adaptive_time_solver.h"
#include "diff_system.h"
#include "euler_solver.h"
#include "numeric_vector.h"



AdaptiveTimeSolver::AdaptiveTimeSolver (sys_type& s)
 : TimeSolver(s),
   core_time_solver(new EulerSolver(s)),
   target_tolerance(1.e-3), upper_tolerance(0.0),
   max_deltat(0.),
   min_deltat(0.),
   max_growth(0.),
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

  // We may have to repeat timesteps entirely if our error is bad
  // enough
  bool max_tolerance_met = false;

  // Calcualting error values each time
  Real single_norm(0.), double_norm(0.), error_norm(0.),
       relative_error(0.);
 
  while (!max_tolerance_met)
    {
      // If we've been asked to reduce deltat if necessary, make sure
      // the core timesolver does so
      core_time_solver->reduce_deltat_on_diffsolver_failure =
        this->reduce_deltat_on_diffsolver_failure;

      if (!quiet)
        {
          std::cout << "\n === Computing adaptive timestep === " 
                    << std::endl;
        }

      // Use the double-length timestep first (so the
      // old_nonlinear_solution won't have to change)
      core_time_solver->solve();

      // Save a copy of the double-length nonlinear solution
      // and the old nonlinear solution
      AutoPtr<NumericVector<Number> > double_solution =
        _system.solution->clone();
      AutoPtr<NumericVector<Number> > old_solution =
        _system.get_vector("_old_nonlinear_solution").clone();

      double_norm = calculate_norm(_system, *double_solution);
      if (!quiet)
        {
	  std::cout << "Double norm = " << double_norm << std::endl;
        }

      // Then reset the initial guess for our single-length calcs
      *(_system.solution) = _system.get_vector("_old_nonlinear_solution");

      // Call two single-length timesteps
      // Be sure that the core_time_solver does not change the
      // timestep here.  (This is unlikely because it just succeeded
      // with a timestep twice as large!)
      // FIXME: even if diffsolver failure is unlikely, we ought to
      // do *something* if it happens
      core_time_solver->reduce_deltat_on_diffsolver_failure = 0;
  
      Real old_time = _system.time;
      Real old_deltat = _system.deltat;
      _system.deltat *= 0.5;
      core_time_solver->solve();
      core_time_solver->advance_timestep();
      core_time_solver->solve();

      single_norm = calculate_norm(_system, *_system.solution);
      if (!quiet)
        {
          std::cout << "Single norm = " << single_norm << std::endl;
        }

      // Reset the core_time_solver's reduce_deltat... value.
      core_time_solver->reduce_deltat_on_diffsolver_failure =
        this->reduce_deltat_on_diffsolver_failure;
  
      // But then back off just in case our advance_timestep() isn't
      // called.
      // FIXME: this probably doesn't work with multistep methods
      _system.get_vector("_old_nonlinear_solution") = *old_solution;
      _system.time = old_time;
      _system.deltat = old_deltat;

      // Find the relative error
      *double_solution -= *(_system.solution);
      error_norm  = calculate_norm(_system, *double_solution);
      relative_error = error_norm / _system.deltat /
        std::max(double_norm, single_norm);

      // If the relative error makes no sense, we're done
      if (!double_norm && !single_norm)
        return;

      if (!quiet)
        {
          std::cout << "Error norm = " << error_norm << std::endl;
          std::cout << "Local relative error = "
		    << (error_norm /
                        std::max(double_norm, single_norm))
                    << std::endl;
          std::cout << "Global relative error = "
		    << (error_norm / _system.deltat / 
                        std::max(double_norm, single_norm)) 
                    << std::endl;
          std::cout << "old delta t = " << _system.deltat << std::endl;
        }

      // If we haven't met our upper error tolerance, we'll have to
      // repeat this timestep entirely
      if (this->upper_tolerance && relative_error > this->upper_tolerance)
        {
	  // Reset the initial guess for our next try
	  *(_system.solution) =
            _system.get_vector("_old_nonlinear_solution");

          // Chop delta t in half
          _system.deltat /= 2.;

          if (!quiet)
            {
              std::cout << "Failed to meet upper error tolerance" 
                        << std::endl;
              std::cout << "Retrying with delta t = "
                        << _system.deltat << std::endl;
            }
        }
      else
        max_tolerance_met = true;
    }

  
  // Otherwise, compare the relative error to the tolerance
  // and adjust deltat
  last_deltat = _system.deltat;

  const Real global_shrink_or_growth_factor =
    std::pow(this->target_tolerance / relative_error,
	     1. / core_time_solver->error_order());

  const Real local_shrink_or_growth_factor =
    std::pow(this->target_tolerance /
	     (error_norm/std::max(double_norm, single_norm)),
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
 
  Real shrink_or_growth_factor =
    this->global_tolerance ? global_shrink_or_growth_factor :
                             local_shrink_or_growth_factor;

  if (this->max_growth && this->max_growth < shrink_or_growth_factor)
    {
      if (!quiet && this->global_tolerance)
        {
	  std::cout << "delta t is constrained by max_growth" << std::endl;
	}
        shrink_or_growth_factor = this->max_growth;
    }

  _system.deltat *= shrink_or_growth_factor;
  
  // Restrict deltat to max-allowable value if necessary
  if ((this->max_deltat != 0.0) && (_system.deltat > this->max_deltat))
    {
      if (!quiet)
	{
	  std::cout << "delta t is constrained by maximum-allowable delta t."
                    << std::endl;
	}
      _system.deltat = this->max_deltat;
    }

  // Restrict deltat to min-allowable value if necessary
  if ((this->min_deltat != 0.0) && (_system.deltat < this->min_deltat))
    {
      if (!quiet)
	{
	  std::cout << "delta t is constrained by minimum-allowable delta t."
                    << std::endl;
	}
      _system.deltat = this->min_deltat;
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
