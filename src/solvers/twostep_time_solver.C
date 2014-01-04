// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


#include "libmesh/twostep_time_solver.h"
#include "libmesh/diff_system.h"
#include "libmesh/euler_solver.h"
#include "libmesh/numeric_vector.h"

namespace libMesh
{



TwostepTimeSolver::TwostepTimeSolver (sys_type& s)
  : AdaptiveTimeSolver(s)

{
  // We start with a reasonable time solver: implicit Euler
  core_time_solver.reset(new EulerSolver(s));
}



TwostepTimeSolver::~TwostepTimeSolver ()
{
}



void TwostepTimeSolver::solve()
{
  libmesh_assert(core_time_solver.get());

  // The core_time_solver will handle any first_solve actions
  first_solve = false;

  // We may have to repeat timesteps entirely if our error is bad
  // enough
  bool max_tolerance_met = false;

  // Calculating error values each time
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
          libMesh::out << "\n === Computing adaptive timestep === "
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
	  libMesh::out << "Double norm = " << double_norm << std::endl;
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
          libMesh::out << "Single norm = " << single_norm << std::endl;
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
          libMesh::out << "Error norm = " << error_norm << std::endl;
          libMesh::out << "Local relative error = "
		        << (error_norm /
                            std::max(double_norm, single_norm))
                        << std::endl;
          libMesh::out << "Global relative error = "
		        << (error_norm / _system.deltat /
                            std::max(double_norm, single_norm))
                        << std::endl;
          libMesh::out << "old delta t = " << _system.deltat << std::endl;
        }

      // If our upper tolerance is negative, that means we want to set
      // it based on the first successful time step
      if (this->upper_tolerance < 0)
        this->upper_tolerance = -this->upper_tolerance * relative_error;

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
              libMesh::out << "Failed to meet upper error tolerance"
                            << std::endl;
              libMesh::out << "Retrying with delta t = "
                            << _system.deltat << std::endl;
            }
        }
      else
        max_tolerance_met = true;
    }


  // Otherwise, compare the relative error to the tolerance
  // and adjust deltat
  last_deltat = _system.deltat;

  // If our target tolerance is negative, that means we want to set
  // it based on the first successful time step
  if (this->target_tolerance < 0)
    this->target_tolerance = -this->target_tolerance * relative_error;

  const Real global_shrink_or_growth_factor =
    std::pow(this->target_tolerance / relative_error,
	     static_cast<Real>(1. / core_time_solver->error_order()));

  const Real local_shrink_or_growth_factor =
    std::pow(this->target_tolerance /
	     (error_norm/std::max(double_norm, single_norm)),
	     static_cast<Real>(1. / (core_time_solver->error_order()+1.)));

  if (!quiet)
    {
      libMesh::out << "The global growth/shrink factor is: "
		    << global_shrink_or_growth_factor << std::endl;
      libMesh::out << "The local growth/shrink factor is: "
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
	  libMesh::out << "delta t is constrained by max_growth" << std::endl;
	}
        shrink_or_growth_factor = this->max_growth;
    }

  _system.deltat *= shrink_or_growth_factor;

  // Restrict deltat to max-allowable value if necessary
  if ((this->max_deltat != 0.0) && (_system.deltat > this->max_deltat))
    {
      if (!quiet)
	{
	  libMesh::out << "delta t is constrained by maximum-allowable delta t."
                        << std::endl;
	}
      _system.deltat = this->max_deltat;
    }

  // Restrict deltat to min-allowable value if necessary
  if ((this->min_deltat != 0.0) && (_system.deltat < this->min_deltat))
    {
      if (!quiet)
	{
	  libMesh::out << "delta t is constrained by minimum-allowable delta t."
                        << std::endl;
	}
      _system.deltat = this->min_deltat;
    }

  if (!quiet)
    {
      libMesh::out << "new delta t = " << _system.deltat << std::endl;
    }
}

} // namespace libMesh
