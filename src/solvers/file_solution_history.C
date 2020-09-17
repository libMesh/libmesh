// The libMesh Finite Element Library.
// Copyright (C) 2002-2020 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

// C++ include files that we need
#include <iostream>
// Local includes
#include "libmesh/file_solution_history.h"

#include "libmesh/diff_system.h"

#include <cmath>
#include <iterator>

namespace libMesh
{

/**
   * Constructor, reference to system to be passed by user, set the
   * stored_sols iterator to some initial value
   */
  FileSolutionHistory::FileSolutionHistory(System & system_)
  : stored_sols(stored_solutions.end()),
  _system(system_), localTimestamp(0),
  timeTotimestamp()
  {
    dual_solution_copies.resize(system_.n_qois());

    libmesh_experimental();
  }


FileSolutionHistory::~FileSolutionHistory ()
{
}

// This function finds, if it can, the entry where we're supposed to
// be storing data, leaves stored_sols unchanged if it cant find an entry
// with the key corresponding to time.
void FileSolutionHistory::find_stored_entry(Real time, bool storing)
{
  if (stored_solutions.begin() == stored_solutions.end())
    return;

  // We will use the map::lower_bound operation to find the key which
  // is the least upper bound among all existing keys for time.
  // (key before map::lower_bound) < time < map::lower_bound, one of these
  // should be within TOLERANCE of time (unless we are creating a new map entry)
  // If the lower bound iterator points to:
  // begin -> we are looking for the solution at the initial time
  // end -> we are creating a new entry
  // anything else, we are looking for an existing entry
  stored_solutions_iterator lower_bound_it = stored_solutions.lower_bound(time);

  // For the key right before the lower bound
  stored_solutions_iterator lower_bound_it_decremented;

  // If we are at end, we could be creating a new entry (depends on the storing bool), return
  // Otherwise, get a decremented iterator for the sandwich test
  if(lower_bound_it == stored_solutions.end())
  {
    // If we are storing and lower_bound_it points to stored_solutions.end(), we assume
    // that this is a brand new entry in the map. We leave stored_sols unchanged.
    if(storing)
    {
      return;
    }
    else
    {
      // We are trying to retrieve and none of the keys was an upper bound.
      // We could have a situation in which the time is greatest key + FPE.
      // So we can check the key before the end and see if it matches time, else we have an error.
      lower_bound_it = std::prev(lower_bound_it);
    }
  }
  else if(lower_bound_it == stored_solutions.begin()) // At the beginning, so we cant go back any further
  {
    stored_sols = stored_solutions.begin();
    return;
  }
  else // A decremented iterator, to perform the sandwich test for the key closest to time
  {
    lower_bound_it_decremented = std::prev(lower_bound_it);
  }

  // Set the stored sols iterator as per the key which is within TOLERANCE of time
  if(std::abs(lower_bound_it->first - time) < TOLERANCE)
  {
    stored_sols = lower_bound_it;
  }
  else if(std::abs(lower_bound_it_decremented->first - time) < TOLERANCE)
  {
    stored_sols = lower_bound_it_decremented;
  }
  else // Neither of the two candidate keys matched our time
  {
    if(storing) // If we are storing, this is fine, we need to create a new entry, so just return
    {
      return;
    }
    else // If we are not storing, then we expected to find something but didnt, so we have a problem
    {
      libmesh_error_msg("Failed to set stored solutions iterator to a valid value.");
    }
  }


}

// This functions writes the solution at the current system time to disk
void FileSolutionHistory::store(bool is_adjoint_solve, Real time)
{
  // This will map the stored_sols iterator to the current time
  this->find_stored_entry(time, true);

  // In an empty history we create the first entry
  if (stored_solutions.begin() == stored_solutions.end())
    {
      stored_solutions[time] = std::string();
      stored_sols = stored_solutions.begin();
    }

  // If we're past the end we can create a new entry
  if (time - stored_sols->first > TOLERANCE )
    {
#ifndef NDEBUG
      ++stored_sols;
      libmesh_assert (stored_sols == stored_solutions.end());
#endif
      stored_solutions[time] = std::string();
      stored_sols = stored_solutions.end();
      --stored_sols;
    }

  // If we're before the beginning we can create a new entry
  else if (stored_sols->first - time > TOLERANCE)
    {
      libmesh_assert (stored_sols == stored_solutions.begin());
      stored_solutions[time] = std::string();
      stored_sols = stored_solutions.begin();
    }

  // We don't support inserting entries elsewhere
  libmesh_assert(std::abs(stored_sols->first - time) < TOLERANCE);

  // The name of the file to in which we store the solution from the current timestep
  std::string & solution_filename = stored_sols->second;

  // This iterator will be used to check if we have already assigned a timestamp for this time key
  // If we have, we are in the adjoint loop, if we have not, we are in the primal loop
  timeTotimestamp_iterator = timeTotimestamp.find(time);

  // Associate the localTimestamp to the current time, if we are in the primal solve loop
  // Then increment the localTimestamp
  if(!is_adjoint_solve)
  {
    // Point solution_filename to the filename generated by the libMesh I/O object
    solution_filename = "primal.out.xda.";
    solution_filename += std::to_string(localTimestamp);

    // Write the current primal solution out to file
    _system.get_equation_systems().write (solution_filename, WRITE, EquationSystems::WRITE_DATA | EquationSystems::WRITE_ADDITIONAL_DATA);

    timeTotimestamp.insert( std::pair<Real, unsigned int>(time, localTimestamp) );

    ++localTimestamp;
  }
  else // We are in the adjoint time stepping loop
  {
    --localTimestamp;

    // For the adjoint solution, we reuse the timestamps generated earlier during the primal time march
    solution_filename = "adjoint.out.xda.";
    solution_filename += std::to_string(localTimestamp);

    // Write the current primal solution out to file
    _system.get_equation_systems().write (solution_filename, WRITE, EquationSystems::WRITE_DATA | EquationSystems::WRITE_ADDITIONAL_DATA);

  }


}

void FileSolutionHistory::retrieve(bool is_adjoint_solve, Real time)
{
  this->find_stored_entry(time, false);

  // To set the deltat while using adaptive timestepping, we will utilize
  // consecutive time entries in the stored solutions iterator
  Real _current_time = stored_sols->first;

  // If we are solving the adjoint, we are moving backwards, so decrement time
  // else we are moving forwards, so increment time
  if(is_adjoint_solve)
  {
    stored_solutions_iterator stored_sols_decrement_time = stored_sols;

    // Recovering deltats needs two different entries from the the
    // stored solutions map
    if(stored_sols_decrement_time != stored_solutions.begin())
    {
      stored_sols_decrement_time--;

      Real _decremented_time = stored_sols_decrement_time->first;

      try
      {
        dynamic_cast<DifferentiableSystem &>(_system).deltat = _current_time - _decremented_time;
      }
      catch(const std::bad_cast& e)
      {
        // For a non-diff system, only fixed time step sizes are supported as of now.
      }
    }
  }
  else
  {
    stored_solutions_iterator stored_sols_increment_time = stored_sols;

    // Recovering deltats needs two different entries from the the
    // stored solutions map
    if(stored_sols_increment_time != std::prev(stored_solutions.end()) )
    {
      stored_sols_increment_time++;

      Real _incremented_time = stored_sols_increment_time->first;

      try
      {
        dynamic_cast<DifferentiableSystem &>(_system).deltat = _incremented_time - _current_time;
      }
      catch(const std::bad_cast& e)
      {
        // For a non-diff system, only fixed time step sizes are supported as of now.
      }
    }
  }

  // Get the time at which we are recovering the solution vectors
  Real recovery_time = stored_sols->first;

  // Do we not have a solution for this time?  Then
  // there's nothing to do.
  if (stored_sols == stored_solutions.end() ||
      std::abs(recovery_time - time) > TOLERANCE)
    {
      //libMesh::out << "No more solutions to recover ! We are at time t = " <<
      //                     _system.time << std::endl;
      return;
    }


  // If we are doing an adjoint solve, we read in the primal solution,
  // but this overwrites the adjoint solution with zero, so we swap
  // the last adjoint solution out to prevent this zeroing
  if(is_adjoint_solve)
  {
    // Reading in the primal xdas overwrites the adjoint solution with zero
    // So swap to retain the old adjoint solution
    for (auto j : make_range(_system.n_qois()))
    {
      dual_solution_copies[j] = _system.get_adjoint_solution(j).clone();
    }

    // Read in the primal solution stored at the current recovery time from the disk
    _system.get_equation_systems().read (stored_sols->second, READ, EquationSystems::READ_DATA | EquationSystems::READ_ADDITIONAL_DATA);

    // Swap back the copy of the last adjoint solution back in place
    for (auto j : make_range(_system.n_qois()))
    {
      (_system.get_adjoint_solution(j)).swap(*dual_solution_copies[j]);
    }
  }
  else
  {
    // Read in the primal solution stored at the current recovery time from the disk
    _system.get_equation_systems().read (stored_sols->second, READ, EquationSystems::READ_DATA | EquationSystems::READ_ADDITIONAL_DATA);
  }

  // We need to call update to put system in a consistent state
  // with the solution that was read in
  _system.update();

}

void FileSolutionHistory::erase(Real time)
{
  // We cant erase the stored_sols iterator which is used in other places
  // So save its current value for the future
  stored_solutions_iterator stored_sols_last = stored_sols;

  // This will map the stored_sols iterator to the current time
  this->find_stored_entry(time, false);

  // map::erase behaviour is undefined if the iterator is pointing
  // to a non-existent element.
  libmesh_assert(stored_sols != stored_solutions.end());

  // We want to keep using the stored_sols iterator, so we have to create
  // a new one to erase the concerned entry
  stored_solutions_iterator stored_sols_copy = stored_sols;

  // If we're asking to erase the entry at stored_sols, then move stored_sols somewhere safer first
  if(stored_sols == stored_sols_last)
    stored_sols--;

  stored_solutions.erase(stored_sols_copy);

}

}
// End namespace libMesh
