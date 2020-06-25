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

// Local includes
#include "libmesh/memory_solution_history.h"

#include "libmesh/diff_system.h"

#include <cmath>

namespace libMesh
{

MemorySolutionHistory::~MemorySolutionHistory ()
{
}

// This function finds, if it can, the entry where we're supposed to
// be storing data
void MemorySolutionHistory::find_stored_entry(Real time)
{
  if (stored_solutions.begin() == stored_solutions.end())
    return;

  libmesh_assert (stored_sols != stored_solutions.end());

  // We will use the map::lower_bound operation to find the least key, which
  // is the least upper bound among all the keys for the time of our interest.
  // (key before map::lower_bound) < time < map::lower_bound, one of these
  // will be within TOLERANCE of time (unless we are creating a new map entry)
  stored_solutions_iterator lower_bound_it = stored_solutions.lower_bound(time);

  // Set the stored sols iterator to whichever key is within TOLERANCE of time
  if(std::abs(lower_bound_it->first - time) < TOLERANCE)
  {
    stored_sols = stored_solutions.find(time);
  }
  else if(std::abs((--lower_bound_it)->first - time) < TOLERANCE)
  {
    stored_sols = stored_solutions.find(time);
  }

}

// This functions saves all the 'projection-worthy' system vectors for
// future use
void MemorySolutionHistory::store(bool /* is_adjoint_solve */, Real time)
{
  this->find_stored_entry(time);

  // In an empty history we create the first entry
  if (stored_solutions.begin() == stored_solutions.end())
    {
      stored_solutions[time] = map_type();
      stored_sols = stored_solutions.begin();
    }

  // If we're past the end we can create a new entry
  if (time - stored_sols->first > TOLERANCE )
    {
#ifndef NDEBUG
      ++stored_sols;
      libmesh_assert (stored_sols == stored_solutions.end());
#endif
      stored_solutions[time] = map_type();
      stored_sols = stored_solutions.end();
      --stored_sols;
    }

  // If we're before the beginning we can create a new entry
  else if (stored_sols->first - time > TOLERANCE)
    {
      libmesh_assert (stored_sols == stored_solutions.begin());
      stored_solutions[time] = map_type();
      stored_sols = stored_solutions.begin();
    }

  // We don't support inserting entries elsewhere
  libmesh_assert(std::abs(stored_sols->first - time) < TOLERANCE);

  // Map of stored vectors for this solution step
  std::map<std::string, std::unique_ptr<NumericVector<Number>>> & saved_vectors = stored_sols->second;

  // Loop over all the system vectors
  for (System::vectors_iterator vec     = _system.vectors_begin(),
                                vec_end = _system.vectors_end();
       vec != vec_end; ++vec)
    {
      // The name of this vector
      const std::string & vec_name = vec->first;

      // If we haven't seen this vector before or if we have and
      // want to overwrite it
      if ((overwrite_previously_stored || !saved_vectors.count(vec_name)) &&
          // and if we think it's worth preserving
          _system.vector_preservation(vec_name))
        {
          // Then we save it.
          saved_vectors[vec_name] = vec->second->clone();
        }
    }

  // Of course, we will usually save the actual solution
  std::string _solution("_solution");
  if ((overwrite_previously_stored || !saved_vectors.count(_solution)) &&
      // and if we think it's worth preserving
      _system.project_solution_on_reinit())
    saved_vectors[_solution] = _system.solution->clone();
}

void MemorySolutionHistory::retrieve(bool is_adjoint_solve, Real time)
{
  this->find_stored_entry(time);

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

  // Print out what time we are recovering vectors at
  //    libMesh::out << "Recovering solution vectors at time: " <<
  //                 recovery_time << std::endl;

  // Do we not have a solution for this time?  Then
  // there's nothing to do.
  if (stored_sols == stored_solutions.end() ||
      std::abs(recovery_time - time) > TOLERANCE)
    {
      //libMesh::out << "No more solutions to recover ! We are at time t = " <<
      //                     _system.time << std::endl;
      return;
    }

  // Get the saved vectors at this timestep
  map_type & saved_vectors = stored_sols->second;

  map_type::iterator vec = saved_vectors.begin();
  map_type::iterator vec_end = saved_vectors.end();

  // Loop over all the saved vectors
  for (; vec != vec_end; ++vec)
    {
      // The name of this vector
      const std::string & vec_name = vec->first;

      // Get the vec_name entry in the saved vectors map and set the
      // current system vec[vec_name] entry to it
      if (vec_name != "_solution")
        _system.get_vector(vec_name) = *(vec->second);
    }

  // Of course, we will *always* have to get the actual solution
  std::string _solution("_solution");
  *(_system.solution) = *(saved_vectors[_solution]);
}

}
