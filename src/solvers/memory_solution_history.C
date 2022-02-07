// The libMesh Finite Element Library.
// Copyright (C) 2002-2022 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/memory_history_data.h"

#include "libmesh/diff_system.h"

#include <cmath>
#include <iterator>

namespace libMesh
{

MemorySolutionHistory::~MemorySolutionHistory ()
{
}

// This functions saves all the projected system vectors for
// future use
void MemorySolutionHistory::store(bool /* is_adjoint_solve */, Real time)
{
  this->find_stored_entry(time, true);

  // In an empty history we create the first entry
  if (stored_data.begin() == stored_data.end())
    {
      stored_data[time] = libmesh_make_unique<MemoryHistoryData>(_system);
      stored_datum = stored_data.begin();
    }

  // If we're past the end we can create a new entry
  if (time - stored_datum->first > TOLERANCE )
    {
#ifndef NDEBUG
      ++stored_datum;
      libmesh_assert (stored_datum == stored_data.end());
#endif
      stored_data[time] = libmesh_make_unique<MemoryHistoryData>(_system);
      stored_datum = stored_data.end();
      --stored_datum;
    }

  // If we're before the beginning we can create a new entry
  else if (stored_datum->first - time > TOLERANCE)
    {
      libmesh_assert (stored_datum == stored_data.begin());
      stored_data[time] = libmesh_make_unique<MemoryHistoryData>(_system);
      stored_datum = stored_data.begin();
    }

  // We don't support inserting entries elsewhere
  libmesh_assert(std::abs(stored_datum->first - time) < TOLERANCE);

  // First we handle the case of the initial data, this is the only case in which
  // stored_data will have size one
  if(stored_data.size() == 1)
  {
    // The initial data should only be stored once.
    (stored_datum->second)->store_initial_solution();
  }
  else if((stored_datum->second)->get_previously_stored() == false) // If we are not at the initial time, we are either creating a new entry or overwriting an existing one
  {
    (stored_datum->second)->store_primal_solution(stored_datum);
  }
  else // We are overwriting an existing history data
  {
    (stored_datum->second)->rewrite_stored_solution();
  }
}

void MemorySolutionHistory::retrieve(bool is_adjoint_solve, Real time)
{
  this->find_stored_entry(time, false);

  // If we are solving the adjoint, the timestep we need to move to the past step
  // is the one taken at that step to get to the current time.
  // At the initial time, be ready for the primal time march again.
  if(is_adjoint_solve)
  {
    if( stored_datum != stored_data.begin() )
    {
      stored_data_iterator stored_datum_past = stored_datum;
      stored_datum_past--;

      _system.deltat = (stored_datum_past->second)->get_deltat_at();
    }
    else
    {
      _system.deltat = (stored_datum->second)->get_deltat_at();
    }
  }
  else
  {
    if( stored_datum != std::prev(stored_data.end()) )
      _system.deltat = (stored_datum->second)->get_deltat_at();
    else
    {
      stored_data_iterator stored_datum_past = stored_datum;
      stored_datum_past--;

      _system.deltat = (stored_datum_past->second)->get_deltat_at();
    }

  }

  // Get the time at which we are recovering the solution vectors
  Real recovery_time = stored_datum->first;

  // Print out what time we are recovering vectors at
  //    libMesh::out << "Recovering solution vectors at time: " <<
  //                 recovery_time << std::endl;

  // Do we not have a solution for this time?  Then
  // there's nothing to do.
  if (stored_datum == stored_data.end() ||
      std::abs(recovery_time - time) > TOLERANCE)
    {
      //libMesh::out << "No more solutions to recover ! We are at time t = " <<
      //                     _system.time << std::endl;
      return;
    }

  (stored_datum->second)->retrieve_primal_solution();

  // We need to call update to put system in a consistent state
  // with the solution that was read in
  _system.update();
}

}
