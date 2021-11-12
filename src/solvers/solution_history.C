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

#include "libmesh/solution_history.h"
#include <cmath>
#include <iterator>

namespace libMesh
{
 // This function finds, if it can, the entry where we're supposed to
 // be storing data, leaves stored_datum unchanged if it cant find an entry
 // with the key corresponding to time.
 void SolutionHistory::find_stored_entry(Real time, bool storing)
 {
  if (stored_data.begin() == stored_data.end())
    return;

  // We will use the map::lower_bound operation to find the key which
  // is the least upper bound among all existing keys for time.
  // (key before map::lower_bound) < time < map::lower_bound, one of these
  // should be within TOLERANCE of time (unless we are creating a new map entry)
  // If the lower bound iterator points to:
  // begin -> we are looking for the solution at the initial time
  // end -> we are creating a new entry
  // anything else, we are looking for an existing entry
  stored_data_iterator lower_bound_it = stored_data.lower_bound(time);

  // For the key right before the lower bound
  stored_data_iterator lower_bound_it_decremented;

  // If we are at end, we could be creating a new entry (depends on the storing bool), return
  // Otherwise, get a decremented iterator for the sandwich test
  if(lower_bound_it == stored_data.end())
  {
    // If we are storing and lower_bound_it points to stored_data.end(), we assume
    // that this is a brand new entry in the map. We leave stored_datum unchanged.
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
  else if(lower_bound_it == stored_data.begin()) // At the beginning, so we cant go back any further
  {
    stored_datum = stored_data.begin();
    return;
  }
  else // A decremented iterator, to perform the sandwich test for the key closest to time
  {
    lower_bound_it_decremented = std::prev(lower_bound_it);
  }

  // Set the stored sols iterator as per the key which is within TOLERANCE of time
  if(std::abs(lower_bound_it->first - time) < TOLERANCE)
  {
    stored_datum = lower_bound_it;
  }
  else if(std::abs(lower_bound_it_decremented->first - time) < TOLERANCE)
  {
    stored_datum = lower_bound_it_decremented;
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

 void SolutionHistory::erase(Real time)
 {
  // We cant erase the stored_datum iterator which is used in other places
  // So save its current value for the future
  stored_data_iterator stored_datum_last = stored_datum;
  //std::map<Real, unsigned int>::iterator timeTotimestamp_iterator_last = timeTotimestamp_iterator;

  // This will map the stored_datum iterator to the current time
  this->find_stored_entry(time, false);

  // map::erase behaviour is undefined if the iterator is pointing
  // to a non-existent element.
  libmesh_assert(stored_datum != stored_data.end());

  // We want to keep using the stored_datum iterator, so we have to create
  // a new one to erase the concerned entry
  stored_data_iterator stored_datum_copy = stored_datum;

  // If we're asking to erase the entry at stored_datum, then move stored_datum somewhere safer first
  if(stored_datum == stored_datum_last)
    stored_datum--;

  stored_data.erase(stored_datum_copy);
 }

}
// End namespace libMesh