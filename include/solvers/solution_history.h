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

#ifndef LIBMESH_SOLUTION_HISTORY_H
#define LIBMESH_SOLUTION_HISTORY_H

// Local Includes
#include "libmesh/system.h"
#include "libmesh/history_data.h"

namespace libMesh
{

/**
 * A SolutionHistory class that enables the storage and retrieval of
 * timesteps and (in the future) adaptive steps.
 * SolutionHistory interfaces between the time solver and HistoryData.
 * SolutionHistory organizes and manages the overall history record as a map,
 * while HistoryData manages individual I/O, memory or preprocessing operations
 * for the history data at a particular time.
 *
 * \author Vikram Garg
 * \date 2012
 * \brief For storing and retrieving timestep data.
 */
class SolutionHistory
{
public:

  /**
   * Constructor
   */
  SolutionHistory() : overwrite_previously_stored(false),
   stored_datum(stored_data.end()) {}

  /**
   * Destructor
   */
  virtual ~SolutionHistory () {}

  /**
   * Function to store a solution, pure virtual
   */
  virtual void store(bool is_adjoint_solve, Real time) = 0;

  /**
   * Function to retrieve a solution, pure virtual
   */
  virtual void retrieve(bool is_adjoint_solve, Real time) = 0;

  /**
   * Erase stored_data entry at time
   */
  void erase(Real time);

  /**
   * Cloning function for a std::unique_ptr, pure virtual, used in the
   * setter function in time_solver.C
   */
  virtual std::unique_ptr<SolutionHistory > clone() const = 0;

  /**
   * Turn on overwrite_previously_stored to overwrite any
   * already-saved data encountered during subsequent store() calls
   */
  void set_overwrite_previously_stored (bool val)
  { overwrite_previously_stored = val; }

protected:

  // Flag to specify whether we want to overwrite previously stored
  // vectors at a given time or not
  bool overwrite_previously_stored;

  // The abstract data structure that indexes and stores history data
  // Any type of history, solution, mesh or any future type will use
  // a realization of this map to store history information. This way,
  // we avoid multiple histories scatterred across the code.
  typedef std::map<Real, std::unique_ptr<HistoryData>> map_type;
  map_type stored_data;
  typedef map_type::iterator stored_data_iterator;
  stored_data_iterator stored_datum;

  // Function to locate entries at a given time
  // Behaviour depends on whether we are calling this function
  // while storing or retrieving/erasing entries.
  // While storing, if no entry in our map matches our time key,
  // we will create a new entry in the map. If we are not storing,
  // not matching a given time key implies an error.
  void find_stored_entry(Real time, bool storing = false);

};

} // end namespace libMesh

#endif // LIBMESH_SOLUTION_HISTORY_H
