// The libMesh Finite Element Library.
// Copyright (C) 2002-2016 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

namespace libMesh
{

/**
 * A SolutionHistory class that enables the storage and retrieval of timesteps
 * and (in the future) adaptive steps
 */
class SolutionHistory
{
public:

  /**
   * Constructor
   */
  SolutionHistory() :
    overwrite_previously_stored(false) {}

  /**
   * Destructor
   */
  virtual ~SolutionHistory () {}

  /**
   * Function to store a solution, pure virtual
   */
  virtual void store() = 0;

  /**
   * Function to retrieve a solution, pure virtual
   */
  virtual void retrieve() = 0;

  /**
   * Cloning function for a UniquePtr, pure virtual, used in the
   * setter function in time_solver.C
   */
  virtual UniquePtr<SolutionHistory > clone() const = 0;

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
};

} // end namespace libMesh

#endif // LIBMESH_SOLUTION_HISTORY_H
