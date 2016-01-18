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



#ifndef LIBMESH_MEMORY_SOLUTION_HISTORY_H
#define LIBMESH_MEMORY_SOLUTION_HISTORY_H

// Local includes
#include "libmesh/numeric_vector.h"
#include "libmesh/solution_history.h"
#include <list>

namespace libMesh
{
/**
 * Subclass of Solution History that stores the solutions
 * and other important vectors in memory.
 */
class MemorySolutionHistory : public SolutionHistory
{
public:

  /**
   * Constructor, reference to system to be passed by user, set the
   * stored_sols iterator to some initial value
   */
  MemorySolutionHistory(System & system_) : stored_sols(stored_solutions.end()), _system(system_)
  { libmesh_experimental(); }

  /**
   * Destructor
   */
  ~MemorySolutionHistory();

  /**
   * Virtual function store which we will be overriding to store timesteps
   */
  virtual void store() libmesh_override;

  /**
   * Virtual function retrieve which we will be overriding to retrieve timesteps
   */
  virtual void retrieve() libmesh_override;

  /**
   * Typedef for Stored Solutions iterator, a list of pairs of the current
   * system time, map of strings and saved vectors
   */
  typedef std::list<std::pair<Real, std::map<std::string, NumericVector<Number> *> > >::iterator stored_solutions_iterator;

  /**
   * Definition of the clone function needed for the setter function
   */
  virtual UniquePtr<SolutionHistory > clone() const libmesh_override
  {
    return UniquePtr<SolutionHistory >(new MemorySolutionHistory(_system));
  }

private:

  // This list of pairs will hold the current time and stored vectors
  // from each timestep
  std::list<std::pair<Real, std::map<std::string, NumericVector<Number> *> > > stored_solutions;

  // The stored solutions iterator
  stored_solutions_iterator stored_sols;

  // A helper function to locate entries at a given time
  void find_stored_entry();

  // A system reference
  System & _system ;
};

} // end namespace libMesh

#endif // LIBMESH_MEMORY_SOLUTION_HISTORY_H
