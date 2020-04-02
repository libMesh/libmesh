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



#ifndef LIBMESH_FILE_SOLUTION_HISTORY_H
#define LIBMESH_FILE_SOLUTION_HISTORY_H

// Local includes
#include "libmesh/numeric_vector.h"
#include "libmesh/solution_history.h"
#include "libmesh/enum_xdr_mode.h"
#include "libmesh/equation_systems.h"
#include "libmesh/libmesh.h"
#include "libmesh/auto_ptr.h" // libmesh_make_unique

// C++ includes
#include <list>

namespace libMesh
{

/**
 * Subclass of Solution History that stores the solutions
 * and other important vectors onto disk.
 *
 * \author Vikram Garg
 * \date 2020
 * \brief Stores past solutions onto disk.
 */
class FileSolutionHistory : public SolutionHistory
{
public:

  /**
   * Constructor, reference to system to be passed by user, set the
   * stored_sols iterator to some initial value
   */
  FileSolutionHistory(System & system_) : stored_sols(stored_solutions.end()), _system(system_), localTimestamp(0), timeTotimestamp()
  { libmesh_experimental(); }

  /**
   * Destructor
   */
  ~FileSolutionHistory();

  /**
   * Virtual function store which we will be overriding to store timesteps
   */
  virtual void store(bool is_adjoint_solve) override;

  /**
   * Virtual function retrieve which we will be overriding to retrieve timesteps
   */
  virtual void retrieve(bool is_adjoint_solve) override;

  /**
   * Typedef for Stored Solutions iterator, a list of pairs of the 
   * system time and filenames of the stored solutions
   */  
  typedef std::pair<Real, std::string> pair_type;
  typedef std::list<pair_type> list_type;  
  typedef list_type::iterator stored_solutions_iterator;

  /**
   * Definition of the clone function needed for the setter function
   */
  virtual std::unique_ptr<SolutionHistory > clone() const override
  {
    return libmesh_make_unique<FileSolutionHistory>(_system);
  }

private:

  // This list of pairs will hold the timestamp and filename of each stored solution
  list_type stored_solutions;

  // The stored solutions iterator
  stored_solutions_iterator stored_sols;

  // A helper function to locate entries at a given time
  void find_stored_entry();

  // A system reference
  System & _system ;

  // A 'timestamp' that belongs specifically to FSH, this will be used to generate filenames
  unsigned int localTimestamp;

  // To assign filenames a timestamp, we will maintain a datastructure within
  // FileSolutionHistory which will map system.time to 'timestamps'
  std::map<Real, unsigned int> timeTotimestamp;

  // An iterator for the timeTotimestamp map, will help member functions distinguish
  // between primal and adjoint time loops
  std::map<Real, unsigned int>::iterator timeTotimestamp_iterator;
};

} // end namespace libMesh

#endif // LIBMESH_FILE_SOLUTION_HISTORY_H
