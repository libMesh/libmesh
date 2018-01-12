// The libMesh Finite Element Library.
// Copyright (C) 2002-2018 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#ifndef LIBMESH_NO_SOLUTION_HISTORY_H
#define LIBMESH_NO_SOLUTION_HISTORY_H

// Local includes
#include "libmesh/solution_history.h"
#include "libmesh/auto_ptr.h" // libmesh_make_unique

namespace libMesh
{

/**
 * 'Save nothing' subclass of Solution History, this is the default.
 *
 * \author Vikram Garg
 * \date 2012
 * \brief For storing and retrieving timestep data.
 */
class NoSolutionHistory : public SolutionHistory
{
public:

  /**
   * Constructor
   */
  NoSolutionHistory() : SolutionHistory() {}

  /**
   * Destructor
   */
  virtual ~NoSolutionHistory() {}

  /**
   * Virtual function store which we will be overriding
   */
  virtual void store() libmesh_override;

  /**
   * Virtual function retrieve which we will be overriding
   */
  virtual void retrieve() libmesh_override;

  /**
   * Definition of the clone function needed for the setter function
   */
  virtual std::unique_ptr<SolutionHistory > clone() const libmesh_override
  {
    return libmesh_make_unique<NoSolutionHistory>();
  }
};

} // end namespace libMesh

#endif // LIBMESH_NO_SOLUTION_HISTORY_H
