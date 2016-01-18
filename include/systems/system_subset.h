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



#ifndef LIBMESH_SYSTEM_SUBSET_H
#define LIBMESH_SYSTEM_SUBSET_H

// Local Includes
#include "libmesh/reference_counted_object.h"

// C++ includes
#include <vector>

namespace libMesh
{

// Forward Declarations
class System;

/**
 * This is a base class for classes which represent subsets of the
 * dofs of a \p System.
 *
 * \author Tim Kroeger
 * \date 2010
 */
class SystemSubset : public ReferenceCountedObject<SystemSubset>
{
public:

  /**
   * Constructor.
   */
  explicit
  SystemSubset (const System & system);

  /**
   * Destructor.
   */
  virtual ~SystemSubset ();

  /**
   * Method that returns the actual set of dofs that the subset
   * consists of.  The result must contain local dofs on each
   * processor only and must not contain duplictates.
   */
  virtual const std::vector<unsigned int> & dof_ids () const = 0;

  /**
   * Returns the \p System to which we belong.
   */
  const System & get_system () const;

protected:

  /**
   * A reference to the \p System we belong to.
   */
  const System & _system;

private:
  /**
   * This isn't a copyable object, so let's make sure nobody tries.
   *
   * We won't even bother implementing this; we'll just make sure that
   * the compiler doesn't implement a default.
   */
  SystemSubset(const SystemSubset &);

  /**
   * This isn't a copyable object, so let's make sure nobody tries.
   *
   * We won't even bother implementing this; we'll just make sure that
   * the compiler doesn't implement a default.
   */
  SystemSubset & operator= (const SystemSubset &);
}; // class SystemSubset

} // namespace libMesh

#endif // LIBMESH_SYSTEM_SUBSET_H
