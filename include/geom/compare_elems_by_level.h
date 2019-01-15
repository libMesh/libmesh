// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_COMPARE_ELEMS_BY_LEVEL_H
#define LIBMESH_COMPARE_ELEMS_BY_LEVEL_H


// Local Includes
#include "libmesh/elem.h"

// C++ Includes
#include <numeric>
#include <set>




namespace libMesh {

/**
 * Specific weak ordering for Elem *'s to be used in a set.
 * We use the id, but first sort by level.  This guarantees
 * when traversing the set from beginning to end the lower
 * level (parent) elements are encountered first.
 */
struct CompareElemIdsByLevel
{
  bool operator()(const Elem * a,
                  const Elem * b) const
  {
    libmesh_assert (a);
    libmesh_assert (b);
    const unsigned int
      al = a->level(), bl = b->level();
    const dof_id_type
      aid = a->id(),   bid = b->id();

    return (al == bl) ? aid < bid : al < bl;
  }
};


} // namespace libMesh

#endif // LIBMESH_COMPARE_ELEMS_BY_LEVEL_H
