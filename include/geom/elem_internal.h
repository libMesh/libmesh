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

#ifndef LIBMESH_ELEM_INTERNAL_H
#define LIBMESH_ELEM_INTERNAL_H

// C++ includes
#include <vector>

namespace libMesh
{

/**
 * The ElemInternal namespace holds helper functions that are used
 * internally by the Elem class. These should not be called directly,
 * call the appropriate member functions on the Elem class instead.
 */
namespace ElemInternal
{

template<class T>
void
total_family_tree(T elem,
                  std::vector<T> & family,
                  bool reset)
{
  // Clear the vector if the flag reset tells us to.
  if (reset)
    family.clear();

  // Add this element to the family tree.
  family.push_back(elem);

  // Recurse into the elements children, if it has them.
  // Do not clear the vector any more.
  if (elem->has_children())
    for (auto & c : elem->child_ref_range())
      if (!c.is_remote())
        internal_total_family_tree (&c, family, false);
}

} // namespace ElemInternal
} // namespace libMesh

#endif
