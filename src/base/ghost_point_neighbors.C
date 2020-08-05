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


// Local Includes
#include "libmesh/ghost_point_neighbors.h"

#include "libmesh/elem.h"
#include "libmesh/remote_elem.h"

// C++ Includes
#include <unordered_set>

namespace libMesh
{

void GhostPointNeighbors::operator()
  (const MeshBase::const_element_iterator & range_begin,
   const MeshBase::const_element_iterator & range_end,
   processor_id_type p,
   GhostPointNeighbors::map_type & coupled_elements)
{
  libmesh_assert(_mesh);

  // Using a connected_nodes set rather than point_neighbors() would
  // give us correct results even in corner cases, such as where two
  // elements meet only at a corner.  ;-)
  //
  // std::unordered_set<const Node *> connected_nodes;
  //
  // Unfortunately, it's not a *fast* enough test to use on a
  // ReplicatedMesh (where it's O(Nprocs) times slower), or as a
  // coupling functor (where it's O(Nelem/Nprocs) times slower), or as
  // a coupling functor on a ReplicatedMesh (where you might as well
  // just give up now)

  // Links between boundary and interior elements on mixed
  // dimensional meshes also give us correct ghosting in this way.
  std::unordered_set<const Elem *> interior_parents;

  // We also preserve neighbors and their neighboring children for
  // active local elements - in most cases this is redundant with the
  // node check, but for non-conforming Tet4 meshes and
  // non-level-one-conforming 2D+3D meshes it is possible for an
  // element and its coarse neighbor to not share any vertices.
  //
  // We also preserve interior parents for active pid elements


  // This code is just for geometric coupling, so we use a null
  // CouplingMatrix pointer.  We'll declare that here so as to avoid
  // confusing the insert() calls later.
  CouplingMatrix * nullcm = nullptr;

  for (const auto & elem : as_range(range_begin, range_end))
    {
      libmesh_assert(_mesh->query_elem_ptr(elem->id()) == elem);

      if (elem->processor_id() != p)
        coupled_elements.emplace(elem, nullcm);

      std::set<const Elem *> elem_point_neighbors;
      elem->find_point_neighbors(elem_point_neighbors);

      for (const auto & neigh : elem_point_neighbors)
        coupled_elements.emplace(neigh, nullcm);

      // An interior_parent isn't on the same manifold so won't be
      // found as a point neighbor, and it may not share nodes so we
      // can't use a connected_nodes test.
      //
      // Trying to preserve interior_parent() only works if it's on
      // the same Mesh, which is *not* guaranteed!  So we'll
      // double-check.
      const Elem * ip = elem->interior_parent();
      if (ip && ip->processor_id() != p &&
          _mesh->query_elem_ptr(ip->id()) == ip)
        coupled_elements.emplace(ip, nullcm);
    }
}

} // namespace libMesh
