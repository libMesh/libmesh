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
  // Using the connected_nodes set rather than point_neighbors() gives
  // us correct results even in corner cases, such as where two
  // elements meet only at a corner.  ;-)

  std::unordered_set<const Node *> connected_nodes;

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
      if (elem->processor_id() != p)
        coupled_elements.emplace(elem, nullcm);

      for (auto neigh : elem->neighbor_ptr_range())
        {
          if (neigh && neigh != remote_elem)
            {
#ifdef LIBMESH_ENABLE_AMR
              if (!neigh->active())
                {
                  std::vector<const Elem*> family;
                  neigh->active_family_tree_by_neighbor(family, elem);

                  for (const Elem * f : family)
                    if (f->processor_id() != p)
                      coupled_elements.emplace(f, nullcm);
                }
              else
#endif
                if (neigh->processor_id() != p)
                  coupled_elements.emplace(neigh, nullcm);
            }
        }

      // It is possible that a refined boundary element will not
      // touch any nodes of its interior_parent, in TRI3/TET4 and in
      // non-level-one rule cases.  So we can't just rely on node
      // connections to preserve interior_parent().  However, trying
      // to preserve interior_parent() manually only works if it's on
      // the same Mesh, which is *not* guaranteed!  So we'll
      // double-check later to make sure our interior parents are in
      // the mesh before we connect them.
      if (elem->dim() < LIBMESH_DIM &&
          elem->interior_parent() &&
          elem->interior_parent()->processor_id() != p)
        interior_parents.insert (elem->interior_parent());

      // Add nodes connected to active local elements
      for (auto n : elem->node_index_range())
        connected_nodes.insert (elem->node_ptr(n));
    }

  // Connect any interior_parents who are really in our mesh
  for (const auto & elem : _mesh->element_ptr_range())
    {
      std::unordered_set<const Elem *>::iterator ip_it =
        interior_parents.find(elem);

      if (ip_it != interior_parents.end())
        {
          coupled_elements.emplace(elem, nullcm);

          // Shrink the set ASAP to speed up subsequent searches
          interior_parents.erase(ip_it);
        }
    }

  // Connect any active elements which are connected to our range's
  // elements' nodes by addin elements connected to nodes on active
  // local elements.
  for (const auto & elem : _mesh->active_element_ptr_range())
    if (elem->processor_id() != p)
      for (auto & n : elem->node_ref_range())
        if (connected_nodes.count(&n))
          coupled_elements.emplace(elem, nullcm);
}

} // namespace libMesh
