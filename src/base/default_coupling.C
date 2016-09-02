// The libMesh Finite Element Library.
// Copyright (C) 2002-2016 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public  License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA



// Local Includes -----------------------------------
#include "libmesh/default_coupling.h"

#include "libmesh/elem.h"
#include "libmesh/periodic_boundaries.h"
#include "libmesh/remote_elem.h"

namespace libMesh
{

void DefaultCoupling::mesh_reinit()
{
  // Unless we have periodic boundary conditions, we don't need
  // to cache anything.
  if (!_periodic_bcs || _periodic_bcs->empty())
    return;

  // If we do have periodic boundary conditions, we'll need a master
  // point locator, so we'd better have a mesh to build it on.
  libmesh_assert(_mesh);

  // Make sure a master point locator has been constructed; we'll need
  // to grab sub-locators soon.
  _mesh->sub_point_locator();
}



void DefaultCoupling::operator()
  (const MeshBase::const_element_iterator & range_begin,
   const MeshBase::const_element_iterator & range_end,
   processor_id_type p,
   map_type & coupled_elements)
{
  LOG_SCOPE("operator()", "DefaultCoupling");

  UniquePtr<PointLocatorBase> point_locator;
  if (_periodic_bcs)
    {
      libmesh_assert(_mesh);
      point_locator = _mesh->sub_point_locator();
    }

  for (MeshBase::const_element_iterator elem_it = range_begin;
       elem_it != range_end; ++elem_it)
    {
      std::vector<const Elem *> active_neighbors;

      const Elem * const elem = *elem_it;

      if (elem->processor_id() != p)
        coupled_elements.insert (std::make_pair(elem,_dof_coupling));

      if (_couple_neighbor_dofs)
        for (unsigned int s=0; s<elem->n_sides(); s++)
          {
            const Elem *neigh = elem->neighbor_ptr(s);

            // If we have a neighbor here
            if (neigh)
              {
                // Mesh ghosting might ask us about what we want to
                // distribute along with non-local elements, and those
                // non-local elements might have remote neighbors, and
                // if they do then we can't say anything about them.
                if (neigh == remote_elem)
                  continue;
              }
#ifdef LIBMESH_ENABLE_PERIODIC
            // We might still have a periodic neighbor here
            else if (_periodic_bcs)
              {
                libmesh_assert(_mesh);

                neigh = elem->topological_neighbor
                  (s, *_mesh, *point_locator, _periodic_bcs);
              }
#endif

            // With no regular *or* periodic neighbors we have nothing
            // to do.
            if (!neigh)
              continue;

	    // With any kind of neighbor, we need to couple to all the
	    // active descendants on our side.
#ifdef LIBMESH_ENABLE_AMR
            neigh->active_family_tree_by_neighbor(active_neighbors,elem);
#else
            active_neighbors.clear();
            active_neighbors.push_back(neigh);
#endif

            for (std::size_t a=0; a != active_neighbors.size(); ++a)
              {
                const Elem * neighbor = active_neighbors[a];

                if (neighbor->processor_id() != p)
                  coupled_elements.insert
                    (std::make_pair(neighbor, _dof_coupling));
              }
          }
    }
}


} // namespace libMesh
