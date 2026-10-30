// The libMesh Finite Element Library.
// Copyright (C) 2002-2026 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/periodic_point_neighbors.h"

#ifdef LIBMESH_ENABLE_PERIODIC

#include "libmesh/boundary_info.h"
#include "libmesh/elem.h"
#include "libmesh/periodic_boundaries.h"
#include "libmesh/remote_elem.h"

// C++ Includes
#include <set>
#include <vector>

namespace libMesh
{

void PeriodicPointNeighbors::operator()
  (const MeshBase::const_element_iterator & range_begin,
   const MeshBase::const_element_iterator & range_end,
   processor_id_type p,
   PeriodicPointNeighbors::map_type & coupled_elements)
{
  if (!_periodic_bcs || _periodic_bcs->empty())
    return;

  libmesh_assert(_mesh);

  std::unique_ptr<PointLocatorBase> point_locator = _mesh->sub_point_locator();

  std::set<const Elem *> periodic_elems_examined;
  const BoundaryInfo & binfo = _mesh->get_boundary_info();
  std::vector<boundary_id_type> appn_bcids;
  std::vector<const Elem *> active_periodic_neighbors;

  // This code is just for geometric coupling, so we use a null
  // CouplingMatrix pointer.
  CouplingMatrix * nullcm = nullptr;

  for (const auto & elem : as_range(range_begin, range_end))
    {
      for (const auto s : elem->side_index_range())
        {
          if (elem->neighbor_ptr(s))
            continue;

          const Elem * const equal_level_periodic_neigh = elem->topological_neighbor
            (s, *_mesh, *point_locator, _periodic_bcs);

          if (!equal_level_periodic_neigh || equal_level_periodic_neigh == remote_elem)
            continue;

#ifdef LIBMESH_ENABLE_AMR
          equal_level_periodic_neigh->active_family_tree_by_topological_neighbor(
            active_periodic_neighbors,
            elem,
            *_mesh,
            *point_locator,
            _periodic_bcs,
            /*reset=*/true);
#else
          active_periodic_neighbors = { equal_level_periodic_neigh };
#endif

          for (const Elem * const active_periodic_neigh : active_periodic_neighbors)
            {
              std::set <const Elem *> active_periodic_point_neighbors;

              // This fills point neighbors *including*
              // active_periodic_neigh. The documentation for this method
              // states that this will return *active* point neighbors
              active_periodic_neigh->find_point_neighbors(active_periodic_point_neighbors);

              for (const Elem * const appn : active_periodic_point_neighbors)
                {
                  // Don't need to ghost RemoteElem or an element we already own or an
                  // element we've already examined
                  if (appn == remote_elem || appn->processor_id() == p ||
                      periodic_elems_examined.count(appn))
                    continue;

                  // We only need to keep point neighbors that are along the periodic boundaries
                  bool on_periodic_boundary = false;
                  for (const auto appn_s : appn->side_index_range())
                    {
                      binfo.boundary_ids(appn, appn_s, appn_bcids);
                      for (const auto appn_bcid : appn_bcids)
                        if (_periodic_bcs->find(appn_bcid) != _periodic_bcs->end())
                          {
                            on_periodic_boundary = true;
                            goto jump;
                          }
                    }
                  jump:
                    if (on_periodic_boundary)
                      coupled_elements.emplace(appn, nullcm);

                    periodic_elems_examined.insert(appn);
                }
            }
        }
    }
}

void PeriodicPointNeighbors::mesh_reinit()
{
  // Unless we have periodic boundary conditions, we don't need
  // anything precomputed.
  if (!_periodic_bcs || _periodic_bcs->empty())
    return;

  // If we do have periodic boundary conditions, we'll need a master
  // point locator, so we'd better have a mesh to build it on.
  libmesh_assert(_mesh);

  // Make sure an up-to-date master point locator has been
  // constructed; we'll need to grab sub-locators soon.
  _mesh->sub_point_locator();
}

} // namespace libMesh

#endif // LIBMESH_ENABLE_PERIODIC
