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

#include "libmesh/fe_base_impl.h"

// Anonymous namespace, for a helper function for periodic boundary
// constraint calculations and builders
namespace
{
using namespace libMesh;

#ifdef LIBMESH_ENABLE_PERIODIC

// Find the "primary" element around a boundary point:
const Elem * primary_boundary_point_neighbor(const Elem * elem,
                                             const Point & p,
                                             const BoundaryInfo & boundary_info,
                                             const std::set<boundary_id_type> & boundary_ids)
{
  // If we don't find a better alternative, the user will have
  // provided the primary element
  const Elem * primary = elem;

  // Container to catch boundary IDs passed back by BoundaryInfo.
  std::vector<boundary_id_type> bc_ids;

  // Pull object out of the loop to reduce heap operations
  std::unique_ptr<const Elem> periodic_side;

  std::set<const Elem *> point_neighbors;
  elem->find_point_neighbors(p, point_neighbors);
  for (const auto & pt_neighbor : point_neighbors)
    {
      // If this point neighbor isn't at least
      // as coarse as the current primary elem, or if it is at
      // the same level but has a lower id, then
      // we won't defer to it.
      if ((pt_neighbor->level() > primary->level()) ||
          (pt_neighbor->level() == primary->level() &&
           pt_neighbor->id() < primary->id()))
        continue;

      // Otherwise, we will defer to the point neighbor, but only if
      // one of its sides is on a relevant boundary and that side
      // contains this vertex
      bool vertex_on_periodic_side = false;
      for (auto ns : pt_neighbor->side_index_range())
        {
          boundary_info.boundary_ids (pt_neighbor, ns, bc_ids);

          bool on_relevant_boundary = false;
          for (const auto & id : boundary_ids)
            if (std::find(bc_ids.begin(), bc_ids.end(), id) != bc_ids.end())
              on_relevant_boundary = true;

          if (!on_relevant_boundary)
            continue;

          pt_neighbor->build_side_ptr(periodic_side, ns);
          if (!periodic_side->contains_point(p))
            continue;

          vertex_on_periodic_side = true;
          break;
        }

      if (vertex_on_periodic_side)
        primary = pt_neighbor;
    }

  return primary;
}

// Find the "primary" element around a boundary edge:
const Elem * primary_boundary_edge_neighbor(const Elem * elem,
                                            const Point & p1,
                                            const Point & p2,
                                            const BoundaryInfo & boundary_info,
                                            const std::set<boundary_id_type> & boundary_ids)
{
  // If we don't find a better alternative, the user will have
  // provided the primary element
  const Elem * primary = elem;

  std::set<const Elem *> edge_neighbors;
  elem->find_edge_neighbors(p1, p2, edge_neighbors);

  // Container to catch boundary IDs handed back by BoundaryInfo
  std::vector<boundary_id_type> bc_ids;

  // Pull object out of the loop to reduce heap operations
  std::unique_ptr<const Elem> periodic_side;

  for (const auto & e_neighbor : edge_neighbors)
    {
      // If this edge neighbor isn't at least
      // as coarse as the current primary elem, or if it is at
      // the same level but has a lower id, then
      // we won't defer to it.
      if ((e_neighbor->level() > primary->level()) ||
          (e_neighbor->level() == primary->level() &&
           e_neighbor->id() < primary->id()))
        continue;

      // Otherwise, we will defer to the edge neighbor, but only if
      // one of its sides is on this periodic boundary and that
      // side contains this edge
      bool vertex_on_periodic_side = false;
      for (auto ns : e_neighbor->side_index_range())
        {
          boundary_info.boundary_ids (e_neighbor, ns, bc_ids);

          bool on_relevant_boundary = false;
          for (const auto & id : boundary_ids)
            if (std::find(bc_ids.begin(), bc_ids.end(), id) != bc_ids.end())
              on_relevant_boundary = true;

          if (!on_relevant_boundary)
            continue;

          e_neighbor->build_side_ptr(periodic_side, ns);
          if (!(periodic_side->contains_point(p1) &&
                periodic_side->contains_point(p2)))
            continue;

          vertex_on_periodic_side = true;
          break;
        }

      if (vertex_on_periodic_side)
        primary = e_neighbor;
    }

  return primary;
}

#endif // LIBMESH_ENABLE_PERIODIC

} // anonymous namespace

namespace libMesh
{
// ------------------------------------------------------------
// Explicit instantiations
template class FEGenericBase<Real, Real>;
template class FEGenericBase<RealGradient, Real>;
}
