// The libMesh Finite Element Library.
// Copyright (C) 2002-2022 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#include "libmesh/intersection_tools.h"

#include "libmesh/point.h"
#include "libmesh/int_range.h"
#include "libmesh/elem.h"
#include "libmesh/elem_corner.h"

namespace libMesh::IntersectionTools
{

WithinSegmentResult within_segment(const Point & s1,
                                   const Point & s2,
                                   const Point & p,
                                   const Real tol)
{
  libmesh_assert(!s1.relative_fuzzy_equals(s2, tol));

  // First, check whether or not the points are collinear
  const auto l = s1 - s2;
  const auto l_norm = l.norm();
  const auto diff1 = p - s1;
  const auto diff1_norm = diff1.norm();
  if ((std::abs(l * diff1) / (l_norm * diff1_norm)) < ((Real)1 - tol))
    return NOT_WITHIN;

  // If the points are collinear, make sure that p is
  // somewhere between [s1, s2]
  const auto diff2_norm = (p - s2).norm();
  const auto tol_scaled = tol * l_norm;
  if (std::abs(diff1_norm + diff2_norm - l_norm) < tol_scaled)
  {
    if (diff1_norm < tol_scaled)
      return AT_BEGINNING;
    if (diff2_norm < tol_scaled)
      return AT_END;
    return BETWEEN;
  }

  return NOT_WITHIN;
}

bool collinear(const Point & p1,
               const Point & p2,
               const Point & p3,
               const Real tol)
{
  // (p1 -> p2) X (p1 - > p3) should be == the zero vector
  const auto p1p2_cross_p1p3 = (p2 - p1).cross(p3 - p1);
  for (const auto i : make_range(LIBMESH_DIM))
    if (std::abs(p1p2_cross_p1p3(i)) > tol)
      return false;
  return true;
}

bool within_edge_on_side(const Elem & elem,
                         const Point & p,
                         const unsigned short s,
                         ElemCorner & corner,
                         const bool linearize,
                         const Real tol)
{
  libmesh_assert_less(s, elem.n_sides());
  libmesh_assert(corner.is_invalid());
  libmesh_assert_equal_to(elem.dim(), 3);

  // For higher order than linear without linearization, make sure
  // that our edges are collinear
  if (elem.default_order() > 1 && !linearize)
    for (const auto e : elem.edge_index_range())
    {
      // we should expect at least 3 nodes/edges for our higher order elems
      libmesh_assert_greater_equal(elem.n_nodes_on_edge(e), 3);

      const unsigned int * edge_nodes_map = elem.nodes_on_edge_ptr(e);
      if (!collinear(elem.point(edge_nodes_map[0]),
                     elem.point(edge_nodes_map[1]),
                     elem.point(edge_nodes_map[2])))
        libmesh_error_msg("Failed to use Cell::without_edge_on_side without linearization "
                          "because an edge was found that is not collinear.");
    }

  const auto vs = elem.n_vertices_on_side(s);
  const unsigned int * side_nodes_map = elem.nodes_on_side_ptr(s);

  // side_nodes_map will point to something like [v0, v1, v2, v3]
  // With the loop below, we're going to check (in this order):
  // [v3 -> v0], [v0 -> v1], [v1 -> v2], [v2 -> v3]
  auto last_v = side_nodes_map[vs - 1];
  for (const auto side_v : make_range(vs))
  {
    const auto other_v = side_nodes_map[side_v];
    const auto within_result = within_segment(elem.point(last_v), elem.point(other_v), p, tol);
    if (within_result == NOT_WITHIN)
    {
      last_v = side_nodes_map[side_v];
      continue;
    }

    if (within_result == BETWEEN)
    {
      corner.set_edge(last_v, other_v);
      libmesh_assert(corner.build_edge(elem)->close_to_point(p, tol));
    }
    else if (within_result == AT_BEGINNING)
    {
      corner.set_vertex(last_v);
      libmesh_assert(elem.point(last_v).absolute_fuzzy_equals(p, tol));
    }
    else
    {
      corner.set_vertex(other_v);
      libmesh_assert(elem.point(other_v).absolute_fuzzy_equals(p, tol));
    }
    return true;
  }

  return false;
}

} // namespace libMesh::IntersectionTools


