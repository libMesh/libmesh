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

#ifndef LIBMESH_FE_REFERENCE_ELEMENT_TRAITS_H
#define LIBMESH_FE_REFERENCE_ELEMENT_TRAITS_H

#include "libmesh/libmesh.h" // invalid_uint

namespace libMesh
{

/**
 * The second-order elements' node maps are not independent facts.
 * Their edge_nodes_map is the first-order equivalent's edge_nodes_map
 * with the mid-edge node appended, and mid-edge nodes are numbered
 * after the vertices in edge order.  Their side_nodes_map row is the
 * first-order equivalent's vertex row, then the mid-edge node of each
 * consecutive vertex pair, then the face node if the side has one.
 *
 * The functions here compute those tables at compile time from the
 * first-order tables, so a second-order element class stores neither:
 * it holds the derived tables as static constexpr members and exposes
 * them through its usual side_nodes_map / edge_nodes_map names as
 * references to the underlying arrays.  Only the first-order tables
 * and the face-node rules are written by hand.
 */

/**
 * The value the element classes pad short side_nodes_map rows with,
 * e.g. the triangular sides of a Prism.
 */
static constexpr unsigned int unused_side_node = 99;

/**
 * A fixed-size table that constexpr functions can return.  Element
 * classes bind their side_nodes_map / edge_nodes_map references to
 * \p values.
 */
template <unsigned int Rows, unsigned int Cols>
struct ReferenceElementTable
{
  unsigned int values[Rows][Cols];
};

/**
 * The face-node rule for elements with no face nodes.
 */
constexpr unsigned int no_face_node (const unsigned int)
{
  return invalid_uint;
}

/**
 * \returns The edge_nodes_map of the second-order element whose
 * first-order equivalent is \p Linear: each edge's two vertices, then
 * its mid-edge node, numbered after the vertices in edge order.
 */
template <class Linear, unsigned int Edges>
constexpr ReferenceElementTable<Edges, 3>
derived_edge_nodes ()
{
  ReferenceElementTable<Edges, 3> t {};
  for (unsigned int e = 0; e != Edges; ++e)
    {
      t.values[e][0] = Linear::edge_nodes_map[e][0];
      t.values[e][1] = Linear::edge_nodes_map[e][1];
      t.values[e][2] = Linear::num_nodes + e;
    }
  return t;
}

/**
 * \returns The number of vertices on side \p s of a first-order
 * element, i.e. the entries of its side_nodes_map row that aren't
 * padding.
 */
template <class Linear>
constexpr unsigned int n_side_vertices (const unsigned int s)
{
  unsigned int n = 0;
  for (unsigned int k = 0; k != Linear::nodes_per_side; ++k)
    if (Linear::side_nodes_map[s][k] != unused_side_node)
      ++n;
  return n;
}

/**
 * \returns The third entry of the row of \p edges joining vertices \p a
 * and \p b, i.e. that edge's mid-edge node.
 */
template <unsigned int Edges>
constexpr unsigned int mid_edge_node (const unsigned int (&edges)[Edges][3],
                                      const unsigned int a,
                                      const unsigned int b)
{
  for (unsigned int e = 0; e != Edges; ++e)
    if ((edges[e][0] == a && edges[e][1] == b) ||
        (edges[e][0] == b && edges[e][1] == a))
      return edges[e][2];
  return invalid_uint;
}

/**
 * \returns The side_nodes_map of a 3D second-order element with
 * \p Sides sides of up to \p Cols nodes each, whose first-order
 * equivalent is \p Linear, whose edge_nodes_map is \p edges, and whose
 * \p face_node(s) is the node at the center of side \p s (or
 * \p invalid_uint if there is none).
 */
template <class Linear, unsigned int Sides, unsigned int Cols,
          unsigned int Edges, class FaceNode>
constexpr ReferenceElementTable<Sides, Cols>
derived_side_nodes (const unsigned int (&edges)[Edges][3],
                    FaceNode face_node)
{
  ReferenceElementTable<Sides, Cols> t {};
  for (unsigned int s = 0; s != Sides; ++s)
    {
      const unsigned int nv = n_side_vertices<Linear>(s);
      unsigned int n = 0;
      for (unsigned int k = 0; k != nv; ++k)
        t.values[s][n++] = Linear::side_nodes_map[s][k];
      for (unsigned int k = 0; k != nv; ++k)
        t.values[s][n++] = mid_edge_node(edges,
                                         Linear::side_nodes_map[s][k],
                                         Linear::side_nodes_map[s][(k+1) % nv]);
      if (face_node(s) != invalid_uint)
        t.values[s][n++] = face_node(s);
      for (; n != Cols; ++n)
        t.values[s][n] = unused_side_node;
    }
  return t;
}

/**
 * \returns The side_nodes_map of a 2D second-order element whose
 * first-order equivalent is \p Linear: each side's two vertices, then
 * its mid-side node, numbered after the vertices in side order.
 */
template <class Linear, unsigned int Sides>
constexpr ReferenceElementTable<Sides, 3>
derived_side_nodes ()
{
  ReferenceElementTable<Sides, 3> t {};
  for (unsigned int s = 0; s != Sides; ++s)
    {
      t.values[s][0] = Linear::side_nodes_map[s][0];
      t.values[s][1] = Linear::side_nodes_map[s][1];
      t.values[s][2] = Linear::num_nodes + s;
    }
  return t;
}

} // namespace libMesh

#endif // LIBMESH_FE_REFERENCE_ELEMENT_TRAITS_H
