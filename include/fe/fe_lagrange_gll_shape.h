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


#ifndef LIBMESH_FE_LAGRANGE_GLL_SHAPE_H
#define LIBMESH_FE_LAGRANGE_GLL_SHAPE_H

// Local includes
#include "libmesh/elem.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_order.h"
#include "libmesh/fe_lagrange_gll_shape_1D.h"

// C++ includes
#include <cmath>

// The tensor product that the Gauss-Lobatto nodal families on the tensor product element
// types build their shape functions from. The families differ in which tensor index a shape
// function index names, and in nothing else, so the evaluation lives here once.

namespace libMesh
{

/**
 * The pairs of coordinates that the second derivative index runs over, in the order libMesh
 * numbers them.
 *
 * The order is such that a given dimension uses the leading \f$d(d+1)/2\f$ of them, so one
 * table serves all three: an edge reads the first entry, a face the first three, and a cell
 * all six.
 */
constexpr unsigned int fe_lagrange_gll_second_deriv_pairs[6][2] =
  {{0, 0},   // d^2 phi / dxi^2
   {0, 1},   // d^2 phi / dxi deta
   {1, 1},   // d^2 phi / deta^2
   {0, 2},   // d^2 phi / dxi dzeta
   {1, 2},   // d^2 phi / deta dzeta
   {2, 2}};  // d^2 phi / dzeta^2


/**
 * \returns The number of second derivatives that a shape function of dimension \p Dim has.
 */
template <unsigned int Dim>
constexpr unsigned int fe_lagrange_gll_n_second_derivs()
{
  return Dim * (Dim + 1) / 2;
}

/**
 * Fills the leading \p dim entries of \p indices with the tensor index that shape function
 * \p i of order \p order names, the entry for a coordinate giving the one-dimensional shape
 * function that coordinate contributes.
 *
 * This is the whole of what separates the Gauss-Lobatto nodal families, so it is declared
 * here and specialized once per family. A family that has no specialization does not link,
 * which is what should happen.
 */
template <FEFamily T>
void fe_lagrange_gll_tensor_index(const Elem * elem,
                                  const Order order,
                                  const unsigned int dim,
                                  const unsigned int i,
                                  unsigned int * indices);


/**
 * The L2_LAGRANGE_GLL degrees of freedom all belong to the element, so their indices run over
 * the tensor grid in lexicographic order with the first coordinate varying fastest, which is
 * also the order in which QBase lays out a tensor product rule. The element carries no index
 * of its own here, so it goes unread.
 */
template <>
inline void fe_lagrange_gll_tensor_index<L2_LAGRANGE_GLL>(const Elem *,
                                                          const Order order,
                                                          const unsigned int dim,
                                                          const unsigned int i,
                                                          unsigned int * indices)
{
  const unsigned int n = static_cast<unsigned int>(order) + 1;

  unsigned int stride = 1;

  for (unsigned int d = 0; d != dim; ++d)
    {
      indices[d] = (i / stride) % n;
      stride *= n;
    }

  libmesh_assert_less (i, stride);
}

/**
 * Fills the leading \p dim entries of \p indices with the tensor index of vertex \p v.
 *
 * A vertex of a tensor product element sits at -1 or 1 in each coordinate of the reference
 * element, which is the first or the last of the interpolation points along it.
 */
inline void fe_lagrange_gll_vertex_index(const Elem & elem,
                                         const unsigned int degree,
                                         const unsigned int dim,
                                         const unsigned int v,
                                         unsigned int * indices)
{
  const Point q = elem.master_point(v);

  for (unsigned int d = 0; d != dim; ++d)
    {
      libmesh_assert_equal_to (std::abs(q(d)), 1.);
      indices[d] = (q(d) > 0.) ? degree : 0;
    }
}


/**
 * Fills \p indices with the tensor index reached by stepping \p step interpolation points from
 * the vertex whose index is \p from toward the vertex whose index is \p to.
 *
 * The two vertices are the ends of an edge, so they differ in one coordinate, and stepping
 * along it counts up from \p from or down from it according to which end it is. Deriving the
 * step from the vertices rather than from a table per edge is what keeps one expression
 * standing for every edge of every tensor product element.
 */
inline void fe_lagrange_gll_step_index(const unsigned int degree,
                                       const unsigned int dim,
                                       const unsigned int * from,
                                       const unsigned int * to,
                                       const unsigned int step,
                                       unsigned int * indices)
{
  unsigned int n_moved = 0;

  for (unsigned int d = 0; d != dim; ++d)
    if (from[d] != to[d])
      {
        indices[d] = (to[d] == degree) ? step : (degree - step);
        ++n_moved;
      }

  // The ends of an edge of a tensor product element differ in one coordinate, and the caller
  // relies on that to lay independent steps over each other
  libmesh_assert_equal_to (n_moved, 1u);
}


/**
 * The LAGRANGE_GLL degrees of freedom are distributed as HIERARCHIC distributes its own: one
 * on each vertex, degree-1 on each edge, (degree-1)^2 on each face, and the rest on the
 * element. libMesh numbers them by node, in node order, and then the element's own, so a
 * shape function index names an entity and a position within it.
 *
 * Two elements sharing an entity have to agree on which interpolation point each of its
 * degrees of freedom sits at, and each numbers that entity's vertices its own way. The
 * agreement comes from Elem::edge_orientation and Elem::face_orientation, which both elements
 * compute from the positions of the same vertices: an edge is walked from its lesser end, and
 * a face from its least vertex along the lesser of the two edges leaving it.
 */
template <>
inline void fe_lagrange_gll_tensor_index<LAGRANGE_GLL>(const Elem * elem,
                                                       const Order order,
                                                       const unsigned int dim,
                                                       const unsigned int i,
                                                       unsigned int * indices)
{
  libmesh_assert(elem);

  const unsigned int degree = static_cast<unsigned int>(order);
  const unsigned int n_v = elem->n_vertices();

  if (i < n_v)
    {
      fe_lagrange_gll_vertex_index(*elem, degree, dim, i, indices);
      return;
    }

  // Below a vertex apiece there is nothing else to hold
  libmesh_assert_greater (degree, 1);

  const unsigned int per_edge = degree - 1;

  // An edge element has no edges of its own, so what follows its vertices is its interior
  if (dim == 1)
    {
      indices[0] = i - n_v + 1;
      libmesh_assert_less (indices[0], degree);
      return;
    }

  const unsigned int n_e = elem->n_edges();

  if (i < n_v + n_e * per_edge)
    {
      const unsigned int e = (i - n_v) / per_edge;
      const unsigned int k = (i - n_v) % per_edge;

      // The node carrying this edge's degrees of freedom follows the vertices in node order
      libmesh_assert (elem->is_node_on_edge(n_v + e, e));

      // edge_orientation reports whether the edge's own first vertex is the greater of the
      // two, so walking from the lesser means starting at the second when it does
      const bool reversed = elem->edge_orientation(e);

      unsigned int from[3], to[3];
      fe_lagrange_gll_vertex_index(*elem, degree, dim,
                                   elem->local_edge_node(e, reversed ? 1 : 0), from);
      fe_lagrange_gll_vertex_index(*elem, degree, dim,
                                   elem->local_edge_node(e, reversed ? 0 : 1), to);

      for (unsigned int d = 0; d != dim; ++d)
        indices[d] = from[d];

      fe_lagrange_gll_step_index(degree, dim, from, to, k + 1, indices);
      return;
    }

  const unsigned int per_face = per_edge * per_edge;
  const unsigned int after_edges = n_v + n_e * per_edge;
  const unsigned int n_f = (dim > 2) ? elem->n_faces() : 0;

  if (i < after_edges + n_f * per_face)
    {
      const unsigned int f = (i - after_edges) / per_face;
      const unsigned int m = (i - after_edges) % per_face;

      // The node carrying this face's degrees of freedom follows the edge nodes
      libmesh_assert (elem->is_node_on_side(n_v + n_e + f, f));

      const std::vector<unsigned int> face_nodes = elem->nodes_on_side(f);
      const unsigned int n_fv = Elem::type_to_n_sides_map[elem->side_type(f)];

      // face_orientation carries twice the position of the face's least vertex plus one when
      // the vertex after it is the lesser of its two neighbors
      const unsigned int orientation = elem->face_orientation(f);
      const unsigned int least = orientation / 2;
      const bool forward = orientation % 2;

      const unsigned int origin = face_nodes[least];
      const unsigned int after  = face_nodes[(least + 1) % n_fv];
      const unsigned int before = face_nodes[(least + n_fv - 1) % n_fv];

      unsigned int from[3], first[3], second[3];
      fe_lagrange_gll_vertex_index(*elem, degree, dim, origin, from);
      fe_lagrange_gll_vertex_index(*elem, degree, dim, forward ? before : after, first);
      fe_lagrange_gll_vertex_index(*elem, degree, dim, forward ? after : before, second);

      for (unsigned int d = 0; d != dim; ++d)
        indices[d] = from[d];

      // The face's two directions move different coordinates, and the third holds the value
      // the whole face shares, so the two steps do not tread on each other
      fe_lagrange_gll_step_index(degree, dim, from, first,  m % per_edge + 1, indices);
      fe_lagrange_gll_step_index(degree, dim, from, second, m / per_edge + 1, indices);
      return;
    }

  // What is left belongs to the element, and runs over the interior of the tensor grid in
  // lexicographic order with the first coordinate varying fastest
  const unsigned int m = i - (after_edges + n_f * per_face);
  unsigned int stride = 1;

  for (unsigned int d = 0; d != dim; ++d)
    {
      indices[d] = (m / stride) % per_edge + 1;
      stride *= per_edge;
    }

  libmesh_assert_less (m, stride);
}


/**
 * \returns The \p Dim dimensional tensor product of one-dimensional Gauss-Lobatto nodal shape
 * functions that shape function \p i of family \p T names, differentiated once in each
 * coordinate that the leading \p n_derivs entries of \p derivs name.
 *
 * A value passes no coordinates, a first derivative one, and a second derivative two, so the
 * three shape routines of every family and dimension reach the basis through this one
 * function.
 */
template <FEFamily T, unsigned int Dim>
inline Real fe_lagrange_gll_shape(const Elem * elem,
                                  const Order order,
                                  const unsigned int i,
                                  const Point & p,
                                  const unsigned int * const derivs,
                                  const unsigned int n_derivs)
{
  unsigned int indices[3];

  fe_lagrange_gll_tensor_index<T>(elem, order, Dim, i, indices);

  Real value = 1.;

  for (unsigned int d = 0; d != Dim; ++d)
    {
      unsigned int n_d = 0;
      for (unsigned int k = 0; k != n_derivs; ++k)
        if (derivs[k] == d)
          ++n_d;

      switch (n_d)
        {
        case 0:
          value *= fe_lagrange_gll_1D_shape(order, indices[d], p(d));
          break;

        case 1:
          value *= fe_lagrange_gll_1D_shape_deriv(order, indices[d], 0, p(d));
          break;

        default:
          value *= fe_lagrange_gll_1D_shape_second_deriv(order, indices[d], 0, p(d));
          break;
        }
    }

  return value;
}

} // namespace libMesh

#endif // LIBMESH_FE_LAGRANGE_GLL_SHAPE_H
