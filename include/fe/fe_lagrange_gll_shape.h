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
#include "libmesh/int_range.h"

// C++ includes
#include <cmath>

// Tensor product evaluation shared by the Gauss-Lobatto nodal families. The families differ
// only in how a shape function index maps to a tensor index, so each family specializes
// fe_lagrange_gll_tensor_index and shares the rest.

namespace libMesh
{

/**
 * The pair of coordinates differentiated by each second derivative index \p j of
 * FE::shape_second_deriv().
 *
 * In dimension d only the first \f$d(d+1)/2\f$ entries are used: one on an edge, three on
 * a face, and all six on a cell.
 */
constexpr unsigned int fe_lagrange_gll_second_deriv_pairs[6][2] =
  {{0, 0},   // d^2 phi / dxi^2
   {0, 1},   // d^2 phi / dxi deta
   {1, 1},   // d^2 phi / deta^2
   {0, 2},   // d^2 phi / dxi dzeta
   {1, 2},   // d^2 phi / deta dzeta
   {2, 2}};  // d^2 phi / dzeta^2


/**
 * \returns The number of distinct second derivatives in dimension \p Dim.
 */
template <unsigned int Dim>
constexpr unsigned int fe_lagrange_gll_n_second_derivs()
{
  return Dim * (Dim + 1) / 2;
}

/**
 * Writes the tensor index of shape function \p i of order \p order into the first \p dim
 * entries of \p indices. Entry d is the index of the one-dimensional shape function in
 * coordinate d.
 *
 * Each Gauss-Lobatto nodal family specializes this function. Using a family without a
 * specialization is a link error.
 */
template <FEFamily T>
void fe_lagrange_gll_tensor_index(const Elem * elem,
                                  const Order order,
                                  const unsigned int dim,
                                  const unsigned int i,
                                  unsigned int * indices);


/**
 * All L2_LAGRANGE_GLL degrees of freedom belong to the element. They are numbered
 * lexicographically over the tensor grid, with the first coordinate varying fastest. QBase
 * orders the points of a tensor product rule the same way, so shape function i is the one
 * that is one at quadrature point i. The numbering doesn't depend on the element, so the
 * element argument is unused.
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

  for (const auto d : make_range(dim))
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
 * agreement comes from Elem::positive_edge_orientation and Elem::oriented_face_vertices, which
 * both elements compute from the positions of the same vertices: an edge is walked from its
 * lesser end, and a face from its least vertex along the lesser of the two edges leaving it.
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

      // positive_edge_orientation reports whether the edge's own first vertex is the greater
      // of the two, so walking from the lesser means starting at the second when it does
      const bool reversed = elem->positive_edge_orientation(e);

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

      // The face's vertices from its least vertex toward the lesser of that vertex's neighbors,
      // so the last is its other neighbor
      const std::vector<unsigned int> face_vertices = elem->oriented_face_vertices(f);

      unsigned int from[3], first[3], second[3];
      fe_lagrange_gll_vertex_index(*elem, degree, dim, face_vertices.front(), from);
      fe_lagrange_gll_vertex_index(*elem, degree, dim, face_vertices[1], first);
      fe_lagrange_gll_vertex_index(*elem, degree, dim, face_vertices.back(), second);

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
 * \returns Shape function \p i of family \p T in dimension \p Dim, differentiated once
 * with respect to each coordinate listed in the first \p n_deriv_coords entries of
 * \p deriv_coords.
 *
 * Each factor of the tensor product is a one-dimensional shape function, differentiated as
 * many times as its coordinate is listed. FE::shape() passes zero coordinates,
 * FE::shape_deriv() passes one, and FE::shape_second_deriv() passes two.
 */
template <FEFamily T, unsigned int Dim>
inline Real fe_lagrange_gll_tensor_product(const Elem * elem,
                                           const Order order,
                                           const unsigned int i,
                                           const Point & p,
                                           const unsigned int * const deriv_coords,
                                           const unsigned int n_deriv_coords)
{
  unsigned int indices[3];

  fe_lagrange_gll_tensor_index<T>(elem, order, Dim, i, indices);

  Real value = 1.;

  for (const auto d : make_range(Dim))
    {
      unsigned int deriv_order = 0;
      for (const auto k : make_range(n_deriv_coords))
        if (deriv_coords[k] == d)
          ++deriv_order;

      switch (deriv_order)
        {
        case 0:
          value *= fe_lagrange_gll_1D_shape(order, indices[d], p(d));
          break;

        case 1:
          value *= fe_lagrange_gll_1D_shape_deriv(order, indices[d], 0, p(d));
          break;

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
        case 2:
          value *= fe_lagrange_gll_1D_shape_second_deriv(order, indices[d], 0, p(d));
          break;
#endif

        default:
          libmesh_error_msg("Invalid derivative order " << deriv_order << ".");
        }
    }

  return value;
}

} // namespace libMesh

#endif // LIBMESH_FE_LAGRANGE_GLL_SHAPE_H
