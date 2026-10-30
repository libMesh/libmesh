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
#include <array>
#include <cmath>
#include <vector>

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
 * The tensor index of a shape function. Entry d is the index of the one-dimensional shape
 * function in coordinate d. Entries at and above the element dimension are zero.
 */
typedef std::array<unsigned int, 3> GLLTensorIndex;


/**
 * \returns The tensor index of shape function \p i of order \p order on an element of
 * dimension \p dim.
 *
 * Each Gauss-Lobatto nodal family specializes this function. Using a family without a
 * specialization is a link error.
 */
template <FEFamily T>
GLLTensorIndex fe_lagrange_gll_tensor_index(const Elem * elem,
                                            const Order order,
                                            const unsigned int dim,
                                            const unsigned int i);


/**
 * All L2_LAGRANGE_GLL degrees of freedom belong to the element. They are numbered
 * lexicographically over the tensor grid, with the first coordinate varying fastest. QBase
 * orders the points of a tensor product rule the same way, so shape function i is the one
 * that is one at quadrature point i. The numbering doesn't depend on the element, so the
 * element argument is unused.
 */
template <>
inline GLLTensorIndex fe_lagrange_gll_tensor_index<L2_LAGRANGE_GLL>(const Elem *,
                                                                    const Order order,
                                                                    const unsigned int dim,
                                                                    const unsigned int i)
{
  const unsigned int n = static_cast<unsigned int>(order) + 1;

  GLLTensorIndex indices = {};
  unsigned int stride = 1;

  for (const auto d : make_range(dim))
    {
      indices[d] = (i / stride) % n;
      stride *= n;
    }

  libmesh_assert_less (i, stride);

  return indices;
}


/**
 * \returns The index of the one-dimensional interpolation point of LAGRANGE_GLL shape
 * function \p i of order \p order on an edge element.
 *
 * The first two shape functions belong to the vertices at -1 and 1, which are the first and
 * the last points. The interior points sit in-between in increasing order. The numbering is the same
 * on every edge element, so it needs no Elem.
 */
inline unsigned int fe_lagrange_gll_1D_index(const Order order,
                                             const unsigned int i)
{
  const unsigned int degree = static_cast<unsigned int>(order);

  libmesh_assert_less_equal (i, degree);

  if (i == 0)
    return 0;
  if (i == 1)
    return degree;
  return i - 1;
}


/**
 * \returns The tensor index of vertex \p v of \p elem.
 *
 * Each reference coordinate of a vertex of a tensor product element is -1 or 1. These are the
 * first and the last interpolation points in that coordinate.
 */
inline GLLTensorIndex fe_lagrange_gll_vertex_index(const Elem & elem,
                                                   const unsigned int degree,
                                                   const unsigned int v)
{
  const Point q = elem.master_point(v);

  GLLTensorIndex indices = {};

  for (const auto d : make_range(elem.dim()))
    {
      libmesh_assert_equal_to (std::abs(q(d)), 1.);
      indices[d] = (q(d) > 0.) ? degree : 0;
    }

  return indices;
}


/**
 * \returns The tensor index that lies \p k interpolation points from \p index, in the
 * direction that leads from vertex \p from to vertex \p to.
 *
 * The two vertices are the ends of an edge. Their tensor indices differ in one coordinate,
 * where one is 0 and the other is the degree. The result adds \p k to that coordinate of
 * \p index, or subtracts it if \p to is the end at 0. The other coordinates stay the same.
 *
 * An edge degree of freedom is one step from the start vertex of its edge. A face degree of
 * freedom is two steps from the start vertex of its face, one along each of two edges.
 */
inline GLLTensorIndex fe_lagrange_gll_advance_index(GLLTensorIndex index,
                                                    const GLLTensorIndex from,
                                                    const GLLTensorIndex to,
                                                    const unsigned int k)
{
  unsigned int n_moved = 0;

  for (const auto d : index_range(index))
    if (from[d] != to[d])
      {
        index[d] = (to[d] > from[d]) ? index[d] + k : index[d] - k;
        ++n_moved;
      }

  libmesh_assert_equal_to (n_moved, 1u);
  libmesh_ignore(n_moved);

  return index;
}


/**
 * An entity is a vertex, an edge, or a face of the element, or its interior. LAGRANGE_GLL
 * places its degrees of freedom on the entities as HIERARCHIC does:
 * - one on each vertex,
 * - degree-1 on each edge,
 * - (degree-1)^2 on each face,
 * - the rest on the interior.
 *
 * We number the vertex, edge, and face degrees of freedom by the node that owns each one,
 * in node order, then the interior degrees of freedom. A shape function index thus
 * identifies an entity and a position in that entity.
 *
 * Two elements that share an edge or a face must put each of its degrees of freedom at the
 * same interpolation point. Each element numbers the corner vertices of that edge or face in
 * its own local order. Elem::positive_edge_orientation and Elem::oriented_face_vertices give
 * the same result on both elements, because they compare the positions of the shared
 * vertices. An edge starts at its lesser vertex. A face starts at its least vertex and goes
 * along the lesser of the two edges at it.
 */
template <>
inline GLLTensorIndex fe_lagrange_gll_tensor_index<LAGRANGE_GLL>(const Elem * elem,
                                                                 const Order order,
                                                                 const unsigned int dim,
                                                                 const unsigned int i)
{
  libmesh_assert(elem);
  libmesh_assert_equal_to (dim, elem->dim());

  if (dim == 1)
    return {fe_lagrange_gll_1D_index(order, i), 0, 0};

  const unsigned int degree = static_cast<unsigned int>(order);
  const unsigned int n_v = elem->n_vertices();

  if (i < n_v)
    return fe_lagrange_gll_vertex_index(*elem, degree, i);

  // At first order, the vertices hold all the degrees of freedom
  libmesh_assert_greater (degree, 1);

  const unsigned int per_edge = degree - 1;
  const unsigned int n_e = elem->n_edges();

  if (i < n_v + n_e * per_edge)
    {
      const unsigned int e = (i - n_v) / per_edge;
      const unsigned int k = (i - n_v) % per_edge;

      // The node carrying this edge's degrees of freedom follows the vertices in node order
      libmesh_assert (elem->is_node_on_edge(n_v + e, e));

      // The edge runs from its lesser vertex to its greater one. positive_edge_orientation is
      // true when the zeroth vertex of the edge is the greater one, in which case the edge
      // starts at the first vertex.
      const unsigned int start = elem->positive_edge_orientation(e) ? 1 : 0;

      const GLLTensorIndex from =
        fe_lagrange_gll_vertex_index(*elem, degree, elem->local_edge_node(e, start));
      const GLLTensorIndex to =
        fe_lagrange_gll_vertex_index(*elem, degree, elem->local_edge_node(e, 1 - start));

      return fe_lagrange_gll_advance_index(from, from, to, k + 1);
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

      // The vertices of the face, from its least vertex toward the lesser neighbor of that
      // vertex. The last entry is the other neighbor.
      const std::vector<unsigned int> face_vertices = elem->oriented_face_vertices(f);

      const GLLTensorIndex from =
        fe_lagrange_gll_vertex_index(*elem, degree, face_vertices.front());
      const GLLTensorIndex first =
        fe_lagrange_gll_vertex_index(*elem, degree, face_vertices[1]);
      const GLLTensorIndex second =
        fe_lagrange_gll_vertex_index(*elem, degree, face_vertices.back());

      // The two edges at the start vertex change different coordinates, so the second step
      // keeps the coordinate that the first step set
      const GLLTensorIndex along_first =
        fe_lagrange_gll_advance_index(from, from, first, m % per_edge + 1);
      return fe_lagrange_gll_advance_index(along_first, from, second, m / per_edge + 1);
    }

  // The remaining degrees of freedom belong to the interior. They cover the interior of the
  // tensor grid in lexicographic order, with the first coordinate varying fastest.
  const unsigned int m = i - (after_edges + n_f * per_face);

  GLLTensorIndex indices = {};
  unsigned int stride = 1;

  for (const auto d : make_range(dim))
    {
      indices[d] = (m / stride) % per_edge + 1;
      stride *= per_edge;
    }

  libmesh_assert_less (m, stride);

  return indices;
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
  const GLLTensorIndex indices = fe_lagrange_gll_tensor_index<T>(elem, order, Dim, i);

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
