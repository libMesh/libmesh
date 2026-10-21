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
#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_order.h"
#include "libmesh/fe_lagrange_gll_shape_1D.h"
#include "libmesh/int_range.h"

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
