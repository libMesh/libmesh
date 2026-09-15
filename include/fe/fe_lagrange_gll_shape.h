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
