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


// Local includes
#include "libmesh/fe.h"
#include "libmesh/elem.h"
#include "libmesh/fe_lagrange_gll_shape_1D.h"
#include "libmesh/int_range.h"


namespace libMesh
{

namespace {

/**
 * \returns The 2 dimensional tensor product of one-dimensional Gauss-Lobatto nodal
 * shape functions, differentiated once with respect to each coordinate listed in
 * \p deriv_coords.
 *
 * The degrees of freedom of the L2_LAGRANGE_GLL basis all belong to the element, so the
 * \p i'th runs over the points of the tensor grid in lexicographic order, the first
 * coordinate varying fastest.
 */
Real tensor_product(const Order order,
                    const unsigned int i,
                    const Point & p,
                    const unsigned int * const deriv_coords,
                    const unsigned int n_deriv_coords)
{
  const unsigned int n = static_cast<unsigned int>(order) + 1;

  libmesh_assert_less (i, n*n);

  Real value = 1.;
  unsigned int stride = 1;

  for (const auto d : make_range(2u))
    {
      const unsigned int i_d = (i / stride) % n;
      stride *= n;

      unsigned int deriv_order = 0;
      for (const auto k : make_range(n_deriv_coords))
        if (deriv_coords[k] == d)
          ++deriv_order;

      switch (deriv_order)
        {
        case 0:
          value *= fe_lagrange_gll_1D_shape(order, i_d, p(d));
          break;

        case 1:
          value *= fe_lagrange_gll_1D_shape_deriv(order, i_d, 0, p(d));
          break;

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
        case 2:
          value *= fe_lagrange_gll_1D_shape_second_deriv(order, i_d, 0, p(d));
          break;
#endif

        default:
          libmesh_error_msg("Invalid derivative order " << deriv_order << ".");
        }
    }

  return value;
}

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
/// The coordinate pairs that the second derivative index runs over
const unsigned int second_deriv_pairs[3][2] = {{0, 0}, {0, 1}, {1, 1}};
#endif // LIBMESH_ENABLE_SECOND_DERIVATIVES

} // anonymous namespace


LIBMESH_DEFAULT_VECTORIZED_FE(2,L2_LAGRANGE_GLL)


template <>
Real FE<2,L2_LAGRANGE_GLL>::shape(const ElemType,
                                  const Order order,
                                  const unsigned int i,
                                  const Point & p)
{
  return tensor_product(order, i, p, nullptr, 0);
}


template <>
Real FE<2,L2_LAGRANGE_GLL>::shape(const Elem * elem,
                                  const Order order,
                                  const unsigned int i,
                                  const Point & p,
                                  const bool add_p_level)
{
  libmesh_assert(elem);

  return tensor_product(order + add_p_level*elem->p_level(), i, p, nullptr, 0);
}


template <>
Real FE<2,L2_LAGRANGE_GLL>::shape(const FEType fet,
                                  const Elem * elem,
                                  const unsigned int i,
                                  const Point & p,
                                  const bool add_p_level)
{
  libmesh_assert(elem);

  return tensor_product(fet.order + add_p_level*elem->p_level(), i, p, nullptr, 0);
}


template <>
Real FE<2,L2_LAGRANGE_GLL>::shape_deriv(const ElemType,
                                        const Order order,
                                        const unsigned int i,
                                        const unsigned int j,
                                        const Point & p)
{
  libmesh_assert_less (j, 2);

  return tensor_product(order, i, p, &j, 1);
}


template <>
Real FE<2,L2_LAGRANGE_GLL>::shape_deriv(const Elem * elem,
                                        const Order order,
                                        const unsigned int i,
                                        const unsigned int j,
                                        const Point & p,
                                        const bool add_p_level)
{
  libmesh_assert(elem);
  libmesh_assert_less (j, 2);

  return tensor_product(order + add_p_level*elem->p_level(), i, p, &j, 1);
}


template <>
Real FE<2,L2_LAGRANGE_GLL>::shape_deriv(const FEType fet,
                                        const Elem * elem,
                                        const unsigned int i,
                                        const unsigned int j,
                                        const Point & p,
                                        const bool add_p_level)
{
  libmesh_assert(elem);
  libmesh_assert_less (j, 2);

  return tensor_product(fet.order + add_p_level*elem->p_level(), i, p, &j, 1);
}


#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

template <>
Real FE<2,L2_LAGRANGE_GLL>::shape_second_deriv(const ElemType,
                                               const Order order,
                                               const unsigned int i,
                                               const unsigned int j,
                                               const Point & p)
{
  libmesh_assert_less (j, 3);

  return tensor_product(order, i, p, second_deriv_pairs[j], 2);
}


template <>
Real FE<2,L2_LAGRANGE_GLL>::shape_second_deriv(const Elem * elem,
                                               const Order order,
                                               const unsigned int i,
                                               const unsigned int j,
                                               const Point & p,
                                               const bool add_p_level)
{
  libmesh_assert(elem);
  libmesh_assert_less (j, 3);

  return tensor_product(order + add_p_level*elem->p_level(), i, p, second_deriv_pairs[j], 2);
}


template <>
Real FE<2,L2_LAGRANGE_GLL>::shape_second_deriv(const FEType fet,
                                               const Elem * elem,
                                               const unsigned int i,
                                               const unsigned int j,
                                               const Point & p,
                                               const bool add_p_level)
{
  libmesh_assert(elem);
  libmesh_assert_less (j, 3);

  return tensor_product(fet.order + add_p_level*elem->p_level(), i, p, second_deriv_pairs[j], 2);
}


#endif // LIBMESH_ENABLE_SECOND_DERIVATIVES

} // namespace libMesh
