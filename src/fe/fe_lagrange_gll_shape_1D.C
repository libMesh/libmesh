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
#include "libmesh/fe_lagrange_gll_shape.h"


namespace libMesh
{


LIBMESH_DEFAULT_VECTORIZED_FE(1,L2_LAGRANGE_GLL)


template <>
Real FE<1,L2_LAGRANGE_GLL>::shape(const ElemType,
                                  const Order order,
                                  const unsigned int i,
                                  const Point & p)
{
  return fe_lagrange_gll_shape<L2_LAGRANGE_GLL, 1>
    (nullptr, order, i, p, nullptr, 0);
}


template <>
Real FE<1,L2_LAGRANGE_GLL>::shape(const Elem * elem,
                                  const Order order,
                                  const unsigned int i,
                                  const Point & p,
                                  const bool add_p_level)
{
  libmesh_assert(elem);

  return fe_lagrange_gll_shape<L2_LAGRANGE_GLL, 1>
    (elem, order + add_p_level*elem->p_level(), i, p, nullptr, 0);
}


template <>
Real FE<1,L2_LAGRANGE_GLL>::shape(const FEType fet,
                                  const Elem * elem,
                                  const unsigned int i,
                                  const Point & p,
                                  const bool add_p_level)
{
  libmesh_assert(elem);

  return fe_lagrange_gll_shape<L2_LAGRANGE_GLL, 1>
    (elem, fet.order + add_p_level*elem->p_level(), i, p, nullptr, 0);
}


template <>
Real FE<1,L2_LAGRANGE_GLL>::shape_deriv(const ElemType,
                                        const Order order,
                                        const unsigned int i,
                                        const unsigned int j,
                                        const Point & p)
{
  libmesh_assert_less (j, 1);

  return fe_lagrange_gll_shape<L2_LAGRANGE_GLL, 1>
    (nullptr, order, i, p, &j, 1);
}


template <>
Real FE<1,L2_LAGRANGE_GLL>::shape_deriv(const Elem * elem,
                                        const Order order,
                                        const unsigned int i,
                                        const unsigned int j,
                                        const Point & p,
                                        const bool add_p_level)
{
  libmesh_assert(elem);
  libmesh_assert_less (j, 1);

  return fe_lagrange_gll_shape<L2_LAGRANGE_GLL, 1>
    (elem, order + add_p_level*elem->p_level(), i, p, &j, 1);
}


template <>
Real FE<1,L2_LAGRANGE_GLL>::shape_deriv(const FEType fet,
                                        const Elem * elem,
                                        const unsigned int i,
                                        const unsigned int j,
                                        const Point & p,
                                        const bool add_p_level)
{
  libmesh_assert(elem);
  libmesh_assert_less (j, 1);

  return fe_lagrange_gll_shape<L2_LAGRANGE_GLL, 1>
    (elem, fet.order + add_p_level*elem->p_level(), i, p, &j, 1);
}


#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

template <>
Real FE<1,L2_LAGRANGE_GLL>::shape_second_deriv(const ElemType,
                                               const Order order,
                                               const unsigned int i,
                                               const unsigned int j,
                                               const Point & p)
{
  libmesh_assert_less (j, fe_lagrange_gll_n_second_derivs<1>());

  return fe_lagrange_gll_shape<L2_LAGRANGE_GLL, 1>
    (nullptr, order, i, p, fe_lagrange_gll_second_deriv_pairs[j], 2);
}


template <>
Real FE<1,L2_LAGRANGE_GLL>::shape_second_deriv(const Elem * elem,
                                               const Order order,
                                               const unsigned int i,
                                               const unsigned int j,
                                               const Point & p,
                                               const bool add_p_level)
{
  libmesh_assert(elem);
  libmesh_assert_less (j, fe_lagrange_gll_n_second_derivs<1>());

  return fe_lagrange_gll_shape<L2_LAGRANGE_GLL, 1>
    (elem, order + add_p_level*elem->p_level(), i, p, fe_lagrange_gll_second_deriv_pairs[j], 2);
}


template <>
Real FE<1,L2_LAGRANGE_GLL>::shape_second_deriv(const FEType fet,
                                               const Elem * elem,
                                               const unsigned int i,
                                               const unsigned int j,
                                               const Point & p,
                                               const bool add_p_level)
{
  libmesh_assert(elem);
  libmesh_assert_less (j, fe_lagrange_gll_n_second_derivs<1>());

  return fe_lagrange_gll_shape<L2_LAGRANGE_GLL, 1>
    (elem, fet.order + add_p_level*elem->p_level(), i, p, fe_lagrange_gll_second_deriv_pairs[j], 2);
}


#endif // LIBMESH_ENABLE_SECOND_DERIVATIVES

} // namespace libMesh
