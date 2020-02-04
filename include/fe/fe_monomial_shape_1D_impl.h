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

#ifndef LIBMESH_FE_MONOMIAL_SHAPE_1D_IMPL_H
#define LIBMESH_FE_MONOMIAL_SHAPE_1D_IMPL_H

// C++ includes

// Local includes
#include "libmesh/fe.h"
#include "libmesh/elem.h"

namespace libMesh
{




template <typename RealType>
typename FEShim<1,MONOMIAL,RealType>::OutputShape FEShim<1,MONOMIAL,RealType>::shape(const ElemType,
                           const Order libmesh_dbg_var(order),
                           const unsigned int i,
                           const Point & p)
{
  const Real xi = p(0);

  libmesh_assert_less_equal (i, static_cast<unsigned int>(order));

  // monomials. since they are hierarchic we only need one case block.
  switch (i)
    {
    case 0:
      return 1.;

    case 1:
      return xi;

    case 2:
      return xi*xi;

    case 3:
      return xi*xi*xi;

    case 4:
      return xi*xi*xi*xi;

    default:
      Real val = 1.;
      for (unsigned int index = 0; index != i; ++index)
        val *= xi;
      return val;
    }
}



template <typename RealType>
typename FEShim<1,MONOMIAL,RealType>::OutputShape FEShim<1,MONOMIAL,RealType>::shape(const Elem * elem,
                           const Order order,
                           const unsigned int i,
                           const Point & p,
                           const bool add_p_level)
{
  libmesh_assert(elem);

  return FEShim<1,MONOMIAL,RealType>::shape(elem->type(), static_cast<Order>(order + add_p_level * elem->p_level()), i, p);
}



template <typename RealType>
typename FEShim<1,MONOMIAL,RealType>::OutputShape FEShim<1,MONOMIAL,RealType>::shape_deriv(const ElemType,
                                 const Order libmesh_dbg_var(order),
                                 const unsigned int i,
                                 const unsigned int libmesh_dbg_var(j),
                                 const Point & p)
{
  // only d()/dxi in 1D!

  libmesh_assert_equal_to (j, 0);

  const Real xi = p(0);

  libmesh_assert_less_equal (i, static_cast<unsigned int>(order));

  // monomials. since they are hierarchic we only need one case block.
  switch (i)
    {
    case 0:
      return 0.;

    case 1:
      return 1.;

    case 2:
      return 2.*xi;

    case 3:
      return 3.*xi*xi;

    case 4:
      return 4.*xi*xi*xi;

    default:
      Real val = i;
      for (unsigned int index = 1; index != i; ++index)
        val *= xi;
      return val;
    }
}



template <typename RealType>
typename FEShim<1,MONOMIAL,RealType>::OutputShape FEShim<1,MONOMIAL,RealType>::shape_deriv(const Elem * elem,
                                 const Order order,
                                 const unsigned int i,
                                 const unsigned int j,
                                 const Point & p,
                                 const bool add_p_level)
{
  libmesh_assert(elem);

  return FEShim<1,MONOMIAL,RealType>::shape_deriv(elem->type(),
                                     static_cast<Order>(order + add_p_level * elem->p_level()), i, j, p);
}


#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

template <typename RealType>
typename FEShim<1,MONOMIAL,RealType>::OutputShape FEShim<1,MONOMIAL,RealType>::shape_second_deriv(const ElemType,
                                        const Order libmesh_dbg_var(order),
                                        const unsigned int i,
                                        const unsigned int libmesh_dbg_var(j),
                                        const Point & p)
{
  // only d()/dxi in 1D!

  libmesh_assert_equal_to (j, 0);

  const Real xi = p(0);

  libmesh_assert_less_equal (i, static_cast<unsigned int>(order));

  switch (i)
    {
    case 0:
    case 1:
      return 0.;

    case 2:
      return 2.;

    case 3:
      return 6.*xi;

    case 4:
      return 12.*xi*xi;

    default:
      Real val = 2.;
      for (unsigned int index = 2; index != i; ++index)
        val *= (index+1) * xi;
      return val;
    }
}



template <typename RealType>
typename FEShim<1,MONOMIAL,RealType>::OutputShape FEShim<1,MONOMIAL,RealType>::shape_second_deriv(const Elem * elem,
                                        const Order order,
                                        const unsigned int i,
                                        const unsigned int j,
                                        const Point & p,
                                        const bool add_p_level)
{
  libmesh_assert(elem);

  return FEShim<1,MONOMIAL,RealType>::shape_second_deriv(elem->type(),
                                            static_cast<Order>(order + add_p_level * elem->p_level()), i, j, p);
}

#endif

} // namespace libMesh

#endif // LIBMESH_FE_MONOMIAL_SHAPE_1D_IMPL_H
