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

#ifndef LIBMESH_FE_LAGRANGE_SHAPE_1D_IMPL_H
#define LIBMESH_FE_LAGRANGE_SHAPE_1D_IMPL_H

// Local includes
#include "libmesh/fe.h"
#include "libmesh/elem.h"
#include "libmesh/fe_lagrange_shape_1D.h"

namespace libMesh
{

template <typename RealType>
typename FEShim<1,LAGRANGE,RealType>::OutputShape FEShim<1,LAGRANGE,RealType>::shape(const ElemType,
                           const Order order,
                           const unsigned int i,
                           const Point & p)
{
  return fe_lagrange_1D_shape(order, i, p(0));
}

template <typename RealType>
typename FEShim<1,L2_LAGRANGE,RealType>::OutputShape FEShim<1,L2_LAGRANGE,RealType>::shape(const ElemType,
                              const Order order,
                              const unsigned int i,
                              const Point & p)
{
  return fe_lagrange_1D_shape(order, i, p(0));
}


template <typename RealType>
typename FEShim<1,LAGRANGE,RealType>::OutputShape FEShim<1,LAGRANGE,RealType>::shape(const Elem * elem,
                           const Order order,
                           const unsigned int i,
                           const Point & p,
                           const bool add_p_level)
{
  libmesh_assert(elem);

  return fe_lagrange_1D_shape(static_cast<Order>(order + add_p_level * elem->p_level()), i, p(0));
}



template <typename RealType>
typename FEShim<1,L2_LAGRANGE,RealType>::OutputShape FEShim<1,L2_LAGRANGE,RealType>::shape(const Elem * elem,
                              const Order order,
                              const unsigned int i,
                              const Point & p,
                              const bool add_p_level)
{
  libmesh_assert(elem);

  return fe_lagrange_1D_shape(static_cast<Order>(order + add_p_level * elem->p_level()), i, p(0));
}



template <typename RealType>
typename FEShim<1,LAGRANGE,RealType>::OutputShape FEShim<1,LAGRANGE,RealType>::shape_deriv(const ElemType,
                                 const Order order,
                                 const unsigned int i,
                                 const unsigned int j,
                                 const Point & p)
{
  return fe_lagrange_1D_shape_deriv(order, i, j, p(0));
}



template <typename RealType>
typename FEShim<1,L2_LAGRANGE,RealType>::OutputShape FEShim<1,L2_LAGRANGE,RealType>::shape_deriv(const ElemType,
                                    const Order order,
                                    const unsigned int i,
                                    const unsigned int j,
                                    const Point & p)
{
  return fe_lagrange_1D_shape_deriv(order, i, j, p(0));
}



template <typename RealType>
typename FEShim<1,LAGRANGE,RealType>::OutputShape FEShim<1,LAGRANGE,RealType>::shape_deriv(const Elem * elem,
                                 const Order order,
                                 const unsigned int i,
                                 const unsigned int j,
                                 const Point & p,
                                 const bool add_p_level)
{
  libmesh_assert(elem);

  return fe_lagrange_1D_shape_deriv(static_cast<Order>(order + add_p_level * elem->p_level()), i, j, p(0));
}



template <typename RealType>
typename FEShim<1,L2_LAGRANGE,RealType>::OutputShape FEShim<1,L2_LAGRANGE,RealType>::shape_deriv(const Elem * elem,
                                    const Order order,
                                    const unsigned int i,
                                    const unsigned int j,
                                    const Point & p,
                                    const bool add_p_level)
{
  libmesh_assert(elem);

  return fe_lagrange_1D_shape_deriv(static_cast<Order>(order + add_p_level * elem->p_level()), i, j, p(0));
}



#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

template <typename RealType>
typename FEShim<1,LAGRANGE,RealType>::OutputShape FEShim<1,LAGRANGE,RealType>::shape_second_deriv(const ElemType,
                                        const Order order,
                                        const unsigned int i,
                                        const unsigned int j,
                                        const Point & p)
{
  return fe_lagrange_1D_shape_second_deriv(order, i, j, p(0));
}



template <typename RealType>
typename FEShim<1,L2_LAGRANGE,RealType>::OutputShape FEShim<1,L2_LAGRANGE,RealType>::shape_second_deriv(const ElemType,
                                           const Order order,
                                           const unsigned int i,
                                           const unsigned int j,
                                           const Point & p)
{
  return fe_lagrange_1D_shape_second_deriv(order, i, j, p(0));
}



template <typename RealType>
typename FEShim<1,LAGRANGE,RealType>::OutputShape FEShim<1,LAGRANGE,RealType>::shape_second_deriv(const Elem * elem,
                                        const Order order,
                                        const unsigned int i,
                                        const unsigned int j,
                                        const Point & p,
                                        const bool add_p_level)
{
  libmesh_assert(elem);

  return fe_lagrange_1D_shape_second_deriv(static_cast<Order>(order + add_p_level * elem->p_level()), i, j, p(0));
}



template <typename RealType>
typename FEShim<1,L2_LAGRANGE,RealType>::OutputShape FEShim<1,L2_LAGRANGE,RealType>::shape_second_deriv(const Elem * elem,
                                           const Order order,
                                           const unsigned int i,
                                           const unsigned int j,
                                           const Point & p,
                                           const bool add_p_level)
{
  libmesh_assert(elem);

  return fe_lagrange_1D_shape_second_deriv(static_cast<Order>(order + add_p_level * elem->p_level()), i, j, p(0));
}

#endif // LIBMESH_ENABLE_SECOND_DERIVATIVES

} // namespace libMesh

#endif // LIBMESH_FE_LAGRANGE_SHAPE_1D_IMPL_H
