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

#ifndef LIBMESH_FE_MONOMIAL_SHAPE_0D_IMPL_H
#define LIBMESH_FE_MONOMIAL_SHAPE_0D_IMPL_H

// C++ includes

// Local includes
#include "libmesh/fe.h"
#include "libmesh/elem.h"

namespace libMesh
{




template <typename RealType>
typename FEShim<0,MONOMIAL,RealType>::OutputShape FEShim<0,MONOMIAL,RealType>::shape(const ElemType,
                           const Order,
                           const unsigned int libmesh_dbg_var(i),
                           const Point &)
{
  libmesh_assert_less (i, 1);
  return 1.;
}



template <typename RealType>
typename FEShim<0,MONOMIAL,RealType>::OutputShape FEShim<0,MONOMIAL,RealType>::shape(const Elem *,
                           const Order,
                           const unsigned int libmesh_dbg_var(i),
                           const Point &,
                           const bool)
{
  libmesh_assert_less (i, 1);
  return 1.;
}



template <typename RealType>
typename FEShim<0,MONOMIAL,RealType>::OutputShape FEShim<0,MONOMIAL,RealType>::shape_deriv(const ElemType,
                                 const Order,
                                 const unsigned int,
                                 const unsigned int,
                                 const Point &)
{
  libmesh_error_msg("No spatial derivatives in 0D!");
  return 0.;
}



template <typename RealType>
typename FEShim<0,MONOMIAL,RealType>::OutputShape FEShim<0,MONOMIAL,RealType>::shape_deriv(const Elem *,
                                 const Order,
                                 const unsigned int,
                                 const unsigned int,
                                 const Point &,
                                 const bool)
{
  libmesh_error_msg("No spatial derivatives in 0D!");
  return 0.;
}



#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

template <typename RealType>
typename FEShim<0,MONOMIAL,RealType>::OutputShape FEShim<0,MONOMIAL,RealType>::shape_second_deriv(const ElemType,
                                        const Order,
                                        const unsigned int,
                                        const unsigned int,
                                        const Point &)
{
  libmesh_error_msg("No spatial derivatives in 0D!");
  return 0.;
}



template <typename RealType>
typename FEShim<0,MONOMIAL,RealType>::OutputShape FEShim<0,MONOMIAL,RealType>::shape_second_deriv(const Elem *,
                                        const Order,
                                        const unsigned int,
                                        const unsigned int,
                                        const Point &,
                                        const bool)
{
  libmesh_error_msg("No spatial derivatives in 0D!");
  return 0.;
}

#endif

} // namespace libMesh

#endif // LIBMESH_FE_MONOMIAL_SHAPE_0D_IMPL_H
