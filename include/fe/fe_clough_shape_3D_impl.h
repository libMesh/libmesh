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

#ifndef LIBMESH_FE_CLOUGH_IMPL_SHAPE_3D_H
#define LIBMESH_FE_CLOUGH_IMPL_SHAPE_3D_H

// C++ includes

// Local includes
#include "libmesh/fe.h"
#include "libmesh/elem.h"


// FIXME: 3D C1 finite elements are still a work in progress

namespace libMesh
{


template <typename RealType>
typename FEShim<3,CLOUGH,RealType>::OutputShape FEShim<3,CLOUGH,RealType>::shape(const ElemType,
                         const Order,
                         const unsigned int,
                         const Point &)
{
  libmesh_error_msg("Clough-Tocher elements require the real element \nto construct gradient-based degrees of freedom.");
  return 0.;
}



template <typename RealType>
typename FEShim<3,CLOUGH,RealType>::OutputShape FEShim<3,CLOUGH,RealType>::shape(const Elem * libmesh_dbg_var(elem),
                         const Order,
                         const unsigned int,
                         const Point &,
                         const bool)
{
  libmesh_assert(elem);

  libmesh_not_implemented();
  return 0.;
}



template <typename RealType>
typename FEShim<3,CLOUGH,RealType>::OutputShape FEShim<3,CLOUGH,RealType>::shape_deriv(const ElemType,
                               const Order,
                               const unsigned int,
                               const unsigned int,
                               const Point &)
{
  libmesh_error_msg("Clough-Tocher elements require the real element \nto construct gradient-based degrees of freedom.");
  return 0.;
}



template <typename RealType>
typename FEShim<3,CLOUGH,RealType>::OutputShape FEShim<3,CLOUGH,RealType>::shape_deriv(const Elem * libmesh_dbg_var(elem),
                               const Order,
                               const unsigned int,
                               const unsigned int,
                               const Point &,
                               const bool)
{
  libmesh_assert(elem);
  libmesh_not_implemented();
  return 0.;
}



#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

template <typename RealType>
typename FEShim<3,CLOUGH,RealType>::OutputShape FEShim<3,CLOUGH,RealType>::shape_second_deriv(const ElemType,
                                      const Order,
                                      const unsigned int,
                                      const unsigned int,
                                      const Point &)
{
  libmesh_error_msg("Clough-Tocher elements require the real element \nto construct gradient-based degrees of freedom.");
  return 0.;
}


template <typename RealType>
typename FEShim<3,CLOUGH,RealType>::OutputShape FEShim<3,CLOUGH,RealType>::shape_second_deriv(const Elem * libmesh_dbg_var(elem),
                                      const Order,
                                      const unsigned int,
                                      const unsigned int,
                                      const Point &,
                                      const bool)
{
  libmesh_assert(elem);
  libmesh_not_implemented();
  return 0.;
}

#endif

} // namespace libMesh

#endif // LIBMESH_FE_CLOUGH_IMPL_SHAPE_3D_H
