// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/enum_to_string.h"

namespace libMesh
{

// The 2d Raviart-Thomas shape functions are just the 2d Nedelec (first kind)
// shape functions after a 90-degree counterclockwise rotation.
template <>
RealGradient FE<2,RAVIART_THOMAS>::shape(const Elem * elem,
                                         const Order order,
                                         const unsigned int i,
                                         const Point & p,
                                         const bool add_p_level)
{
  RealGradient ND1 = FE<2,NEDELEC_ONE>::shape(elem, order, i, p, add_p_level);
  return RealGradient(-ND1(1), ND1(0));
}

template <>
RealGradient FE<2,L2_RAVIART_THOMAS>::shape(const Elem * elem,
                                            const Order order,
                                            const unsigned int i,
                                            const Point & p,
                                            const bool add_p_level)
{
  return FE<2,RAVIART_THOMAS>::shape(elem, order, i, p, add_p_level);
}

template <>
RealGradient FE<2,RAVIART_THOMAS>::shape(const ElemType,
                                         const Order,
                                         const unsigned int,
                                         const Point &)
{
  libmesh_error_msg("Raviart-Thomas elements require the element type \nbecause edge orientation is needed.");
  return RealGradient();
}

template <>
RealGradient FE<2,L2_RAVIART_THOMAS>::shape(const ElemType,
                                            const Order,
                                            const unsigned int,
                                            const Point &)
{
  libmesh_error_msg("Raviart-Thomas elements require the element type \nbecause edge orientation is needed.");
  return RealGradient();
}

template <>
RealGradient FE<2,RAVIART_THOMAS>::shape(const FEType fet,
                                         const Elem * elem,
                                         const unsigned int i,
                                         const Point & p,
                                         const bool add_p_level)
{
  return FE<2,RAVIART_THOMAS>::shape(elem, fet.order, i, p, add_p_level);
}

template <>
RealGradient FE<2,L2_RAVIART_THOMAS>::shape(const FEType fet,
                                            const Elem * elem,
                                            const unsigned int i,
                                            const Point & p,
                                            const bool add_p_level)
{
  return FE<2,L2_RAVIART_THOMAS>::shape(elem, fet.order, i, p, add_p_level);
}

template <>
RealGradient FE<2,RAVIART_THOMAS>::shape_deriv(const Elem * elem,
                                               const Order order,
                                               const unsigned int i,
                                               const unsigned int j,
                                               const Point & p,
                                               const bool add_p_level)
{
  RealGradient ND1 = FE<2,NEDELEC_ONE>::shape_deriv(elem, order, i, j, p, add_p_level);
  return RealGradient(-ND1(1), ND1(0));
}

template <>
RealGradient FE<2,L2_RAVIART_THOMAS>::shape_deriv(const Elem * elem,
                                                  const Order order,
                                                  const unsigned int i,
                                                  const unsigned int j,
                                                  const Point & p,
                                                  const bool add_p_level)
{
  return FE<2,RAVIART_THOMAS>::shape_deriv(elem, order, i, j, p, add_p_level);
}

template <>
RealGradient FE<2,RAVIART_THOMAS>::shape_deriv(const ElemType,
                                               const Order,
                                               const unsigned int,
                                               const unsigned int,
                                               const Point &)
{
  libmesh_error_msg("Raviart-Thomas elements require the element type \nbecause edge orientation is needed.");
  return RealGradient();
}

template <>
RealGradient FE<2,L2_RAVIART_THOMAS>::shape_deriv(const ElemType,
                                                  const Order,
                                                  const unsigned int,
                                                  const unsigned int,
                                                  const Point &)
{
  libmesh_error_msg("Raviart-Thomas elements require the element type \nbecause edge orientation is needed.");
  return RealGradient();
}

template <>
RealGradient FE<2,RAVIART_THOMAS>::shape_deriv(const FEType fet,
                                               const Elem * elem,
                                               const unsigned int i,
                                               const unsigned int j,
                                               const Point & p,
                                               const bool add_p_level)
{
  return FE<2,RAVIART_THOMAS>::shape_deriv(elem, fet.order, i, j, p, add_p_level);
}

template <>
RealGradient FE<2,L2_RAVIART_THOMAS>::shape_deriv(const FEType fet,
                                                  const Elem * elem,
                                                  const unsigned int i,
                                                  const unsigned int j,
                                                  const Point & p,
                                                  const bool add_p_level)
{
  return FE<2,L2_RAVIART_THOMAS>::shape_deriv(elem, fet.order, i, j, p, add_p_level);
}

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

template <>
RealGradient FE<2,RAVIART_THOMAS>::shape_second_deriv(const Elem * elem,
                                                      const Order order,
                                                      const unsigned int i,
                                                      const unsigned int j,
                                                      const Point & p,
                                                      const bool add_p_level)
{
  RealGradient ND1 = FE<2,NEDELEC_ONE>::shape_second_deriv(elem, order, i, j, p, add_p_level);
  return RealGradient(-ND1(1), ND1(0));
}

template <>
RealGradient FE<2,L2_RAVIART_THOMAS>::shape_second_deriv(const Elem * elem,
                                                         const Order order,
                                                         const unsigned int i,
                                                         const unsigned int j,
                                                         const Point & p,
                                                         const bool add_p_level)
{
  return FE<2,RAVIART_THOMAS>::shape_second_deriv(elem, order, i, j, p, add_p_level);
}

template <>
RealGradient FE<2,RAVIART_THOMAS>::shape_second_deriv(const ElemType,
                                                      const Order,
                                                      const unsigned int,
                                                      const unsigned int,
                                                      const Point &)
{
  libmesh_error_msg("Raviart-Thomas elements require the element type \nbecause edge orientation is needed.");
  return RealGradient();
}

template <>
RealGradient FE<2,L2_RAVIART_THOMAS>::shape_second_deriv(const ElemType,
                                                         const Order,
                                                         const unsigned int,
                                                         const unsigned int,
                                                         const Point &)
{
  libmesh_error_msg("Raviart-Thomas elements require the element type \nbecause edge orientation is needed.");
  return RealGradient();
}

template <>
RealGradient FE<2,RAVIART_THOMAS>::shape_second_deriv(const FEType fet,
                                                      const Elem * elem,
                                                      const unsigned int i,
                                                      const unsigned int j,
                                                      const Point & p,
                                                      const bool add_p_level)
{
  return FE<2,RAVIART_THOMAS>::shape_second_deriv(elem, fet.order, i, j, p, add_p_level);
}

template <>
RealGradient FE<2,L2_RAVIART_THOMAS>::shape_second_deriv(const FEType fet,
                                                         const Elem * elem,
                                                         const unsigned int i,
                                                         const unsigned int j,
                                                         const Point & p,
                                                         const bool add_p_level)
{
  return FE<2,L2_RAVIART_THOMAS>::shape_second_deriv(elem, fet.order, i, j, p, add_p_level);
}

#endif

} // namespace libMesh
