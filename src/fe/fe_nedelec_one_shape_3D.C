// The libMesh Finite Element Library.
// Copyright (C) 2002-2024 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

template <>
RealGradient FE<3,NEDELEC_ONE>::shape(const Elem * elem,
                                      const Order order,
                                      const unsigned int i,
                                      const Point & p,
                                      const bool add_p_level)
{
#if LIBMESH_DIM == 3
  libmesh_assert(elem);

  const Order totalorder = order + add_p_level*elem->p_level();
  libmesh_assert_less(i, n_dofs(elem->type(), totalorder));

  const char sign = elem->edge_orientation(i) ? 1 : -1;

  const Real xi   = p(0);
  const Real eta  = p(1);
  const Real zeta = p(2);

  switch (totalorder)
    {
      // linear Nedelec (first kind) shape functions
    case FIRST:
      {
        switch (elem->type())
          {
          case HEX20:
          case HEX27:
            {
              // Even with a loose inverse_map tolerance we ought to
              // be nearly on the element interior in master
              // coordinates
              libmesh_assert_less_equal ( std::fabs(xi),   1.0+10*TOLERANCE );
              libmesh_assert_less_equal ( std::fabs(eta),  1.0+10*TOLERANCE );
              libmesh_assert_less_equal ( std::fabs(zeta), 1.0+10*TOLERANCE );

              switch(i)
                {
                case 0:
                  return sign * RealGradient( -0.125*(1.0-eta-zeta+eta*zeta), 0.0, 0.0 );
                case 1:
                  return sign * RealGradient( 0.0, -0.125*(1.0+xi-zeta-xi*zeta), 0.0 );
                case 2:
                  return sign * RealGradient( 0.125*(1.0+eta-zeta-eta*zeta), 0.0, 0.0 );
                case 3:
                  return sign * RealGradient( 0.0, -0.125*(1.0-xi-zeta+xi*zeta), 0.0 );
                case 4:
                  return sign * RealGradient( 0.0, 0.0, -0.125*(1.0-xi-eta+xi*eta) );
                case 5:
                  return sign * RealGradient( 0.0, 0.0, -0.125*(1.0+xi-eta-xi*eta) );
                case 6:
                  return sign * RealGradient( 0.0, 0.0, -0.125*(1.0+xi+eta+xi*eta) );
                case 7:
                  return sign * RealGradient( 0.0, 0.0, -0.125*(1.0-xi+eta-xi*eta) );
                case 8:
                  return sign * RealGradient( -0.125*(1.0-eta+zeta-eta*zeta), 0.0, 0.0 );
                case 9:
                  return sign * RealGradient( 0.0, -0.125*(1.0+xi+zeta+xi*zeta), 0.0 );
                case 10:
                  return sign * RealGradient( 0.125*(1.0+eta+zeta+eta*zeta), 0.0, 0.0 );
                case 11:
                  return sign * RealGradient( 0.0, -0.125*(1.0-xi+zeta-xi*zeta), 0.0 );

                default:
                  libmesh_error_msg("Invalid i = " << i);
                }
            }

          case TET10:
          case TET14:
            {
              switch(i)
                {
                case 0:
                  return sign * RealGradient( -1.0+eta+zeta, -xi, -xi );
                case 1:
                  return sign * RealGradient( eta, -xi, 0.0 );
                case 2:
                  return sign * RealGradient( -eta, -1.0+xi+zeta, -eta );
                case 3:
                  return sign * RealGradient( -zeta, -zeta, -1.0+xi+eta );
                case 4:
                  return sign * RealGradient( zeta, 0.0, -xi );
                case 5:
                  return sign * RealGradient( 0.0, zeta, -eta );

                default:
                  libmesh_error_msg("Invalid i = " << i);
                }
            }

          default:
            libmesh_error_msg("ERROR: Unsupported 3D element type!: " << Utility::enum_to_string(elem->type()));
          }
      }

      // unsupported order
    default:
      libmesh_error_msg("ERROR: Unsupported 3D FE order!: " << totalorder);
    }

#else // LIBMESH_DIM != 3
  libmesh_ignore(elem, order, i, p, add_p_level);
  libmesh_not_implemented();
#endif
}


template <>
RealGradient FE<3,NEDELEC_ONE>::shape(const ElemType,
                                      const Order,
                                      const unsigned int,
                                      const Point &)
{
  libmesh_error_msg("Nedelec elements require the element type \nbecause edge orientation is needed.");
  return RealGradient();
}


template <>
RealGradient FE<3,NEDELEC_ONE>::shape(const FEType fet,
                                      const Elem * elem,
                                      const unsigned int i,
                                      const Point & p,
                                      const bool add_p_level)
{
  return FE<3,NEDELEC_ONE>::shape(elem, fet.order, i, p, add_p_level);
}


template <>
RealGradient FE<3,NEDELEC_ONE>::shape_deriv(const Elem * elem,
                                            const Order order,
                                            const unsigned int i,
                                            const unsigned int j,
                                            const Point & p,
                                            const bool add_p_level)
{
#if LIBMESH_DIM == 3
  libmesh_assert(elem);
  libmesh_assert_less (j, 3);

  const Order totalorder = order + add_p_level*elem->p_level();
  libmesh_assert_less(i, n_dofs(elem->type(), totalorder));

  const char sign = elem->edge_orientation(i) ? 1 : -1;

  const Real xi   = p(0);
  const Real eta  = p(1);
  const Real zeta = p(2);

  switch (totalorder)
    {
      // linear Nedelec (first kind) shape function first derivatives
    case FIRST:
      {
        switch (elem->type())
          {
          case HEX20:
          case HEX27:
            {
              // Even with a loose inverse_map tolerance we ought to
              // be nearly on the element interior in master
              // coordinates
              libmesh_assert_less_equal ( std::fabs(xi),   1.0+10*TOLERANCE );
              libmesh_assert_less_equal ( std::fabs(eta),  1.0+10*TOLERANCE );
              libmesh_assert_less_equal ( std::fabs(zeta), 1.0+10*TOLERANCE );

              switch (j)
                {
                  // d()/dxi
                case 0:
                  {
                    switch(i)
                      {
                      case 0:
                      case 2:
                      case 8:
                      case 10:
                        return RealGradient();
                      case 1:
                        return sign * RealGradient( 0.0, -0.125*(1.0-zeta) );
                      case 3:
                        return sign * RealGradient( 0.0, -0.125*(-1.0+zeta) );
                      case 4:
                        return sign * RealGradient( 0.0, 0.0, -0.125*(-1.0+eta) );
                      case 5:
                        return sign * RealGradient( 0.0, 0.0, -0.125*(1.0-eta) );
                      case 6:
                        return sign * RealGradient( 0.0, 0.0, -0.125*(1.0+eta) );
                      case 7:
                        return sign * RealGradient( 0.0, 0.0, -0.125*(-1.0-eta) );
                      case 9:
                        return sign * RealGradient( 0.0, -0.125*(1.0+zeta), 0.0 );
                      case 11:
                        return sign * RealGradient( 0.0, -0.125*(-1.0-zeta), 0.0 );

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      } // switch(i)

                  } // j = 0

                  // d()/deta
                case 1:
                  {
                    switch(i)
                      {
                      case 1:
                      case 3:
                      case 9:
                      case 11:
                        return RealGradient();
                      case 0:
                        return sign * RealGradient( -0.125*(-1.0+zeta), 0.0, 0.0 );
                      case 2:
                        return sign * RealGradient( 0.125*(1.0-zeta), 0.0, 0.0 );
                      case 4:
                        return sign * RealGradient( 0.0, 0.0, -0.125*(-1.0+xi) );
                      case 5:
                        return sign * RealGradient( 0.0, 0.0, -0.125*(-1.0-xi) );
                      case 6:
                        return sign * RealGradient( 0.0, 0.0, -0.125*(1.0+xi) );
                      case 7:
                        return sign * RealGradient( 0.0, 0.0, -0.125*(1.0-xi) );
                      case 8:
                        return sign * RealGradient( -0.125*(-1.0-zeta), 0.0, 0.0 );
                      case 10:
                        return sign * RealGradient( 0.125*(1.0+zeta), 0.0, 0.0 );

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      } // switch(i)

                  } // j = 1

                  // d()/dzeta
                case 2:
                  {
                    switch(i)
                      {
                      case 4:
                      case 5:
                      case 6:
                      case 7:
                        return RealGradient();
                      case 0:
                        return sign * RealGradient( -0.125*(-1.0+eta), 0.0, 0.0 );
                      case 1:
                        return sign * RealGradient( 0.0, -0.125*(-1.0-xi), 0.0 );
                      case 2:
                        return sign * RealGradient( 0.125*(-1.0-eta), 0.0, 0.0 );
                      case 3:
                        return sign * RealGradient( 0.0, -0.125*(-1.0+xi), 0.0 );
                      case 8:
                        return sign * RealGradient( -0.125*(1.0-eta), 0.0, 0.0 );
                      case 9:
                        return sign * RealGradient( 0.0, -0.125*(1.0+xi), 0.0 );
                      case 10:
                        return sign * RealGradient( 0.125*(1.0+eta), 0.0, 0.0 );
                      case 11:
                        return sign * RealGradient( 0.0, -0.125*(1.0-xi), 0.0 );

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      } // switch(i)

                  } // j = 2

                default:
                  libmesh_error_msg("Invalid j = " << j);
                }
            }

          case TET10:
          case TET14:
            {
              switch (j)
                {
                  // d()/dxi
                case 0:
                  {
                    switch(i)
                      {
                      case 0:
                        return sign * RealGradient( 0.0, -1.0, -1.0 );
                      case 1:
                        return sign * RealGradient( 0.0, -1.0, 0.0 );
                      case 2:
                        return sign * RealGradient( 0.0, 1.0, 0.0 );
                      case 3:
                        return sign * RealGradient( 0.0, 0.0, 1.0 );
                      case 4:
                        return sign * RealGradient( 0.0, 0.0, -1.0 );
                      case 5:
                        return RealGradient();

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      } // switch(i)

                  } // j = 0

                  // d()/deta
                case 1:
                  {
                    switch(i)
                      {
                      case 0:
                      case 1:
                        return sign * RealGradient( 1.0, 0.0, 0.0 );
                      case 2:
                        return sign * RealGradient( -1.0, 0.0, -1.0 );
                      case 3:
                        return sign * RealGradient( 0.0, 0.0, 1.0 );
                      case 4:
                        return RealGradient();
                      case 5:
                        return sign * RealGradient( 0.0, 0.0, -1.0 );

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      } // switch(i)

                  } // j = 1

                  // d()/dzeta
                case 2:
                  {
                    switch(i)
                      {
                        case 0:
                          return sign * RealGradient( 1.0, 0.0, 0.0 );
                        case 1:
                          return RealGradient();
                        case 2:
                        case 5:
                          return sign * RealGradient( 0.0, 1.0, 0.0 );
                        case 3:
                          return sign * RealGradient( -1.0, -1.0, 0.0 );
                        case 4:
                          return sign * RealGradient( 1.0, 0.0, 0.0 );

                        default:
                          libmesh_error_msg("Invalid i = " << i);
                      } // switch(i)

                  } // j = 2

                default:
                  libmesh_error_msg("Invalid j = " << j);
                }
            }

          default:
            libmesh_error_msg("ERROR: Unsupported 3D element type!: " << Utility::enum_to_string(elem->type()));
          }
      }
      // unsupported order
    default:
      libmesh_error_msg("ERROR: Unsupported 3D FE order!: " << totalorder);
    }

#else // LIBMESH_DIM != 3
  libmesh_ignore(elem, order, i, j, p, add_p_level);
  libmesh_not_implemented();
#endif
}


template <>
RealGradient FE<3,NEDELEC_ONE>::shape_deriv(const ElemType,
                                            const Order,
                                            const unsigned int,
                                            const unsigned int,
                                            const Point &)
{
  libmesh_error_msg("Nedelec elements require the element type \nbecause edge orientation is needed.");
  return RealGradient();
}


template <>
RealGradient FE<3,NEDELEC_ONE>::shape_deriv(const FEType fet,
                                            const Elem * elem,
                                            const unsigned int i,
                                            const unsigned int j,
                                            const Point & p,
                                            const bool add_p_level)
{
  return FE<3,NEDELEC_ONE>::shape_deriv(elem, fet.order, i, j, p, add_p_level);
}


#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

template <>
RealGradient FE<3,NEDELEC_ONE>::shape_second_deriv(const Elem * elem,
                                                   const Order order,
                                                   const unsigned int i,
                                                   const unsigned int j,
                                                   const Point &,
                                                   const bool add_p_level)
{
#if LIBMESH_DIM == 3

  libmesh_assert(elem);

  // j = 0 ==> d^2 phi / dxi^2
  // j = 1 ==> d^2 phi / dxi deta
  // j = 2 ==> d^2 phi / deta^2
  // j = 3 ==> d^2 phi / dxi dzeta
  // j = 4 ==> d^2 phi / deta dzeta
  // j = 5 ==> d^2 phi / dzeta^2
  libmesh_assert_less (j, 6);

  const Order totalorder = order + add_p_level*elem->p_level();
  libmesh_assert_less(i, n_dofs(elem->type(), totalorder));

  const char sign = elem->edge_orientation(i) ? 1 : -1;

  switch (totalorder)
    {
      // linear Nedelec (first kind) shape function second derivatives
    case FIRST:
      {
        switch (elem->type())
          {
          case HEX20:
          case HEX27:
            {
              switch (j)
                {
                  // d^2()/dxi^2, d^2()/deta^2, d^2()/dzeta^2
                case 0:
                case 2:
                case 5:
                  return RealGradient();

                  // d^2()/dxideta
                case 1:
                  {
                    switch(i)
                      {
                      case 0:
                      case 1:
                      case 2:
                      case 3:
                      case 8:
                      case 9:
                      case 10:
                      case 11:
                        return RealGradient();
                      case 4:
                      case 6:
                        return sign * RealGradient( 0.0, 0.0, -0.125 );
                      case 5:
                      case 7:
                        return sign * RealGradient( 0.0, 0.0, 0.125 );

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      } // switch(i)

                  } // j = 1

                  // d^2()/dxidzeta
                case 3:
                  {
                    switch(i)
                      {
                      case 0:
                      case 2:
                      case 4:
                      case 5:
                      case 6:
                      case 7:
                      case 8:
                      case 10:
                        return RealGradient();
                      case 1:
                      case 3:
                      case 11:
                        return sign * RealGradient( 0.0, 0.125 );
                      case 9:
                        return sign * RealGradient( 0.0, -0.125, 0.0 );

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      } // switch(i)

                  } // j = 3

                  // d^2()/detadzeta
                case 4:
                  {
                    switch(i)
                      {
                      case 0:
                        return sign * RealGradient( -0.125, 0.0, 0.0 );
                      case 1:
                      case 3:
                      case 4:
                      case 5:
                      case 6:
                      case 7:
                      case 9:
                      case 11:
                        return RealGradient();
                      case 2:
                      case 8:
                      case 10:
                        return sign * RealGradient( 0.125, 0.0, 0.0 );

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      } // switch(i)

                  } // j = 4

                default:
                  libmesh_error_msg("Invalid j = " << j);
                }
            }

            // All second derivatives for linear tets are zero.
          case TET10:
          case TET14:
            return RealGradient();

          default:
            libmesh_error_msg("ERROR: Unsupported 3D element type!: " << Utility::enum_to_string(elem->type()));

          } //switch(type)

      }

      // unsupported order
    default:
      libmesh_error_msg("ERROR: Unsupported 3D FE order!: " << totalorder);
    }

#else // LIBMESH_DIM != 3
  libmesh_assert(true || p(0));
  libmesh_ignore(elem, order, i, j, add_p_level);
  libmesh_not_implemented();
#endif
}


template <>
RealGradient FE<3,NEDELEC_ONE>::shape_second_deriv(const ElemType,
                                                   const Order,
                                                   const unsigned int,
                                                   const unsigned int,
                                                   const Point &)
{
  libmesh_error_msg("Nedelec elements require the element type \nbecause edge orientation is needed.");
  return RealGradient();
}


template <>
RealGradient FE<3,NEDELEC_ONE>::shape_second_deriv(const FEType fet,
                                                   const Elem * elem,
                                                   const unsigned int i,
                                                   const unsigned int j,
                                                   const Point & p,
                                                   const bool add_p_level)
{
  return FE<3,NEDELEC_ONE>::shape_second_deriv(elem, fet.order, i, j, p, add_p_level);
}

#endif

} // namespace libMesh
