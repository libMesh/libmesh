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

// An excellent discussion of Nedelec shape functions is given in
// https://www.dealii.org/reports/nedelec/nedelec.pdf
// An excellent summary of Nedelec shape functions is also given in
// https://defelement.com/elements/nedelec1.html
template <>
RealGradient FE<2,NEDELEC_ONE>::shape(const Elem * elem,
                                      const Order order,
                                      const unsigned int i,
                                      const Point & p,
                                      const bool add_p_level)
{
#if LIBMESH_DIM > 1
  libmesh_assert(elem);

  const Order total_order = static_cast<Order>(order + add_p_level * elem->p_level());
  libmesh_assert_less(i, n_dofs(elem->type(), total_order));

  const char sign = i >= total_order * elem->n_edges() || elem->point(i / total_order) > elem->point((i / total_order + 1) % elem->n_vertices()) ? 1 : -1;
  const unsigned int ii = sign > 0 ? i : (i / total_order * 2 + 1) * total_order - 1 - i;

  const Real xi  = p(0);
  const Real eta = p(1);

  switch (total_order)
    {
      // linear Nedelec (first kind) shape functions
    case FIRST:
      {
        switch (elem->type())
          {
          case QUAD8:
          case QUAD9:
            {
              // Even with a loose inverse_map tolerance we ought to
              // be nearly on the element interior in master
              // coordinates
              libmesh_assert_less_equal ( std::fabs(xi), 1.0+10*TOLERANCE );
              libmesh_assert_less_equal ( std::fabs(eta), 1.0+10*TOLERANCE );

              switch(ii)
                {
                case 0:
                  return sign * RealGradient( -0.25*(1.0-eta), 0.0 );
                case 1:
                  return sign * RealGradient( 0.0, -0.25*(1.0+xi) );
                case 2:
                  return sign * RealGradient( 0.25*(1.0+eta), 0.0 );
                case 3:
                  return sign * RealGradient( 0.0, -0.25*(xi-1.0) );

                default:
                  libmesh_error_msg("Invalid i = " << i);
                }
            }

          case TRI6:
          case TRI7:
            {
              switch(ii)
                {
                case 0:
                  return sign * RealGradient( -1.0+eta, -xi );
                case 1:
                  return sign * RealGradient( eta, -xi );
                case 2:
                  return sign * RealGradient( eta, -xi+1.0 );

                default:
                  libmesh_error_msg("Invalid i = " << i);
                }
            }

          default:
            libmesh_error_msg("ERROR: Unsupported 2D element type!: " << Utility::enum_to_string(elem->type()));
          }
      }

      // quadratic Nedelec (first kind) shape functions
    case SECOND:
      {
        switch (elem->type())
          {
          case QUAD8:
          case QUAD9:
            {
              // Even with a loose inverse_map tolerance we ought to
              // be nearly on the element interior in master
              // coordinates
              libmesh_assert_less_equal ( std::fabs(xi), 1.0+10*TOLERANCE );
              libmesh_assert_less_equal ( std::fabs(eta), 1.0+10*TOLERANCE );

              const Real x = 0.5 * (xi + 1.0);
              const Real y = 0.5 * (eta + 1.0);

              switch(ii)
                {
                case 0:
                  return sign * RealGradient( 0.5*(-18.0*x*y*y+24.0*x*y-6.0*x+12.0*y*y-16.0*y+4.0), 0.0 );
                case 1:
                  return sign * RealGradient( 0.5*( 18.0*x*y*y-24.0*x*y+6.0*x-6.0*y*y+8.0*y-2.0), 0.0 );
                case 2:
                  return sign * RealGradient( 0.0, x*(-9.0*x*y+6.0*x+6.0*y-4.0) );
                case 3:
                  return sign * RealGradient( 0.0, x*( 9.0*x*y-3.0*x-6.0*y+2.0) );
                case 4:
                  return sign * RealGradient( y*(-9.0*x*y+6.0*x+3.0*y-2.0), 0.0 );
                case 5:
                  return sign * RealGradient( y*( 9.0*x*y-6.0*x-6.0*y+4.0), 0.0 );
                case 6:
                  return sign * RealGradient( 0.0, 0.5*(-18.0*x*x*y+6.0*x*x+24.0*x*y-8.0*x-6.0*y+2.0) );
                case 7:
                  return sign * RealGradient( 0.0, 0.5*( 18.0*x*x*y-12.0*x*x-24.0*x*y+16.0*x+6.0*y-4.0) );
                case 8:
                  return RealGradient( 0.0, 3.0*x*(3*x*y-2.0*x-3.0*y+2.0) );
                case 9:
                  return RealGradient( 3.0*y*(-3.0*x*y+3.0*x+2.0*y-2.0), 0.0 );
                case 10:
                  return RealGradient( 3.0*y*(3.0*x*y-3.0*x-y+1.0), 0.0 );
                case 11:
                  return RealGradient( 0.0, 3.0*x*(-3.0*x*y+x+3.0*y-1.0) );

                default:
                  libmesh_error_msg("Invalid i = " << i);
                }
            }

          case TRI6:
          case TRI7:
            {
              switch(ii)
                {
                case 0:
                  return sign * RealGradient( 8.0*xi*eta-6.0*xi+8.0*eta*eta-12.0*eta+4.0, 2.0*xi*(-4.0*xi-4.0*eta+3.0) );
                case 1:
                  return sign * RealGradient( -8.0*xi*eta+6.0*xi+2.0*eta-2.0, 4.0*xi*(2.0*xi-1.0) );
                case 2:
                  return sign * RealGradient( 2.0*eta*(1.0-4.0*xi), 4.0*xi*(2.0*xi-1.0) );
                case 3:
                  return sign * RealGradient( 4.0*eta*(1.0-2.0*eta), 2.0*xi*(4.0*eta-1.0) );
                case 4:
                  return sign * RealGradient( 4.0*eta*(1.0-2.0*eta), 8.0*xi*eta-2.0*xi-6.0*eta+2.0 );
                case 5:
                  return sign * RealGradient( 2.0*eta*(4.0*xi+4*eta-3.0), -8.0*xi*xi-8.0*xi*eta+12.0*xi+6.0*eta-4.0 );
                case 6:
                  return RealGradient( 8.0*eta*(-xi-2.0*eta+2.0), 8.0*xi*(xi+2.0*eta-1.0) );
                case 7:
                  return RealGradient( 8.0*eta*(2.0*xi+eta-1.0), 8.0*xi*(-2.0*xi-eta+2.0) );

                default:
                  libmesh_error_msg("Invalid i = " << i);
                }
            }

          default:
            libmesh_error_msg("ERROR: Unsupported 2D element type!: " << Utility::enum_to_string(elem->type()));
          }
      }

      // unsupported order
    default:
      libmesh_error_msg("ERROR: Unsupported 2D FE order!: " << total_order);
    }
#else // LIBMESH_DIM > 1
  libmesh_ignore(elem, order, i, p, add_p_level);
  libmesh_not_implemented();
#endif
}




template <>
RealGradient FE<2,NEDELEC_ONE>::shape(const ElemType,
                                      const Order,
                                      const unsigned int,
                                      const Point &)
{
  libmesh_error_msg("Nedelec elements require the element type \nbecause edge orientation is needed.");
  return RealGradient();
}


template <>
RealGradient FE<2,NEDELEC_ONE>::shape(const FEType fet,
                                      const Elem * elem,
                                      const unsigned int i,
                                      const Point & p,
                                      const bool add_p_level)
{
  return FE<2,NEDELEC_ONE>::shape(elem, fet.order, i, p, add_p_level);
}

template <>
RealGradient FE<2,NEDELEC_ONE>::shape_deriv(const Elem * elem,
                                            const Order order,
                                            const unsigned int i,
                                            const unsigned int j,
                                            const Point & p,
                                            const bool add_p_level)
{
#if LIBMESH_DIM > 1
  libmesh_assert(elem);
  libmesh_assert_less (j, 2);

  const Order total_order = static_cast<Order>(order + add_p_level * elem->p_level());
  libmesh_assert_less(i, n_dofs(elem->type(), total_order));

  const char sign = i >= total_order * elem->n_edges() || elem->point(i / total_order) > elem->point((i / total_order + 1) % elem->n_vertices()) ? 1 : -1;
  const unsigned int ii = sign > 0 ? i : (i / total_order * 2 + 1) * total_order - 1 - i;

  const Real xi  = p(0);
  const Real eta = p(1);

  switch (total_order)
    {
      // linear Nedelec (first kind) shape function first derivatives
    case FIRST:
      {
        switch (elem->type())
          {
          case QUAD8:
          case QUAD9:
            {
              switch (j)
                {
                  // d()/dxi
                case 0:
                  {
                    switch(ii)
                      {
                      case 0:
                      case 2:
                        return RealGradient();
                      case 1:
                      case 3:
                        return sign * RealGradient( 0.0, -0.25 );

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  } // j = 0

                  // d()/deta
                case 1:
                  {
                    switch(ii)
                      {
                      case 1:
                      case 3:
                        return RealGradient();
                      case 0:
                      case 2:
                        return sign * RealGradient( 0.25 );

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  } // j = 1

                default:
                  libmesh_error_msg("Invalid j = " << j);
                }
            }

          case TRI6:
          case TRI7:
            {
              switch (j)
                {
                  // d()/dxi
                case 0:
                  return sign * RealGradient( 0.0, -1.0 );

                  // d()/deta
                case 1:
                  return sign * RealGradient( 1.0 );

                default:
                  libmesh_error_msg("Invalid j = " << j);
                }
            }

          default:
            libmesh_error_msg("ERROR: Unsupported 2D element type!: " << Utility::enum_to_string(elem->type()));
          }
      }

      // quadratic Nedelec (first kind) shape function first derivatives
    case SECOND:
      {
        switch (elem->type())
          {
          case QUAD8:
          case QUAD9:
            {
              // Even with a loose inverse_map tolerance we ought to
              // be nearly on the element interior in master
              // coordinates
              libmesh_assert_less_equal ( std::fabs(xi), 1.0+10*TOLERANCE );
              libmesh_assert_less_equal ( std::fabs(eta), 1.0+10*TOLERANCE );

              const Real x = 0.5 * (xi + 1.0);
              const Real y = 0.5 * (eta + 1.0);

              switch (j)
                {
                  // d()/dxi
                case 0:
                  {
                    switch(ii)
                      {
                      case 0:
                        return sign * RealGradient( 0.25*(-18.0*y*y+24.0*y-6.0), 0.0 );
                      case 1:
                        return sign * RealGradient( 0.25*( 18.0*y*y-24.0*y+6.0), 0.0 );
                      case 2:
                        return sign * RealGradient( 0.0, 0.25*(-36.0*x*y+24.0*x+12.0*y-8.0) );
                      case 3:
                        return sign * RealGradient( 0.0, 0.25*( 36.0*x*y-12.0*x-12.0*y+4.0) );
                      case 4:
                        return sign * RealGradient( 0.25*(-18.0*y*y+12.0*y), 0.0 );
                      case 5:
                        return sign * RealGradient( 0.25*( 18.0*y*y-12.0*y), 0.0 );
                      case 6:
                        return sign * RealGradient( 0.0, 0.25*(-36.0*x*y+12.0*x+24.0*y-8.0) );
                      case 7:
                        return sign * RealGradient( 0.0, 0.25*( 36.0*x*y-24.0*x-24.0*y+16.0) );
                      case 8:
                        return RealGradient( 0.0, 1.5*(6.0*x*y-4.0*x-3.0*y+2.0) );
                      case 9:
                        return RealGradient( 1.5*y*(-3.0*y+3.0), 0.0 );
                      case 10:
                        return RealGradient( 1.5*y*(3.0*y-3.0), 0.0 );
                      case 11:
                        return RealGradient( 0.0, 1.5*(-6.0*x*y+2.0*x+3.0*y-1.0) );

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  } // j = 0

                  // d()/deta
                case 1:
                  {
                    switch(ii)
                      {
                      case 0:
                        return sign * RealGradient( 0.25*(-36.0*x*y+24.0*x+24.0*y-16.0), 0.0 );
                      case 1:
                        return sign * RealGradient( 0.25*( 36.0*x*y-24.0*x-12.0*y+8.0), 0.0 );
                      case 2:
                        return sign * RealGradient( 0.0, 0.25*x*(-18.0*x+12.0) );
                      case 3:
                        return sign * RealGradient( 0.0, 0.25*x*( 18.0*x-12.0) );
                      case 4:
                        return sign * RealGradient( 0.25*(-36.0*x*y+12.0*x+12.0*y-4.0), 0.0 );
                      case 5:
                        return sign * RealGradient( 0.25*( 36.0*x*y-12.0*x-24.0*y+8.0), 0.0 );
                      case 6:
                        return sign * RealGradient( 0.0, 0.25*(-18.0*x*x+24.0*x-6.0) );
                      case 7:
                        return sign * RealGradient( 0.0, 0.25*( 18.0*x*x-24.0*x+6.0) );
                      case 8:
                        return RealGradient( 0.0, 1.5*x*(3.0*x-3.0) );
                      case 9:
                        return RealGradient( 1.5*(-6.0*x*y+3.0*x+4.0*y-2.0), 0.0 );
                      case 10:
                        return RealGradient( 1.5*(6.0*x*y-3.0*x-2.0*y+1.0), 0.0 );
                      case 11:
                        return RealGradient( 0.0, 1.5*x*(-3.0*x+3.0) );

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  } // j = 1

                default:
                  libmesh_error_msg("Invalid j = " << j);
                }
            }

          case TRI6:
          case TRI7:
            {
              switch (j)
                {
                  // d()/dxi
                case 0:
                  {
                    switch(ii)
                      {
                      case 0:
                        return sign * RealGradient( 8.0*eta-6.0, -16.0*xi-8.0*eta+6.0 );
                      case 1:
                        return sign * RealGradient( -8.0*eta+6.0, 16.0*xi-4.0 );
                      case 2:
                        return sign * RealGradient( -8.0*eta, 16.0*xi-4.0 );
                      case 3:
                        return sign * RealGradient( 0.0, 8.0*eta-2.0 );
                      case 4:
                        return sign * RealGradient( 0.0, 8.0*eta-2.0 );
                      case 5:
                        return sign * RealGradient( 8.0*eta, -16.0*xi-8.0*eta+12.0 );
                      case 6:
                        return RealGradient( -8.0*eta, 16.0*xi+16.0*eta-8.0 );
                      case 7:
                        return RealGradient( 16.0*eta, -32.0*xi-8.0*eta+16.0 );

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  } // j = 0

                  // d()/deta
                case 1:
                  {
                    switch(ii)
                      {
                      case 0:
                        return sign * RealGradient( 8.0*xi+16.0*eta-12.0, -8.0*xi );
                      case 1:
                        return sign * RealGradient( -8.0*xi+2.0, 0.0 );
                      case 2:
                        return sign * RealGradient( 2.0-8.0*xi, 0.0 );
                      case 3:
                        return sign * RealGradient( 4.0-16.0*eta, 8.0*xi );
                      case 4:
                        return sign * RealGradient( 4.0-16.0*eta,8.0*xi-6.0 );
                      case 5:
                        return sign * RealGradient( 8.0*xi+16*eta-6.0,-8.0*xi+6.0 );
                      case 6:
                        return RealGradient( -8.0*xi-32.0*eta+16.0, 16.0*xi );
                      case 7:
                        return RealGradient( 16.0*xi+16.0*eta-8.0,-8.0*xi );

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  } // j = 1

                default:
                  libmesh_error_msg("Invalid j = " << j);
                }
            }

          default:
            libmesh_error_msg("ERROR: Unsupported 2D element type!: " << Utility::enum_to_string(elem->type()));
          }
      }

      // unsupported order
    default:
      libmesh_error_msg("ERROR: Unsupported 2D FE order!: " << total_order);
    }
#else // LIBMESH_DIM > 1
  libmesh_ignore(elem, order, i, j, add_p_level);
  libmesh_not_implemented();
#endif
}



template <>
RealGradient FE<2,NEDELEC_ONE>::shape_deriv(const ElemType,
                                            const Order,
                                            const unsigned int,
                                            const unsigned int,
                                            const Point &)
{
  libmesh_error_msg("Nedelec elements require the element type \nbecause edge orientation is needed.");
  return RealGradient();
}

template <>
RealGradient FE<2,NEDELEC_ONE>::shape_deriv(const FEType fet,
                                            const Elem * elem,
                                            const unsigned int i,
                                            const unsigned int j,
                                            const Point & p,
                                            const bool add_p_level)
{
  return FE<2,NEDELEC_ONE>::shape_deriv(elem, fet.order, i, j, p, add_p_level);
}





#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

template <>
RealGradient FE<2,NEDELEC_ONE>::shape_second_deriv(const Elem * elem,
                                                   const Order order,
                                                   const unsigned int i,
                                                   const unsigned int j,
                                                   const Point & p,
                                                   const bool add_p_level)
{
#if LIBMESH_DIM > 1
  libmesh_assert(elem);

  // j = 0 ==> d^2 phi / dxi^2
  // j = 1 ==> d^2 phi / dxi deta
  // j = 2 ==> d^2 phi / deta^2
  libmesh_assert_less (j, 3);

  const Order total_order = static_cast<Order>(order + add_p_level * elem->p_level());
  libmesh_assert_less(i, n_dofs(elem->type(), total_order));

  const char sign = i >= total_order * elem->n_edges() || elem->point(i / total_order) > elem->point((i / total_order + 1) % elem->n_vertices()) ? 1 : -1;
  const unsigned int ii = sign > 0 ? i : (i / total_order * 2 + 1) * total_order - 1 - i;

  const Real xi  = p(0);
  const Real eta = p(1);

  switch (total_order)
    {
      // linear Nedelec (first kind) shape function second derivatives
    case FIRST:
      {
        switch (elem->type())
          {
          case QUAD8:
          case QUAD9:
          case TRI6:
          case TRI7:
            // All second derivatives for linear quads and triangles are zero.
            return RealGradient();

          default:
            libmesh_error_msg("ERROR: Unsupported 2D element type!: " << Utility::enum_to_string(elem->type()));

          } // end switch (type)
      } // end case FIRST

      // quadratic Nedelec (first kind) shape function second derivatives
    case SECOND:
      {
        switch (elem->type())
          {
          case QUAD8:
          case QUAD9:
            {
              // Even with a loose inverse_map tolerance we ought to
              // be nearly on the element interior in master
              // coordinates
              libmesh_assert_less_equal ( std::fabs(xi), 1.0+10*TOLERANCE );
              libmesh_assert_less_equal ( std::fabs(eta), 1.0+10*TOLERANCE );

              const Real x = 0.5 * (xi + 1.0);
              const Real y = 0.5 * (eta + 1.0);

              switch (j)
                {
                  //d^2 () / dxi^2
                case 0:
                  {
                    switch(ii)
                      {
                      case 0:
                      case 1:
                      case 4:
                      case 5:
                      case 9:
                      case 10:
                        return RealGradient();
                      case 2:
                        return sign * RealGradient( 0.0, 0.125*(-36.0*y+24.0) );
                      case 3:
                        return sign * RealGradient( 0.0, 0.125*( 36.0*y-12.0) );
                      case 6:
                        return sign * RealGradient( 0.0, 0.125*(-36.0*y+12.0) );
                      case 7:
                        return sign * RealGradient( 0.0, 0.125*( 36.0*y-24.0) );
                      case 8:
                        return RealGradient( 0.0,  0.75*(6.0*y-4.0) );
                      case 11:
                        return RealGradient( 0.0,  0.75*(-6.0*y+2.0) );

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  } // j = 0

                  // ^2()/dxi deta
                case 1:
                  {
                    switch(ii)
                      {
                      case 0:
                        return sign * RealGradient( 0.125*(-36.0*y+24.0), 0.0 );
                      case 1:
                        return sign * RealGradient( 0.125*( 36.0*y-24.0), 0.0 );
                      case 2:
                       return sign * RealGradient( 0.0, 0.125*(-36.0*x+12.0) );
                      case 3:
                       return sign * RealGradient( 0.0, 0.125*( 36.0*x-12.0) );
                      case 4:
                        return sign * RealGradient( 0.125*(-36.0*y+12.0), 0.0 );
                      case 5:
                        return sign * RealGradient( 0.125*( 36.0*y-12.0), 0.0 );
                      case 6:
                        return sign * RealGradient( 0.0, 0.125*(-36.0*x+24.0) );
                      case 7:
                        return sign * RealGradient( 0.0, 0.125*( 36.0*x-24.0) );
                      case 8:
                        return RealGradient( 0.0,  0.75*(6.0*x-3.0) );
                      case 9:
                        return RealGradient( 0.75*(-6.0*y), 0.0 );
                      case 10:
                        return RealGradient( 0.75*(6.0*y), 0.0 );
                      case 11:
                        return RealGradient( 0.0, 0.75*(-6.0*x+3.0) );

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  } // j = 1

                  //d^2 () / deta^2
                case 2:
                  {
                    switch(ii)
                      {
                      case 2:
                      case 3:
                      case 6:
                      case 7:
                      case 8:
                      case 11:
                        return RealGradient();
                      case 0:
                        return sign * RealGradient( 0.125*(-36.0*x+24.0), 0.0 );
                      case 1:
                        return sign * RealGradient( 0.125*( 36.0*x-12.0), 0.0 );
                      case 4:
                        return sign * RealGradient( 0.125*(-36.0*x+12.0), 0.0 );
                      case 5:
                        return sign * RealGradient( 0.125*( 36.0*x-24.0), 0.0 );
                      case 9:
                        return RealGradient( 0.75*(-6.0*x+4.0), 0.0 );
                      case 10:
                        return RealGradient( 0.75*( 6.0*x-2.0), 0.0 );

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  } // j = 2

                default:
                  libmesh_error_msg("Invalid j = " << j);
                }
            }

          case TRI6:
          case TRI7:
            {
              switch (j)
                {
                  //d^2 () / dxi^2
                case 0:
                  {
                    switch(ii)
                      {
                      case 3:
                      case 4:
                        return RealGradient();
                      case 0:
                        return sign * RealGradient( 0.0, -16.0 );
                      case 1:
                        return sign * RealGradient( 0.0, 16.0 );
                      case 2:
                        return sign * RealGradient( 0.0, 16.0 );
                      case 5:
                        return sign * RealGradient( 0.0, -16.0 );
                      case 6:
                        return RealGradient( 0.0, 16.0 );
                      case 7:
                        return RealGradient( 0.0,-32.0 );
                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  } // j = 0

                  // ^2()/dxi deta
                case 1:
                  {
                    switch(ii)
                      {
                      case 0:
                        return sign * RealGradient( 8.0, -8.0 );
                      case 1:
                        return sign * RealGradient( -8.0, 0.0 );
                      case 2:
                        return sign * RealGradient( -8.0, 0.0 );
                      case 3:
                        return sign * RealGradient( 0.0, 8.0 );
                      case 4:
                        return sign * RealGradient( 0.0, 8.0 );
                      case 5:
                        return sign * RealGradient( 8.0, -8.0 );
                      case 6:
                        return RealGradient( -8.0, 16.0 );
                      case 7:
                        return RealGradient( 16.0, -8.0 );
                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  } // j = 1

                  //d^2 () / deta^2
                case 2:
                  {
                    switch(ii)
                      {
                      case 1:
                      case 2:
                        return RealGradient();
                      case 0:
                        return sign * RealGradient( 16.0, 0.0 );
                      case 3:
                        return sign * RealGradient( -16.0, 0.0 );
                      case 4:
                        return sign * RealGradient( -16.0, 0.0 );
                      case 5:
                        return sign * RealGradient( 16.0, 0.0 );
                      case 6:
                        return RealGradient( -32.0, 0.0 );
                      case 7:
                        return RealGradient( 16.0, 0.0 );
                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  } // j = 2

                default:
                  libmesh_error_msg("Invalid j = " << j);
                }
            }

          default:
            libmesh_error_msg("ERROR: Unsupported 2D element type!: " << Utility::enum_to_string(elem->type()));

          } // end switch (type)
      } // end case SECOND

      // unsupported order
    default:
      libmesh_error_msg("ERROR: Unsupported 2D FE order!: " << total_order);

    } // end switch (order)

#else // LIBMESH_DIM > 1
  libmesh_assert(true || i || j);
  libmesh_ignore(elem, order, add_p_level);
  libmesh_not_implemented();
#endif
}



template <>
RealGradient FE<2,NEDELEC_ONE>::shape_second_deriv(const ElemType,
                                                   const Order,
                                                   const unsigned int,
                                                   const unsigned int,
                                                   const Point &)
{
  libmesh_error_msg("Nedelec elements require the element type \nbecause edge orientation is needed.");
  return RealGradient();
}


template <>
RealGradient FE<2,NEDELEC_ONE>::shape_second_deriv(const FEType fet,
                                                   const Elem * elem,
                                                   const unsigned int i,
                                                   const unsigned int j,
                                                   const Point & p,
                                                   const bool add_p_level)
{
  return FE<2,NEDELEC_ONE>::shape_second_deriv(elem, fet.order, i, j, p, add_p_level);
}



#endif

} // namespace libMesh
