// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


// C++ includes

// Local includes
#include "libmesh/fe.h"
#include "libmesh/elem.h"

namespace libMesh
{




template <>
Real FE<3,L2_LAGRANGE>::shape(const ElemType type,
                              const Order order,
                              const unsigned int i,
                              const Point & p)
{
#if LIBMESH_DIM == 3


  switch (order)
    {
      // linear Lagrange shape functions
    case FIRST:
      {
        switch (type)
          {
            // trilinear hexahedral shape functions
          case HEX8:
          case HEX20:
          case HEX27:
            {
              libmesh_assert_less (i, 8);

              // Compute hex shape functions as a tensor-product
              const Real xi   = p(0);
              const Real eta  = p(1);
              const Real zeta = p(2);

              //                                0  1  2  3  4  5  6  7
              static const unsigned int i0[] = {0, 1, 1, 0, 0, 1, 1, 0};
              static const unsigned int i1[] = {0, 0, 1, 1, 0, 0, 1, 1};
              static const unsigned int i2[] = {0, 0, 0, 0, 1, 1, 1, 1};

              return (FE<1,L2_LAGRANGE>::shape(EDGE2, FIRST, i0[i], xi)*
                      FE<1,L2_LAGRANGE>::shape(EDGE2, FIRST, i1[i], eta)*
                      FE<1,L2_LAGRANGE>::shape(EDGE2, FIRST, i2[i], zeta));
            }

            // linear tetrahedral shape functions
          case TET4:
          case TET10:
            {
              libmesh_assert_less (i, 4);

              // Area coordinates, pg. 205, Vol. I, Carey, Oden, Becker FEM
              const Real zeta1 = p(0);
              const Real zeta2 = p(1);
              const Real zeta3 = p(2);
              const Real zeta0 = 1. - zeta1 - zeta2 - zeta3;

              switch(i)
                {
                case 0:
                  return zeta0;

                case 1:
                  return zeta1;

                case 2:
                  return zeta2;

                case 3:
                  return zeta3;

                default:
                  libmesh_error_msg("Invalid i = " << i);
                }
            }

            // linear prism shape functions
          case PRISM6:
          case PRISM15:
          case PRISM18:
            {
              libmesh_assert_less (i, 6);

              // Compute prism shape functions as a tensor-product
              // of a triangle and an edge

              Point p2d(p(0),p(1));
              Point p1d(p(2));

              //                                0  1  2  3  4  5
              static const unsigned int i0[] = {0, 0, 0, 1, 1, 1};
              static const unsigned int i1[] = {0, 1, 2, 0, 1, 2};

              return (FE<2,L2_LAGRANGE>::shape(TRI3,  FIRST, i1[i], p2d)*
                      FE<1,L2_LAGRANGE>::shape(EDGE2, FIRST, i0[i], p1d));
            }

            // linear pyramid shape functions
          case PYRAMID5:
            {
              libmesh_assert_less (i, 5);

              const Real xi   = p(0);
              const Real eta  = p(1);
              const Real zeta = p(2);
              const Real eps  = 1.e-35;

              switch(i)
                {
                case 0:
                  return .25*(zeta + xi - 1.)*(zeta + eta - 1.)/((1. - zeta) + eps);

                case 1:
                  return .25*(zeta - xi - 1.)*(zeta + eta - 1.)/((1. - zeta) + eps);

                case 2:
                  return .25*(zeta - xi - 1.)*(zeta - eta - 1.)/((1. - zeta) + eps);

                case 3:
                  return .25*(zeta + xi - 1.)*(zeta - eta - 1.)/((1. - zeta) + eps);

                case 4:
                  return zeta;

                default:
                  libmesh_error_msg("Invalid i = " << i);
                }
            }


          default:
            libmesh_error_msg("ERROR: Unsupported 3D element type!: " << type);
          }
      }


      // quadratic Lagrange shape functions
    case SECOND:
      {
        switch (type)
          {

            // serendipity hexahedral quadratic shape functions
          case HEX20:
            {
              libmesh_assert_less (i, 20);

              const Real xi   = p(0);
              const Real eta  = p(1);
              const Real zeta = p(2);

              // these functions are defined for (x,y,z) in [0,1]^3
              // so transform the locations
              const Real x = .5*(xi   + 1.);
              const Real y = .5*(eta  + 1.);
              const Real z = .5*(zeta + 1.);

              switch (i)
                {
                case 0:
                  return (1. - x)*(1. - y)*(1. - z)*(1. - 2.*x - 2.*y - 2.*z);

                case 1:
                  return x*(1. - y)*(1. - z)*(2.*x - 2.*y - 2.*z - 1.);

                case 2:
                  return x*y*(1. - z)*(2.*x + 2.*y - 2.*z - 3.);

                case 3:
                  return (1. - x)*y*(1. - z)*(2.*y - 2.*x - 2.*z - 1.);

                case 4:
                  return (1. - x)*(1. - y)*z*(2.*z - 2.*x - 2.*y - 1.);

                case 5:
                  return x*(1. - y)*z*(2.*x - 2.*y + 2.*z - 3.);

                case 6:
                  return x*y*z*(2.*x + 2.*y + 2.*z - 5.);

                case 7:
                  return (1. - x)*y*z*(2.*y - 2.*x + 2.*z - 3.);

                case 8:
                  return 4.*x*(1. - x)*(1. - y)*(1. - z);

                case 9:
                  return 4.*x*y*(1. - y)*(1. - z);

                case 10:
                  return 4.*x*(1. - x)*y*(1. - z);

                case 11:
                  return 4.*(1. - x)*y*(1. - y)*(1. - z);

                case 12:
                  return 4.*(1. - x)*(1. - y)*z*(1. - z);

                case 13:
                  return 4.*x*(1. - y)*z*(1. - z);

                case 14:
                  return 4.*x*y*z*(1. - z);

                case 15:
                  return 4.*(1. - x)*y*z*(1. - z);

                case 16:
                  return 4.*x*(1. - x)*(1. - y)*z;

                case 17:
                  return 4.*x*y*(1. - y)*z;

                case 18:
                  return 4.*x*(1. - x)*y*z;

                case 19:
                  return 4.*(1. - x)*y*(1. - y)*z;

                default:
                  libmesh_error_msg("Invalid i = " << i);
                }
            }

            // triquadratic hexahedral shape functions
          case HEX27:
            {
              libmesh_assert_less (i, 27);

              // Compute hex shape functions as a tensor-product
              const Real xi   = p(0);
              const Real eta  = p(1);
              const Real zeta = p(2);

              // The only way to make any sense of this
              // is to look at the mgflo/mg2/mgf documentation
              // and make the cut-out cube!
              //                                0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26
              static const unsigned int i0[] = {0, 1, 1, 0, 0, 1, 1, 0, 2, 1, 2, 0, 0, 1, 1, 0, 2, 1, 2, 0, 2, 2, 1, 2, 0, 2, 2};
              static const unsigned int i1[] = {0, 0, 1, 1, 0, 0, 1, 1, 0, 2, 1, 2, 0, 0, 1, 1, 0, 2, 1, 2, 2, 0, 2, 1, 2, 2, 2};
              static const unsigned int i2[] = {0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 2, 2, 2, 2, 1, 1, 1, 1, 0, 2, 2, 2, 2, 1, 2};

              return (FE<1,L2_LAGRANGE>::shape(EDGE3, SECOND, i0[i], xi)*
                      FE<1,L2_LAGRANGE>::shape(EDGE3, SECOND, i1[i], eta)*
                      FE<1,L2_LAGRANGE>::shape(EDGE3, SECOND, i2[i], zeta));
            }

            // quadratic tetrahedral shape functions
          case TET10:
            {
              libmesh_assert_less (i, 10);

              // Area coordinates, pg. 205, Vol. I, Carey, Oden, Becker FEM
              const Real zeta1 = p(0);
              const Real zeta2 = p(1);
              const Real zeta3 = p(2);
              const Real zeta0 = 1. - zeta1 - zeta2 - zeta3;

              switch(i)
                {
                case 0:
                  return zeta0*(2.*zeta0 - 1.);

                case 1:
                  return zeta1*(2.*zeta1 - 1.);

                case 2:
                  return zeta2*(2.*zeta2 - 1.);

                case 3:
                  return zeta3*(2.*zeta3 - 1.);

                case 4:
                  return 4.*zeta0*zeta1;

                case 5:
                  return 4.*zeta1*zeta2;

                case 6:
                  return 4.*zeta2*zeta0;

                case 7:
                  return 4.*zeta0*zeta3;

                case 8:
                  return 4.*zeta1*zeta3;

                case 9:
                  return 4.*zeta2*zeta3;

                default:
                  libmesh_error_msg("Invalid i = " << i);
                }
            }

            // quadratic prism shape functions
          case PRISM18:
            {
              libmesh_assert_less (i, 18);

              // Compute prism shape functions as a tensor-product
              // of a triangle and an edge

              Point p2d(p(0),p(1));
              Point p1d(p(2));

              //                                0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17
              static const unsigned int i0[] = {0, 0, 0, 1, 1, 1, 0, 0, 0, 2, 2, 2, 1, 1, 1, 2, 2, 2};
              static const unsigned int i1[] = {0, 1, 2, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 3, 4, 5};

              return (FE<2,L2_LAGRANGE>::shape(TRI6,  SECOND, i1[i], p2d)*
                      FE<1,L2_LAGRANGE>::shape(EDGE3, SECOND, i0[i], p1d));
            }


          default:
            libmesh_error_msg("ERROR: Unsupported 3D element type!: " << type);
          }
      }


      // unsupported order
    default:
      libmesh_error_msg("ERROR: Unsupported 3D FE order!: " << order);
    }

#endif
}



template <>
Real FE<3,L2_LAGRANGE>::shape(const Elem * elem,
                              const Order order,
                              const unsigned int i,
                              const Point & p)
{
  libmesh_assert(elem);

  // call the orientation-independent shape functions
  return FE<3,L2_LAGRANGE>::shape(elem->type(), static_cast<Order>(order + elem->p_level()), i, p);
}




template <>
Real FE<3,L2_LAGRANGE>::shape_deriv(const ElemType type,
                                    const Order order,
                                    const unsigned int i,
                                    const unsigned int j,
                                    const Point & p)
{
#if LIBMESH_DIM == 3

  libmesh_assert_less (j, 3);

  switch (order)
    {
      // linear Lagrange shape functions
    case FIRST:
      {
        switch (type)
          {
            // trilinear hexahedral shape functions
          case HEX8:
          case HEX20:
          case HEX27:
            {
              libmesh_assert_less (i, 8);

              // Compute hex shape functions as a tensor-product
              const Real xi   = p(0);
              const Real eta  = p(1);
              const Real zeta = p(2);

              static const unsigned int i0[] = {0, 1, 1, 0, 0, 1, 1, 0};
              static const unsigned int i1[] = {0, 0, 1, 1, 0, 0, 1, 1};
              static const unsigned int i2[] = {0, 0, 0, 0, 1, 1, 1, 1};

              switch(j)
                {
                case 0:
                  return (FE<1,L2_LAGRANGE>::shape_deriv(EDGE2, FIRST, i0[i], 0, xi)*
                          FE<1,L2_LAGRANGE>::shape      (EDGE2, FIRST, i1[i], eta)*
                          FE<1,L2_LAGRANGE>::shape      (EDGE2, FIRST, i2[i], zeta));

                case 1:
                  return (FE<1,L2_LAGRANGE>::shape      (EDGE2, FIRST, i0[i], xi)*
                          FE<1,L2_LAGRANGE>::shape_deriv(EDGE2, FIRST, i1[i], 0, eta)*
                          FE<1,L2_LAGRANGE>::shape      (EDGE2, FIRST, i2[i], zeta));

                case 2:
                  return (FE<1,L2_LAGRANGE>::shape      (EDGE2, FIRST, i0[i], xi)*
                          FE<1,L2_LAGRANGE>::shape      (EDGE2, FIRST, i1[i], eta)*
                          FE<1,L2_LAGRANGE>::shape_deriv(EDGE2, FIRST, i2[i], 0, zeta));

                default:
                  libmesh_error_msg("Invalid j = " << j);
                }
            }

            // linear tetrahedral shape functions
          case TET4:
          case TET10:
            {
              libmesh_assert_less (i, 4);

              // Area coordinates, pg. 205, Vol. I, Carey, Oden, Becker FEM
              const Real dzeta0dxi = -1.;
              const Real dzeta1dxi =  1.;
              const Real dzeta2dxi =  0.;
              const Real dzeta3dxi =  0.;

              const Real dzeta0deta = -1.;
              const Real dzeta1deta =  0.;
              const Real dzeta2deta =  1.;
              const Real dzeta3deta =  0.;

              const Real dzeta0dzeta = -1.;
              const Real dzeta1dzeta =  0.;
              const Real dzeta2dzeta =  0.;
              const Real dzeta3dzeta =  1.;

              switch (j)
                {
                  // d()/dxi
                case 0:
                  {
                    switch(i)
                      {
                      case 0:
                        return dzeta0dxi;

                      case 1:
                        return dzeta1dxi;

                      case 2:
                        return dzeta2dxi;

                      case 3:
                        return dzeta3dxi;

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  }

                  // d()/deta
                case 1:
                  {
                    switch(i)
                      {
                      case 0:
                        return dzeta0deta;

                      case 1:
                        return dzeta1deta;

                      case 2:
                        return dzeta2deta;

                      case 3:
                        return dzeta3deta;

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  }

                  // d()/dzeta
                case 2:
                  {
                    switch(i)
                      {
                      case 0:
                        return dzeta0dzeta;

                      case 1:
                        return dzeta1dzeta;

                      case 2:
                        return dzeta2dzeta;

                      case 3:
                        return dzeta3dzeta;

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  }

                default:
                  libmesh_error_msg("Invalid shape function derivative j = " << j);
                }
            }

            // linear prism shape functions
          case PRISM6:
          case PRISM15:
          case PRISM18:
            {
              libmesh_assert_less (i, 6);

              // Compute prism shape functions as a tensor-product
              // of a triangle and an edge

              Point p2d(p(0),p(1));
              Point p1d(p(2));

              //                                0  1  2  3  4  5
              static const unsigned int i0[] = {0, 0, 0, 1, 1, 1};
              static const unsigned int i1[] = {0, 1, 2, 0, 1, 2};

              switch (j)
                {
                  // d()/dxi
                case 0:
                  return (FE<2,L2_LAGRANGE>::shape_deriv(TRI3,  FIRST, i1[i], 0, p2d)*
                          FE<1,L2_LAGRANGE>::shape(EDGE2, FIRST, i0[i], p1d));

                  // d()/deta
                case 1:
                  return (FE<2,L2_LAGRANGE>::shape_deriv(TRI3,  FIRST, i1[i], 1, p2d)*
                          FE<1,L2_LAGRANGE>::shape(EDGE2, FIRST, i0[i], p1d));

                  // d()/dzeta
                case 2:
                  return (FE<2,L2_LAGRANGE>::shape(TRI3,  FIRST, i1[i], p2d)*
                          FE<1,L2_LAGRANGE>::shape_deriv(EDGE2, FIRST, i0[i], 0, p1d));

                default:
                  libmesh_error_msg("Invalid shape function derivative j = " << j);
                }
            }

            // linear pyramid shape functions
          case PYRAMID5:
            {
              libmesh_assert_less (i, 5);

              const Real xi   = p(0);
              const Real eta  = p(1);
              const Real zeta = p(2);
              const Real eps  = 1.e-35;

              switch (j)
                {
                  // d/dxi
                case 0:
                  switch(i)
                    {
                    case 0:
                      return  .25*(zeta + eta - 1.)/((1. - zeta) + eps);

                    case 1:
                      return -.25*(zeta + eta - 1.)/((1. - zeta) + eps);

                    case 2:
                      return -.25*(zeta - eta - 1.)/((1. - zeta) + eps);

                    case 3:
                      return  .25*(zeta - eta - 1.)/((1. - zeta) + eps);

                    case 4:
                      return 0;

                    default:
                      libmesh_error_msg("Invalid i = " << i);
                    }


                  // d/deta
                case 1:
                  switch(i)
                    {
                    case 0:
                      return  .25*(zeta + xi - 1.)/((1. - zeta) + eps);

                    case 1:
                      return  .25*(zeta - xi - 1.)/((1. - zeta) + eps);

                    case 2:
                      return -.25*(zeta - xi - 1.)/((1. - zeta) + eps);

                    case 3:
                      return -.25*(zeta + xi - 1.)/((1. - zeta) + eps);

                    case 4:
                      return 0;

                    default:
                      libmesh_error_msg("Invalid i = " << i);
                    }


                  // d/dzeta
                case 2:
                  switch(i)
                    {
                    case 0:
                      {
                        const Real a=1.;
                        const Real b=1.;

                        return .25*(((zeta + a*xi - 1.)*(zeta + b*eta - 1.) +
                                     (1. - zeta)*((zeta + a*xi -1.) + (zeta + b*eta - 1.)))/
                                    ((1. - zeta)*(1. - zeta) + eps));
                      }

                    case 1:
                      {
                        const Real a=-1.;
                        const Real b=1.;

                        return .25*(((zeta + a*xi - 1.)*(zeta + b*eta - 1.) +
                                     (1. - zeta)*((zeta + a*xi -1.) + (zeta + b*eta - 1.)))/
                                    ((1. - zeta)*(1. - zeta) + eps));
                      }

                    case 2:
                      {
                        const Real a=-1.;
                        const Real b=-1.;

                        return .25*(((zeta + a*xi - 1.)*(zeta + b*eta - 1.) +
                                     (1. - zeta)*((zeta + a*xi -1.) + (zeta + b*eta - 1.)))/
                                    ((1. - zeta)*(1. - zeta) + eps));
                      }

                    case 3:
                      {
                        const Real a=1.;
                        const Real b=-1.;

                        return .25*(((zeta + a*xi - 1.)*(zeta + b*eta - 1.) +
                                     (1. - zeta)*((zeta + a*xi -1.) + (zeta + b*eta - 1.)))/
                                    ((1. - zeta)*(1. - zeta) + eps));
                      }

                    case 4:
                      return 1.;

                    default:
                      libmesh_error_msg("Invalid i = " << i);
                    }


                default:
                  libmesh_error_msg("Invalid j = " << j);
                }
            }


          default:
            libmesh_error_msg("ERROR: Unsupported 3D element type!: " << type);
          }
      }


      // quadratic Lagrange shape functions
    case SECOND:
      {
        switch (type)
          {

            // serendipity hexahedral quadratic shape functions
          case HEX20:
            {
              libmesh_assert_less (i, 20);

              const Real xi   = p(0);
              const Real eta  = p(1);
              const Real zeta = p(2);

              // these functions are defined for (x,y,z) in [0,1]^3
              // so transform the locations
              const Real x = .5*(xi   + 1.);
              const Real y = .5*(eta  + 1.);
              const Real z = .5*(zeta + 1.);

              // and don't forget the chain rule!

              switch (j)
                {

                  // d/dx*dx/dxi
                case 0:
                  switch (i)
                    {
                    case 0:
                      return .5*(1. - y)*(1. - z)*((1. - x)*(-2.) +
                                                   (-1.)*(1. - 2.*x - 2.*y - 2.*z));

                    case 1:
                      return .5*(1. - y)*(1. - z)*(x*(2.) +
                                                   (1.)*(2.*x - 2.*y - 2.*z - 1.));

                    case 2:
                      return .5*y*(1. - z)*(x*(2.) +
                                            (1.)*(2.*x + 2.*y - 2.*z - 3.));

                    case 3:
                      return .5*y*(1. - z)*((1. - x)*(-2.) +
                                            (-1.)*(2.*y - 2.*x - 2.*z - 1.));

                    case 4:
                      return .5*(1. - y)*z*((1. - x)*(-2.) +
                                            (-1.)*(2.*z - 2.*x - 2.*y - 1.));

                    case 5:
                      return .5*(1. - y)*z*(x*(2.) +
                                            (1.)*(2.*x - 2.*y + 2.*z - 3.));

                    case 6:
                      return .5*y*z*(x*(2.) +
                                     (1.)*(2.*x + 2.*y + 2.*z - 5.));

                    case 7:
                      return .5*y*z*((1. - x)*(-2.) +
                                     (-1.)*(2.*y - 2.*x + 2.*z - 3.));

                    case 8:
                      return 2.*(1. - y)*(1. - z)*(1. - 2.*x);

                    case 9:
                      return 2.*y*(1. - y)*(1. - z);

                    case 10:
                      return 2.*y*(1. - z)*(1. - 2.*x);

                    case 11:
                      return 2.*y*(1. - y)*(1. - z)*(-1.);

                    case 12:
                      return 2.*(1. - y)*z*(1. - z)*(-1.);

                    case 13:
                      return 2.*(1. - y)*z*(1. - z);

                    case 14:
                      return 2.*y*z*(1. - z);

                    case 15:
                      return 2.*y*z*(1. - z)*(-1.);

                    case 16:
                      return 2.*(1. - y)*z*(1. - 2.*x);

                    case 17:
                      return 2.*y*(1. - y)*z;

                    case 18:
                      return 2.*y*z*(1. - 2.*x);

                    case 19:
                      return 2.*y*(1. - y)*z*(-1.);

                    default:
                      libmesh_error_msg("Invalid i = " << i);
                    }


                  // d/dy*dy/deta
                case 1:
                  switch (i)
                    {
                    case 0:
                      return .5*(1. - x)*(1. - z)*((1. - y)*(-2.) +
                                                   (-1.)*(1. - 2.*x - 2.*y - 2.*z));

                    case 1:
                      return .5*x*(1. - z)*((1. - y)*(-2.) +
                                            (-1.)*(2.*x - 2.*y - 2.*z - 1.));

                    case 2:
                      return .5*x*(1. - z)*(y*(2.) +
                                            (1.)*(2.*x + 2.*y - 2.*z - 3.));

                    case 3:
                      return .5*(1. - x)*(1. - z)*(y*(2.) +
                                                   (1.)*(2.*y - 2.*x - 2.*z - 1.));

                    case 4:
                      return .5*(1. - x)*z*((1. - y)*(-2.) +
                                            (-1.)*(2.*z - 2.*x - 2.*y - 1.));

                    case 5:
                      return .5*x*z*((1. - y)*(-2.) +
                                     (-1.)*(2.*x - 2.*y + 2.*z - 3.));

                    case 6:
                      return .5*x*z*(y*(2.) +
                                     (1.)*(2.*x + 2.*y + 2.*z - 5.));

                    case 7:
                      return .5*(1. - x)*z*(y*(2.) +
                                            (1.)*(2.*y - 2.*x + 2.*z - 3.));

                    case 8:
                      return 2.*x*(1. - x)*(1. - z)*(-1.);

                    case 9:
                      return 2.*x*(1. - z)*(1. - 2.*y);

                    case 10:
                      return 2.*x*(1. - x)*(1. - z);

                    case 11:
                      return 2.*(1. - x)*(1. - z)*(1. - 2.*y);

                    case 12:
                      return 2.*(1. - x)*z*(1. - z)*(-1.);

                    case 13:
                      return 2.*x*z*(1. - z)*(-1.);

                    case 14:
                      return 2.*x*z*(1. - z);

                    case 15:
                      return 2.*(1. - x)*z*(1. - z);

                    case 16:
                      return 2.*x*(1. - x)*z*(-1.);

                    case 17:
                      return 2.*x*z*(1. - 2.*y);

                    case 18:
                      return 2.*x*(1. - x)*z;

                    case 19:
                      return 2.*(1. - x)*z*(1. - 2.*y);

                    default:
                      libmesh_error_msg("Invalid i = " << i);
                    }


                  // d/dz*dz/dzeta
                case 2:
                  switch (i)
                    {
                    case 0:
                      return .5*(1. - x)*(1. - y)*((1. - z)*(-2.) +
                                                   (-1.)*(1. - 2.*x - 2.*y - 2.*z));

                    case 1:
                      return .5*x*(1. - y)*((1. - z)*(-2.) +
                                            (-1.)*(2.*x - 2.*y - 2.*z - 1.));

                    case 2:
                      return .5*x*y*((1. - z)*(-2.) +
                                     (-1.)*(2.*x + 2.*y - 2.*z - 3.));

                    case 3:
                      return .5*(1. - x)*y*((1. - z)*(-2.) +
                                            (-1.)*(2.*y - 2.*x - 2.*z - 1.));

                    case 4:
                      return .5*(1. - x)*(1. - y)*(z*(2.) +
                                                   (1.)*(2.*z - 2.*x - 2.*y - 1.));

                    case 5:
                      return .5*x*(1. - y)*(z*(2.) +
                                            (1.)*(2.*x - 2.*y + 2.*z - 3.));

                    case 6:
                      return .5*x*y*(z*(2.) +
                                     (1.)*(2.*x + 2.*y + 2.*z - 5.));

                    case 7:
                      return .5*(1. - x)*y*(z*(2.) +
                                            (1.)*(2.*y - 2.*x + 2.*z - 3.));

                    case 8:
                      return 2.*x*(1. - x)*(1. - y)*(-1.);

                    case 9:
                      return 2.*x*y*(1. - y)*(-1.);

                    case 10:
                      return 2.*x*(1. - x)*y*(-1.);

                    case 11:
                      return 2.*(1. - x)*y*(1. - y)*(-1.);

                    case 12:
                      return 2.*(1. - x)*(1. - y)*(1. - 2.*z);

                    case 13:
                      return 2.*x*(1. - y)*(1. - 2.*z);

                    case 14:
                      return 2.*x*y*(1. - 2.*z);

                    case 15:
                      return 2.*(1. - x)*y*(1. - 2.*z);

                    case 16:
                      return 2.*x*(1. - x)*(1. - y);

                    case 17:
                      return 2.*x*y*(1. - y);

                    case 18:
                      return 2.*x*(1. - x)*y;

                    case 19:
                      return 2.*(1. - x)*y*(1. - y);

                    default:
                      libmesh_error_msg("Invalid i = " << i);
                    }

                default:
                  libmesh_error_msg("Invalid shape function derivative j = " << j);
                }
            }

            // triquadratic hexahedral shape functions
          case HEX27:
            {
              libmesh_assert_less (i, 27);

              // Compute hex shape functions as a tensor-product
              const Real xi   = p(0);
              const Real eta  = p(1);
              const Real zeta = p(2);

              // The only way to make any sense of this
              // is to look at the mgflo/mg2/mgf documentation
              // and make the cut-out cube!
              //                                0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26
              static const unsigned int i0[] = {0, 1, 1, 0, 0, 1, 1, 0, 2, 1, 2, 0, 0, 1, 1, 0, 2, 1, 2, 0, 2, 2, 1, 2, 0, 2, 2};
              static const unsigned int i1[] = {0, 0, 1, 1, 0, 0, 1, 1, 0, 2, 1, 2, 0, 0, 1, 1, 0, 2, 1, 2, 2, 0, 2, 1, 2, 2, 2};
              static const unsigned int i2[] = {0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 2, 2, 2, 2, 1, 1, 1, 1, 0, 2, 2, 2, 2, 1, 2};

              switch(j)
                {
                case 0:
                  return (FE<1,L2_LAGRANGE>::shape_deriv(EDGE3, SECOND, i0[i], 0, xi)*
                          FE<1,L2_LAGRANGE>::shape      (EDGE3, SECOND, i1[i], eta)*
                          FE<1,L2_LAGRANGE>::shape      (EDGE3, SECOND, i2[i], zeta));

                case 1:
                  return (FE<1,L2_LAGRANGE>::shape      (EDGE3, SECOND, i0[i], xi)*
                          FE<1,L2_LAGRANGE>::shape_deriv(EDGE3, SECOND, i1[i], 0, eta)*
                          FE<1,L2_LAGRANGE>::shape      (EDGE3, SECOND, i2[i], zeta));

                case 2:
                  return (FE<1,L2_LAGRANGE>::shape      (EDGE3, SECOND, i0[i], xi)*
                          FE<1,L2_LAGRANGE>::shape      (EDGE3, SECOND, i1[i], eta)*
                          FE<1,L2_LAGRANGE>::shape_deriv(EDGE3, SECOND, i2[i], 0, zeta));

                default:
                  libmesh_error_msg("Invalid j = " << j);
                }
            }

            // quadratic tetrahedral shape functions
          case TET10:
            {
              libmesh_assert_less (i, 10);

              // Area coordinates, pg. 205, Vol. I, Carey, Oden, Becker FEM
              const Real zeta1 = p(0);
              const Real zeta2 = p(1);
              const Real zeta3 = p(2);
              const Real zeta0 = 1. - zeta1 - zeta2 - zeta3;

              const Real dzeta0dxi = -1.;
              const Real dzeta1dxi =  1.;
              const Real dzeta2dxi =  0.;
              const Real dzeta3dxi =  0.;

              const Real dzeta0deta = -1.;
              const Real dzeta1deta =  0.;
              const Real dzeta2deta =  1.;
              const Real dzeta3deta =  0.;

              const Real dzeta0dzeta = -1.;
              const Real dzeta1dzeta =  0.;
              const Real dzeta2dzeta =  0.;
              const Real dzeta3dzeta =  1.;

              switch (j)
                {
                  // d()/dxi
                case 0:
                  {
                    switch(i)
                      {
                      case 0:
                        return (4.*zeta0 - 1.)*dzeta0dxi;

                      case 1:
                        return (4.*zeta1 - 1.)*dzeta1dxi;

                      case 2:
                        return (4.*zeta2 - 1.)*dzeta2dxi;

                      case 3:
                        return (4.*zeta3 - 1.)*dzeta3dxi;

                      case 4:
                        return 4.*(zeta0*dzeta1dxi + dzeta0dxi*zeta1);

                      case 5:
                        return 4.*(zeta1*dzeta2dxi + dzeta1dxi*zeta2);

                      case 6:
                        return 4.*(zeta0*dzeta2dxi + dzeta0dxi*zeta2);

                      case 7:
                        return 4.*(zeta0*dzeta3dxi + dzeta0dxi*zeta3);

                      case 8:
                        return 4.*(zeta1*dzeta3dxi + dzeta1dxi*zeta3);

                      case 9:
                        return 4.*(zeta2*dzeta3dxi + dzeta2dxi*zeta3);

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  }

                  // d()/deta
                case 1:
                  {
                    switch(i)
                      {
                      case 0:
                        return (4.*zeta0 - 1.)*dzeta0deta;

                      case 1:
                        return (4.*zeta1 - 1.)*dzeta1deta;

                      case 2:
                        return (4.*zeta2 - 1.)*dzeta2deta;

                      case 3:
                        return (4.*zeta3 - 1.)*dzeta3deta;

                      case 4:
                        return 4.*(zeta0*dzeta1deta + dzeta0deta*zeta1);

                      case 5:
                        return 4.*(zeta1*dzeta2deta + dzeta1deta*zeta2);

                      case 6:
                        return 4.*(zeta0*dzeta2deta + dzeta0deta*zeta2);

                      case 7:
                        return 4.*(zeta0*dzeta3deta + dzeta0deta*zeta3);

                      case 8:
                        return 4.*(zeta1*dzeta3deta + dzeta1deta*zeta3);

                      case 9:
                        return 4.*(zeta2*dzeta3deta + dzeta2deta*zeta3);

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  }

                  // d()/dzeta
                case 2:
                  {
                    switch(i)
                      {
                      case 0:
                        return (4.*zeta0 - 1.)*dzeta0dzeta;

                      case 1:
                        return (4.*zeta1 - 1.)*dzeta1dzeta;

                      case 2:
                        return (4.*zeta2 - 1.)*dzeta2dzeta;

                      case 3:
                        return (4.*zeta3 - 1.)*dzeta3dzeta;

                      case 4:
                        return 4.*(zeta0*dzeta1dzeta + dzeta0dzeta*zeta1);

                      case 5:
                        return 4.*(zeta1*dzeta2dzeta + dzeta1dzeta*zeta2);

                      case 6:
                        return 4.*(zeta0*dzeta2dzeta + dzeta0dzeta*zeta2);

                      case 7:
                        return 4.*(zeta0*dzeta3dzeta + dzeta0dzeta*zeta3);

                      case 8:
                        return 4.*(zeta1*dzeta3dzeta + dzeta1dzeta*zeta3);

                      case 9:
                        return 4.*(zeta2*dzeta3dzeta + dzeta2dzeta*zeta3);

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  }

                default:
                  libmesh_error_msg("Invalid j = " << j);
                }
            }



            // quadratic prism shape functions
          case PRISM18:
            {
              libmesh_assert_less (i, 18);

              // Compute prism shape functions as a tensor-product
              // of a triangle and an edge

              Point p2d(p(0),p(1));
              Point p1d(p(2));

              //                                0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17
              static const unsigned int i0[] = {0, 0, 0, 1, 1, 1, 0, 0, 0, 2, 2, 2, 1, 1, 1, 2, 2, 2};
              static const unsigned int i1[] = {0, 1, 2, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 3, 4, 5};

              switch (j)
                {
                  // d()/dxi
                case 0:
                  return (FE<2,L2_LAGRANGE>::shape_deriv(TRI6,  SECOND, i1[i], 0, p2d)*
                          FE<1,L2_LAGRANGE>::shape(EDGE3, SECOND, i0[i], p1d));

                  // d()/deta
                case 1:
                  return (FE<2,L2_LAGRANGE>::shape_deriv(TRI6,  SECOND, i1[i], 1, p2d)*
                          FE<1,L2_LAGRANGE>::shape(EDGE3, SECOND, i0[i], p1d));

                  // d()/dzeta
                case 2:
                  return (FE<2,L2_LAGRANGE>::shape(TRI6,  SECOND, i1[i], p2d)*
                          FE<1,L2_LAGRANGE>::shape_deriv(EDGE3, SECOND, i0[i], 0, p1d));

                default:
                  libmesh_error_msg("Invalid shape function derivative j = " << j);
                }
            }


          default:
            libmesh_error_msg("ERROR: Unsupported 3D element type!: " << type);
          }
      }


      // unsupported order
    default:
      libmesh_error_msg("ERROR: Unsupported 3D FE order!: " << order);
    }

#endif
}



template <>
Real FE<3,L2_LAGRANGE>::shape_deriv(const Elem * elem,
                                    const Order order,
                                    const unsigned int i,
                                    const unsigned int j,
                                    const Point & p)
{
  libmesh_assert(elem);

  // call the orientation-independent shape function derivatives
  return FE<3,L2_LAGRANGE>::shape_deriv(elem->type(), static_cast<Order>(order + elem->p_level()), i, j, p);
}



template <>
Real FE<3,L2_LAGRANGE>::shape_second_deriv(const ElemType type,
                                           const Order order,
                                           const unsigned int i,
                                           const unsigned int j,
                                           const Point & p)
{
#if LIBMESH_DIM == 3

  libmesh_assert_less (j, 6);

  switch (order)
    {
      // linear Lagrange shape functions
    case FIRST:
      {
        return 0.;
      }

      // quadratic Lagrange shape functions
    case SECOND:
      {
        switch (type)
          {

            // serendipity hexahedral quadratic shape functions
          case HEX20:
            {
              static bool warning_given_HEX20 = false;

              if (!warning_given_HEX20)
                libMesh::err << "Second derivatives for 3D Lagrangian HEX20"
                             << " elements are not yet implemented!"
                             << std::endl;
              warning_given_HEX20 = true;
            }
            libmesh_fallthrough();

          case HEX27:
            // triquadratic hexahedral shape functions
            {
              libmesh_assert_less (i, 27);

              // Compute hex shape functions as a tensor-product
              const Real xi   = p(0);
              const Real eta  = p(1);
              const Real zeta = p(2);

              // The only way to make any sense of this
              // is to look at the mgflo/mg2/mgf documentation
              // and make the cut-out cube!
              //                                0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26
              static const unsigned int i0[] = {0, 1, 1, 0, 0, 1, 1, 0, 2, 1, 2, 0, 0, 1, 1, 0, 2, 1, 2, 0, 2, 2, 1, 2, 0, 2, 2};
              static const unsigned int i1[] = {0, 0, 1, 1, 0, 0, 1, 1, 0, 2, 1, 2, 0, 0, 1, 1, 0, 2, 1, 2, 2, 0, 2, 1, 2, 2, 2};
              static const unsigned int i2[] = {0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 2, 2, 2, 2, 1, 1, 1, 1, 0, 2, 2, 2, 2, 1, 2};

              switch(j)
                {
                  // d^2()/dxi^2
                case 0:
                  return (FE<1,L2_LAGRANGE>::shape_second_deriv(EDGE3, SECOND, i0[i], 0, xi)*
                          FE<1,L2_LAGRANGE>::shape      (EDGE3, SECOND, i1[i], eta)*
                          FE<1,L2_LAGRANGE>::shape      (EDGE3, SECOND, i2[i], zeta));

                  // d^2()/dxideta
                case 1:
                  return (FE<1,L2_LAGRANGE>::shape_deriv(EDGE3, SECOND, i0[i], 0, xi)*
                          FE<1,L2_LAGRANGE>::shape_deriv(EDGE3, SECOND, i1[i], 0, eta)*
                          FE<1,L2_LAGRANGE>::shape      (EDGE3, SECOND, i2[i], zeta));

                  // d^2()/deta^2
                case 2:
                  return (FE<1,L2_LAGRANGE>::shape      (EDGE3, SECOND, i0[i], xi)*
                          FE<1,L2_LAGRANGE>::shape_second_deriv(EDGE3, SECOND, i1[i], 0, eta)*
                          FE<1,L2_LAGRANGE>::shape      (EDGE3, SECOND, i2[i], zeta));

                  // d^2()/dxidzeta
                case 3:
                  return (FE<1,L2_LAGRANGE>::shape_deriv(EDGE3, SECOND, i0[i], 0, xi)*
                          FE<1,L2_LAGRANGE>::shape      (EDGE3, SECOND, i1[i], eta)*
                          FE<1,L2_LAGRANGE>::shape_deriv(EDGE3, SECOND, i2[i], 0, zeta));

                  // d^2()/detadzeta
                case 4:
                  return (FE<1,L2_LAGRANGE>::shape      (EDGE3, SECOND, i0[i], xi)*
                          FE<1,L2_LAGRANGE>::shape_deriv(EDGE3, SECOND, i1[i], 0, eta)*
                          FE<1,L2_LAGRANGE>::shape_deriv(EDGE3, SECOND, i2[i], 0, zeta));

                  // d^2()/dzeta^2
                case 5:
                  return (FE<1,L2_LAGRANGE>::shape      (EDGE3, SECOND, i0[i], xi)*
                          FE<1,L2_LAGRANGE>::shape      (EDGE3, SECOND, i1[i], eta)*
                          FE<1,L2_LAGRANGE>::shape_second_deriv(EDGE3, SECOND, i2[i], 0, zeta));

                default:
                  libmesh_error_msg("Invalid j = " << j);
                }
            }

            // quadratic tetrahedral shape functions
          case TET10:
            {
              // The area coordinates are the same as used for the
              // shape() and shape_deriv() functions.
              // const Real zeta0 = 1. - zeta1 - zeta2 - zeta3;
              // const Real zeta1 = p(0);
              // const Real zeta2 = p(1);
              // const Real zeta3 = p(2);
              static const Real dzetadxi[4][3] =
                {
                  {-1., -1., -1.},
                  {1.,   0.,  0.},
                  {0.,   1.,  0.},
                  {0.,   0.,  1.}
                };

              // Convert from j -> (j,k) indices for independent variable
              // (0=xi, 1=eta, 2=zeta)
              static const unsigned short int independent_var_indices[6][2] =
                {
                  {0, 0}, // d^2 phi / dxi^2
                  {0, 1}, // d^2 phi / dxi deta
                  {1, 1}, // d^2 phi / deta^2
                  {0, 2}, // d^2 phi / dxi dzeta
                  {1, 2}, // d^2 phi / deta dzeta
                  {2, 2}  // d^2 phi / dzeta^2
                };

              // Convert from i -> zeta indices.  Each quadratic shape
              // function for the Tet10 depends on up to two of the zeta
              // area coordinate functions (see the shape() function above).
              // This table just tells which two area coords it uses.
              static const unsigned short int zeta_indices[10][2] =
                {
                  {0, 0},
                  {1, 1},
                  {2, 2},
                  {3, 3},
                  {0, 1},
                  {1, 2},
                  {2, 0},
                  {0, 3},
                  {1, 3},
                  {2, 3},
                };

              // Look up the independent variable indices for this value of j.
              const unsigned int my_j = independent_var_indices[j][0];
              const unsigned int my_k = independent_var_indices[j][1];

              if (i<4)
                {
                  return 4.*dzetadxi[i][my_j]*dzetadxi[i][my_k];
                }

              else if (i<10)
                {
                  const unsigned short int my_m = zeta_indices[i][0];
                  const unsigned short int my_n = zeta_indices[i][1];

                  return 4.*(dzetadxi[my_n][my_j]*dzetadxi[my_m][my_k] +
                             dzetadxi[my_m][my_j]*dzetadxi[my_n][my_k] );
                }
              else
                libmesh_error_msg("Invalid shape function index " << i);
            }


            // quadratic prism shape functions
          case PRISM18:
            {
              libmesh_assert_less (i, 18);

              // Compute prism shape functions as a tensor-product
              // of a triangle and an edge

              Point p2d(p(0),p(1));
              Point p1d(p(2));

              //                                0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17
              static const unsigned int i0[] = {0, 0, 0, 1, 1, 1, 0, 0, 0, 2, 2, 2, 1, 1, 1, 2, 2, 2};
              static const unsigned int i1[] = {0, 1, 2, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 3, 4, 5};

              switch (j)
                {
                  // d^2()/dxi^2
                case 0:
                  return (FE<2,L2_LAGRANGE>::shape_second_deriv(TRI6, SECOND, i1[i], 0, p2d)*
                          FE<1,L2_LAGRANGE>::shape(EDGE3, SECOND, i0[i], p1d));

                  // d^2()/dxideta
                case 1:
                  return (FE<2,L2_LAGRANGE>::shape_second_deriv(TRI6, SECOND, i1[i], 1, p2d)*
                          FE<1,L2_LAGRANGE>::shape(EDGE3, SECOND, i0[i], p1d));

                  // d^2()/deta^2
                case 2:
                  return (FE<2,L2_LAGRANGE>::shape_second_deriv(TRI6, SECOND, i1[i], 2, p2d)*
                          FE<1,L2_LAGRANGE>::shape(EDGE3, SECOND, i0[i], p1d));

                  // d^2()/dxidzeta
                case 3:
                  return (FE<2,L2_LAGRANGE>::shape_deriv(TRI6,  SECOND, i1[i], 0, p2d)*
                          FE<1,L2_LAGRANGE>::shape_deriv(EDGE3, SECOND, i0[i], 0, p1d));

                  // d^2()/detadzeta
                case 4:
                  return (FE<2,L2_LAGRANGE>::shape_deriv(TRI6,  SECOND, i1[i], 1, p2d)*
                          FE<1,L2_LAGRANGE>::shape_deriv(EDGE3, SECOND, i0[i], 0, p1d));

                  // d^2()/dzeta^2
                case 5:
                  return (FE<2,L2_LAGRANGE>::shape(TRI6,  SECOND, i1[i], p2d)*
                          FE<1,L2_LAGRANGE>::shape_second_deriv(EDGE3, SECOND, i0[i], 0, p1d));

                default:
                  libmesh_error_msg("Invalid shape function derivative j = " << j);
                }
            }



          default:
            libmesh_error_msg("ERROR: Unsupported 3D element type!: " << type);
          }
      }


      // unsupported order
    default:
      libmesh_error_msg("ERROR: Unsupported 3D FE order!: " << order);
    }

#endif
}



template <>
Real FE<3,L2_LAGRANGE>::shape_second_deriv(const Elem * elem,
                                           const Order order,
                                           const unsigned int i,
                                           const unsigned int j,
                                           const Point & p)
{
  libmesh_assert(elem);

  // call the orientation-independent shape function derivatives
  return FE<3,L2_LAGRANGE>::shape_second_deriv(elem->type(), static_cast<Order>(order + elem->p_level()), i, j, p);
}

} // namespace libMesh
