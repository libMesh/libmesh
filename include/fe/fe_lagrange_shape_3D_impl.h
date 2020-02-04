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

#ifndef LIBMESH_FE_LAGRANGE_SHAPE_3D_IMPL_H
#define LIBMESH_FE_LAGRANGE_SHAPE_3D_IMPL_H

// Local includes
#include "libmesh/fe.h"
#include "libmesh/elem.h"
#include "libmesh/fe_lagrange_shape_1D.h"

// Anonymous namespace for functions shared by LAGRANGE and
// L2_LAGRANGE implementations. Implementations appear at the bottom
// of this file.
namespace
{
using namespace libMesh;

template <typename RealType>
RealType fe_lagrange_3D_shape(const ElemType,
                              const Order order,
                              const unsigned int i,
                              const PointTempl<RealType> & p);

template <typename RealType>
RealType fe_lagrange_3D_shape_deriv(const ElemType type,
                                    const Order order,
                                    const unsigned int i,
                                    const unsigned int j,
                                    const PointTempl<RealType> & p);

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

template <typename RealType>
RealType fe_lagrange_3D_shape_second_deriv(const ElemType type,
                                           const Order order,
                                           const unsigned int i,
                                           const unsigned int j,
                                           const PointTempl<RealType> & p);

#endif // LIBMESH_ENABLE_SECOND_DERIVATIVES

} // anonymous namespace

namespace libMesh
{

template <typename RealType>
typename FEShim<3,LAGRANGE,RealType>::OutputShape FEShim<3,LAGRANGE,RealType>::shape(const ElemType type,
                           const Order order,
                           const unsigned int i,
                           const Point & p)
{
  return fe_lagrange_3D_shape(type, order, i, p);
}



template <typename RealType>
typename FEShim<3,L2_LAGRANGE,RealType>::OutputShape FEShim<3,L2_LAGRANGE,RealType>::shape(const ElemType type,
                              const Order order,
                              const unsigned int i,
                              const Point & p)
{
  return fe_lagrange_3D_shape(type, order, i, p);
}



template <typename RealType>
typename FEShim<3,LAGRANGE,RealType>::OutputShape FEShim<3,LAGRANGE,RealType>::shape(const Elem * elem,
                           const Order order,
                           const unsigned int i,
                           const Point & p,
                           const bool add_p_level)
{
  libmesh_assert(elem);

  // call the orientation-independent shape functions
  return fe_lagrange_3D_shape(elem->type(), static_cast<Order>(order + add_p_level * elem->p_level()), i, p);
}



template <typename RealType>
typename FEShim<3,L2_LAGRANGE,RealType>::OutputShape FEShim<3,L2_LAGRANGE,RealType>::shape(const Elem * elem,
                              const Order order,
                              const unsigned int i,
                              const Point & p,
                              const bool add_p_level)
{
  libmesh_assert(elem);

  // call the orientation-independent shape functions
  return fe_lagrange_3D_shape(elem->type(), static_cast<Order>(order + add_p_level * elem->p_level()), i, p);
}



template <typename RealType>
typename FEShim<3,LAGRANGE,RealType>::OutputShape FEShim<3,LAGRANGE,RealType>::shape_deriv(const ElemType type,
                                 const Order order,
                                 const unsigned int i,
                                 const unsigned int j,
                                 const Point & p)
{
  return fe_lagrange_3D_shape_deriv(type, order, i, j, p);
}



template <typename RealType>
typename FEShim<3,L2_LAGRANGE,RealType>::OutputShape FEShim<3,L2_LAGRANGE,RealType>::shape_deriv(const ElemType type,
                                    const Order order,
                                    const unsigned int i,
                                    const unsigned int j,
                                    const Point & p)
{
  return fe_lagrange_3D_shape_deriv(type, order, i, j, p);
}



template <typename RealType>
typename FEShim<3,LAGRANGE,RealType>::OutputShape FEShim<3,LAGRANGE,RealType>::shape_deriv(const Elem * elem,
                                 const Order order,
                                 const unsigned int i,
                                 const unsigned int j,
                                 const Point & p,
                                 const bool add_p_level)
{
  libmesh_assert(elem);

  // call the orientation-independent shape function derivatives
  return fe_lagrange_3D_shape_deriv(elem->type(), static_cast<Order>(order + add_p_level * elem->p_level()), i, j, p);
}


template <typename RealType>
typename FEShim<3,L2_LAGRANGE,RealType>::OutputShape FEShim<3,L2_LAGRANGE,RealType>::shape_deriv(const Elem * elem,
                                    const Order order,
                                    const unsigned int i,
                                    const unsigned int j,
                                    const Point & p,
                                    const bool add_p_level)
{
  libmesh_assert(elem);

  // call the orientation-independent shape function derivatives
  return fe_lagrange_3D_shape_deriv(elem->type(), static_cast<Order>(order + add_p_level * elem->p_level()), i, j, p);
}



#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

template <typename RealType>
typename FEShim<3,LAGRANGE,RealType>::OutputShape FEShim<3,LAGRANGE,RealType>::shape_second_deriv(const ElemType type,
                                        const Order order,
                                        const unsigned int i,
                                        const unsigned int j,
                                        const Point & p)
{
  return fe_lagrange_3D_shape_second_deriv(type, order, i, j, p);
}



template <typename RealType>
typename FEShim<3,L2_LAGRANGE,RealType>::OutputShape FEShim<3,L2_LAGRANGE,RealType>::shape_second_deriv(const ElemType type,
                                           const Order order,
                                           const unsigned int i,
                                           const unsigned int j,
                                           const Point & p)
{
  return fe_lagrange_3D_shape_second_deriv(type, order, i, j, p);
}



template <typename RealType>
typename FEShim<3,LAGRANGE,RealType>::OutputShape FEShim<3,LAGRANGE,RealType>::shape_second_deriv(const Elem * elem,
                                        const Order order,
                                        const unsigned int i,
                                        const unsigned int j,
                                        const Point & p,
                                        const bool add_p_level)
{
  libmesh_assert(elem);

  // call the orientation-independent shape function derivatives
  return fe_lagrange_3D_shape_second_deriv
    (elem->type(), static_cast<Order>(order + add_p_level * elem->p_level()), i, j, p);
}



template <typename RealType>
typename FEShim<3,L2_LAGRANGE,RealType>::OutputShape FEShim<3,L2_LAGRANGE,RealType>::shape_second_deriv(const Elem * elem,
                                           const Order order,
                                           const unsigned int i,
                                           const unsigned int j,
                                           const Point & p,
                                           const bool add_p_level)
{
  libmesh_assert(elem);

  // call the orientation-independent shape function derivatives
  return fe_lagrange_3D_shape_second_deriv
    (elem->type(), static_cast<Order>(order + add_p_level * elem->p_level()), i, j, p);
}

#endif // LIBMESH_ENABLE_SECOND_DERIVATIVES

} // namespace libMesh



namespace
{
using namespace libMesh;

template <typename RealType>
RealType fe_lagrange_3D_shape(const ElemType type,
                              const Order order,
                              const unsigned int i,
                              const PointTempl<RealType> & p)
{
  typedef PointTempl<RealType> Point;

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
              const auto xi   = p(0);
              const auto eta  = p(1);
              const auto zeta = p(2);

              //                                0  1  2  3  4  5  6  7
              static const unsigned int i0[] = {0, 1, 1, 0, 0, 1, 1, 0};
              static const unsigned int i1[] = {0, 0, 1, 1, 0, 0, 1, 1};
              static const unsigned int i2[] = {0, 0, 0, 0, 1, 1, 1, 1};

              return (fe_lagrange_1D_linear_shape(i0[i], xi)*
                      fe_lagrange_1D_linear_shape(i1[i], eta)*
                      fe_lagrange_1D_linear_shape(i2[i], zeta));
            }

            // linear tetrahedral shape functions
          case TET4:
          case TET10:
            {
              libmesh_assert_less (i, 4);

              // Area coordinates, pg. 205, Vol. I, Carey, Oden, Becker FEM
              const auto zeta1 = p(0);
              const auto zeta2 = p(1);
              const auto zeta3 = p(2);
              const auto zeta0 = 1. - zeta1 - zeta2 - zeta3;

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
              auto p1d = p(2);

              //                                0  1  2  3  4  5
              static const unsigned int i0[] = {0, 0, 0, 1, 1, 1};
              static const unsigned int i1[] = {0, 1, 2, 0, 1, 2};

              return (FEShim<2,LAGRANGE,RealType>::shape(TRI3,  FIRST, i1[i], p2d)*
                      fe_lagrange_1D_linear_shape(i0[i], p1d));
            }

            // linear pyramid shape functions
          case PYRAMID5:
          case PYRAMID13:
          case PYRAMID14:
            {
              libmesh_assert_less (i, 5);

              const auto xi   = p(0);
              const auto eta  = p(1);
              const auto zeta = p(2);
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

              const auto xi   = p(0);
              const auto eta  = p(1);
              const auto zeta = p(2);

              // these functions are defined for (x,y,z) in [0,1]^3
              // so transform the locations
              const auto x = .5*(xi   + 1.);
              const auto y = .5*(eta  + 1.);
              const auto z = .5*(zeta + 1.);

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
              const auto xi   = p(0);
              const auto eta  = p(1);
              const auto zeta = p(2);

              // The only way to make any sense of this
              // is to look at the mgflo/mg2/mgf documentation
              // and make the cut-out cube!
              //                                0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26
              static const unsigned int i0[] = {0, 1, 1, 0, 0, 1, 1, 0, 2, 1, 2, 0, 0, 1, 1, 0, 2, 1, 2, 0, 2, 2, 1, 2, 0, 2, 2};
              static const unsigned int i1[] = {0, 0, 1, 1, 0, 0, 1, 1, 0, 2, 1, 2, 0, 0, 1, 1, 0, 2, 1, 2, 2, 0, 2, 1, 2, 2, 2};
              static const unsigned int i2[] = {0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 2, 2, 2, 2, 1, 1, 1, 1, 0, 2, 2, 2, 2, 1, 2};

              return (fe_lagrange_1D_quadratic_shape(i0[i], xi)*
                      fe_lagrange_1D_quadratic_shape(i1[i], eta)*
                      fe_lagrange_1D_quadratic_shape(i2[i], zeta));
            }

            // quadratic tetrahedral shape functions
          case TET10:
            {
              libmesh_assert_less (i, 10);

              // Area coordinates, pg. 205, Vol. I, Carey, Oden, Becker FEM
              const auto zeta1 = p(0);
              const auto zeta2 = p(1);
              const auto zeta3 = p(2);
              const auto zeta0 = 1. - zeta1 - zeta2 - zeta3;

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

            // "serendipity" prism
          case PRISM15:
            {
              libmesh_assert_less (i, 15);

              const auto xi   = p(0);
              const auto eta  = p(1);
              const auto zeta = p(2);

              switch(i)
                {
                case 0:
                  return (1. - zeta)*(xi + eta - 1.)*(xi + eta + 0.5*zeta);

                case 1:
                  return (1. - zeta)*xi*(xi - 1. - 0.5*zeta);

                case 2: // phi1 with xi <- eta
                  return (1. - zeta)*eta*(eta - 1. - 0.5*zeta);

                case 3: // phi0 with zeta <- (-zeta)
                  return (1. + zeta)*(xi + eta - 1.)*(xi + eta - 0.5*zeta);

                case 4: // phi1 with zeta <- (-zeta)
                  return (1. + zeta)*xi*(xi - 1. + 0.5*zeta);

                case 5: // phi4 with xi <- eta
                  return (1. + zeta)*eta*(eta - 1. + 0.5*zeta);

                case 6:
                  return 2.*(1. - zeta)*xi*(1. - xi - eta);

                case 7:
                  return 2.*(1. - zeta)*xi*eta;

                case 8:
                  return 2.*(1. - zeta)*eta*(1. - xi - eta);

                case 9:
                  return (1. - zeta)*(1. + zeta)*(1. - xi - eta);

                case 10:
                  return (1. - zeta)*(1. + zeta)*xi;

                case 11: // phi10 with xi <-> eta
                  return (1. - zeta)*(1. + zeta)*eta;

                case 12: // phi6 with zeta <- (-zeta)
                  return 2.*(1. + zeta)*xi*(1. - xi - eta);

                case 13: // phi7 with zeta <- (-zeta)
                  return 2.*(1. + zeta)*xi*eta;

                case 14: // phi8 with zeta <- (-zeta)
                  return 2.*(1. + zeta)*eta*(1. - xi - eta);

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
              auto p1d = p(2);

              //                                0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17
              static const unsigned int i0[] = {0, 0, 0, 1, 1, 1, 0, 0, 0, 2, 2, 2, 1, 1, 1, 2, 2, 2};
              static const unsigned int i1[] = {0, 1, 2, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 3, 4, 5};

              return (FEShim<2,LAGRANGE,RealType>::shape(TRI6,  SECOND, i1[i], p2d)*
                      fe_lagrange_1D_quadratic_shape(i0[i], p1d));
            }

            // G. Bedrosian, "Shape functions and integration formulas for
            // three-dimensional finite element analysis", Int. J. Numerical
            // Methods Engineering, vol 35, p. 95-108, 1992.
          case PYRAMID13:
            {
              libmesh_assert_less (i, 13);

              const auto xi   = p(0);
              const auto eta  = p(1);
              const auto zeta = p(2);
              const Real eps  = 1.e-35;

              // Denominators are perturbed by epsilon to avoid
              // divide-by-zero issues.
              auto den = (1. - zeta + eps);

              switch(i)
                {
                case 0:
                  return 0.25*(-xi - eta - 1.)*((1. - xi)*(1. - eta) - zeta + xi*eta*zeta/den);

                case 1:
                  return 0.25*(-eta + xi - 1.)*((1. + xi)*(1. - eta) - zeta - xi*eta*zeta/den);

                case 2:
                  return 0.25*(xi + eta - 1.)*((1. + xi)*(1. + eta) - zeta + xi*eta*zeta/den);

                case 3:
                  return 0.25*(eta - xi - 1.)*((1. - xi)*(1. + eta) - zeta - xi*eta*zeta/den);

                case 4:
                  return zeta*(2.*zeta - 1.);

                case 5:
                  return 0.5*(1. + xi - zeta)*(1. - xi - zeta)*(1. - eta - zeta)/den;

                case 6:
                  return 0.5*(1. + eta - zeta)*(1. - eta - zeta)*(1. + xi - zeta)/den;

                case 7:
                  return 0.5*(1. + xi - zeta)*(1. - xi - zeta)*(1. + eta - zeta)/den;

                case 8:
                  return 0.5*(1. + eta - zeta)*(1. - eta - zeta)*(1. - xi - zeta)/den;

                case 9:
                  return zeta*(1. - xi - zeta)*(1. - eta - zeta)/den;

                case 10:
                  return zeta*(1. + xi - zeta)*(1. - eta - zeta)/den;

                case 11:
                  return zeta*(1. + eta - zeta)*(1. + xi - zeta)/den;

                case 12:
                  return zeta*(1. - xi - zeta)*(1. + eta - zeta)/den;

                default:
                  libmesh_error_msg("Invalid i = " << i);
                }
            }

            // Quadratic shape functions, as defined in R. Graglia, "Higher order
            // bases on pyramidal elements", IEEE Trans Antennas and Propagation,
            // vol 47, no 5, May 1999.
          case PYRAMID14:
            {
              libmesh_assert_less (i, 14);

              const auto xi   = p(0);
              const auto eta  = p(1);
              const auto zeta = p(2);
              const Real eps  = 1.e-35;

              // The "normalized coordinates" defined by Graglia.  These are
              // the planes which define the faces of the pyramid.
              auto
                p1 = 0.5*(1. - eta - zeta), // back
                p2 = 0.5*(1. + xi  - zeta), // left
                p3 = 0.5*(1. + eta - zeta), // front
                p4 = 0.5*(1. - xi  - zeta); // right

              // Denominators are perturbed by epsilon to avoid
              // divide-by-zero issues.
              auto
                den = (-1. + zeta + eps),
                den2 = den*den;

              switch(i)
                {
                case 0:
                  return p4*p1*(xi*eta - zeta + zeta*zeta)/den2;

                case 1:
                  return -p1*p2*(xi*eta + zeta - zeta*zeta)/den2;

                case 2:
                  return p2*p3*(xi*eta - zeta + zeta*zeta)/den2;

                case 3:
                  return -p3*p4*(xi*eta + zeta - zeta*zeta)/den2;

                case 4:
                  return zeta*(2.*zeta - 1.);

                case 5:
                  return -4.*p2*p1*p4*eta/den2;

                case 6:
                  return 4.*p1*p2*p3*xi/den2;

                case 7:
                  return 4.*p2*p3*p4*eta/den2;

                case 8:
                  return -4.*p3*p4*p1*xi/den2;

                case 9:
                  return -4.*p1*p4*zeta/den;

                case 10:
                  return -4.*p2*p1*zeta/den;

                case 11:
                  return -4.*p3*p2*zeta/den;

                case 12:
                  return -4.*p4*p3*zeta/den;

                case 13:
                  return 16.*p1*p2*p3*p4/den2;

                default:
                  libmesh_error_msg("Invalid i = " << i);
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

#else // LIBMESH_DIM != 3
  libmesh_ignore(type, order, i, p);
  libmesh_not_implemented();
#endif
}


template <typename RealType>
RealType fe_lagrange_3D_shape_deriv(const ElemType type,
                                    const Order order,
                                    const unsigned int i,
                                    const unsigned int j,
                                    const PointTempl<RealType> & p)
{
  typedef PointTempl<RealType> Point;

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
              const auto xi   = p(0);
              const auto eta  = p(1);
              const auto zeta = p(2);

              static const unsigned int i0[] = {0, 1, 1, 0, 0, 1, 1, 0};
              static const unsigned int i1[] = {0, 0, 1, 1, 0, 0, 1, 1};
              static const unsigned int i2[] = {0, 0, 0, 0, 1, 1, 1, 1};

              switch(j)
                {
                case 0:
                  return (fe_lagrange_1D_linear_shape_deriv(i0[i], 0, xi)*
                          fe_lagrange_1D_linear_shape      (i1[i], eta)*
                          fe_lagrange_1D_linear_shape      (i2[i], zeta));

                case 1:
                  return (fe_lagrange_1D_linear_shape      (i0[i], xi)*
                          fe_lagrange_1D_linear_shape_deriv(i1[i], 0, eta)*
                          fe_lagrange_1D_linear_shape      (i2[i], zeta));

                case 2:
                  return (fe_lagrange_1D_linear_shape      (i0[i], xi)*
                          fe_lagrange_1D_linear_shape      (i1[i], eta)*
                          fe_lagrange_1D_linear_shape_deriv(i2[i], 0, zeta));

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
              auto p1d = p(2);

              //                                0  1  2  3  4  5
              static const unsigned int i0[] = {0, 0, 0, 1, 1, 1};
              static const unsigned int i1[] = {0, 1, 2, 0, 1, 2};

              switch (j)
                {
                  // d()/dxi
                case 0:
                  return (FEShim<2,LAGRANGE,RealType>::shape_deriv(TRI3,  FIRST, i1[i], 0, p2d)*
                          fe_lagrange_1D_linear_shape(i0[i], p1d));

                  // d()/deta
                case 1:
                  return (FEShim<2,LAGRANGE,RealType>::shape_deriv(TRI3,  FIRST, i1[i], 1, p2d)*
                          fe_lagrange_1D_linear_shape(i0[i], p1d));

                  // d()/dzeta
                case 2:
                  return (FEShim<2,LAGRANGE,RealType>::shape(TRI3,  FIRST, i1[i], p2d)*
                          fe_lagrange_1D_linear_shape_deriv(i0[i], 0, p1d));

                default:
                  libmesh_error_msg("Invalid shape function derivative j = " << j);
                }
            }

            // linear pyramid shape functions
          case PYRAMID5:
          case PYRAMID13:
          case PYRAMID14:
            {
              libmesh_assert_less (i, 5);

              const auto xi   = p(0);
              const auto eta  = p(1);
              const auto zeta = p(2);
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
                  {
                    // We computed the derivatives with general eps and
                    // then let eps tend to zero in the numerators...
                    auto
                      num = zeta*(2. - zeta) - 1.,
                      den = (1. - zeta + eps)*(1. - zeta + eps);

                    switch(i)
                      {
                      case 0:
                      case 2:
                        return .25*(num + xi*eta)/den;

                      case 1:
                      case 3:
                        return .25*(num - xi*eta)/den;

                      case 4:
                        return 1.;

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
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

              const auto xi   = p(0);
              const auto eta  = p(1);
              const auto zeta = p(2);

              // these functions are defined for (x,y,z) in [0,1]^3
              // so transform the locations
              const auto x = .5*(xi   + 1.);
              const auto y = .5*(eta  + 1.);
              const auto z = .5*(zeta + 1.);

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
              const auto xi   = p(0);
              const auto eta  = p(1);
              const auto zeta = p(2);

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
                  return (fe_lagrange_1D_quadratic_shape_deriv(i0[i], 0, xi)*
                          fe_lagrange_1D_quadratic_shape      (i1[i], eta)*
                          fe_lagrange_1D_quadratic_shape      (i2[i], zeta));

                case 1:
                  return (fe_lagrange_1D_quadratic_shape      (i0[i], xi)*
                          fe_lagrange_1D_quadratic_shape_deriv(i1[i], 0, eta)*
                          fe_lagrange_1D_quadratic_shape      (i2[i], zeta));

                case 2:
                  return (fe_lagrange_1D_quadratic_shape      (i0[i], xi)*
                          fe_lagrange_1D_quadratic_shape      (i1[i], eta)*
                          fe_lagrange_1D_quadratic_shape_deriv(i2[i], 0, zeta));

                default:
                  libmesh_error_msg("Invalid j = " << j);
                }
            }

            // quadratic tetrahedral shape functions
          case TET10:
            {
              libmesh_assert_less (i, 10);

              // Area coordinates, pg. 205, Vol. I, Carey, Oden, Becker FEM
              const auto zeta1 = p(0);
              const auto zeta2 = p(1);
              const auto zeta3 = p(2);
              const auto zeta0 = 1. - zeta1 - zeta2 - zeta3;

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


            // "serendipity" prism
          case PRISM15:
            {
              libmesh_assert_less (i, 15);

              const auto xi   = p(0);
              const auto eta  = p(1);
              const auto zeta = p(2);

              switch (j)
                {
                  // d()/dxi
                case 0:
                  {
                    switch(i)
                      {
                      case 0:
                        return (2.*xi + 2.*eta + 0.5*zeta - 1.)*(1. - zeta);
                      case 1:
                        return (2.*xi - 1. - 0.5*zeta)*(1. - zeta);
                      case 2:
                        return 0.;
                      case 3:
                        return (2.*xi + 2.*eta - 0.5*zeta - 1.)*(1. + zeta);
                      case 4:
                        return (2.*xi - 1. + 0.5*zeta)*(1. + zeta);
                      case 5:
                        return 0.;
                      case 6:
                        return (4.*xi + 2.*eta - 2.)*(zeta - 1.);
                      case 7:
                        return -2.*(zeta - 1.)*eta;
                      case 8:
                        return 2.*(zeta - 1.)*eta;
                      case 9:
                        return (zeta - 1.)*(1. + zeta);
                      case 10:
                        return (1. - zeta)*(1. + zeta);
                      case 11:
                        return 0.;
                      case 12:
                        return (-4.*xi - 2.*eta + 2.)*(1. + zeta);
                      case 13:
                        return 2.*(1. + zeta)*eta;
                      case 14:
                        return -2.*(1. + zeta)*eta;
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
                        return (2.*xi + 2.*eta + 0.5*zeta - 1.)*(1. - zeta);
                      case 1:
                        return 0.;
                      case 2:
                        return (2.*eta - 1. - 0.5*zeta)*(1. - zeta);
                      case 3:
                        return (2.*xi + 2.*eta - 0.5*zeta - 1.)*(1. + zeta);
                      case 4:
                        return 0.;
                      case 5:
                        return (2.*eta - 1. + 0.5*zeta)*(1. + zeta);
                      case 6:
                        return 2.*(zeta - 1.)*xi;
                      case 7:
                        return 2.*(1. - zeta)*xi;
                      case 8:
                        return (2.*xi + 4.*eta - 2.)*(zeta - 1.);
                      case 9:
                        return (zeta - 1.)*(1. + zeta);
                      case 10:
                        return 0.;
                      case 11:
                        return (1. - zeta)*(1. + zeta);
                      case 12:
                        return -2.*(1. + zeta)*xi;
                      case 13:
                        return 2.*(1. + zeta)*xi;
                      case 14:
                        return (-2.*xi - 4.*eta + 2.)*(1. + zeta);

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
                        return (-xi - eta - zeta + 0.5)*(xi + eta - 1.);
                      case 1:
                        return -0.5*xi*(2.*xi - 1. - 2.*zeta);
                      case 2:
                        return -0.5*eta*(2.*eta - 1. - 2.*zeta);
                      case 3:
                        return (xi + eta - zeta - 0.5)*(xi + eta - 1.);
                      case 4:
                        return 0.5*xi*(2.*xi - 1. + 2.*zeta);
                      case 5:
                        return 0.5*eta*(2.*eta - 1. + 2.*zeta);
                      case 6:
                        return 2.*xi*(xi + eta - 1.);
                      case 7:
                        return -2.*xi*eta;
                      case 8:
                        return 2.*eta*(xi + eta - 1.);
                      case 9:
                        return 2.*zeta*(xi + eta - 1.);
                      case 10:
                        return -2.*xi*zeta;
                      case 11:
                        return -2.*eta*zeta;
                      case 12:
                        return 2.*xi*(1. - xi - eta);
                      case 13:
                        return 2.*xi*eta;
                      case 14:
                        return 2.*eta*(1. - xi - eta);

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
              auto p1d = p(2);

              //                                0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17
              static const unsigned int i0[] = {0, 0, 0, 1, 1, 1, 0, 0, 0, 2, 2, 2, 1, 1, 1, 2, 2, 2};
              static const unsigned int i1[] = {0, 1, 2, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 3, 4, 5};

              switch (j)
                {
                  // d()/dxi
                case 0:
                  return (FEShim<2,LAGRANGE,RealType>::shape_deriv(TRI6,  SECOND, i1[i], 0, p2d)*
                          fe_lagrange_1D_quadratic_shape(i0[i], p1d));

                  // d()/deta
                case 1:
                  return (FEShim<2,LAGRANGE,RealType>::shape_deriv(TRI6,  SECOND, i1[i], 1, p2d)*
                          fe_lagrange_1D_quadratic_shape(i0[i], p1d));

                  // d()/dzeta
                case 2:
                  return (FEShim<2,LAGRANGE,RealType>::shape(TRI6,  SECOND, i1[i], p2d)*
                          fe_lagrange_1D_quadratic_shape_deriv(i0[i], 0, p1d));

                default:
                  libmesh_error_msg("Invalid shape function derivative j = " << j);
                }
            }

            // G. Bedrosian, "Shape functions and integration formulas for
            // three-dimensional finite element analysis", Int. J. Numerical
            // Methods Engineering, vol 35, p. 95-108, 1992.
          case PYRAMID13:
            {
              libmesh_assert_less (i, 13);

              const auto xi   = p(0);
              const auto eta  = p(1);
              const auto zeta = p(2);
              const Real eps  = 1.e-35;

              // Denominators are perturbed by epsilon to avoid
              // divide-by-zero issues.
              RealType
                den = (-1. + zeta + eps),
                den2 = den*den,
                xi2 = xi*xi,
                eta2 = eta*eta,
                zeta2 = zeta*zeta,
                zeta3 = zeta2*zeta;

              switch (j)
                {
                  // d/dxi
                case 0:
                  switch(i)
                    {
                    case 0:
                      return 0.25*(-zeta - eta + 2.*eta*zeta - 2.*xi + 2.*zeta*xi + 2.*eta*xi + zeta2 + eta2)/den;

                    case 1:
                      return -0.25*(-zeta - eta + 2.*eta*zeta + 2.*xi - 2.*zeta*xi - 2.*eta*xi + zeta2 + eta2)/den;

                    case 2:
                      return -0.25*(-zeta + eta - 2.*eta*zeta + 2.*xi - 2.*zeta*xi + 2.*eta*xi + zeta2 + eta2)/den;

                    case 3:
                      return 0.25*(-zeta + eta - 2.*eta*zeta - 2.*xi + 2.*zeta*xi - 2.*eta*xi + zeta2 + eta2)/den;

                    case 4:
                      return 0.;

                    case 5:
                      return -(-1. + eta + zeta)*xi/den;

                    case 6:
                      return 0.5*(-1. + eta + zeta)*(1. + eta - zeta)/den;

                    case 7:
                      return (1. + eta - zeta)*xi/den;

                    case 8:
                      return -0.5*(-1. + eta + zeta)*(1. + eta - zeta)/den;

                    case 9:
                      return -(-1. + eta + zeta)*zeta/den;

                    case 10:
                      return (-1. + eta + zeta)*zeta/den;

                    case 11:
                      return -(1. + eta - zeta)*zeta/den;

                    case 12:
                      return (1. + eta - zeta)*zeta/den;

                    default:
                      libmesh_error_msg("Invalid i = " << i);
                    }

                  // d/deta
                case 1:
                  switch(i)
                    {
                    case 0:
                      return 0.25*(-zeta - 2.*eta + 2.*eta*zeta - xi + 2.*zeta*xi + 2.*eta*xi + zeta2 + xi2)/den;

                    case 1:
                      return -0.25*(zeta + 2.*eta - 2.*eta*zeta - xi + 2.*zeta*xi + 2.*eta*xi - zeta2 - xi2)/den;

                    case 2:
                      return -0.25*(-zeta + 2.*eta - 2.*eta*zeta + xi - 2.*zeta*xi + 2.*eta*xi + zeta2 + xi2)/den;

                    case 3:
                      return 0.25*(zeta - 2.*eta + 2.*eta*zeta + xi - 2.*zeta*xi + 2.*eta*xi - zeta2 - xi2)/den;

                    case 4:
                      return 0.;

                    case 5:
                      return -0.5*(-1. + xi + zeta)*(1. + xi - zeta)/den;

                    case 6:
                      return (1. + xi - zeta)*eta/den;

                    case 7:
                      return 0.5*(-1. + xi + zeta)*(1. + xi - zeta)/den;

                    case 8:
                      return -(-1. + xi + zeta)*eta/den;

                    case 9:
                      return -(-1. + xi + zeta)*zeta/den;

                    case 10:
                      return (1. + xi - zeta)*zeta/den;

                    case 11:
                      return -(1. + xi - zeta)*zeta/den;

                    case 12:
                      return (-1. + xi + zeta)*zeta/den;

                    default:
                      libmesh_error_msg("Invalid i = " << i);
                    }

                  // d/dzeta
                case 2:
                  {
                    switch(i)
                      {
                      case 0:
                        return -0.25*(xi + eta + 1.)*(-1. + 2.*zeta - zeta2 + eta*xi)/den2;

                      case 1:
                        return 0.25*(eta - xi + 1.)*(1. - 2.*zeta + zeta2 + eta*xi)/den2;

                      case 2:
                        return 0.25*(xi + eta - 1.)*(-1. + 2.*zeta - zeta2 + eta*xi)/den2;

                      case 3:
                        return -0.25*(eta - xi - 1.)*(1. - 2.*zeta + zeta2 + eta*xi)/den2;

                      case 4:
                        return 4.*zeta - 1.;

                      case 5:
                        return 0.5*(-2 + eta + 6.*zeta + eta*xi2 + eta*zeta2 - 6.*zeta2 + 2.*zeta3 - 2.*eta*zeta)/den2;

                      case 6:
                        return -0.5*(2 - 6.*zeta + xi + xi*zeta2 + eta2*xi + 6.*zeta2 - 2.*zeta3 - 2.*zeta*xi)/den2;

                      case 7:
                        return -0.5*(2 + eta - 6.*zeta + eta*xi2 + eta*zeta2 + 6.*zeta2 - 2.*zeta3 - 2.*eta*zeta)/den2;

                      case 8:
                        return 0.5*(-2 + 6.*zeta + xi + xi*zeta2 + eta2*xi - 6.*zeta2 + 2.*zeta3 - 2.*zeta*xi)/den2;

                      case 9:
                        return (1. - eta - 4.*zeta - xi - xi*zeta2 - eta*zeta2 + eta*xi + 5.*zeta2 - 2.*zeta3 + 2.*eta*zeta + 2.*zeta*xi)/den2;

                      case 10:
                        return -(-1. + eta + 4.*zeta - xi - xi*zeta2 + eta*zeta2 + eta*xi - 5.*zeta2 + 2.*zeta3 - 2.*eta*zeta + 2.*zeta*xi)/den2;

                      case 11:
                        return (1. + eta - 4.*zeta + xi + xi*zeta2 + eta*zeta2 + eta*xi + 5.*zeta2 - 2.*zeta3 - 2.*eta*zeta - 2.*zeta*xi)/den2;

                      case 12:
                        return -(-1. - eta + 4.*zeta + xi + xi*zeta2 - eta*zeta2 + eta*xi - 5.*zeta2 + 2.*zeta3 + 2.*eta*zeta - 2.*zeta*xi)/den2;

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  }

                default:
                  libmesh_error_msg("Invalid j = " << j);
                }
            }

            // Quadratic shape functions, as defined in R. Graglia, "Higher order
            // bases on pyramidal elements", IEEE Trans Antennas and Propagation,
            // vol 47, no 5, May 1999.
          case PYRAMID14:
            {
              libmesh_assert_less (i, 14);

              const auto xi   = p(0);
              const auto eta  = p(1);
              const auto zeta = p(2);
              const Real eps  = 1.e-35;

              // The "normalized coordinates" defined by Graglia.  These are
              // the planes which define the faces of the pyramid.
              RealType
                p1 = 0.5*(1. - eta - zeta), // back
                p2 = 0.5*(1. + xi  - zeta), // left
                p3 = 0.5*(1. + eta - zeta), // front
                p4 = 0.5*(1. - xi  - zeta); // right

              // Denominators are perturbed by epsilon to avoid
              // divide-by-zero issues.
              RealType
                den = (-1. + zeta + eps),
                den2 = den*den,
                den3 = den2*den;

              switch (j)
                {
                  // d/dxi
                case 0:
                  switch(i)
                    {
                    case 0:
                      return 0.5*p1*(-xi*eta + zeta - zeta*zeta + 2.*p4*eta)/den2;

                    case 1:
                      return -0.5*p1*(xi*eta + zeta - zeta*zeta + 2.*p2*eta)/den2;

                    case 2:
                      return 0.5*p3*(xi*eta - zeta + zeta*zeta + 2.*p2*eta)/den2;

                    case 3:
                      return -0.5*p3*(-xi*eta - zeta + zeta*zeta + 2.*p4*eta)/den2;

                    case 4:
                      return 0.;

                    case 5:
                      return 2.*p1*eta*xi/den2;

                    case 6:
                      return 2.*p1*p3*(xi + 2.*p2)/den2;

                    case 7:
                      return -2.*p3*eta*xi/den2;

                    case 8:
                      return -2.*p1*p3*(-xi + 2.*p4)/den2;

                    case 9:
                      return 2.*p1*zeta/den;

                    case 10:
                      return -2.*p1*zeta/den;

                    case 11:
                      return -2.*p3*zeta/den;

                    case 12:
                      return 2.*p3*zeta/den;

                    case 13:
                      return -8.*p1*p3*xi/den2;

                    default:
                      libmesh_error_msg("Invalid i = " << i);
                    }

                  // d/deta
                case 1:
                  switch(i)
                    {
                    case 0:
                      return -0.5*p4*(xi*eta - zeta + zeta*zeta - 2.*p1*xi)/den2;

                    case 1:
                      return 0.5*p2*(xi*eta + zeta - zeta*zeta - 2.*p1*xi)/den2;

                    case 2:
                      return 0.5*p2*(xi*eta - zeta + zeta*zeta + 2.*p3*xi)/den2;

                    case 3:
                      return -0.5*p4*(xi*eta + zeta - zeta*zeta + 2.*p3*xi)/den2;

                    case 4:
                      return 0.;

                    case 5:
                      return 2.*p2*p4*(eta - 2.*p1)/den2;

                    case 6:
                      return -2.*p2*xi*eta/den2;

                    case 7:
                      return 2.*p2*p4*(eta + 2.*p3)/den2;

                    case 8:
                      return 2.*p4*xi*eta/den2;

                    case 9:
                      return 2.*p4*zeta/den;

                    case 10:
                      return 2.*p2*zeta/den;

                    case 11:
                      return -2.*p2*zeta/den;

                    case 12:
                      return -2.*p4*zeta/den;

                    case 13:
                      return -8.*p2*p4*eta/den2;

                    default:
                      libmesh_error_msg("Invalid i = " << i);
                    }


                  // d/dzeta
                case 2:
                  {
                    switch(i)
                      {
                      case 0:
                        return -0.5*p1*(xi*eta - zeta + zeta*zeta)/den2
                          - 0.5*p4*(xi*eta - zeta + zeta*zeta)/den2
                          + p4*p1*(2.*zeta - 1)/den2
                          - 2.*p4*p1*(xi*eta - zeta + zeta*zeta)/den3;

                      case 1:
                        return 0.5*p2*(xi*eta + zeta - zeta*zeta)/den2
                          + 0.5*p1*(xi*eta + zeta - zeta*zeta)/den2
                          - p1*p2*(1 - 2.*zeta)/den2
                          + 2.*p1*p2*(xi*eta + zeta - zeta*zeta)/den3;

                      case 2:
                        return -0.5*p3*(xi*eta - zeta + zeta*zeta)/den2
                          - 0.5*p2*(xi*eta - zeta + zeta*zeta)/den2
                          + p2*p3*(2.*zeta - 1)/den2
                          - 2.*p2*p3*(xi*eta - zeta + zeta*zeta)/den3;

                      case 3:
                        return 0.5*p4*(xi*eta + zeta - zeta*zeta)/den2
                          + 0.5*p3*(xi*eta + zeta - zeta*zeta)/den2
                          - p3*p4*(1 - 2.*zeta)/den2
                          + 2.*p3*p4*(xi*eta + zeta - zeta*zeta)/den3;

                      case 4:
                        return 4.*zeta - 1.;

                      case 5:
                        return 2.*p4*p1*eta/den2
                          + 2.*p2*p4*eta/den2
                          + 2.*p1*p2*eta/den2
                          + 8.*p2*p1*p4*eta/den3;

                      case 6:
                        return -2.*p2*p3*xi/den2
                          - 2.*p1*p3*xi/den2
                          - 2.*p1*p2*xi/den2
                          - 8.*p1*p2*p3*xi/den3;

                      case 7:
                        return -2.*p3*p4*eta/den2
                          - 2.*p2*p4*eta/den2
                          - 2.*p2*p3*eta/den2
                          - 8.*p2*p3*p4*eta/den3;

                      case 8:
                        return 2.*p4*p1*xi/den2
                          + 2.*p1*p3*xi/den2
                          + 2.*p3*p4*xi/den2
                          + 8.*p3*p4*p1*xi/den3;

                      case 9:
                        return 2.*p4*zeta/den
                          + 2.*p1*zeta/den
                          - 4.*p1*p4/den
                          + 4.*p1*p4*zeta/den2;

                      case 10:
                        return 2.*p1*zeta/den
                          + 2.*p2*zeta/den
                          - 4.*p2*p1/den
                          + 4.*p2*p1*zeta/den2;

                      case 11:
                        return 2.*p2*zeta/den
                          + 2.*p3*zeta/den
                          - 4.*p3*p2/den
                          + 4.*p3*p2*zeta/den2;

                      case 12:
                        return 2.*p3*zeta/den
                          + 2.*p4*zeta/den
                          - 4.*p4*p3/den
                          + 4.*p4*p3*zeta/den2;

                      case 13:
                        return -8.*p2*p3*p4/den2
                          - 8.*p3*p4*p1/den2
                          - 8.*p2*p1*p4/den2
                          - 8.*p1*p2*p3/den2
                          - 32.*p1*p2*p3*p4/den3;

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  }

                default:
                  libmesh_error_msg("Invalid j = " << j);
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

#else // LIBMESH_DIM != 3
  libmesh_ignore(type, order, i, j, p);
  libmesh_not_implemented();
#endif
}



#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

template <typename RealType>
RealType fe_lagrange_3D_shape_second_deriv(const ElemType type,
                                           const Order order,
                                           const unsigned int i,
                                           const unsigned int j,
                                           const PointTempl<RealType> & p)
{
  typedef PointTempl<RealType> Point;

#if LIBMESH_DIM == 3

  libmesh_assert_less (j, 6);

  switch (order)
    {
      // linear Lagrange shape functions
    case FIRST:
      {
        switch (type)
          {
            // Linear tets have all second derivatives = 0
          case TET4:
          case TET10:
            {
              return 0.;
            }

            // The following elements use either tensor product or
            // rational basis functions, and therefore probably have
            // second derivatives, but we have not implemented them
            // yet...
          case PRISM6:
          case PRISM15:
          case PRISM18:
            {
              libmesh_assert_less (i, 6);

              // Compute prism shape functions as a tensor-product
              // of a triangle and an edge

              Point p2d(p(0),p(1));
              auto p1d = p(2);

              //                                0  1  2  3  4  5
              static const unsigned int i0[] = {0, 0, 0, 1, 1, 1};
              static const unsigned int i1[] = {0, 1, 2, 0, 1, 2};

              switch (j)
                {
                  // All repeated second derivatives and the xi-eta derivative are zero on PRISMs
                case 0: // d^2()/dxi^2
                case 1: // d^2()/dxideta
                case 2: // d^2()/deta^2
                case 5: // d^2()/dzeta^2
                  {
                    return 0.;
                  }

                case 3: // d^2()/dxidzeta
                  return (FEShim<2,LAGRANGE,RealType>::shape_deriv(TRI3,  FIRST, i1[i], 0, p2d)*
                          fe_lagrange_1D_linear_shape_deriv(i0[i], 0, p1d));

                case 4: // d^2()/detadzeta
                  return (FEShim<2,LAGRANGE,RealType>::shape_deriv(TRI3,  FIRST, i1[i], 1, p2d)*
                          fe_lagrange_1D_linear_shape_deriv(i0[i], 0, p1d));

                default:
                  libmesh_error_msg("Invalid j = " << j);
                }
            }

          case PYRAMID5:
          case PYRAMID13:
          case PYRAMID14:
            {
              libmesh_assert_less (i, 5);

              const auto xi   = p(0);
              const auto eta  = p(1);
              const auto zeta = p(2);
              const Real eps  = 1.e-35;

              switch (j)
                {
                  // xi-xi and eta-eta derivatives are all zero for PYRAMID5.
                case 0: // d^2()/dxi^2
                case 2: // d^2()/deta^2
                  return 0.;

                case 1: // d^2()/dxideta
                  {
                    switch (i)
                      {
                      case 0:
                      case 2:
                        return 0.25/(1. - zeta + eps);
                      case 1:
                      case 3:
                        return -0.25/(1. - zeta + eps);
                      case 4:
                        return 0.;
                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  }

                case 3: // d^2()/dxidzeta
                  {
                    auto den = (1. - zeta + eps)*(1. - zeta + eps);

                    switch (i)
                      {
                      case 0:
                      case 2:
                        return 0.25*eta/den;
                      case 1:
                      case 3:
                        return -0.25*eta/den;
                      case 4:
                        return 0.;
                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  }

                case 4: // d^2()/detadzeta
                  {
                    auto den = (1. - zeta + eps)*(1. - zeta + eps);

                    switch (i)
                      {
                      case 0:
                      case 2:
                        return 0.25*xi/den;
                      case 1:
                      case 3:
                        return -0.25*xi/den;
                      case 4:
                        return 0.;
                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  }

                case 5: // d^2()/dzeta^2
                  {
                    auto den = (1. - zeta + eps)*(1. - zeta + eps)*(1. - zeta + eps);

                    switch (i)
                      {
                      case 0:
                      case 2:
                        return 0.5*xi*eta/den;
                      case 1:
                      case 3:
                        return -0.5*xi*eta/den;
                      case 4:
                        return 0.;
                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  }

                default:
                  libmesh_error_msg("Invalid j = " << j);
                }
            }

            // Trilinear shape functions on HEX8s have nonzero mixed second derivatives
          case HEX8:
          case HEX20:
          case HEX27:
            {
              libmesh_assert_less (i, 8);

              // Compute hex shape functions as a tensor-product
              const auto xi   = p(0);
              const auto eta  = p(1);
              const auto zeta = p(2);

              static const unsigned int i0[] = {0, 1, 1, 0, 0, 1, 1, 0};
              static const unsigned int i1[] = {0, 0, 1, 1, 0, 0, 1, 1};
              static const unsigned int i2[] = {0, 0, 0, 0, 1, 1, 1, 1};

              switch (j)
                {
                  // All repeated second derivatives are zero on HEX8
                case 0: // d^2()/dxi^2
                case 2: // d^2()/deta^2
                case 5: // d^2()/dzeta^2
                  {
                    return 0.;
                  }

                case 1: // d^2()/dxideta
                  return (fe_lagrange_1D_linear_shape_deriv(i0[i], 0, xi)*
                          fe_lagrange_1D_linear_shape_deriv(i1[i], 0, eta)*
                          fe_lagrange_1D_linear_shape      (i2[i], zeta));

                case 3: // d^2()/dxidzeta
                  return (fe_lagrange_1D_linear_shape_deriv(i0[i], 0, xi)*
                          fe_lagrange_1D_linear_shape      (i1[i], eta)*
                          fe_lagrange_1D_linear_shape_deriv(i2[i], 0, zeta));

                case 4: // d^2()/detadzeta
                  return (fe_lagrange_1D_linear_shape      (i0[i], xi)*
                          fe_lagrange_1D_linear_shape_deriv(i1[i], 0, eta)*
                          fe_lagrange_1D_linear_shape_deriv(i2[i], 0, zeta));

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

              const auto xi   = p(0);
              const auto eta  = p(1);
              const auto zeta = p(2);

              // these functions are defined for (x,y,z) in [0,1]^3
              // so transform the locations
              const auto x = .5*(xi   + 1.);
              const auto y = .5*(eta  + 1.);
              const auto z = .5*(zeta + 1.);

              switch(j)
                {
                case 0: // d^2()/dxi^2
                  {
                    switch(i)
                      {
                      case 0:
                      case 1:
                        return (1. - y) * (1. - z);
                      case 2:
                      case 3:
                        return y * (1. - z);
                      case 4:
                      case 5:
                        return (1. - y) * z;
                      case 6:
                      case 7:
                        return y * z;
                      case 8:
                        return -2. * (1. - y) * (1. - z);
                      case 10:
                        return -2. * y * (1. - z);
                      case 16:
                        return -2. * (1. - y) * z;
                      case 18:
                        return -2. * y * z;
                      case 9:
                      case 11:
                      case 12:
                      case 13:
                      case 14:
                      case 15:
                      case 17:
                      case 19:
                        return 0;
                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  }
                case 1: // d^2()/dxideta
                  {
                    switch(i)
                      {
                      case 0:
                        return (1.25 - x - y - .5*z) * (1. - z);
                      case 1:
                        return (-x + y + .5*z - .25) * (1. - z);
                      case 2:
                        return (x + y - .5*z - .75) * (1. - z);
                      case 3:
                        return (-y + x + .5*z - .25) * (1. - z);
                      case 4:
                        return -.25*z * (4.*x + 4.*y - 2.*z - 3);
                      case 5:
                        return -.25*z * (-4.*y + 4.*x + 2.*z - 1.);
                      case 6:
                        return .25*z * (-5 + 4.*x + 4.*y + 2.*z);
                      case 7:
                        return .25*z * (4.*x - 4.*y - 2.*z + 1.);
                      case 8:
                        return (-1. + 2.*x) * (1. - z);
                      case 9:
                        return (1. - 2.*y) * (1. - z);
                      case 10:
                        return (1. - 2.*x) * (1. - z);
                      case 11:
                        return (-1. + 2.*y) * (1. - z);
                      case 12:
                        return z * (1. - z);
                      case 13:
                        return -z * (1. - z);
                      case 14:
                        return z * (1. - z);
                      case 15:
                        return -z * (1. - z);
                      case 16:
                        return (-1. + 2.*x) * z;
                      case 17:
                        return (1. - 2.*y) * z;
                      case 18:
                        return (1. - 2.*x) * z;
                      case 19:
                        return (-1. + 2.*y) * z;
                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  }
                case 2: // d^2()/deta^2
                  switch(i)
                    {
                    case 0:
                    case 3:
                      return (1. - x) * (1. - z);
                    case 1:
                    case 2:
                      return x * (1. - z);
                    case 4:
                    case 7:
                      return (1. - x) * z;
                    case 5:
                    case 6:
                      return x * z;
                    case 9:
                      return -2. * x * (1. - z);
                    case 11:
                      return -2. * (1. - x) * (1. - z);
                    case 17:
                      return -2. * x * z;
                    case 19:
                      return -2. * (1. - x) * z;
                    case 8:
                    case 10:
                    case 12:
                    case 13:
                    case 14:
                    case 15:
                    case 16:
                    case 18:
                      return 0.;
                    default:
                      libmesh_error_msg("Invalid i = " << i);
                    }
                case 3: // d^2()/dxidzeta
                  switch(i)
                    {
                    case 0:
                      return (1.25 - x - .5*y - z) * (1. - y);
                    case 1:
                      return (-x + .5*y + z - .25) * (1. - y);
                    case 2:
                      return -.25*y * (2.*y + 4.*x - 4.*z - 1.);
                    case 3:
                      return -.25*y * (-2.*y + 4.*x + 4.*z - 3);
                    case 4:
                      return (-z + x + .5*y - .25) * (1. - y);
                    case 5:
                      return (x - .5*y + z - .75) * (1. - y);
                    case 6:
                      return .25*y * (2.*y + 4.*x + 4.*z - 5);
                    case 7:
                      return .25*y * (-2.*y + 4.*x - 4.*z + 1.);
                    case 8:
                      return (-1. + 2.*x) * (1. - y);
                    case 9:
                      return -y * (1. - y);
                    case 10:
                      return (-1. + 2.*x) * y;
                    case 11:
                      return y * (1. - y);
                    case 12:
                      return (-1. + 2.*z) * (1. - y);
                    case 13:
                      return (1. - 2.*z) * (1. - y);
                    case 14:
                      return (1. - 2.*z) * y;
                    case 15:
                      return (-1. + 2.*z) * y;
                    case 16:
                      return (1. - 2.*x) * (1. - y);
                    case 17:
                      return y * (1. - y);
                    case 18:
                      return (1. - 2.*x) * y;
                    case 19:
                      return -y * (1. - y);
                    default:
                      libmesh_error_msg("Invalid i = " << i);
                    }
                case 4: // d^2()/detadzeta
                  switch(i)
                    {
                    case 0:
                      return (1.25 - .5*x - y - z) * (1. - x);
                    case 1:
                      return .25*x * (2.*x - 4.*y - 4.*z + 3.);
                    case 2:
                      return -.25*x * (2.*x + 4.*y - 4.*z - 1.);
                    case 3:
                      return (-y + .5*x + z - .25) * (1. - x);
                    case 4:
                      return (-z + .5*x + y - .25) * (1. - x);
                    case 5:
                      return -.25*x * (2.*x - 4.*y + 4.*z - 1.);
                    case 6:
                      return .25*x * (2.*x + 4.*y + 4.*z - 5);
                    case 7:
                      return (y - .5*x + z - .75) * (1. - x);
                    case 8:
                      return x * (1. - x);
                    case 9:
                      return (-1. + 2.*y) * x;
                    case 10:
                      return -x * (1. - x);
                    case 11:
                      return (-1. + 2.*y) * (1. - x);
                    case 12:
                      return (-1. + 2.*z) * (1. - x);
                    case 13:
                      return (-1. + 2.*z) * x;
                    case 14:
                      return (1. - 2.*z) * x;
                    case 15:
                      return (1. - 2.*z) * (1. - x);
                    case 16:
                      return -x * (1. - x);
                    case 17:
                      return (1. - 2.*y) * x;
                    case 18:
                      return x * (1. - x);
                    case 19:
                      return (1. - 2.*y) * (1. - x);
                    default:
                      libmesh_error_msg("Invalid i = " << i);
                    }
                case 5: // d^2()/dzeta^2
                  switch(i)
                    {
                    case 0:
                    case 4:
                      return (1. - x) * (1. - y);
                    case 1:
                    case 5:
                      return x * (1. - y);
                    case 2:
                    case 6:
                      return x * y;
                    case 3:
                    case 7:
                      return (1. - x) * y;
                    case 12:
                      return -2. * (1. - x) * (1. - y);
                    case 13:
                      return -2. * x * (1. - y);
                    case 14:
                      return -2. * x * y;
                    case 15:
                      return -2. * (1. - x) * y;
                    case 8:
                    case 9:
                    case 10:
                    case 11:
                    case 16:
                    case 17:
                    case 18:
                    case 19:
                      return 0.;
                    default:
                      libmesh_error_msg("Invalid i = " << i);
                    }
                default:
                  libmesh_error_msg("Invalid j = " << j);
                }
            }

            // triquadratic hexahedral shape functions
          case HEX27:
            {
              libmesh_assert_less (i, 27);

              // Compute hex shape functions as a tensor-product
              const auto xi   = p(0);
              const auto eta  = p(1);
              const auto zeta = p(2);

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
                  return (fe_lagrange_1D_quadratic_shape_second_deriv(i0[i], 0, xi)*
                          fe_lagrange_1D_quadratic_shape             (i1[i], eta)*
                          fe_lagrange_1D_quadratic_shape             (i2[i], zeta));

                  // d^2()/dxideta
                case 1:
                  return (fe_lagrange_1D_quadratic_shape_deriv(i0[i], 0, xi)*
                          fe_lagrange_1D_quadratic_shape_deriv(i1[i], 0, eta)*
                          fe_lagrange_1D_quadratic_shape      (i2[i], zeta));

                  // d^2()/deta^2
                case 2:
                  return (fe_lagrange_1D_quadratic_shape             (i0[i], xi)*
                          fe_lagrange_1D_quadratic_shape_second_deriv(i1[i], 0, eta)*
                          fe_lagrange_1D_quadratic_shape             (i2[i], zeta));

                  // d^2()/dxidzeta
                case 3:
                  return (fe_lagrange_1D_quadratic_shape_deriv(i0[i], 0, xi)*
                          fe_lagrange_1D_quadratic_shape      (i1[i], eta)*
                          fe_lagrange_1D_quadratic_shape_deriv(i2[i], 0, zeta));

                  // d^2()/detadzeta
                case 4:
                  return (fe_lagrange_1D_quadratic_shape      (i0[i], xi)*
                          fe_lagrange_1D_quadratic_shape_deriv(i1[i], 0, eta)*
                          fe_lagrange_1D_quadratic_shape_deriv(i2[i], 0, zeta));

                  // d^2()/dzeta^2
                case 5:
                  return (fe_lagrange_1D_quadratic_shape             (i0[i], xi)*
                          fe_lagrange_1D_quadratic_shape             (i1[i], eta)*
                          fe_lagrange_1D_quadratic_shape_second_deriv(i2[i], 0, zeta));

                default:
                  libmesh_error_msg("Invalid j = " << j);
                }
            }

            // quadratic tetrahedral shape functions
          case TET10:
            {
              // The area coordinates are the same as used for the
              // shape() and shape_deriv() functions.
              // const auto zeta0 = 1. - zeta1 - zeta2 - zeta3;
              // const auto zeta1 = p(0);
              // const auto zeta2 = p(1);
              // const auto zeta3 = p(2);
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



            // "serendipity" prism
          case PRISM15:
            {
              libmesh_assert_less (i, 15);

              const auto xi   = p(0);
              const auto eta  = p(1);
              const auto zeta = p(2);

              switch (j)
                {
                  // d^2()/dxi^2
                case 0:
                  {
                    switch(i)
                      {
                      case 0:
                      case 1:
                        return 2.*(1. - zeta);
                      case 2:
                      case 5:
                      case 7:
                      case 8:
                      case 9:
                      case 10:
                      case 11:
                      case 13:
                      case 14:
                        return 0.;
                      case 3:
                      case 4:
                        return 2.*(1. + zeta);
                      case 6:
                        return 4.*(zeta - 1);
                      case 12:
                        return -4.*(1. + zeta);
                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  }

                  // d^2()/dxideta
                case 1:
                  {
                    switch(i)
                      {
                      case 0:
                      case 7:
                        return 2.*(1. - zeta);
                      case 1:
                      case 2:
                      case 4:
                      case 5:
                      case 9:
                      case 10:
                      case 11:
                        return 0.;
                      case 3:
                      case 13:
                        return 2.*(1. + zeta);
                      case 6:
                      case 8:
                        return 2.*(zeta - 1.);
                      case 12:
                      case 14:
                        return -2.*(1. + zeta);
                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  }

                  // d^2()/deta^2
                case 2:
                  {
                    switch(i)
                      {
                      case 0:
                      case 2:
                        return 2.*(1. - zeta);
                      case 1:
                      case 4:
                      case 6:
                      case 7:
                      case 9:
                      case 10:
                      case 11:
                      case 12:
                      case 13:
                        return 0.;
                      case 3:
                      case 5:
                        return 2.*(1. + zeta);
                      case 8:
                        return 4.*(zeta - 1.);
                      case 14:
                        return -4.*(1. + zeta);
                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  }

                  // d^2()/dxidzeta
                case 3:
                  {
                    switch(i)
                      {
                      case 0:
                        return 1.5 - zeta - 2.*xi - 2.*eta;
                      case 1:
                        return 0.5 + zeta - 2.*xi;
                      case 2:
                      case 5:
                      case 11:
                        return 0.;
                      case 3:
                        return -1.5 - zeta + 2.*xi + 2.*eta;
                      case 4:
                        return -0.5 + zeta + 2.*xi;
                      case 6:
                        return 4.*xi + 2.*eta - 2.;
                      case 7:
                        return -2.*eta;
                      case 8:
                        return 2.*eta;
                      case 9:
                        return 2.*zeta;
                      case 10:
                        return -2.*zeta;
                      case 12:
                        return -4.*xi - 2.*eta + 2.;
                      case 13:
                        return 2.*eta;
                      case 14:
                        return -2.*eta;
                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  }

                  // d^2()/detadzeta
                case 4:
                  {
                    switch(i)
                      {
                      case 0:
                        return 1.5 - zeta - 2.*xi - 2.*eta;
                      case 1:
                      case 4:
                      case 10:
                        return 0.;
                      case 2:
                        return .5 + zeta - 2.*eta;
                      case 3:
                        return -1.5 - zeta + 2.*xi + 2.*eta;
                      case 5:
                        return -.5 + zeta + 2.*eta;
                      case 6:
                        return 2.*xi;
                      case 7:
                        return -2.*xi;
                      case 8:
                        return 2.*xi + 4.*eta - 2.;
                      case 9:
                        return 2.*zeta;
                      case 11:
                        return -2.*zeta;
                      case 12:
                        return -2.*xi;
                      case 13:
                        return 2.*xi;
                      case 14:
                        return -2.*xi - 4.*eta + 2.;
                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  }

                  // d^2()/dzeta^2
                case 5:
                  {
                    switch(i)
                      {
                      case 0:
                      case 3:
                        return 1. - xi - eta;
                      case 1:
                      case 4:
                        return xi;
                      case 2:
                      case 5:
                        return eta;
                      case 6:
                      case 7:
                      case 8:
                      case 12:
                      case 13:
                      case 14:
                        return 0.;
                      case 9:
                        return 2.*xi + 2.*eta - 2.;
                      case 10:
                        return -2.*xi;
                      case 11:
                        return -2.*eta;
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
              auto p1d = p(2);

              //                                0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17
              static const unsigned int i0[] = {0, 0, 0, 1, 1, 1, 0, 0, 0, 2, 2, 2, 1, 1, 1, 2, 2, 2};
              static const unsigned int i1[] = {0, 1, 2, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 3, 4, 5};

              switch (j)
                {
                  // d^2()/dxi^2
                case 0:
                  return (FEShim<2,LAGRANGE,RealType>::shape_second_deriv(TRI6, SECOND, i1[i], 0, p2d)*
                          fe_lagrange_1D_quadratic_shape(i0[i], p1d));

                  // d^2()/dxideta
                case 1:
                  return (FEShim<2,LAGRANGE,RealType>::shape_second_deriv(TRI6, SECOND, i1[i], 1, p2d)*
                          fe_lagrange_1D_quadratic_shape(i0[i], p1d));

                  // d^2()/deta^2
                case 2:
                  return (FEShim<2,LAGRANGE,RealType>::shape_second_deriv(TRI6, SECOND, i1[i], 2, p2d)*
                          fe_lagrange_1D_quadratic_shape(i0[i], p1d));

                  // d^2()/dxidzeta
                case 3:
                  return (FEShim<2,LAGRANGE,RealType>::shape_deriv(TRI6,  SECOND, i1[i], 0, p2d)*
                          fe_lagrange_1D_quadratic_shape_deriv(i0[i], 0, p1d));

                  // d^2()/detadzeta
                case 4:
                  return (FEShim<2,LAGRANGE,RealType>::shape_deriv(TRI6,  SECOND, i1[i], 1, p2d)*
                          fe_lagrange_1D_quadratic_shape_deriv(i0[i], 0, p1d));

                  // d^2()/dzeta^2
                case 5:
                  return (FEShim<2,LAGRANGE,RealType>::shape(TRI6,  SECOND, i1[i], p2d)*
                          fe_lagrange_1D_quadratic_shape_second_deriv(i0[i], 0, p1d));

                default:
                  libmesh_error_msg("Invalid shape function derivative j = " << j);
                }
            }


            // Quadratic shape functions, as defined in R. Graglia, "Higher order
            // bases on pyramidal elements", IEEE Trans Antennas and Propagation,
            // vol 47, no 5, May 1999.
          case PYRAMID14:
            {
              libmesh_assert_less (i, 14);

              const auto xi   = p(0);
              const auto eta  = p(1);
              const auto zeta = p(2);
              const Real eps  = 1.e-35;

              // The "normalized coordinates" defined by Graglia.  These are
              // the planes which define the faces of the pyramid.
              auto
                p1 = 0.5*(1. - eta - zeta), // back
                p2 = 0.5*(1. + xi  - zeta), // left
                p3 = 0.5*(1. + eta - zeta), // front
                p4 = 0.5*(1. - xi  - zeta); // right

              // Denominators are perturbed by epsilon to avoid
              // divide-by-zero issues.
              auto
                den = (-1. + zeta + eps),
                den2 = den*den,
                den3 = den2*den,
                den4 = den2*den2;

              // These terms are used in several of the derivatives
              auto
                numer_mp = xi*eta - zeta + zeta*zeta,
                numer_pm = xi*eta + zeta - zeta*zeta;

              switch (j)
                {
                case 0: // d^2()/dxi^2
                  {
                    switch(i)
                      {
                      case 0:
                      case 1:
                        return -p1*eta/den2;
                      case 2:
                      case 3:
                        return p3*eta/den2;
                      case 4:
                      case 9:
                      case 10:
                      case 11:
                      case 12:
                        return 0.;
                      case 5:
                        return 2.*p1*eta/den2;
                      case 6:
                      case 8:
                        return 4.*p1*p3/den2;
                      case 7:
                        return -2.*p3*eta/den2;
                      case 13:
                        return -8.*p1*p3/den2;
                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  }

                case 1: // d^2()/dxideta
                  {
                    switch(i)
                      {
                      case 0:
                        return 0.25*numer_mp/den2
                          - 0.5*p1*xi/den2
                          - 0.5*p4*eta/den2
                          + p4*p1/den2;

                      case 1:
                        return 0.25*numer_pm/den2
                          - 0.5*p1*xi/den2
                          + 0.5*p2*eta/den2
                          - p1*p2/den2;

                      case 2:
                        return 0.25*numer_mp/den2
                          + 0.5*p3*xi/den2
                          + 0.5*p2*eta/den2
                          + p2*p3/den2;

                      case 3:
                        return 0.25*numer_pm/den2
                          + 0.5*p3*xi/den2
                          - 0.5*p4*eta/den2
                          - p3*p4/den2;

                      case 4:
                        return 0.;

                      case 5:
                        return p4*eta/den2
                          - 2.*p4*p1/den2
                          - p2*eta/den2
                          + 2.*p1*p2/den2;

                      case 6:
                        return -p3*xi/den2
                          + p1*xi/den2
                          - 2.*p2*p3/den2
                          + 2.*p1*p2/den2;

                      case 7:
                        return p4*eta/den2
                          + 2.*p3*p4/den2
                          - p2*eta/den2
                          - 2.*p2*p3/den2;

                      case 8:
                        return -p3*xi/den2
                          + p1*xi/den2
                          - 2.*p4*p1/den2
                          + 2.*p3*p4/den2;

                      case 9:
                      case 11:
                        return -zeta/den;

                      case 10:
                      case 12:
                        return zeta/den;

                      case 13:
                        return 4.*p4*p1/den2
                          - 4.*p3*p4/den2
                          + 4.*p2*p3/den2
                          - 4.*p1*p2/den2;

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  }


                case 2: // d^2()/deta^2
                  {
                    switch(i)
                      {
                      case 0:
                      case 3:
                        return -p4*xi/den2;
                      case 1:
                      case 2:
                        return p2*xi/den2;
                      case 4:
                      case 9:
                      case 10:
                      case 11:
                      case 12:
                        return 0.;
                      case 5:
                      case 7:
                        return 4.*p2*p4/den2;
                      case 6:
                        return -2.*p2*xi/den2;
                      case 8:
                        return 2.*p4*xi/den2;
                      case 13:
                        return -8.*p2*p4/den2;
                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  }


                case 3: // d^2()/dxidzeta
                  {
                    switch(i)
                      {
                      case 0:
                        return 0.25*numer_mp/den2
                          - 0.5*p1*(2.*zeta - 1.)/den2
                          + p1*numer_mp/den3
                          - 0.5*p1*eta/den2
                          - 0.5*p4*eta/den2
                          - 2.*p4*p1*eta/den3;

                      case 1:
                        return 0.25*numer_pm/den2
                          - 0.5*p1*(1 - 2.*zeta)/den2
                          + p1*numer_pm/den3
                          + 0.5*p2*eta/den2
                          + 0.5*p1*eta/den2
                          + 2.*p1*p2*eta/den3;

                      case 2:
                        return -0.25*numer_mp/den2
                          + 0.5*p3*(2.*zeta - 1.)/den2
                          - p3*numer_mp/den3
                          - 0.5*p3*eta/den2
                          - 0.5*p2*eta/den2
                          - 2.*p2*p3*eta/den3;

                      case 3:
                        return -0.25*numer_pm/den2
                          + 0.5*p3*(1 - 2.*zeta)/den2
                          - p3*numer_pm/den3
                          + 0.5*p4*eta/den2
                          + 0.5*p3*eta/den2
                          + 2.*p3*p4*eta/den3;

                      case 4:
                        return 0.;

                      case 5:
                        return p4*eta/den2
                          + 4.*p4*p1*eta/den3
                          - p2*eta/den2
                          - 4.*p1*p2*eta/den3;

                      case 6:
                        return -p3*xi/den2
                          - p1*xi/den2
                          - 4.*p1*p3*xi/den3
                          - 2.*p2*p3/den2
                          - 2.*p1*p3/den2
                          - 2.*p1*p2/den2
                          - 8.*p1*p2*p3/den3;

                      case 7:
                        return -p4*eta/den2
                          - 4.*p3*p4*eta/den3
                          + p2*eta/den2
                          + 4.*p2*p3*eta/den3;

                      case 8:
                        return -p3*xi/den2
                          - p1*xi/den2
                          - 4.*p1*p3*xi/den3
                          + 2.*p4*p1/den2
                          + 2.*p1*p3/den2
                          + 2.*p3*p4/den2
                          + 8.*p3*p4*p1/den3;

                      case 9:
                        return -zeta/den
                          + 2.*p1/den
                          - 2.*p1*zeta/den2;

                      case 10:
                        return zeta/den
                          - 2.*p1/den
                          + 2.*p1*zeta/den2;

                      case 11:
                        return zeta/den
                          - 2.*p3/den
                          + 2.*p3*zeta/den2;

                      case 12:
                        return -zeta/den
                          + 2.*p3/den
                          - 2.*p3*zeta/den2;

                      case 13:
                        return -4.*p4*p1/den2
                          - 4.*p3*p4/den2
                          - 16.*p3*p4*p1/den3
                          + 4.*p2*p3/den2
                          + 4.*p1*p2/den2
                          + 16.*p1*p2*p3/den3;

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  }

                case 4: // d^2()/detadzeta
                  {
                    switch(i)
                      {
                      case 0:
                        return 0.25*numer_mp/den2
                          - 0.5*p4*(2.*zeta - 1.)/den2
                          + p4*numer_mp/den3
                          - 0.5*p1*xi/den2
                          - 0.5*p4*xi/den2
                          - 2.*p4*p1*xi/den3;

                      case 1:
                        return -0.25*numer_pm/den2
                          + 0.5*p2*(1. - 2.*zeta)/den2
                          - p2*numer_pm/den3
                          + 0.5*p2*xi/den2
                          + 0.5*p1*xi/den2
                          + 2.*p1*p2*xi/den3;

                      case 2:
                        return -0.25*numer_mp/den2
                          + 0.5*p2*(2.*zeta - 1.)/den2
                          - p2*numer_mp/den3
                          - 0.5*p3*xi/den2
                          - 0.5*p2*xi/den2
                          - 2.*p2*p3*xi/den3;

                      case 3:
                        return 0.25*numer_pm/den2
                          - 0.5*p4*(1. - 2.*zeta)/den2
                          + p4*numer_pm/den3
                          + 0.5*p4*xi/den2
                          + 0.5*p3*xi/den2
                          + 2.*p3*p4*xi/den3;

                      case 4:
                        return 0.;

                      case 5:
                        return -p4*eta/den2
                          - p2*eta/den2
                          - 4.*p2*p4*eta/den3
                          + 2.*p4*p1/den2
                          + 2.*p2*p4/den2
                          + 2.*p1*p2/den2
                          + 8.*p2*p1*p4/den3;

                      case 6:
                        return p3*xi/den2
                          + 4.*p2*p3*xi/den3
                          - p1*xi/den2
                          - 4.*p1*p2*xi/den3;

                      case 7:
                        return -p4*eta/den2
                          - p2*eta/den2
                          - 4.*p2*p4*eta/den3
                          - 2.*p3*p4/den2
                          - 2.*p2*p4/den2
                          - 2.*p2*p3/den2
                          - 8.*p2*p3*p4/den3;

                      case 8:
                        return p1*xi/den2
                          + 4.*p4*p1*xi/den3
                          - p3*xi/den2
                          - 4.*p3*p4*xi/den3;

                      case 9:
                        return -zeta/den
                          + 2.*p4/den
                          - 2.*p4*zeta/den2;

                      case 10:
                        return -zeta/den
                          + 2.*p2/den
                          - 2.*p2*zeta/den2;

                      case 11:
                        return zeta/den
                          - 2.*p2/den
                          + 2.*p2*zeta/den2;

                      case 12:
                        return zeta/den
                          - 2.*p4/den
                          + 2.*p4*zeta/den2;

                      case 13:
                        return 4.*p3*p4/den2
                          + 4.*p2*p3/den2
                          + 16.*p2*p3*p4/den3
                          - 4.*p4*p1/den2
                          - 4.*p1*p2/den2
                          - 16.*p2*p1*p4/den3;

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  }

                case 5: // d^2()/dzeta^2
                  {
                    switch(i)
                      {
                      case 0:
                        return 0.5*numer_mp/den2
                          - p1*(2.*zeta - 1.)/den2
                          + 2.*p1*numer_mp/den3
                          - p4*(2.*zeta - 1.)/den2
                          + 2.*p4*numer_mp/den3
                          + 2.*p4*p1/den2
                          - 4.*p4*p1*(2.*zeta - 1.)/den3
                          + 6.*p4*p1*numer_mp/den4;

                      case 1:
                        return -0.5*numer_pm/den2
                          + p2*(1 - 2.*zeta)/den2
                          - 2.*p2*numer_pm/den3
                          + p1*(1 - 2.*zeta)/den2
                          - 2.*p1*numer_pm/den3
                          + 2.*p1*p2/den2
                          + 4.*p1*p2*(1 - 2.*zeta)/den3
                          - 6.*p1*p2*numer_pm/den4;

                      case 2:
                        return 0.5*numer_mp/den2
                          - p3*(2.*zeta - 1.)/den2
                          + 2.*p3*numer_mp/den3
                          - p2*(2.*zeta - 1.)/den2
                          + 2.*p2*numer_mp/den3
                          + 2.*p2*p3/den2
                          - 4.*p2*p3*(2.*zeta - 1.)/den3
                          + 6.*p2*p3*numer_mp/den4;

                      case 3:
                        return -0.5*numer_pm/den2
                          + p4*(1 - 2.*zeta)/den2
                          - 2.*p4*numer_pm/den3
                          + p3*(1 - 2.*zeta)/den2
                          - 2.*p3*numer_pm/den3
                          + 2.*p3*p4/den2
                          + 4.*p3*p4*(1 - 2.*zeta)/den3
                          - 6.*p3*p4*numer_pm/den4;

                      case 4:
                        return 4.;

                      case 5:
                        return -2.*p1*eta/den2
                          - 2.*p4*eta/den2
                          - 8.*p4*p1*eta/den3
                          - 2.*p2*eta/den2
                          - 8.*p2*p4*eta/den3
                          - 8.*p1*p2*eta/den3
                          - 24.*p2*p1*p4*eta/den4;

                      case 6:
                        return 2.*p3*xi/den2
                          + 2.*p2*xi/den2
                          + 8.*p2*p3*xi/den3
                          + 2.*p1*xi/den2
                          + 8.*p1*p3*xi/den3
                          + 8.*p1*p2*xi/den3
                          + 24.*p1*p2*p3*xi/den4;

                      case 7:
                        return 2.*p4*eta/den2
                          + 2.*p3*eta/den2
                          + 8.*p3*p4*eta/den3
                          + 2.*p2*eta/den2
                          + 8.*p2*p4*eta/den3
                          + 8.*p2*p3*eta/den3
                          + 24.*p2*p3*p4*eta/den4;

                      case 8:
                        return -2.*p1*xi/den2
                          - 2.*p4*xi/den2
                          - 8.*p4*p1*xi/den3
                          - 2.*p3*xi/den2
                          - 8.*p1*p3*xi/den3
                          - 8.*p3*p4*xi/den3
                          - 24.*p3*p4*p1*xi/den4;

                      case 9:
                        return -2.*zeta/den
                          + 4.*p4/den
                          - 4.*p4*zeta/den2
                          + 4.*p1/den
                          - 4.*p1*zeta/den2
                          + 8.*p4*p1/den2
                          - 8.*p1*p4*zeta/den3;

                      case 10:
                        return -2.*zeta/den
                          + 4.*p1/den
                          - 4.*p1*zeta/den2
                          + 4.*p2/den
                          - 4.*p2*zeta/den2
                          + 8.*p1*p2/den2
                          - 8.*p2*p1*zeta/den3;

                      case 11:
                        return -2.*zeta/den
                          + 4.*p2/den
                          - 4.*p2*zeta/den2
                          + 4.*p3/den
                          - 4.*p3*zeta/den2
                          + 8.*p2*p3/den2
                          - 8.*p3*p2*zeta/den3;

                      case 12:
                        return -2.*zeta/den
                          + 4.*p3/den
                          - 4.*p3*zeta/den2
                          + 4.*p4/den
                          - 4.*p4*zeta/den2
                          + 8.*p3*p4/den2
                          - 8.*p4*p3*zeta/den3;

                      case 13:
                        return 8.*p3*p4/den2
                          + 8.*p2*p4/den2
                          + 8.*p2*p3/den2
                          + 32.*p2*p3*p4/den3
                          + 8.*p4*p1/den2
                          + 8.*p1*p3/den2
                          + 32.*p3*p4*p1/den3
                          + 8.*p1*p2/den2
                          + 32.*p2*p1*p4/den3
                          + 32.*p1*p2*p3/den3
                          + 96.*p1*p2*p3*p4/den4;

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  }

                default:
                  libmesh_error_msg("Invalid j = " << j);
                }
            }

            // G. Bedrosian, "Shape functions and integration formulas for
            // three-dimensional finite element analysis", Int. J. Numerical
            // Methods Engineering, vol 35, p. 95-108, 1992.
          case PYRAMID13:
            {
              libmesh_assert_less (i, 13);

              const auto xi   = p(0);
              const auto eta  = p(1);
              const auto zeta = p(2);
              const Real eps  = 1.e-35;

              // Denominators are perturbed by epsilon to avoid
              // divide-by-zero issues.
              auto
                den = (-1. + zeta + eps),
                den2 = den*den,
                den3 = den2*den,
                xi2 = xi*xi,
                eta2 = eta*eta,
                zeta2 = zeta*zeta,
                zeta3 = zeta2*zeta;

              switch (j)
                {
                case 0: // d^2()/dxi^2
                  {
                    switch(i)
                      {
                      case 0:
                      case 1:
                        return 0.5*(-1. + zeta + eta)/den;

                      case 2:
                      case 3:
                        return 0.5*(-1. + zeta - eta)/den;

                      case 4:
                      case 6:
                      case 8:
                      case 9:
                      case 10:
                      case 11:
                      case 12:
                        return 0.;

                      case 5:
                        return (1. - eta - zeta)/den;

                      case 7:
                        return (1. + eta - zeta)/den;

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  }

                case 1: // d^2()/dxideta
                  {
                    switch(i)
                      {
                      case 0:
                        return  0.25*(-1. + 2.*zeta + 2.*xi + 2.*eta)/den;

                      case 1:
                        return -0.25*(-1. + 2.*zeta - 2.*xi + 2.*eta)/den;

                      case 2:
                        return -0.25*(1. - 2.*zeta + 2.*xi + 2.*eta)/den;

                      case 3:
                        return  0.25*(1. - 2.*zeta - 2.*xi + 2.*eta)/den;

                      case 4:
                        return 0.;

                      case 5:
                        return -xi/den;

                      case 6:
                        return eta/den;

                      case 7:
                        return xi/den;

                      case 8:
                        return -eta/den;

                      case 9:
                        return -zeta/den;

                      case 10:
                        return zeta/den;

                      case 11:
                        return -zeta/den;

                      case 12:
                        return zeta/den;

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  }


                case 2: // d^2()/deta^2
                  {
                    switch(i)
                      {
                      case 0:
                      case 3:
                        return 0.5*(-1. + zeta + xi)/den;

                      case 1:
                      case 2:
                        return 0.5*(-1. + zeta - xi)/den;

                      case 4:
                      case 5:
                      case 7:
                      case 9:
                      case 10:
                      case 11:
                      case 12:
                        return 0.;

                      case 6:
                        return (1. + xi - zeta)/den;

                      case 8:
                        return (1. - xi - zeta)/den;

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  }


                case 3: // d^2()/dxidzeta
                  {
                    switch(i)
                      {
                      case 0:
                        return -0.25*(-1. + 2.*zeta - zeta2 + eta + 2.*eta*xi + eta2)/den2;

                      case 1:
                        return 0.25*(-1. + 2.*zeta - zeta2 + eta - 2.*eta*xi + eta2)/den2;

                      case 2:
                        return 0.25*(-1. + 2.*zeta - zeta2 - eta + 2.*eta*xi + eta2)/den2;

                      case 3:
                        return -0.25*(-1. + 2.*zeta - zeta2 - eta - 2.*eta*xi + eta2)/den2;

                      case 4:
                        return 0.;

                      case 5:
                        return eta*xi/den2;

                      case 6:
                        return -0.5*(1. + zeta2 + eta2 - 2.*zeta)/den2;

                      case 7:
                        return -eta*xi/den2;

                      case 8:
                        return 0.5*(1. + zeta2 + eta2 - 2.*zeta)/den2;

                      case 9:
                        return (-1. - zeta2 + eta + 2.*zeta)/den2;

                      case 10:
                        return -(-1. - zeta2 + eta + 2.*zeta)/den2;

                      case 11:
                        return (1. + zeta2 + eta - 2.*zeta)/den2;

                      case 12:
                        return -(1. + zeta2 + eta - 2.*zeta)/den2;

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  }

                case 4: // d^2()/detadzeta
                  {
                    switch(i)
                      {
                      case 0:
                        return -0.25*(-1. + 2.*zeta - zeta2 + xi + 2.*eta*xi + xi2)/den2;

                      case 1:
                        return 0.25*(1. - 2.*zeta + zeta2 + xi + 2.*eta*xi - xi2)/den2;

                      case 2:
                        return 0.25*(-1. + 2.*zeta - zeta2 - xi + 2.*eta*xi + xi2)/den2;

                      case 3:
                        return -0.25*(1. - 2.*zeta + zeta2 - xi + 2.*eta*xi - xi2)/den2;

                      case 4:
                        return 0.;

                      case 5:
                        return 0.5*(1. + xi2 + zeta2 - 2.*zeta)/den2;

                      case 6:
                        return -eta*xi/den2;

                      case 7:
                        return -0.5*(1. + xi2 + zeta2 - 2.*zeta)/den2;

                      case 8:
                        return eta*xi/den2;

                      case 9:
                        return (-1. - zeta2 + xi + 2.*zeta)/den2;

                      case 10:
                        return -(1. + zeta2 + xi - 2.*zeta)/den2;

                      case 11:
                        return (1. + zeta2 + xi - 2.*zeta)/den2;

                      case 12:
                        return -(-1. - zeta2 + xi + 2.*zeta)/den2;

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  }

                case 5: // d^2()/dzeta^2
                  {
                    switch(i)
                      {
                      case 0:
                        return 0.5*(xi + eta + 1.)*eta*xi/den3;

                      case 1:
                        return -0.5*(eta - xi + 1.)*eta*xi/den3;

                      case 2:
                        return -0.5*(xi + eta - 1.)*eta*xi/den3;

                      case 3:
                        return 0.5*(eta - xi - 1.)*eta*xi/den3;

                      case 4:
                        return 4.;

                      case 5:
                        return -(1. - 3.*zeta + 3.*zeta2 - zeta3 + eta*xi2)/den3;

                      case 6:
                        return (-1. + 3.*zeta - 3.*zeta2 + zeta3 + eta2*xi)/den3;

                      case 7:
                        return (-1. + 3.*zeta - 3.*zeta2 + zeta3 + eta*xi2)/den3;

                      case 8:
                        return -(1. - 3.*zeta + 3.*zeta2 - zeta3 + eta2*xi)/den3;

                      case 9:
                        return -2.*(-1. + 3.*zeta - 3.*zeta2 + zeta3 + eta*xi)/den3;

                      case 10:
                        return 2.*(1. - 3.*zeta + 3.*zeta2 - zeta3 + eta*xi)/den3;

                      case 11:
                        return -2.*(-1. + 3.*zeta - 3.*zeta2 + zeta3 + eta*xi)/den3;

                      case 12:
                        return 2.*(1. - 3.*zeta + 3.*zeta2 - zeta3 + eta*xi)/den3;

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  }

                default:
                  libmesh_error_msg("Invalid j = " << j);
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

#else // LIBMESH_DIM != 3
  libmesh_ignore(type, order, i, j, p);
  libmesh_not_implemented();
#endif
}

#endif // LIBMESH_ENABLE_SECOND_DERIVATIVES


} // anonymous namespace

#endif // LIBMESH_FE_LAGRANGE_SHAPE_3D_IMPL_H
