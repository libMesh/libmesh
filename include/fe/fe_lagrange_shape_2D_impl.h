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

#ifndef LIBMESH_FE_LAGRANGE_SHAPE_2D_IMPL_H
#define LIBMESH_FE_LAGRANGE_SHAPE_2D_IMPL_H

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

Real fe_lagrange_2D_shape(const ElemType,
                          const Order order,
                          const unsigned int i,
                          const Point & p);

Real fe_lagrange_2D_shape_deriv(const ElemType type,
                                const Order order,
                                const unsigned int i,
                                const unsigned int j,
                                const Point & p);

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

Real fe_lagrange_2D_shape_second_deriv(const ElemType type,
                                       const Order order,
                                       const unsigned int i,
                                       const unsigned int j,
                                       const Point & p);

#endif // LIBMESH_ENABLE_SECOND_DERIVATIVES

} // anonymous namespace



namespace libMesh
{

template <typename RealType>
typename FEShim<2,LAGRANGE,RealType>::OutputShape FEShim<2,LAGRANGE,RealType>::shape(const ElemType type,
                           const Order order,
                           const unsigned int i,
                           const Point & p)
{
  return fe_lagrange_2D_shape(type, order, i, p);
}



template <typename RealType>
typename FEShim<2,L2_LAGRANGE,RealType>::OutputShape FEShim<2,L2_LAGRANGE,RealType>::shape(const ElemType type,
                              const Order order,
                              const unsigned int i,
                              const Point & p)
{
  return fe_lagrange_2D_shape(type, order, i, p);
}


template <typename RealType>
typename FEShim<2,LAGRANGE,RealType>::OutputShape FEShim<2,LAGRANGE,RealType>::shape(const Elem * elem,
                           const Order order,
                           const unsigned int i,
                           const Point & p,
                           const bool add_p_level)
{
  libmesh_assert(elem);

  // call the orientation-independent shape functions
  return fe_lagrange_2D_shape(elem->type(), static_cast<Order>(order + add_p_level * elem->p_level()), i, p);
}



template <typename RealType>
typename FEShim<2,L2_LAGRANGE,RealType>::OutputShape FEShim<2,L2_LAGRANGE,RealType>::shape(const Elem * elem,
                              const Order order,
                              const unsigned int i,
                              const Point & p,
                              const bool add_p_level)
{
  libmesh_assert(elem);

  // call the orientation-independent shape functions
  return fe_lagrange_2D_shape(elem->type(), static_cast<Order>(order + add_p_level * elem->p_level()), i, p);
}



template <typename RealType>
typename FEShim<2,LAGRANGE,RealType>::OutputShape FEShim<2,LAGRANGE,RealType>::shape_deriv(const ElemType type,
                                 const Order order,
                                 const unsigned int i,
                                 const unsigned int j,
                                 const Point & p)
{
  return fe_lagrange_2D_shape_deriv(type, order, i, j, p);
}



template <typename RealType>
typename FEShim<2,L2_LAGRANGE,RealType>::OutputShape FEShim<2,L2_LAGRANGE,RealType>::shape_deriv(const ElemType type,
                                    const Order order,
                                    const unsigned int i,
                                    const unsigned int j,
                                    const Point & p)
{
  return fe_lagrange_2D_shape_deriv(type, order, i, j, p);
}



template <typename RealType>
typename FEShim<2,LAGRANGE,RealType>::OutputShape FEShim<2,LAGRANGE,RealType>::shape_deriv(const Elem * elem,
                                 const Order order,
                                 const unsigned int i,
                                 const unsigned int j,
                                 const Point & p,
                                 const bool add_p_level)
{
  libmesh_assert(elem);

  // call the orientation-independent shape functions
  return fe_lagrange_2D_shape_deriv(elem->type(), static_cast<Order>(order + add_p_level * elem->p_level()), i, j, p);
}



template <typename RealType>
typename FEShim<2,L2_LAGRANGE,RealType>::OutputShape FEShim<2,L2_LAGRANGE,RealType>::shape_deriv(const Elem * elem,
                                    const Order order,
                                    const unsigned int i,
                                    const unsigned int j,
                                    const Point & p,
                                    const bool add_p_level)
{
  libmesh_assert(elem);


  // call the orientation-independent shape functions
  return fe_lagrange_2D_shape_deriv(elem->type(), static_cast<Order>(order + add_p_level * elem->p_level()), i, j, p);
}



#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

template <typename RealType>
typename FEShim<2,LAGRANGE,RealType>::OutputShape FEShim<2,LAGRANGE,RealType>::shape_second_deriv(const ElemType type,
                                        const Order order,
                                        const unsigned int i,
                                        const unsigned int j,
                                        const Point & p)
{
  return fe_lagrange_2D_shape_second_deriv(type, order, i, j, p);
}



template <typename RealType>
typename FEShim<2,L2_LAGRANGE,RealType>::OutputShape FEShim<2,L2_LAGRANGE,RealType>::shape_second_deriv(const ElemType type,
                                           const Order order,
                                           const unsigned int i,
                                           const unsigned int j,
                                           const Point & p)
{
  return fe_lagrange_2D_shape_second_deriv(type, order, i, j, p);
}



template <typename RealType>
typename FEShim<2,LAGRANGE,RealType>::OutputShape FEShim<2,LAGRANGE,RealType>::shape_second_deriv(const Elem * elem,
                                        const Order order,
                                        const unsigned int i,
                                        const unsigned int j,
                                        const Point & p,
                                        const bool add_p_level)
{
  libmesh_assert(elem);

  // call the orientation-independent shape functions
  return fe_lagrange_2D_shape_second_deriv(elem->type(), static_cast<Order>(order + add_p_level * elem->p_level()), i, j, p);
}



template <typename RealType>
typename FEShim<2,L2_LAGRANGE,RealType>::OutputShape FEShim<2,L2_LAGRANGE,RealType>::shape_second_deriv(const Elem * elem,
                                           const Order order,
                                           const unsigned int i,
                                           const unsigned int j,
                                           const Point & p,
                                           const bool add_p_level)
{
  libmesh_assert(elem);

  // call the orientation-independent shape functions
  return fe_lagrange_2D_shape_second_deriv(elem->type(), static_cast<Order>(order + add_p_level * elem->p_level()), i, j, p);
}

#endif // LIBMESH_ENABLE_SECOND_DERIVATIVES

} // namespace libMesh



// Anonymous namespace function definitions
namespace
{
using namespace libMesh;

Real fe_lagrange_2D_shape(const ElemType type,
                          const Order order,
                          const unsigned int i,
                          const Point & p)
{
#if LIBMESH_DIM > 1

  switch (order)
    {
      // linear Lagrange shape functions
    case FIRST:
      {
        switch (type)
          {
          case QUAD4:
          case QUADSHELL4:
          case QUAD8:
          case QUADSHELL8:
          case QUAD9:
            {
              // Compute quad shape functions as a tensor-product
              const Real xi  = p(0);
              const Real eta = p(1);

              libmesh_assert_less (i, 4);

              //                                0  1  2  3
              static const unsigned int i0[] = {0, 1, 1, 0};
              static const unsigned int i1[] = {0, 0, 1, 1};

              return (fe_lagrange_1D_linear_shape(i0[i], xi)*
                      fe_lagrange_1D_linear_shape(i1[i], eta));
            }

          case TRI3:
          case TRISHELL3:
          case TRI6:
            {
              const Real zeta1 = p(0);
              const Real zeta2 = p(1);
              const Real zeta0 = 1. - zeta1 - zeta2;

              libmesh_assert_less (i, 3);

              switch(i)
                {
                case 0:
                  return zeta0;

                case 1:
                  return zeta1;

                case 2:
                  return zeta2;

                default:
                  libmesh_error_msg("Invalid shape function index i = " << i);
                }
            }

          default:
            libmesh_error_msg("ERROR: Unsupported 2D element type: " << type);
          }
      }


      // quadratic Lagrange shape functions
    case SECOND:
      {
        switch (type)
          {
          case QUAD8:
          case QUADSHELL8:
            {
              const Real xi  = p(0);
              const Real eta = p(1);

              libmesh_assert_less (i, 8);

              switch (i)
                {
                case 0:
                  return .25*(1. - xi)*(1. - eta)*(-1. - xi - eta);

                case 1:
                  return .25*(1. + xi)*(1. - eta)*(-1. + xi - eta);

                case 2:
                  return .25*(1. + xi)*(1. + eta)*(-1. + xi + eta);

                case 3:
                  return .25*(1. - xi)*(1. + eta)*(-1. - xi + eta);

                case 4:
                  return .5*(1. - xi*xi)*(1. - eta);

                case 5:
                  return .5*(1. + xi)*(1. - eta*eta);

                case 6:
                  return .5*(1. - xi*xi)*(1. + eta);

                case 7:
                  return .5*(1. - xi)*(1. - eta*eta);

                default:
                  libmesh_error_msg("Invalid shape function index i = " << i);
                }
            }

          case QUAD9:
            {
              // Compute quad shape functions as a tensor-product
              const Real xi  = p(0);
              const Real eta = p(1);

              libmesh_assert_less (i, 9);

              //                                0  1  2  3  4  5  6  7  8
              static const unsigned int i0[] = {0, 1, 1, 0, 2, 1, 2, 0, 2};
              static const unsigned int i1[] = {0, 0, 1, 1, 0, 2, 1, 2, 2};

              return (fe_lagrange_1D_quadratic_shape(i0[i], xi)*
                      fe_lagrange_1D_quadratic_shape(i1[i], eta));
            }

          case TRI6:
            {
              const Real zeta1 = p(0);
              const Real zeta2 = p(1);
              const Real zeta0 = 1. - zeta1 - zeta2;

              libmesh_assert_less (i, 6);

              switch(i)
                {
                case 0:
                  return 2.*zeta0*(zeta0-0.5);

                case 1:
                  return 2.*zeta1*(zeta1-0.5);

                case 2:
                  return 2.*zeta2*(zeta2-0.5);

                case 3:
                  return 4.*zeta0*zeta1;

                case 4:
                  return 4.*zeta1*zeta2;

                case 5:
                  return 4.*zeta2*zeta0;

                default:
                  libmesh_error_msg("Invalid shape function index i = " << i);
                }
            }

          default:
            libmesh_error_msg("ERROR: Unsupported 2D element type: " << type);
          }
      }



      // unsupported order
    default:
      libmesh_error_msg("ERROR: Unsupported 2D FE order: " << order);
    }
#else // LIBMESH_DIM > 1
  libmesh_ignore(type, order, i, p);
  libmesh_not_implemented();
#endif
}



Real fe_lagrange_2D_shape_deriv(const ElemType type,
                                const Order order,
                                const unsigned int i,
                                const unsigned int j,
                                const Point & p)
{
#if LIBMESH_DIM > 1

  libmesh_assert_less (j, 2);

  switch (order)
    {
      // linear Lagrange shape functions
    case FIRST:
      {
        switch (type)
          {
          case QUAD4:
          case QUADSHELL4:
          case QUAD8:
          case QUADSHELL8:
          case QUAD9:
            {
              // Compute quad shape functions as a tensor-product
              const Real xi  = p(0);
              const Real eta = p(1);

              libmesh_assert_less (i, 4);

              //                                0  1  2  3
              static const unsigned int i0[] = {0, 1, 1, 0};
              static const unsigned int i1[] = {0, 0, 1, 1};

              switch (j)
                {
                  // d()/dxi
                case 0:
                  return (fe_lagrange_1D_linear_shape_deriv(i0[i], 0, xi)*
                          fe_lagrange_1D_linear_shape      (i1[i], eta));

                  // d()/deta
                case 1:
                  return (fe_lagrange_1D_linear_shape      (i0[i], xi)*
                          fe_lagrange_1D_linear_shape_deriv(i1[i], 0, eta));

                default:
                  libmesh_error_msg("ERROR: Invalid derivative index j = " << j);
                }
            }

          case TRI3:
          case TRISHELL3:
          case TRI6:
            {
              libmesh_assert_less (i, 3);

              const Real dzeta0dxi  = -1.;
              const Real dzeta1dxi  = 1.;
              const Real dzeta2dxi  = 0.;

              const Real dzeta0deta = -1.;
              const Real dzeta1deta = 0.;
              const Real dzeta2deta = 1.;

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

                      default:
                        libmesh_error_msg("Invalid shape function index i = " << i);
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

                      default:
                        libmesh_error_msg("Invalid shape function index i = " << i);
                      }
                  }
                default:
                  libmesh_error_msg("ERROR: Invalid derivative index j = " << j);
                }
            }

          default:
            libmesh_error_msg("ERROR: Unsupported 2D element type: " << type);
          }
      }


      // quadratic Lagrange shape functions
    case SECOND:
      {
        switch (type)
          {
          case QUAD8:
          case QUADSHELL8:
            {
              const Real xi  = p(0);
              const Real eta = p(1);

              libmesh_assert_less (i, 8);

              switch (j)
                {
                  // d/dxi
                case 0:
                  switch (i)
                    {
                    case 0:
                      return .25*(1. - eta)*((1. - xi)*(-1.) +
                                             (-1.)*(-1. - xi - eta));

                    case 1:
                      return .25*(1. - eta)*((1. + xi)*(1.) +
                                             (1.)*(-1. + xi - eta));

                    case 2:
                      return .25*(1. + eta)*((1. + xi)*(1.) +
                                             (1.)*(-1. + xi + eta));

                    case 3:
                      return .25*(1. + eta)*((1. - xi)*(-1.) +
                                             (-1.)*(-1. - xi + eta));

                    case 4:
                      return .5*(-2.*xi)*(1. - eta);

                    case 5:
                      return .5*(1.)*(1. - eta*eta);

                    case 6:
                      return .5*(-2.*xi)*(1. + eta);

                    case 7:
                      return .5*(-1.)*(1. - eta*eta);

                    default:
                      libmesh_error_msg("Invalid shape function index i = " << i);
                    }

                  // d/deta
                case 1:
                  switch (i)
                    {
                    case 0:
                      return .25*(1. - xi)*((1. - eta)*(-1.) +
                                            (-1.)*(-1. - xi - eta));

                    case 1:
                      return .25*(1. + xi)*((1. - eta)*(-1.) +
                                            (-1.)*(-1. + xi - eta));

                    case 2:
                      return .25*(1. + xi)*((1. + eta)*(1.) +
                                            (1.)*(-1. + xi + eta));

                    case 3:
                      return .25*(1. - xi)*((1. + eta)*(1.) +
                                            (1.)*(-1. - xi + eta));

                    case 4:
                      return .5*(1. - xi*xi)*(-1.);

                    case 5:
                      return .5*(1. + xi)*(-2.*eta);

                    case 6:
                      return .5*(1. - xi*xi)*(1.);

                    case 7:
                      return .5*(1. - xi)*(-2.*eta);

                    default:
                      libmesh_error_msg("Invalid shape function index i = " << i);
                    }

                default:
                  libmesh_error_msg("ERROR: Invalid derivative index j = " << j);
                }
            }

          case QUAD9:
            {
              // Compute quad shape functions as a tensor-product
              const Real xi  = p(0);
              const Real eta = p(1);

              libmesh_assert_less (i, 9);

              //                                0  1  2  3  4  5  6  7  8
              static const unsigned int i0[] = {0, 1, 1, 0, 2, 1, 2, 0, 2};
              static const unsigned int i1[] = {0, 0, 1, 1, 0, 2, 1, 2, 2};

              switch (j)
                {
                  // d()/dxi
                case 0:
                  return (fe_lagrange_1D_quadratic_shape_deriv(i0[i], 0, xi)*
                          fe_lagrange_1D_quadratic_shape      (i1[i], eta));

                  // d()/deta
                case 1:
                  return (fe_lagrange_1D_quadratic_shape      (i0[i], xi)*
                          fe_lagrange_1D_quadratic_shape_deriv(i1[i], 0, eta));

                default:
                  libmesh_error_msg("ERROR: Invalid derivative index j = " << j);
                }
            }

          case TRI6:
            {
              libmesh_assert_less (i, 6);

              const Real zeta1 = p(0);
              const Real zeta2 = p(1);
              const Real zeta0 = 1. - zeta1 - zeta2;

              const Real dzeta0dxi  = -1.;
              const Real dzeta1dxi  = 1.;
              const Real dzeta2dxi  = 0.;

              const Real dzeta0deta = -1.;
              const Real dzeta1deta = 0.;
              const Real dzeta2deta = 1.;

              switch(j)
                {
                case 0:
                  {
                    switch(i)
                      {
                      case 0:
                        return (4.*zeta0-1.)*dzeta0dxi;

                      case 1:
                        return (4.*zeta1-1.)*dzeta1dxi;

                      case 2:
                        return (4.*zeta2-1.)*dzeta2dxi;

                      case 3:
                        return 4.*zeta1*dzeta0dxi + 4.*zeta0*dzeta1dxi;

                      case 4:
                        return 4.*zeta2*dzeta1dxi + 4.*zeta1*dzeta2dxi;

                      case 5:
                        return 4.*zeta2*dzeta0dxi + 4*zeta0*dzeta2dxi;

                      default:
                        libmesh_error_msg("Invalid shape function index i = " << i);
                      }
                  }

                case 1:
                  {
                    switch(i)
                      {
                      case 0:
                        return (4.*zeta0-1.)*dzeta0deta;

                      case 1:
                        return (4.*zeta1-1.)*dzeta1deta;

                      case 2:
                        return (4.*zeta2-1.)*dzeta2deta;

                      case 3:
                        return 4.*zeta1*dzeta0deta + 4.*zeta0*dzeta1deta;

                      case 4:
                        return 4.*zeta2*dzeta1deta + 4.*zeta1*dzeta2deta;

                      case 5:
                        return 4.*zeta2*dzeta0deta + 4*zeta0*dzeta2deta;

                      default:
                        libmesh_error_msg("Invalid shape function index i = " << i);
                      }
                  }
                default:
                  libmesh_error_msg("ERROR: Invalid derivative index j = " << j);
                }
            }

          default:
            libmesh_error_msg("ERROR: Unsupported 2D element type: " << type);
          }
      }

      // unsupported order
    default:
      libmesh_error_msg("ERROR: Unsupported 2D FE order: " << order);
    }
#else // LIBMESH_DIM > 1
  libmesh_ignore(type, order, i, j, p);
  libmesh_not_implemented();
#endif
}



#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

Real fe_lagrange_2D_shape_second_deriv(const ElemType type,
                                       const Order order,
                                       const unsigned int i,
                                       const unsigned int j,
                                       const Point & p)
{
#if LIBMESH_DIM > 1

  // j = 0 ==> d^2 phi / dxi^2
  // j = 1 ==> d^2 phi / dxi deta
  // j = 2 ==> d^2 phi / deta^2
  libmesh_assert_less (j, 3);

  switch (order)
    {
      // linear Lagrange shape functions
    case FIRST:
      {
        switch (type)
          {
          case QUAD4:
          case QUADSHELL4:
          case QUAD8:
          case QUADSHELL8:
          case QUAD9:
            {
              // Compute quad shape functions as a tensor-product
              const Real xi  = p(0);
              const Real eta = p(1);

              libmesh_assert_less (i, 4);

              //                                0  1  2  3
              static const unsigned int i0[] = {0, 1, 1, 0};
              static const unsigned int i1[] = {0, 0, 1, 1};

              switch (j)
                {
                  // d^2() / dxi^2
                case 0:
                  return 0.;

                  // d^2() / dxi deta
                case 1:
                  return (fe_lagrange_1D_linear_shape_deriv(i0[i], 0, xi)*
                          fe_lagrange_1D_linear_shape_deriv(i1[i], 0, eta));

                  // d^2() / deta^2
                case 2:
                  return 0.;

                default:
                  libmesh_error_msg("ERROR: Invalid derivative index j = " << j);
                }
            }

          case TRI3:
          case TRISHELL3:
          case TRI6:
            {
              // All second derivatives for linear triangles are zero.
              return 0.;
            }

          default:
            libmesh_error_msg("ERROR: Unsupported 2D element type: " << type);

          } // end switch (type)
      } // end case FIRST


      // quadratic Lagrange shape functions
    case SECOND:
      {
        switch (type)
          {
          case QUAD8:
          case QUADSHELL8:
            {
              const Real xi  = p(0);
              const Real eta = p(1);

              libmesh_assert_less (j, 3);

              switch (j)
                {
                  // d^2() / dxi^2
                case 0:
                  {
                    switch (i)
                      {
                      case 0:
                      case 1:
                        return 0.5*(1.-eta);

                      case 2:
                      case 3:
                        return 0.5*(1.+eta);

                      case 4:
                        return eta - 1.;

                      case 5:
                      case 7:
                        return 0.0;

                      case 6:
                        return -1. - eta;

                      default:
                        libmesh_error_msg("Invalid shape function index i = " << i);
                      }
                  }

                  // d^2() / dxi deta
                case 1:
                  {
                    switch (i)
                      {
                      case 0:
                        return 0.25*( 1. - 2.*xi - 2.*eta);

                      case 1:
                        return 0.25*(-1. - 2.*xi + 2.*eta);

                      case 2:
                        return 0.25*( 1. + 2.*xi + 2.*eta);

                      case 3:
                        return 0.25*(-1. + 2.*xi - 2.*eta);

                      case 4:
                        return xi;

                      case 5:
                        return -eta;

                      case 6:
                        return -xi;

                      case 7:
                        return eta;

                      default:
                        libmesh_error_msg("Invalid shape function index i = " << i);
                      }
                  }

                  // d^2() / deta^2
                case 2:
                  {
                    switch (i)
                      {
                      case 0:
                      case 3:
                        return 0.5*(1.-xi);

                      case 1:
                      case 2:
                        return 0.5*(1.+xi);

                      case 4:
                      case 6:
                        return 0.0;

                      case 5:
                        return -1.0 - xi;

                      case 7:
                        return xi - 1.0;

                      default:
                        libmesh_error_msg("Invalid shape function index i = " << i);
                      }
                  }

                default:
                  libmesh_error_msg("ERROR: Invalid derivative index j = " << j);
                } // end switch (j)
            } // end case QUAD8

          case QUAD9:
            {
              // Compute QUAD9 second derivatives as tensor product
              const Real xi  = p(0);
              const Real eta = p(1);

              libmesh_assert_less (i, 9);

              //                                0  1  2  3  4  5  6  7  8
              static const unsigned int i0[] = {0, 1, 1, 0, 2, 1, 2, 0, 2};
              static const unsigned int i1[] = {0, 0, 1, 1, 0, 2, 1, 2, 2};

              switch (j)
                {
                  // d^2() / dxi^2
                case 0:
                  return (fe_lagrange_1D_quadratic_shape_second_deriv(i0[i], 0, xi)*
                          fe_lagrange_1D_quadratic_shape             (i1[i], eta));

                  // d^2() / dxi deta
                case 1:
                  return (fe_lagrange_1D_quadratic_shape_deriv(i0[i], 0, xi)*
                          fe_lagrange_1D_quadratic_shape_deriv(i1[i], 0, eta));

                  // d^2() / deta^2
                case 2:
                  return (fe_lagrange_1D_quadratic_shape             (i0[i], xi)*
                          fe_lagrange_1D_quadratic_shape_second_deriv(i1[i], 0, eta));

                default:
                  libmesh_error_msg("ERROR: Invalid derivative index j = " << j);
                }  // end switch (j)
            } // end case QUAD9

          case TRI6:
            {
              const Real dzeta0dxi  = -1.;
              const Real dzeta1dxi  = 1.;
              const Real dzeta2dxi  = 0.;

              const Real dzeta0deta = -1.;
              const Real dzeta1deta = 0.;
              const Real dzeta2deta = 1.;

              libmesh_assert_less (j, 3);

              switch (j)
                {
                  // d^2() / dxi^2
                case 0:
                  {
                    switch (i)
                      {
                      case 0:
                        return 4.*dzeta0dxi*dzeta0dxi;

                      case 1:
                        return 4.*dzeta1dxi*dzeta1dxi;

                      case 2:
                        return 4.*dzeta2dxi*dzeta2dxi;

                      case 3:
                        return 8.*dzeta0dxi*dzeta1dxi;

                      case 4:
                        return 8.*dzeta1dxi*dzeta2dxi;

                      case 5:
                        return 8.*dzeta0dxi*dzeta2dxi;

                      default:
                        libmesh_error_msg("Invalid shape function index i = " << i);
                      }
                  }

                  // d^2() / dxi deta
                case 1:
                  {
                    switch (i)
                      {
                      case 0:
                        return 4.*dzeta0dxi*dzeta0deta;

                      case 1:
                        return 4.*dzeta1dxi*dzeta1deta;

                      case 2:
                        return 4.*dzeta2dxi*dzeta2deta;

                      case 3:
                        return 4.*dzeta1deta*dzeta0dxi + 4.*dzeta0deta*dzeta1dxi;

                      case 4:
                        return 4.*dzeta2deta*dzeta1dxi + 4.*dzeta1deta*dzeta2dxi;

                      case 5:
                        return 4.*dzeta2deta*dzeta0dxi + 4.*dzeta0deta*dzeta2dxi;

                      default:
                        libmesh_error_msg("Invalid shape function index i = " << i);
                      }
                  }

                  // d^2() / deta^2
                case 2:
                  {
                    switch (i)
                      {
                      case 0:
                        return 4.*dzeta0deta*dzeta0deta;

                      case 1:
                        return 4.*dzeta1deta*dzeta1deta;

                      case 2:
                        return 4.*dzeta2deta*dzeta2deta;

                      case 3:
                        return 8.*dzeta0deta*dzeta1deta;

                      case 4:
                        return 8.*dzeta1deta*dzeta2deta;

                      case 5:
                        return 8.*dzeta0deta*dzeta2deta;

                      default:
                        libmesh_error_msg("Invalid shape function index i = " << i);
                      }
                  }

                default:
                  libmesh_error_msg("ERROR: Invalid derivative index j = " << j);
                } // end switch (j)
            }  // end case TRI6

          default:
            libmesh_error_msg("ERROR: Unsupported 2D element type: " << type);
          }
      } // end case SECOND



      // unsupported order
    default:
      libmesh_error_msg("ERROR: Unsupported 2D FE order: " << order);

    } // end switch (order)

#else // LIBMESH_DIM > 1
  libmesh_ignore(type, order, i, j, p);
  libmesh_not_implemented();
#endif
}

#endif // LIBMESH_ENABLE_SECOND_DERIVATIVES

} // anonymous namespace

#endif // LIBMESH_FE_LAGRANGE_SHAPE_2D_IMPL_H
