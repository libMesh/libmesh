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

#ifndef LIBMESH_FE_BERNSTEIN_SHAPE_1D_IMPL_H
#define LIBMESH_FE_BERNSTEIN_SHAPE_1D_IMPL_H

// Local includes
#include "libmesh/libmesh_config.h"
#ifdef LIBMESH_ENABLE_HIGHER_ORDER_SHAPES

#include <cstdlib> // *must* precede <cmath> for proper std:abs() on PGI, Sun Studio CC
#include <cmath>

#include "libmesh/libmesh_common.h"
#include "libmesh/fe.h"
#include "libmesh/elem.h"
#include "libmesh/utility.h"


namespace libMesh
{


template <typename RealType>
typename FEShim<1,BERNSTEIN,RealType>::OutputShape FEShim<1,BERNSTEIN,RealType>::shape(const ElemType,
                            const Order order,
                            const unsigned int i,
                            const Point & p)
{
  const auto xi = p(0);
  using Utility::pow;

  switch (order)
    {
    case FIRST:

      switch(i)
        {
        case 0:
          return (1.-xi)/2.;
        case 1:
          return (1.+xi)/2.;
        default:
          libmesh_error_msg("Invalid shape function index i = " << i);
        }

    case SECOND:

      switch(i)
        {
        case 0:
          return (1./4.)*pow<2>(1.-xi);
        case 1:
          return (1./4.)*pow<2>(1.+xi);
        case 2:
          return (1./2.)*(1.-xi)*(1.+xi);
        default:
          libmesh_error_msg("Invalid shape function index i = " << i);
        }

    case THIRD:

      switch(i)
        {
        case 0:
          return (1./8.)*pow<3>(1.-xi);
        case 1:
          return (1./8.)*pow<3>(1.+xi);
        case 2:
          return (3./8.)*(1.+xi)*pow<2>(1.-xi);
        case 3:
          return (3./8.)*pow<2>(1.+xi)*(1.-xi);
        default:
          libmesh_error_msg("Invalid shape function index i = " << i);
        }

    case FOURTH:

      switch(i)
        {
        case 0:
          return (1./16.)*pow<4>(1.-xi);
        case 1:
          return (1./16.)*pow<4>(1.+xi);
        case 2:
          return (1./ 4.)*(1.+xi)*pow<3>(1.-xi);
        case 3:
          return (3./ 8.)*pow<2>(1.+xi)*pow<2>(1.-xi);
        case 4:
          return (1./ 4.)*pow<3>(1.+xi)*(1.-xi);
        default:
          libmesh_error_msg("Invalid shape function index i = " << i);
        }


    case FIFTH:

      switch(i)
        {
        case 0:
          return (1./32.)*pow<5>(1.-xi);
        case 1:
          return (1./32.)*pow<5>(1.+xi);
        case 2:
          return (5./32.)*(1.+xi)*pow<4>(1.-xi);
        case 3:
          return (5./16.)*pow<2>(1.+xi)*pow<3>(1.-xi);
        case 4:
          return (5./16.)*pow<3>(1.+xi)*pow<2>(1.-xi);
        case 5:
          return (5./32.)*pow<4>(1.+xi)*(1.-xi);
        default:
          libmesh_error_msg("Invalid shape function index i = " << i);
        }


    case SIXTH:

      switch (i)
        {
        case 0:
          return ( 1./64.)*pow<6>(1.-xi);
        case 1:
          return ( 1./64.)*pow<6>(1.+xi);
        case 2:
          return ( 3./32.)*(1.+xi)*pow<5>(1.-xi);
        case 3:
          return (15./64.)*pow<2>(1.+xi)*pow<4>(1.-xi);
        case 4:
          return ( 5./16.)*pow<3>(1.+xi)*pow<3>(1.-xi);
        case 5:
          return (15./64.)*pow<4>(1.+xi)*pow<2>(1.-xi);
        case 6:
          return ( 3./32.)*pow<5>(1.+xi)*(1.-xi);
        default:
          libmesh_error_msg("Invalid shape function index i = " << i);
        }

    default:
      {
        libmesh_assert (order>6);

        // Use this for arbitrary orders.
        // Note that this implementation is less efficient.
        const int p_order = static_cast<int>(order);
        const int m       = p_order-i+1;
        const int n       = (i-1);

        Real binomial_p_i = 1;

        // the binomial coefficient (p choose n)
        // Using an unsigned long here will work for any of the orders we support.
        // Explicitly construct a Real to prevent conversion warnings
        if (i>1)
          binomial_p_i = Real(Utility::binomial(static_cast<unsigned long>(p_order),
                                                static_cast<unsigned long>(n)));

        switch(i)
          {
          case 0:
            return binomial_p_i * std::pow((1-xi)/2, p_order);
          case 1:
            return binomial_p_i * std::pow((1+xi)/2, p_order);
          default:
            {
              return binomial_p_i * std::pow((1+xi)/2,n)
                * std::pow((1-xi)/2,m);
            }
          }
      }
    }
}



template <typename RealType>
typename FEShim<1,BERNSTEIN,RealType>::OutputShape FEShim<1,BERNSTEIN,RealType>::shape(const Elem * elem,
                            const Order order,
                            const unsigned int i,
                            const Point & p,
                            const bool add_p_level)
{
  libmesh_assert(elem);

  return FEShim<1,BERNSTEIN,RealType>::shape
    (elem->type(),
     static_cast<Order>(order + add_p_level*elem->p_level()), i, p);
}



template <typename RealType>
typename FEShim<1,BERNSTEIN,RealType>::OutputShape FEShim<1,BERNSTEIN,RealType>::shape_deriv(const ElemType,
                                  const Order order,
                                  const unsigned int i,
                                  const unsigned int libmesh_dbg_var(j),
                                  const Point & p)
{
  // only d()/dxi in 1D!

  libmesh_assert_equal_to (j, 0);

  const auto xi = p(0);

  using Utility::pow;

  switch (order)
    {
    case FIRST:

      switch(i)
        {
        case 0:
          return -.5;
        case 1:
          return .5;
        default:
          libmesh_error_msg("Invalid shape function index i = " << i);
        }

    case SECOND:

      switch(i)
        {
        case 0:
          return (xi-1.)*.5;
        case 1:
          return (xi+1.)*.5;
        case 2:
          return -xi;
        default:
          libmesh_error_msg("Invalid shape function index i = " << i);
        }

    case THIRD:

      switch(i)
        {
        case 0:
          return -0.375*pow<2>(1.-xi);
        case 1:
          return  0.375*pow<2>(1.+xi);
        case 2:
          return -0.375 -.75*xi +1.125*pow<2>(xi);
        case 3:
          return  0.375 -.75*xi -1.125*pow<2>(xi);
        default:
          libmesh_error_msg("Invalid shape function index i = " << i);
        }

    case FOURTH:

      switch(i)
        {
        case 0:
          return -0.25*pow<3>(1.-xi);
        case 1:
          return  0.25*pow<3>(1.+xi);
        case 2:
          return -0.5 +1.5*pow<2>(xi)-pow<3>(xi);
        case 3:
          return  1.5*(pow<3>(xi)-xi);
        case 4:
          return  0.5 -1.5*pow<2>(xi)-pow<3>(xi);
        default:
          libmesh_error_msg("Invalid shape function index i = " << i);
        }

    case FIFTH:

      switch(i)
        {
        case 0:
          return -(5./32.)*pow<4>(xi-1.);
        case 1:
          return  (5./32.)*pow<4>(xi+1.);
        case 2:
          return  (5./32.)*pow<4>(1.-xi)         -(5./8.)*(1.+xi)*pow<3>(1.-xi);
        case 3:
          return  (5./ 8.)*(1.+xi)*pow<3>(1.-xi) -(15./16.)*pow<2>(1.+xi)*pow<2>(1.-xi);
        case 4:
          return -(5./ 8.)*pow<3>(1.+xi)*(1.-xi) +(15./16.)*pow<2>(1.+xi)*pow<2>(1.-xi);
        case 5:
          return  (5./ 8.)*pow<3>(1.+xi)*(1.-xi) -(5./32.)*pow<4>(1.+xi);
        default:
          libmesh_error_msg("Invalid shape function index i = " << i);
        }

    case SIXTH:

      switch(i)
        {
        case 0:
          return -( 3./32.)*pow<5>(1.-xi);
        case 1:
          return  ( 3./32.)*pow<5>(1.+xi);
        case 2:
          return  ( 3./32.)*pow<5>(1.-xi)-(15./32.)*(1.+xi)*pow<4>(1.-xi);
        case 3:
          return  (15./32.)*(1.+xi)*pow<4>(1.-xi)-(15./16.)*pow<2>(1.+xi)*pow<3>(1.-xi);
        case 4:
          return -(15./ 8.)*xi +(15./4.)*pow<3>(xi)-(15./8.)*pow<5>(xi);
        case 5:
          return -(15./32.)*(1.-xi)*pow<4>(1.+xi)+(15./16.)*pow<2>(1.-xi)*pow<3>(1.+xi);
        case 6:
          return  (15./32.)*pow<4>(1.+xi)*(1.-xi)-(3./32.)*pow<5>(1.+xi);
        default:
          libmesh_error_msg("Invalid shape function index i = " << i);
        }


    default:
      {
        libmesh_assert (order>6);

        // Use this for arbitrary orders
        const int p_order = static_cast<int>(order);
        const int m       = p_order-(i-1);
        const int n       = (i-1);

        Real binomial_p_i = 1;

        // the binomial coefficient (p choose n)
        // Using an unsigned long here will work for any of the orders we support.
        // Explicitly construct a Real to prevent conversion warnings
        if (i>1)
          binomial_p_i = Real(Utility::binomial(static_cast<unsigned long>(p_order),
                                                static_cast<unsigned long>(n)));

        switch(i)
          {
          case 0:
            return binomial_p_i * (-1./2.) * p_order * std::pow((1-xi)/2, p_order-1);
          case 1:
            return binomial_p_i * ( 1./2.) * p_order * std::pow((1+xi)/2, p_order-1);

          default:
            {
              return binomial_p_i * (1./2. * n * std::pow((1+xi)/2,n-1) * std::pow((1-xi)/2,m)
                                     - 1./2. * m * std::pow((1+xi)/2,n)   * std::pow((1-xi)/2,m-1));
            }
          }
      }

    }
}



template <typename RealType>
typename FEShim<1,BERNSTEIN,RealType>::OutputShape FEShim<1,BERNSTEIN,RealType>::shape_deriv(const Elem * elem,
                                  const Order order,
                                  const unsigned int i,
                                  const unsigned int j,
                                  const Point & p,
                                  const bool add_p_level)
{
  libmesh_assert(elem);

  return FEShim<1,BERNSTEIN,RealType>::shape_deriv
    (elem->type(),
     static_cast<Order>(order + add_p_level*elem->p_level()), i, j, p);
}


#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

template <typename RealType>
typename FEShim<1,BERNSTEIN,RealType>::OutputShape FEShim<1,BERNSTEIN,RealType>::shape_second_deriv(const ElemType,
                                         const Order order,
                                         const unsigned int i,
                                         const unsigned int libmesh_dbg_var(j),
                                         const Point & p)
{
  // only d^2()/dxi^2 in 1D!

  libmesh_assert_equal_to (j, 0);

  const auto xi = p(0);

  using Utility::pow;

  switch (order)
    {
    case FIRST:

      switch(i)
        {
        case 0:
        case 1:
          return 0;
        default:
          libmesh_error_msg("Invalid shape function index i = " << i);
        }

    case SECOND:

      switch(i)
        {
        case 0:
        case 1:
          return .5;
        case 2:
          return -1;
        default:
          libmesh_error_msg("Invalid shape function index i = " << i);
        }

    case THIRD:

      switch(i)
        {
        case 0:
          return 0.75*(1.-xi);
        case 1:
          return 0.75*(1.+xi);
        case 2:
          return -.75 + 2.25*xi;
        case 3:
          return -.75 - 2.25*xi;
        default:
          libmesh_error_msg("Invalid shape function index i = " << i);
        }

    case FOURTH:

      switch(i)
        {
        case 0:
          return  0.75*pow<2>(1.-xi);
        case 1:
          return  0.75*pow<2>(1.+xi);
        case 2:
          return  3*(xi - pow<2>(xi));
        case 3:
          return  1.5*(3*pow<2>(xi)-1);
        case 4:
          return  -3*xi-3*pow<2>(xi);
        default:
          libmesh_error_msg("Invalid shape function index i = " << i);
        }

    case FIFTH:

      switch(i)
        {
        case 0:
          return -(5./8.)*pow<3>(xi-1.);
        case 1:
          return  (5./8.)*pow<3>(xi+1.);
        case 2:
          return -(5./4.)*pow<3>(1.-xi) + (15./8.)*(1.+xi)*pow<2>(1.-xi);
        case 3:
          return -(15./ 4.)*(1.+xi)*pow<2>(1.-xi) + (5./ 8.)*pow<3>(1.-xi)
          + (15./8.)*pow<2>(1.+xi)*(1.-xi);
        case 4:
          return  (5./ 8.)*pow<3>(1.+xi) - (15./ 4.)*pow<2>(1.+xi)*(1.-xi)
          +(15./8.)*(1.+xi)*pow<2>(1.-xi);
        case 5:
          return -(5./ 8.)*pow<3>(1.+xi) + (15./ 8.)*pow<2>(1.+xi)*(1.-xi)
          -(5./8.)*pow<3>(1.+xi);
        default:
          libmesh_error_msg("Invalid shape function index i = " << i);
        }

    case SIXTH:

      switch(i)
        {
        case 0:
          return  ( 15./32.)*pow<4>(1.-xi);
        case 1:
          return  ( 15./32.)*pow<4>(1.+xi);
        case 2:
          return -( 15./8.)*pow<4>(1.-xi) +
                  ( 15./8.)*(1.+xi)*pow<3>(1.-xi);
        case 3:
          return -(15./4.)*(1.+xi)*pow<3>(1.-xi)
                  + (15./32.)*pow<4>(1.-xi)
                  + (45./16.)*pow<2>(1.+xi)*pow<2>(1.-xi);
        case 4:
          return -(15./ 8.) +(45./4.)*pow<2>(xi) - (75./8.)*pow<4>(xi);
        case 5:
          return -(15./4.)*(1.-xi)*pow<3>(1.+xi)
                  + (15./32.)*pow<4>(1.+xi)
                  + (45./16.)*pow<2>(1.-xi)*pow<2>(1.+xi);
        case 6:
          return -(15./16.)*pow<4>(1.+xi)
                 + (15./8.)*pow<3>(1.+xi)*(1.-xi);
        default:
          libmesh_error_msg("Invalid shape function index i = " << i);
        }


    default:
      {
        libmesh_assert (order>6);

        // Use this for arbitrary orders
        const int p_order = static_cast<int>(order);
        const int m       = p_order-(i-1);
        const int n       = (i-1);

        Real binomial_p_i = 1;

        // the binomial coefficient (p choose n)
        // Using an unsigned long here will work for any of the orders we support.
        // Explicitly construct a Real to prevent conversion warnings
        if (i>1)
          binomial_p_i = Real(Utility::binomial(static_cast<unsigned long>(p_order),
                                                static_cast<unsigned long>(n)));

        switch(i)
          {
          case 0:
            return binomial_p_i * (1./4.) * p_order * (p_order-1) * std::pow((1-xi)/2, p_order-2);
          case 1:
            return binomial_p_i * (1./4.) * p_order * (p_order-1) * std::pow((1+xi)/2, p_order-2);

          default:
            {
              RealType val = 0;

              if (n == 1)
                val +=
                  binomial_p_i * (-1./4. * m * std::pow((1-xi)/2,m-1));
              else
                val +=
                  binomial_p_i * (-1./4. * n * m * std::pow((1+xi)/2,n-1) * std::pow((1-xi)/2,m-1) +
                                  1./4. * n * (n-1) * std::pow((1+xi)/2,n-2) * std::pow((1-xi)/2,m));

              if (m == 1)
                val += binomial_p_i * (-1./4. * n * std::pow((1+xi)/2,n-1));
              else
                val +=
                  binomial_p_i * (1./4. * m * (m-1) * std::pow((1+xi)/2,n)   * std::pow((1-xi)/2,m-2)
                                  - 1./4. * m * n * std::pow((1+xi)/2,n-1)   * std::pow((1-xi)/2,m-1));

              return val;
            }
          }
      }

    }
}



template <typename RealType>
typename FEShim<1,BERNSTEIN,RealType>::OutputShape FEShim<1,BERNSTEIN,RealType>::shape_second_deriv(const Elem * elem,
                                         const Order order,
                                         const unsigned int i,
                                         const unsigned int j,
                                         const Point & p,
                                         const bool add_p_level)
{
  libmesh_assert(elem);

  return FEShim<1,BERNSTEIN,RealType>::shape_second_deriv
    (elem->type(),
     static_cast<Order>(order + add_p_level*elem->p_level()), i, j, p);
}

#endif

} // namespace libMesh


#endif //LIBMESH_ENABLE_HIGHER_ORDER_SHAPES

#endif // LIBMESH_FE_BERNSTEIN_SHAPE_1D_IMPL_H
