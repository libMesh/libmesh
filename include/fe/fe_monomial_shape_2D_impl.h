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

#ifndef LIBMESH_FE_MONOMIAL_SHAPE_2D_IMPL_H
#define LIBMESH_FE_MONOMIAL_SHAPE_2D_IMPL_H

// C++ includes

// Local includes
#include "libmesh/fe.h"
#include "libmesh/elem.h"

namespace libMesh
{




template <typename RealType>
typename FEShim<2,MONOMIAL,RealType>::OutputShape FEShim<2,MONOMIAL,RealType>::shape(const ElemType,
                           const Order libmesh_dbg_var(order),
                           const unsigned int i,
                           const Point & p)
{
#if LIBMESH_DIM > 1

  libmesh_assert_less (i, (static_cast<unsigned int>(order)+1)*
                       (static_cast<unsigned int>(order)+2)/2);

  const Real xi  = p(0);
  const Real eta = p(1);

  switch (i)
    {
      // constant
    case 0:
      return 1.;

      // linear
    case 1:
      return xi;

    case 2:
      return eta;

      // quadratics
    case 3:
      return xi*xi;

    case 4:
      return xi*eta;

    case 5:
      return eta*eta;

      // cubics
    case 6:
      return xi*xi*xi;

    case 7:
      return xi*xi*eta;

    case 8:
      return xi*eta*eta;

    case 9:
      return eta*eta*eta;

      // quartics
    case 10:
      return xi*xi*xi*xi;

    case 11:
      return xi*xi*xi*eta;

    case 12:
      return xi*xi*eta*eta;

    case 13:
      return xi*eta*eta*eta;

    case 14:
      return eta*eta*eta*eta;

    default:
      unsigned int o = 0;
      for (; i >= (o+1)*(o+2)/2; o++) { }
      unsigned int ny = i - (o*(o+1)/2);
      unsigned int nx = o - ny;
      Real val = 1.;
      for (unsigned int index=0; index != nx; index++)
        val *= xi;
      for (unsigned int index=0; index != ny; index++)
        val *= eta;
      return val;
    }

#else // LIBMESH_DIM == 1
  libmesh_ignore(i, p);
  libmesh_assert(order);
  libmesh_not_implemented();
#endif
}



template <typename RealType>
typename FEShim<2,MONOMIAL,RealType>::OutputShape FEShim<2,MONOMIAL,RealType>::shape(const Elem * elem,
                           const Order order,
                           const unsigned int i,
                           const Point & p,
                           const bool add_p_level)
{
  libmesh_assert(elem);

  // by default call the orientation-independent shape functions
  return FEShim<2,MONOMIAL,RealType>::shape(elem->type(), static_cast<Order>(order + add_p_level * elem->p_level()), i, p);
}



template <typename RealType>
typename FEShim<2,MONOMIAL,RealType>::OutputShape FEShim<2,MONOMIAL,RealType>::shape_deriv(const ElemType,
                                 const Order libmesh_dbg_var(order),
                                 const unsigned int i,
                                 const unsigned int j,
                                 const Point & p)
{
#if LIBMESH_DIM > 1


  libmesh_assert_less (j, 2);

  libmesh_assert_less (i, (static_cast<unsigned int>(order)+1)*
                       (static_cast<unsigned int>(order)+2)/2);

  const Real xi  = p(0);
  const Real eta = p(1);

  // monomials. since they are hierarchic we only need one case block.

  switch (j)
    {
      // d()/dxi
    case 0:
      {
        switch (i)
          {
            // constants
          case 0:
            return 0.;

            // linears
          case 1:
            return 1.;

          case 2:
            return 0.;

            // quadratics
          case 3:
            return 2.*xi;

          case 4:
            return eta;

          case 5:
            return 0.;

            // cubics
          case 6:
            return 3.*xi*xi;

          case 7:
            return 2.*xi*eta;

          case 8:
            return eta*eta;

          case 9:
            return 0.;

            // quartics
          case 10:
            return 4.*xi*xi*xi;

          case 11:
            return 3.*xi*xi*eta;

          case 12:
            return 2.*xi*eta*eta;

          case 13:
            return eta*eta*eta;

          case 14:
            return 0.;

          default:
            unsigned int o = 0;
            for (; i >= (o+1)*(o+2)/2; o++) { }
            unsigned int ny = i - (o*(o+1)/2);
            unsigned int nx = o - ny;
            Real val = nx;
            for (unsigned int index=1; index < nx; index++)
              val *= xi;
            for (unsigned int index=0; index != ny; index++)
              val *= eta;
            return val;
          }
      }


      // d()/deta
    case 1:
      {
        switch (i)
          {
            // constants
          case 0:
            return 0.;

            // linears
          case 1:
            return 0.;

          case 2:
            return 1.;

            // quadratics
          case 3:
            return 0.;

          case 4:
            return xi;

          case 5:
            return 2.*eta;

            // cubics
          case 6:
            return 0.;

          case 7:
            return xi*xi;

          case 8:
            return 2.*xi*eta;

          case 9:
            return 3.*eta*eta;

            // quartics
          case 10:
            return 0.;

          case 11:
            return xi*xi*xi;

          case 12:
            return 2.*xi*xi*eta;

          case 13:
            return 3.*xi*eta*eta;

          case 14:
            return 4.*eta*eta*eta;

          default:
            unsigned int o = 0;
            for (; i >= (o+1)*(o+2)/2; o++) { }
            unsigned int ny = i - (o*(o+1)/2);
            unsigned int nx = o - ny;
            Real val = ny;
            for (unsigned int index=0; index != nx; index++)
              val *= xi;
            for (unsigned int index=1; index < ny; index++)
              val *= eta;
            return val;
          }
      }

    default:
      libmesh_error_msg("Invalid shape function derivative j = " << j);
    }

#else // LIBMESH_DIM == 1
  libmesh_ignore(i, j, p);
  libmesh_assert(order);
  libmesh_not_implemented();
#endif
}



template <typename RealType>
typename FEShim<2,MONOMIAL,RealType>::OutputShape FEShim<2,MONOMIAL,RealType>::shape_deriv(const Elem * elem,
                                 const Order order,
                                 const unsigned int i,
                                 const unsigned int j,
                                 const Point & p,
                                 const bool add_p_level)
{
  libmesh_assert(elem);

  // by default call the orientation-independent shape functions
  return FEShim<2,MONOMIAL,RealType>::shape_deriv(elem->type(), static_cast<Order>(order + add_p_level * elem->p_level()), i, j, p);
}


#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

template <typename RealType>
typename FEShim<2,MONOMIAL,RealType>::OutputShape FEShim<2,MONOMIAL,RealType>::shape_second_deriv(const ElemType,
                                        const Order libmesh_dbg_var(order),
                                        const unsigned int i,
                                        const unsigned int j,
                                        const Point & p)
{
#if LIBMESH_DIM > 1


  libmesh_assert_less_equal (j, 2);

  libmesh_assert_less (i, (static_cast<unsigned int>(order)+1)*
                       (static_cast<unsigned int>(order)+2)/2);

  const Real xi  = p(0);
  const Real eta = p(1);

  // monomials. since they are hierarchic we only need one case block.

  switch (j)
    {
      // d^2()/dxi^2
    case 0:
      {
        switch (i)
          {
            // constants
          case 0:
            // linears
          case 1:
          case 2:
            return 0.;

            // quadratics
          case 3:
            return 2.;

          case 4:
          case 5:
            return 0.;

            // cubics
          case 6:
            return 6.*xi;

          case 7:
            return 2.*eta;

          case 8:
          case 9:
            return 0.;

            // quartics
          case 10:
            return 12.*xi*xi;

          case 11:
            return 6.*xi*eta;

          case 12:
            return 2.*eta*eta;

          case 13:
          case 14:
            return 0.;

          default:
            unsigned int o = 0;
            for (; i >= (o+1)*(o+2)/2; o++) { }
            unsigned int ny = i - (o*(o+1)/2);
            unsigned int nx = o - ny;
            Real val = nx * (nx - 1);
            for (unsigned int index=2; index < nx; index++)
              val *= xi;
            for (unsigned int index=0; index != ny; index++)
              val *= eta;
            return val;
          }
      }

      // d^2()/dxideta
    case 1:
      {
        switch (i)
          {
            // constants
          case 0:

            // linears
          case 1:
          case 2:
            return 0.;

            // quadratics
          case 3:
            return 0.;

          case 4:
            return 1.;

          case 5:
            return 0.;

            // cubics
          case 6:
            return 0.;
          case 7:
            return 2.*xi;

          case 8:
            return 2.*eta;

          case 9:
            return 0.;

            // quartics
          case 10:
            return 0.;

          case 11:
            return 3.*xi*xi;

          case 12:
            return 4.*xi*eta;

          case 13:
            return 3.*eta*eta;

          case 14:
            return 0.;

          default:
            unsigned int o = 0;
            for (; i >= (o+1)*(o+2)/2; o++) { }
            unsigned int ny = i - (o*(o+1)/2);
            unsigned int nx = o - ny;
            Real val = nx * ny;
            for (unsigned int index=1; index < nx; index++)
              val *= xi;
            for (unsigned int index=1; index < ny; index++)
              val *= eta;
            return val;
          }
      }

      // d^2()/deta^2
    case 2:
      {
        switch (i)
          {
            // constants
          case 0:

            // linears
          case 1:
          case 2:
            return 0.;

            // quadratics
          case 3:
          case 4:
            return 0.;

          case 5:
            return 2.;

            // cubics
          case 6:
            return 0.;

          case 7:
            return 0.;

          case 8:
            return 2.*xi;

          case 9:
            return 6.*eta;

            // quartics
          case 10:
          case 11:
            return 0.;

          case 12:
            return 2.*xi*xi;

          case 13:
            return 6.*xi*eta;

          case 14:
            return 12.*eta*eta;

          default:
            unsigned int o = 0;
            for (; i >= (o+1)*(o+2)/2; o++) { }
            unsigned int ny = i - (o*(o+1)/2);
            unsigned int nx = o - ny;
            Real val = ny * (ny - 1);
            for (unsigned int index=0; index != nx; index++)
              val *= xi;
            for (unsigned int index=2; index < ny; index++)
              val *= eta;
            return val;
          }
      }

    default:
      libmesh_error_msg("Invalid shape function derivative j = " << j);
    }

#else // LIBMESH_DIM == 1
  libmesh_assert(order);
  libmesh_ignore(i, j, p);
  libmesh_not_implemented();
#endif
}



template <typename RealType>
typename FEShim<2,MONOMIAL,RealType>::OutputShape FEShim<2,MONOMIAL,RealType>::shape_second_deriv(const Elem * elem,
                                        const Order order,
                                        const unsigned int i,
                                        const unsigned int j,
                                        const Point & p,
                                        const bool add_p_level)
{
  libmesh_assert(elem);

  // by default call the orientation-independent shape functions
  return FEShim<2,MONOMIAL,RealType>::shape_second_deriv(elem->type(), static_cast<Order>(order + add_p_level * elem->p_level()), i, j, p);
}

#endif

} // namespace libMesh

#endif // LIBMESH_FE_MONOMIAL_SHAPE_2D_IMPL_H
