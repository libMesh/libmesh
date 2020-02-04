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

#ifndef LIBMESH_FE_MONOMIAL_SHAPE_3D_IMPL_H
#define LIBMESH_FE_MONOMIAL_SHAPE_3D_IMPL_H

// C++ includes

// Local includes
#include "libmesh/fe.h"
#include "libmesh/elem.h"

namespace libMesh
{




template <typename RealType>
typename FEShim<3,MONOMIAL,RealType>::OutputShape FEShim<3,MONOMIAL,RealType>::shape(const ElemType,
                           const Order libmesh_dbg_var(order),
                           const unsigned int i,
                           const Point & p)
{
#if LIBMESH_DIM == 3

  const Real xi   = p(0);
  const Real eta  = p(1);
  const Real zeta = p(2);

  libmesh_assert_less (i, (static_cast<unsigned int>(order)+1)*
                       (static_cast<unsigned int>(order)+2)*
                       (static_cast<unsigned int>(order)+3)/6);

  // monomials. since they are hierarchic we only need one case block.
  switch (i)
    {
      // constant
    case 0:
      return 1.;

      // linears
    case 1:
      return xi;

    case 2:
      return eta;

    case 3:
      return zeta;

      // quadratics
    case 4:
      return xi*xi;

    case 5:
      return xi*eta;

    case 6:
      return eta*eta;

    case 7:
      return xi*zeta;

    case 8:
      return zeta*eta;

    case 9:
      return zeta*zeta;

      // cubics
    case 10:
      return xi*xi*xi;

    case 11:
      return xi*xi*eta;

    case 12:
      return xi*eta*eta;

    case 13:
      return eta*eta*eta;

    case 14:
      return xi*xi*zeta;

    case 15:
      return xi*eta*zeta;

    case 16:
      return eta*eta*zeta;

    case 17:
      return xi*zeta*zeta;

    case 18:
      return eta*zeta*zeta;

    case 19:
      return zeta*zeta*zeta;

      // quartics
    case 20:
      return xi*xi*xi*xi;

    case 21:
      return xi*xi*xi*eta;

    case 22:
      return xi*xi*eta*eta;

    case 23:
      return xi*eta*eta*eta;

    case 24:
      return eta*eta*eta*eta;

    case 25:
      return xi*xi*xi*zeta;

    case 26:
      return xi*xi*eta*zeta;

    case 27:
      return xi*eta*eta*zeta;

    case 28:
      return eta*eta*eta*zeta;

    case 29:
      return xi*xi*zeta*zeta;

    case 30:
      return xi*eta*zeta*zeta;

    case 31:
      return eta*eta*zeta*zeta;

    case 32:
      return xi*zeta*zeta*zeta;

    case 33:
      return eta*zeta*zeta*zeta;

    case 34:
      return zeta*zeta*zeta*zeta;

    default:
      unsigned int o = 0;
      for (; i >= (o+1)*(o+2)*(o+3)/6; o++) { }
      unsigned int i2 = i - (o*(o+1)*(o+2)/6);
      unsigned int block=o, nz = 0;
      for (; block < i2; block += (o-nz+1)) { nz++; }
      const unsigned int nx = block - i2;
      const unsigned int ny = o - nx - nz;
      Real val = 1.;
      for (unsigned int index=0; index != nx; index++)
        val *= xi;
      for (unsigned int index=0; index != ny; index++)
        val *= eta;
      for (unsigned int index=0; index != nz; index++)
        val *= zeta;
      return val;
    }

#else // LIBMESH_DIM != 3
  libmesh_assert(order);
  libmesh_ignore(i, p);
  libmesh_not_implemented();
#endif
}



template <typename RealType>
typename FEShim<3,MONOMIAL,RealType>::OutputShape FEShim<3,MONOMIAL,RealType>::shape(const Elem * elem,
                           const Order order,
                           const unsigned int i,
                           const Point & p,
                           const bool add_p_level)
{
  libmesh_assert(elem);

  // call the orientation-independent shape functions
  return FEShim<3,MONOMIAL,RealType>::shape(elem->type(), static_cast<Order>(order + add_p_level * elem->p_level()), i, p);
}



template <typename RealType>
typename FEShim<3,MONOMIAL,RealType>::OutputShape FEShim<3,MONOMIAL,RealType>::shape_deriv(const ElemType,
                                 const Order libmesh_dbg_var(order),
                                 const unsigned int i,
                                 const unsigned int j,
                                 const Point & p)
{
#if LIBMESH_DIM == 3

  libmesh_assert_less (j, 3);

  libmesh_assert_less (i, (static_cast<unsigned int>(order)+1)*
                       (static_cast<unsigned int>(order)+2)*
                       (static_cast<unsigned int>(order)+3)/6);


  const Real xi   = p(0);
  const Real eta  = p(1);
  const Real zeta = p(2);

  // monomials. since they are hierarchic we only need one case block.
  switch (j)
    {
      // d()/dxi
    case 0:
      {
        switch (i)
          {
            // constant
          case 0:
            return 0.;

            // linear
          case 1:
            return 1.;

          case 2:
            return 0.;

          case 3:
            return 0.;

            // quadratic
          case 4:
            return 2.*xi;

          case 5:
            return eta;

          case 6:
            return 0.;

          case 7:
            return zeta;

          case 8:
            return 0.;

          case 9:
            return 0.;

            // cubic
          case 10:
            return 3.*xi*xi;

          case 11:
            return 2.*xi*eta;

          case 12:
            return eta*eta;

          case 13:
            return 0.;

          case 14:
            return 2.*xi*zeta;

          case 15:
            return eta*zeta;

          case 16:
            return 0.;

          case 17:
            return zeta*zeta;

          case 18:
            return 0.;

          case 19:
            return 0.;

            // quartics
          case 20:
            return 4.*xi*xi*xi;

          case 21:
            return 3.*xi*xi*eta;

          case 22:
            return 2.*xi*eta*eta;

          case 23:
            return eta*eta*eta;

          case 24:
            return 0.;

          case 25:
            return 3.*xi*xi*zeta;

          case 26:
            return 2.*xi*eta*zeta;

          case 27:
            return eta*eta*zeta;

          case 28:
            return 0.;

          case 29:
            return 2.*xi*zeta*zeta;

          case 30:
            return eta*zeta*zeta;

          case 31:
            return 0.;

          case 32:
            return zeta*zeta*zeta;

          case 33:
            return 0.;

          case 34:
            return 0.;

          default:
            unsigned int o = 0;
            for (; i >= (o+1)*(o+2)*(o+3)/6; o++) { }
            unsigned int i2 = i - (o*(o+1)*(o+2)/6);
            unsigned int block=o, nz = 0;
            for (; block < i2; block += (o-nz+1)) { nz++; }
            const unsigned int nx = block - i2;
            const unsigned int ny = o - nx - nz;
            Real val = nx;
            for (unsigned int index=1; index < nx; index++)
              val *= xi;
            for (unsigned int index=0; index != ny; index++)
              val *= eta;
            for (unsigned int index=0; index != nz; index++)
              val *= zeta;
            return val;
          }
      }


      // d()/deta
    case 1:
      {
        switch (i)
          {
            // constant
          case 0:
            return 0.;

            // linear
          case 1:
            return 0.;

          case 2:
            return 1.;

          case 3:
            return 0.;

            // quadratic
          case 4:
            return 0.;

          case 5:
            return xi;

          case 6:
            return 2.*eta;

          case 7:
            return 0.;

          case 8:
            return zeta;

          case 9:
            return 0.;

            // cubic
          case 10:
            return 0.;

          case 11:
            return xi*xi;

          case 12:
            return 2.*xi*eta;

          case 13:
            return 3.*eta*eta;

          case 14:
            return 0.;

          case 15:
            return xi*zeta;

          case 16:
            return 2.*eta*zeta;

          case 17:
            return 0.;

          case 18:
            return zeta*zeta;

          case 19:
            return 0.;

            // quartics
          case 20:
            return 0.;

          case 21:
            return xi*xi*xi;

          case 22:
            return 2.*xi*xi*eta;

          case 23:
            return 3.*xi*eta*eta;

          case 24:
            return 4.*eta*eta*eta;

          case 25:
            return 0.;

          case 26:
            return xi*xi*zeta;

          case 27:
            return 2.*xi*eta*zeta;

          case 28:
            return 3.*eta*eta*zeta;

          case 29:
            return 0.;

          case 30:
            return xi*zeta*zeta;

          case 31:
            return 2.*eta*zeta*zeta;

          case 32:
            return 0.;

          case 33:
            return zeta*zeta*zeta;

          case 34:
            return 0.;

          default:
            unsigned int o = 0;
            for (; i >= (o+1)*(o+2)*(o+3)/6; o++) { }
            unsigned int i2 = i - (o*(o+1)*(o+2)/6);
            unsigned int block=o, nz = 0;
            for (; block < i2; block += (o-nz+1)) { nz++; }
            const unsigned int nx = block - i2;
            const unsigned int ny = o - nx - nz;
            Real val = ny;
            for (unsigned int index=0; index != nx; index++)
              val *= xi;
            for (unsigned int index=1; index < ny; index++)
              val *= eta;
            for (unsigned int index=0; index != nz; index++)
              val *= zeta;
            return val;
          }
      }


      // d()/dzeta
    case 2:
      {
        switch (i)
          {
            // constant
          case 0:
            return 0.;

            // linear
          case 1:
            return 0.;

          case 2:
            return 0.;

          case 3:
            return 1.;

            // quadratic
          case 4:
            return 0.;

          case 5:
            return 0.;

          case 6:
            return 0.;

          case 7:
            return xi;

          case 8:
            return eta;

          case 9:
            return 2.*zeta;

            // cubic
          case 10:
            return 0.;

          case 11:
            return 0.;

          case 12:
            return 0.;

          case 13:
            return 0.;

          case 14:
            return xi*xi;

          case 15:
            return xi*eta;

          case 16:
            return eta*eta;

          case 17:
            return 2.*xi*zeta;

          case 18:
            return 2.*eta*zeta;

          case 19:
            return 3.*zeta*zeta;

            // quartics
          case 20:
            return 0.;

          case 21:
            return 0.;

          case 22:
            return 0.;

          case 23:
            return 0.;

          case 24:
            return 0.;

          case 25:
            return xi*xi*xi;

          case 26:
            return xi*xi*eta;

          case 27:
            return xi*eta*eta;

          case 28:
            return eta*eta*eta;

          case 29:
            return 2.*xi*xi*zeta;

          case 30:
            return 2.*xi*eta*zeta;

          case 31:
            return 2.*eta*eta*zeta;

          case 32:
            return 3.*xi*zeta*zeta;

          case 33:
            return 3.*eta*zeta*zeta;

          case 34:
            return 4.*zeta*zeta*zeta;

          default:
            unsigned int o = 0;
            for (; i >= (o+1)*(o+2)*(o+3)/6; o++) { }
            unsigned int i2 = i - (o*(o+1)*(o+2)/6);
            unsigned int block=o, nz = 0;
            for (; block < i2; block += (o-nz+1)) { nz++; }
            const unsigned int nx = block - i2;
            const unsigned int ny = o - nx - nz;
            Real val = nz;
            for (unsigned int index=0; index != nx; index++)
              val *= xi;
            for (unsigned int index=0; index != ny; index++)
              val *= eta;
            for (unsigned int index=1; index < nz; index++)
              val *= zeta;
            return val;
          }
      }

    default:
      libmesh_error_msg("Invalid shape function derivative j = " << j);
    }

#else // LIBMESH_DIM != 3
  libmesh_assert(order);
  libmesh_ignore(i, j, p);
  libmesh_not_implemented();
#endif
}



template <typename RealType>
typename FEShim<3,MONOMIAL,RealType>::OutputShape FEShim<3,MONOMIAL,RealType>::shape_deriv(const Elem * elem,
                                 const Order order,
                                 const unsigned int i,
                                 const unsigned int j,
                                 const Point & p,
                                 const bool add_p_level)
{
  libmesh_assert(elem);

  // call the orientation-independent shape function derivatives
  return FEShim<3,MONOMIAL,RealType>::shape_deriv(elem->type(), static_cast<Order>(order + add_p_level * elem->p_level()), i, j, p);
}


#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

template <typename RealType>
typename FEShim<3,MONOMIAL,RealType>::OutputShape FEShim<3,MONOMIAL,RealType>::shape_second_deriv(const ElemType,
                                        const Order libmesh_dbg_var(order),
                                        const unsigned int i,
                                        const unsigned int j,
                                        const Point & p)
{
#if LIBMESH_DIM == 3

  libmesh_assert_less (j, 6);

  libmesh_assert_less (i, (static_cast<unsigned int>(order)+1)*
                       (static_cast<unsigned int>(order)+2)*
                       (static_cast<unsigned int>(order)+3)/6);

  const Real xi   = p(0);
  const Real eta  = p(1);
  const Real zeta = p(2);

  // monomials. since they are hierarchic we only need one case block.
  switch (j)
    {
      // d^2()/dxi^2
    case 0:
      {
        switch (i)
          {
            // constant
          case 0:

            // linear
          case 1:
          case 2:
          case 3:
            return 0.;

            // quadratic
          case 4:
            return 2.;

          case 5:
          case 6:
          case 7:
          case 8:
          case 9:
            return 0.;

            // cubic
          case 10:
            return 6.*xi;

          case 11:
            return 2.*eta;

          case 12:
          case 13:
            return 0.;

          case 14:
            return 2.*zeta;

          case 15:
          case 16:
          case 17:
          case 18:
          case 19:
            return 0.;

            // quartics
          case 20:
            return 12.*xi*xi;

          case 21:
            return 6.*xi*eta;

          case 22:
            return 2.*eta*eta;

          case 23:
          case 24:
            return 0.;

          case 25:
            return 6.*xi*zeta;

          case 26:
            return 2.*eta*zeta;

          case 27:
          case 28:
            return 0.;

          case 29:
            return 2.*zeta*zeta;

          case 30:
          case 31:
          case 32:
          case 33:
          case 34:
            return 0.;

          default:
            unsigned int o = 0;
            for (; i >= (o+1)*(o+2)*(o+3)/6; o++) { }
            unsigned int i2 = i - (o*(o+1)*(o+2)/6);
            unsigned int block=o, nz = 0;
            for (; block < i2; block += (o-nz+1)) { nz++; }
            const unsigned int nx = block - i2;
            const unsigned int ny = o - nx - nz;
            Real val = nx * (nx - 1);
            for (unsigned int index=2; index < nx; index++)
              val *= xi;
            for (unsigned int index=0; index != ny; index++)
              val *= eta;
            for (unsigned int index=0; index != nz; index++)
              val *= zeta;
            return val;
          }
      }


      // d^2()/dxideta
    case 1:
      {
        switch (i)
          {
            // constant
          case 0:

            // linear
          case 1:
          case 2:
          case 3:
            return 0.;

            // quadratic
          case 4:
            return 0.;

          case 5:
            return 1.;

          case 6:
          case 7:
          case 8:
          case 9:
            return 0.;

            // cubic
          case 10:
            return 0.;

          case 11:
            return 2.*xi;

          case 12:
            return 2.*eta;

          case 13:
          case 14:
            return 0.;

          case 15:
            return zeta;

          case 16:
          case 17:
          case 18:
          case 19:
            return 0.;

            // quartics
          case 20:
            return 0.;

          case 21:
            return 3.*xi*xi;

          case 22:
            return 4.*xi*eta;

          case 23:
            return 3.*eta*eta;

          case 24:
          case 25:
            return 0.;

          case 26:
            return 2.*xi*zeta;

          case 27:
            return 2.*eta*zeta;

          case 28:
          case 29:
            return 0.;

          case 30:
            return zeta*zeta;

          case 31:
          case 32:
          case 33:
          case 34:
            return 0.;

          default:
            unsigned int o = 0;
            for (; i >= (o+1)*(o+2)*(o+3)/6; o++) { }
            unsigned int i2 = i - (o*(o+1)*(o+2)/6);
            unsigned int block=o, nz = 0;
            for (; block < i2; block += (o-nz+1)) { nz++; }
            const unsigned int nx = block - i2;
            const unsigned int ny = o - nx - nz;
            Real val = nx * ny;
            for (unsigned int index=1; index < nx; index++)
              val *= xi;
            for (unsigned int index=1; index < ny; index++)
              val *= eta;
            for (unsigned int index=0; index != nz; index++)
              val *= zeta;
            return val;
          }
      }


      // d^2()/deta^2
    case 2:
      {
        switch (i)
          {
            // constant
          case 0:

            // linear
          case 1:
          case 2:
          case 3:
            return 0.;

            // quadratic
          case 4:
          case 5:
            return 0.;

          case 6:
            return 2.;

          case 7:
          case 8:
          case 9:
            return 0.;

            // cubic
          case 10:
          case 11:
            return 0.;

          case 12:
            return 2.*xi;
          case 13:
            return 6.*eta;

          case 14:
          case 15:
            return 0.;

          case 16:
            return 2.*zeta;

          case 17:
          case 18:
          case 19:
            return 0.;

            // quartics
          case 20:
          case 21:
            return 0.;

          case 22:
            return 2.*xi*xi;

          case 23:
            return 6.*xi*eta;

          case 24:
            return 12.*eta*eta;

          case 25:
          case 26:
            return 0.;

          case 27:
            return 2.*xi*zeta;

          case 28:
            return 6.*eta*zeta;

          case 29:
          case 30:
            return 0.;

          case 31:
            return 2.*zeta*zeta;

          case 32:
          case 33:
          case 34:
            return 0.;

          default:
            unsigned int o = 0;
            for (; i >= (o+1)*(o+2)*(o+3)/6; o++) { }
            unsigned int i2 = i - (o*(o+1)*(o+2)/6);
            unsigned int block=o, nz = 0;
            for (; block < i2; block += (o-nz+1)) { nz++; }
            const unsigned int nx = block - i2;
            const unsigned int ny = o - nx - nz;
            Real val = ny * (ny - 1);
            for (unsigned int index=0; index != nx; index++)
              val *= xi;
            for (unsigned int index=2; index < ny; index++)
              val *= eta;
            for (unsigned int index=0; index != nz; index++)
              val *= zeta;
            return val;
          }
      }


      // d^2()/dxidzeta
    case 3:
      {
        switch (i)
          {
            // constant
          case 0:

            // linear
          case 1:
          case 2:
          case 3:
            return 0.;

            // quadratic
          case 4:
          case 5:
          case 6:
            return 0.;

          case 7:
            return 1.;

          case 8:
          case 9:
            return 0.;

            // cubic
          case 10:
          case 11:
          case 12:
          case 13:
            return 0.;

          case 14:
            return 2.*xi;

          case 15:
            return eta;

          case 16:
            return 0.;

          case 17:
            return 2.*zeta;

          case 18:
          case 19:
            return 0.;

            // quartics
          case 20:
          case 21:
          case 22:
          case 23:
          case 24:
            return 0.;

          case 25:
            return 3.*xi*xi;

          case 26:
            return 2.*xi*eta;

          case 27:
            return eta*eta;

          case 28:
            return 0.;

          case 29:
            return 4.*xi*zeta;

          case 30:
            return 2.*eta*zeta;

          case 31:
            return 0.;

          case 32:
            return 3.*zeta*zeta;

          case 33:
          case 34:
            return 0.;

          default:
            unsigned int o = 0;
            for (; i >= (o+1)*(o+2)*(o+3)/6; o++) { }
            unsigned int i2 = i - (o*(o+1)*(o+2)/6);
            unsigned int block=o, nz = 0;
            for (; block < i2; block += (o-nz+1)) { nz++; }
            const unsigned int nx = block - i2;
            const unsigned int ny = o - nx - nz;
            Real val = nx * nz;
            for (unsigned int index=1; index < nx; index++)
              val *= xi;
            for (unsigned int index=0; index != ny; index++)
              val *= eta;
            for (unsigned int index=1; index < nz; index++)
              val *= zeta;
            return val;
          }
      }

      // d^2()/detadzeta
    case 4:
      {
        switch (i)
          {
            // constant
          case 0:

            // linear
          case 1:
          case 2:
          case 3:
            return 0.;

            // quadratic
          case 4:
          case 5:
          case 6:
          case 7:
            return 0.;

          case 8:
            return 1.;

          case 9:
            return 0.;

            // cubic
          case 10:
          case 11:
          case 12:
          case 13:
          case 14:
            return 0.;

          case 15:
            return xi;

          case 16:
            return 2.*eta;

          case 17:
            return 0.;

          case 18:
            return 2.*zeta;

          case 19:
            return 0.;

            // quartics
          case 20:
          case 21:
          case 22:
          case 23:
          case 24:
          case 25:
            return 0.;

          case 26:
            return xi*xi;

          case 27:
            return 2.*xi*eta;

          case 28:
            return 3.*eta*eta;

          case 29:
            return 0.;

          case 30:
            return 2.*xi*zeta;

          case 31:
            return 4.*eta*zeta;

          case 32:
            return 0.;

          case 33:
            return 3.*zeta*zeta;

          case 34:
            return 0.;

          default:
            unsigned int o = 0;
            for (; i >= (o+1)*(o+2)*(o+3)/6; o++) { }
            unsigned int i2 = i - (o*(o+1)*(o+2)/6);
            unsigned int block=o, nz = 0;
            for (; block < i2; block += (o-nz+1)) { nz++; }
            const unsigned int nx = block - i2;
            const unsigned int ny = o - nx - nz;
            Real val = ny * nz;
            for (unsigned int index=0; index != nx; index++)
              val *= xi;
            for (unsigned int index=1; index < ny; index++)
              val *= eta;
            for (unsigned int index=1; index < nz; index++)
              val *= zeta;
            return val;
          }
      }


      // d^2()/dzeta^2
    case 5:
      {
        switch (i)
          {
            // constant
          case 0:

            // linear
          case 1:
          case 2:
          case 3:
            return 0.;

            // quadratic
          case 4:
          case 5:
          case 6:
          case 7:
          case 8:
            return 0.;

          case 9:
            return 2.;

            // cubic
          case 10:
          case 11:
          case 12:
          case 13:
          case 14:
          case 15:
          case 16:
            return 0.;

          case 17:
            return 2.*xi;

          case 18:
            return 2.*eta;

          case 19:
            return 6.*zeta;

            // quartics
          case 20:
          case 21:
          case 22:
          case 23:
          case 24:
          case 25:
          case 26:
          case 27:
          case 28:
            return 0.;

          case 29:
            return 2.*xi*xi;

          case 30:
            return 2.*xi*eta;

          case 31:
            return 2.*eta*eta;

          case 32:
            return 6.*xi*zeta;

          case 33:
            return 6.*eta*zeta;

          case 34:
            return 12.*zeta*zeta;

          default:
            unsigned int o = 0;
            for (; i >= (o+1)*(o+2)*(o+3)/6; o++) { }
            unsigned int i2 = i - (o*(o+1)*(o+2)/6);
            unsigned int block=o, nz = 0;
            for (; block < i2; block += (o-nz+1)) { nz++; }
            const unsigned int nx = block - i2;
            const unsigned int ny = o - nx - nz;
            Real val = nz * (nz - 1);
            for (unsigned int index=0; index != nx; index++)
              val *= xi;
            for (unsigned int index=0; index != ny; index++)
              val *= eta;
            for (unsigned int index=2; index < nz; index++)
              val *= zeta;
            return val;
          }
      }

    default:
      libmesh_error_msg("Invalid j = " << j);
    }

#else // LIBMESH_DIM != 3
  libmesh_assert(order);
  libmesh_ignore(i, j, p);
  libmesh_not_implemented();
#endif
}



template <typename RealType>
typename FEShim<3,MONOMIAL,RealType>::OutputShape FEShim<3,MONOMIAL,RealType>::shape_second_deriv(const Elem * elem,
                                        const Order order,
                                        const unsigned int i,
                                        const unsigned int j,
                                        const Point & p,
                                        const bool add_p_level)
{
  libmesh_assert(elem);

  // call the orientation-independent shape function derivatives
  return FEShim<3,MONOMIAL,RealType>::shape_second_deriv(elem->type(), static_cast<Order>(order + add_p_level * elem->p_level()), i, j, p);
}

#endif

} // namespace libMesh

#endif // LIBMESH_FE_MONOMIAL_SHAPE_3D_IMPL_H
