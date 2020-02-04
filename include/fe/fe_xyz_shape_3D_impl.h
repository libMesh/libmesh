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

#ifndef LIBMESH_FE_XYZ_SHAPE_3D_IMPL_H
#define LIBMESH_FE_XYZ_SHAPE_3D_IMPL_H

// C++ includes

// Local includes
#include "libmesh/fe.h"
#include "libmesh/elem.h"



namespace libMesh
{


template <typename RealType>
typename FEShim<3,XYZ,RealType>::OutputShape FEShim<3,XYZ,RealType>::shape(const ElemType,
                      const Order,
                      const unsigned int,
                      const Point &)
{
  libmesh_error_msg("XYZ polynomials require the element because the centroid is needed.");
  return 0.;
}



template <typename RealType>
typename FEShim<3,XYZ,RealType>::OutputShape FEShim<3,XYZ,RealType>::shape(const Elem * elem,
                      const Order libmesh_dbg_var(order),
                      const unsigned int i,
                      const Point & point_in,
                      const bool libmesh_dbg_var(add_p_level))
{
#if LIBMESH_DIM == 3
  libmesh_assert(elem);

  Point centroid = elem->centroid();
  Point max_distance = Point(0.,0.,0.);
  for (unsigned int p = 0; p < elem->n_nodes(); p++)
    for (unsigned int d = 0; d < 3; d++)
      {
        const Real distance = std::abs(centroid(d) - elem->point(p)(d));
        max_distance(d) = std::max(distance, max_distance(d));
      }

  const Real x  = point_in(0);
  const Real y  = point_in(1);
  const Real z  = point_in(2);
  const Real xc = centroid(0);
  const Real yc = centroid(1);
  const Real zc = centroid(2);
  const Real distx = max_distance(0);
  const Real disty = max_distance(1);
  const Real distz = max_distance(2);
  const Real dx = (x - xc)/distx;
  const Real dy = (y - yc)/disty;
  const Real dz = (z - zc)/distz;

#ifndef NDEBUG
  // totalorder is only used in the assertion below, so
  // we avoid declaring it when asserts are not active.
  const unsigned int totalorder = order + add_p_level * elem->p_level();
#endif
  libmesh_assert_less (i, (totalorder+1) * (totalorder+2) *
                       (totalorder+3)/6);

  // monomials. since they are hierarchic we only need one case block.
  switch (i)
    {
      // constant
    case 0:
      return 1.;

      // linears
    case 1:
      return dx;

    case 2:
      return dy;

    case 3:
      return dz;

      // quadratics
    case 4:
      return dx*dx;

    case 5:
      return dx*dy;

    case 6:
      return dy*dy;

    case 7:
      return dx*dz;

    case 8:
      return dz*dy;

    case 9:
      return dz*dz;

      // cubics
    case 10:
      return dx*dx*dx;

    case 11:
      return dx*dx*dy;

    case 12:
      return dx*dy*dy;

    case 13:
      return dy*dy*dy;

    case 14:
      return dx*dx*dz;

    case 15:
      return dx*dy*dz;

    case 16:
      return dy*dy*dz;

    case 17:
      return dx*dz*dz;

    case 18:
      return dy*dz*dz;

    case 19:
      return dz*dz*dz;

      // quartics
    case 20:
      return dx*dx*dx*dx;

    case 21:
      return dx*dx*dx*dy;

    case 22:
      return dx*dx*dy*dy;

    case 23:
      return dx*dy*dy*dy;

    case 24:
      return dy*dy*dy*dy;

    case 25:
      return dx*dx*dx*dz;

    case 26:
      return dx*dx*dy*dz;

    case 27:
      return dx*dy*dy*dz;

    case 28:
      return dy*dy*dy*dz;

    case 29:
      return dx*dx*dz*dz;

    case 30:
      return dx*dy*dz*dz;

    case 31:
      return dy*dy*dz*dz;

    case 32:
      return dx*dz*dz*dz;

    case 33:
      return dy*dz*dz*dz;

    case 34:
      return dz*dz*dz*dz;

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
        val *= dx;
      for (unsigned int index=0; index != ny; index++)
        val *= dy;
      for (unsigned int index=0; index != nz; index++)
        val *= dz;
      return val;
    }

#else // LIBMESH_DIM != 3
  libmesh_assert(true || order || add_p_level);
  libmesh_ignore(elem, i, point_in);
  libmesh_not_implemented();
#endif
}



template <typename RealType>
typename FEShim<3,XYZ,RealType>::OutputShape FEShim<3,XYZ,RealType>::shape_deriv(const ElemType,
                            const Order,
                            const unsigned int,
                            const unsigned int,
                            const Point &)
{
  libmesh_error_msg("XYZ polynomials require the element \nbecause the centroid is needed.");
  return 0.;
}



template <typename RealType>
typename FEShim<3,XYZ,RealType>::OutputShape FEShim<3,XYZ,RealType>::shape_deriv(const Elem * elem,
                            const Order libmesh_dbg_var(order),
                            const unsigned int i,
                            const unsigned int j,
                            const Point & point_in,
                            const bool libmesh_dbg_var(add_p_level))
{
#if LIBMESH_DIM == 3

  libmesh_assert(elem);
  libmesh_assert_less (j, 3);

  Point centroid = elem->centroid();
  Point max_distance = Point(0.,0.,0.);
  for (unsigned int p = 0; p < elem->n_nodes(); p++)
    for (unsigned int d = 0; d < 3; d++)
      {
        const Real distance = std::abs(centroid(d) - elem->point(p)(d));
        max_distance(d) = std::max(distance, max_distance(d));
      }

  const Real x  = point_in(0);
  const Real y  = point_in(1);
  const Real z  = point_in(2);
  const Real xc = centroid(0);
  const Real yc = centroid(1);
  const Real zc = centroid(2);
  const Real distx = max_distance(0);
  const Real disty = max_distance(1);
  const Real distz = max_distance(2);
  const Real dx = (x - xc)/distx;
  const Real dy = (y - yc)/disty;
  const Real dz = (z - zc)/distz;

#ifndef NDEBUG
  // totalorder is only used in the assertion below, so
  // we avoid declaring it when asserts are not active.
  const unsigned int totalorder = static_cast<Order>(order + add_p_level * elem->p_level());
#endif
  libmesh_assert_less (i, (totalorder+1) * (totalorder+2) *
                       (totalorder+3)/6);

  switch (j)
    {
      // d()/dx
    case 0:
      {
        switch (i)
          {
            // constant
          case 0:
            return 0.;

            // linear
          case 1:
            return 1./distx;

          case 2:
            return 0.;

          case 3:
            return 0.;

            // quadratic
          case 4:
            return 2.*dx/distx;

          case 5:
            return dy/distx;

          case 6:
            return 0.;

          case 7:
            return dz/distx;

          case 8:
            return 0.;

          case 9:
            return 0.;

            // cubic
          case 10:
            return 3.*dx*dx/distx;

          case 11:
            return 2.*dx*dy/distx;

          case 12:
            return dy*dy/distx;

          case 13:
            return 0.;

          case 14:
            return 2.*dx*dz/distx;

          case 15:
            return dy*dz/distx;

          case 16:
            return 0.;

          case 17:
            return dz*dz/distx;

          case 18:
            return 0.;

          case 19:
            return 0.;

            // quartics
          case 20:
            return 4.*dx*dx*dx/distx;

          case 21:
            return 3.*dx*dx*dy/distx;

          case 22:
            return 2.*dx*dy*dy/distx;

          case 23:
            return dy*dy*dy/distx;

          case 24:
            return 0.;

          case 25:
            return 3.*dx*dx*dz/distx;

          case 26:
            return 2.*dx*dy*dz/distx;

          case 27:
            return dy*dy*dz/distx;

          case 28:
            return 0.;

          case 29:
            return 2.*dx*dz*dz/distx;

          case 30:
            return dy*dz*dz/distx;

          case 31:
            return 0.;

          case 32:
            return dz*dz*dz/distx;

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
              val *= dx;
            for (unsigned int index=0; index != ny; index++)
              val *= dy;
            for (unsigned int index=0; index != nz; index++)
              val *= dz;
            return val/distx;
          }
      }


      // d()/dy
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
            return 1./disty;

          case 3:
            return 0.;

            // quadratic
          case 4:
            return 0.;

          case 5:
            return dx/disty;

          case 6:
            return 2.*dy/disty;

          case 7:
            return 0.;

          case 8:
            return dz/disty;

          case 9:
            return 0.;

            // cubic
          case 10:
            return 0.;

          case 11:
            return dx*dx/disty;

          case 12:
            return 2.*dx*dy/disty;

          case 13:
            return 3.*dy*dy/disty;

          case 14:
            return 0.;

          case 15:
            return dx*dz/disty;

          case 16:
            return 2.*dy*dz/disty;

          case 17:
            return 0.;

          case 18:
            return dz*dz/disty;

          case 19:
            return 0.;

            // quartics
          case 20:
            return 0.;

          case 21:
            return dx*dx*dx/disty;

          case 22:
            return 2.*dx*dx*dy/disty;

          case 23:
            return 3.*dx*dy*dy/disty;

          case 24:
            return 4.*dy*dy*dy/disty;

          case 25:
            return 0.;

          case 26:
            return dx*dx*dz/disty;

          case 27:
            return 2.*dx*dy*dz/disty;

          case 28:
            return 3.*dy*dy*dz/disty;

          case 29:
            return 0.;

          case 30:
            return dx*dz*dz/disty;

          case 31:
            return 2.*dy*dz*dz/disty;

          case 32:
            return 0.;

          case 33:
            return dz*dz*dz/disty;

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
              val *= dx;
            for (unsigned int index=1; index < ny; index++)
              val *= dy;
            for (unsigned int index=0; index != nz; index++)
              val *= dz;
            return val/disty;
          }
      }


      // d()/dz
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
            return 1./distz;

            // quadratic
          case 4:
            return 0.;

          case 5:
            return 0.;

          case 6:
            return 0.;

          case 7:
            return dx/distz;

          case 8:
            return dy/distz;

          case 9:
            return 2.*dz/distz;

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
            return dx*dx/distz;

          case 15:
            return dx*dy/distz;

          case 16:
            return dy*dy/distz;

          case 17:
            return 2.*dx*dz/distz;

          case 18:
            return 2.*dy*dz/distz;

          case 19:
            return 3.*dz*dz/distz;

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
            return dx*dx*dx/distz;

          case 26:
            return dx*dx*dy/distz;

          case 27:
            return dx*dy*dy/distz;

          case 28:
            return dy*dy*dy/distz;

          case 29:
            return 2.*dx*dx*dz/distz;

          case 30:
            return 2.*dx*dy*dz/distz;

          case 31:
            return 2.*dy*dy*dz/distz;

          case 32:
            return 3.*dx*dz*dz/distz;

          case 33:
            return 3.*dy*dz*dz/distz;

          case 34:
            return 4.*dz*dz*dz/distz;

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
              val *= dx;
            for (unsigned int index=0; index != ny; index++)
              val *= dy;
            for (unsigned int index=1; index < nz; index++)
              val *= dz;
            return val/distz;
          }
      }


    default:
      libmesh_error_msg("Invalid j = " << j);
    }

#else // LIBMESH_DIM != 3
  libmesh_assert(true || order || add_p_level);
  libmesh_ignore(elem, i, j, point_in);
  libmesh_not_implemented();
#endif
}


#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

template <typename RealType>
typename FEShim<3,XYZ,RealType>::OutputShape FEShim<3,XYZ,RealType>::shape_second_deriv(const ElemType,
                                   const Order,
                                   const unsigned int,
                                   const unsigned int,
                                   const Point &)
{
  libmesh_error_msg("XYZ polynomials require the element \nbecause the centroid is needed.");
  return 0.;
}



template <typename RealType>
typename FEShim<3,XYZ,RealType>::OutputShape FEShim<3,XYZ,RealType>::shape_second_deriv(const Elem * elem,
                                   const Order libmesh_dbg_var(order),
                                   const unsigned int i,
                                   const unsigned int j,
                                   const Point & point_in,
                                   const bool libmesh_dbg_var(add_p_level))
{
#if LIBMESH_DIM == 3

  libmesh_assert(elem);
  libmesh_assert_less (j, 6);

  Point centroid = elem->centroid();
  Point max_distance = Point(0.,0.,0.);
  for (const Point & p : elem->node_ref_range())
    for (unsigned int d = 0; d < 3; d++)
      {
        const Real distance = std::abs(centroid(d) - p(d));
        max_distance(d) = std::max(distance, max_distance(d));
      }

  const Real x  = point_in(0);
  const Real y  = point_in(1);
  const Real z  = point_in(2);
  const Real xc = centroid(0);
  const Real yc = centroid(1);
  const Real zc = centroid(2);
  const Real distx = max_distance(0);
  const Real disty = max_distance(1);
  const Real distz = max_distance(2);
  const Real dx = (x - xc)/distx;
  const Real dy = (y - yc)/disty;
  const Real dz = (z - zc)/distz;
  const Real dist2x = pow(distx,2.);
  const Real dist2y = pow(disty,2.);
  const Real dist2z = pow(distz,2.);
  const Real distxy = distx * disty;
  const Real distxz = distx * distz;
  const Real distyz = disty * distz;

#ifndef NDEBUG
  // totalorder is only used in the assertion below, so
  // we avoid declaring it when asserts are not active.
  const unsigned int totalorder = static_cast<Order>(order + add_p_level * elem->p_level());
#endif
  libmesh_assert_less (i, (totalorder+1) * (totalorder+2) *
                       (totalorder+3)/6);

  // monomials. since they are hierarchic we only need one case block.
  switch (j)
    {
      // d^2()/dx^2
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
            return 2./dist2x;

          case 5:
          case 6:
          case 7:
          case 8:
          case 9:
            return 0.;

            // cubic
          case 10:
            return 6.*dx/dist2x;

          case 11:
            return 2.*dy/dist2x;

          case 12:
          case 13:
            return 0.;

          case 14:
            return 2.*dz/dist2x;

          case 15:
          case 16:
          case 17:
          case 18:
          case 19:
            return 0.;

            // quartics
          case 20:
            return 12.*dx*dx/dist2x;

          case 21:
            return 6.*dx*dy/dist2x;

          case 22:
            return 2.*dy*dy/dist2x;

          case 23:
          case 24:
            return 0.;

          case 25:
            return 6.*dx*dz/dist2x;

          case 26:
            return 2.*dy*dz/dist2x;

          case 27:
          case 28:
            return 0.;

          case 29:
            return 2.*dz*dz/dist2x;

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
              val *= dx;
            for (unsigned int index=0; index != ny; index++)
              val *= dy;
            for (unsigned int index=0; index != nz; index++)
              val *= dz;
            return val/dist2x;
          }
      }


      // d^2()/dxdy
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
            return 1./distxy;

          case 6:
          case 7:
          case 8:
          case 9:
            return 0.;

            // cubic
          case 10:
            return 0.;

          case 11:
            return 2.*dx/distxy;

          case 12:
            return 2.*dy/distxy;

          case 13:
          case 14:
            return 0.;

          case 15:
            return dz/distxy;

          case 16:
          case 17:
          case 18:
          case 19:
            return 0.;

            // quartics
          case 20:
            return 0.;

          case 21:
            return 3.*dx*dx/distxy;

          case 22:
            return 4.*dx*dy/distxy;

          case 23:
            return 3.*dy*dy/distxy;

          case 24:
          case 25:
            return 0.;

          case 26:
            return 2.*dx*dz/distxy;

          case 27:
            return 2.*dy*dz/distxy;

          case 28:
          case 29:
            return 0.;

          case 30:
            return dz*dz/distxy;

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
              val *= dx;
            for (unsigned int index=1; index < ny; index++)
              val *= dy;
            for (unsigned int index=0; index != nz; index++)
              val *= dz;
            return val/distxy;
          }
      }


      // d^2()/dy^2
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
            return 2./dist2y;

          case 7:
          case 8:
          case 9:
            return 0.;

            // cubic
          case 10:
          case 11:
            return 0.;

          case 12:
            return 2.*dx/dist2y;
          case 13:
            return 6.*dy/dist2y;

          case 14:
          case 15:
            return 0.;

          case 16:
            return 2.*dz/dist2y;

          case 17:
          case 18:
          case 19:
            return 0.;

            // quartics
          case 20:
          case 21:
            return 0.;

          case 22:
            return 2.*dx*dx/dist2y;

          case 23:
            return 6.*dx*dy/dist2y;

          case 24:
            return 12.*dy*dy/dist2y;

          case 25:
          case 26:
            return 0.;

          case 27:
            return 2.*dx*dz/dist2y;

          case 28:
            return 6.*dy*dz/dist2y;

          case 29:
          case 30:
            return 0.;

          case 31:
            return 2.*dz*dz/dist2y;

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
              val *= dx;
            for (unsigned int index=2; index < ny; index++)
              val *= dy;
            for (unsigned int index=0; index != nz; index++)
              val *= dz;
            return val/dist2y;
          }
      }


      // d^2()/dxdz
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
            return 1./distxz;

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
            return 2.*dx/distxz;

          case 15:
            return dy/distxz;

          case 16:
            return 0.;

          case 17:
            return 2.*dz/distxz;

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
            return 3.*dx*dx/distxz;

          case 26:
            return 2.*dx*dy/distxz;

          case 27:
            return dy*dy/distxz;

          case 28:
            return 0.;

          case 29:
            return 4.*dx*dz/distxz;

          case 30:
            return 2.*dy*dz/distxz;

          case 31:
            return 0.;

          case 32:
            return 3.*dz*dz/distxz;

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
              val *= dx;
            for (unsigned int index=0; index != ny; index++)
              val *= dy;
            for (unsigned int index=1; index < nz; index++)
              val *= dz;
            return val/distxz;
          }
      }

      // d^2()/dydz
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
            return 1./distyz;

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
            return dx/distyz;

          case 16:
            return 2.*dy/distyz;

          case 17:
            return 0.;

          case 18:
            return 2.*dz/distyz;

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
            return dx*dx/distyz;

          case 27:
            return 2.*dx*dy/distyz;

          case 28:
            return 3.*dy*dy/distyz;

          case 29:
            return 0.;

          case 30:
            return 2.*dx*dz/distyz;

          case 31:
            return 4.*dy*dz/distyz;

          case 32:
            return 0.;

          case 33:
            return 3.*dz*dz/distyz;

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
              val *= dx;
            for (unsigned int index=1; index < ny; index++)
              val *= dy;
            for (unsigned int index=1; index < nz; index++)
              val *= dz;
            return val/distyz;
          }
      }


      // d^2()/dz^2
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
            return 2./dist2z;

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
            return 2.*dx/dist2z;

          case 18:
            return 2.*dy/dist2z;

          case 19:
            return 6.*dz/dist2z;

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
            return 2.*dx*dx/dist2z;

          case 30:
            return 2.*dx*dy/dist2z;

          case 31:
            return 2.*dy*dy/dist2z;

          case 32:
            return 6.*dx*dz/dist2z;

          case 33:
            return 6.*dy*dz/dist2z;

          case 34:
            return 12.*dz*dz/dist2z;

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
              val *= dx;
            for (unsigned int index=0; index != ny; index++)
              val *= dy;
            for (unsigned int index=2; index < nz; index++)
              val *= dz;
            return val/dist2z;
          }
      }


    default:
      libmesh_error_msg("Invalid j = " << j);
    }

#else // LIBMESH_DIM != 3
  libmesh_assert(true || order || add_p_level);
  libmesh_ignore(elem, i, j, point_in);
  libmesh_not_implemented();
#endif
}

#endif

} // namespace libMesh

#endif // LIBMESH_FE_XYZ_SHAPE_3D_IMPL_H
