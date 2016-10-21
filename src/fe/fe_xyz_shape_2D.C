// The libMesh Finite Element Library.
// Copyright (C) 2002-2016 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


// C++ inlcludes

// Local includes
#include "libmesh/fe.h"
#include "libmesh/elem.h"



namespace libMesh
{


template <>
Real FE<2,XYZ>::shape(const ElemType,
                      const Order,
                      const unsigned int,
                      const Point &)
{
  libmesh_error_msg("XYZ polynomials require the element \nbecause the centroid is needed.");
  return 0.;
}



template <>
Real FE<2,XYZ>::shape(const Elem * elem,
                      const Order libmesh_dbg_var(order),
                      const unsigned int i,
                      const Point & point_in)
{
#if LIBMESH_DIM > 1

  libmesh_assert(elem);

  Point centroid = elem->centroid();
  Point max_distance = Point(0.,0.,0.);
  for (unsigned int p = 0; p < elem->n_nodes(); p++)
    for (unsigned int d = 0; d < 2; d++)
      {
        const Real distance = std::abs(centroid(d) - elem->point(p)(d));
        max_distance(d) = std::max(distance, max_distance(d));
      }

  const Real x  = point_in(0);
  const Real y  = point_in(1);
  const Real xc = centroid(0);
  const Real yc = centroid(1);
  const Real distx = max_distance(0);
  const Real disty = max_distance(1);
  const Real dx = (x - xc)/distx;
  const Real dy = (y - yc)/disty;

#ifndef NDEBUG
  // totalorder is only used in the assertion below, so
  // we avoid declaring it when asserts are not active.
  const unsigned int totalorder = order + elem->p_level();
#endif
  libmesh_assert_less (i, (totalorder+1)*(totalorder+2)/2);


  // monomials. since they are hierarchic we only need one case block.
  switch (i)
    {
      // constant
    case 0:
      return 1.;

      // linear
    case 1:
      return dx;

    case 2:
      return dy;

      // quadratics
    case 3:
      return dx*dx;

    case 4:
      return dx*dy;

    case 5:
      return dy*dy;

      // cubics
    case 6:
      return dx*dx*dx;

    case 7:
      return dx*dx*dy;

    case 8:
      return dx*dy*dy;

    case 9:
      return dy*dy*dy;

      // quartics
    case 10:
      return dx*dx*dx*dx;

    case 11:
      return dx*dx*dx*dy;

    case 12:
      return dx*dx*dy*dy;

    case 13:
      return dx*dy*dy*dy;

    case 14:
      return dy*dy*dy*dy;

    default:
      unsigned int o = 0;
      for (; i >= (o+1)*(o+2)/2; o++) { }
      unsigned int i2 = i - (o*(o+1)/2);
      Real val = 1.;
      for (unsigned int index=i2; index != o; index++)
        val *= dx;
      for (unsigned int index=0; index != i2; index++)
        val *= dy;
      return val;
    }

  libmesh_error_msg("We'll never get here!");
  return 0.;

#endif
}



template <>
Real FE<2,XYZ>::shape_deriv(const ElemType,
                            const Order,
                            const unsigned int,
                            const unsigned int,
                            const Point &)
{
  libmesh_error_msg("XYZ polynomials require the element \nbecause the centroid is needed.");
  return 0.;
}



template <>
Real FE<2,XYZ>::shape_deriv(const Elem * elem,
                            const Order libmesh_dbg_var(order),
                            const unsigned int i,
                            const unsigned int j,
                            const Point & point_in)
{
#if LIBMESH_DIM > 1


  libmesh_assert_less (j, 2);
  libmesh_assert(elem);

  Point centroid = elem->centroid();
  Point max_distance = Point(0.,0.,0.);
  for (unsigned int p = 0; p < elem->n_nodes(); p++)
    for (unsigned int d = 0; d < 2; d++)
      {
        const Real distance = std::abs(centroid(d) - elem->point(p)(d));
        max_distance(d) = std::max(distance, max_distance(d));
      }

  const Real x  = point_in(0);
  const Real y  = point_in(1);
  const Real xc = centroid(0);
  const Real yc = centroid(1);
  const Real distx = max_distance(0);
  const Real disty = max_distance(1);
  const Real dx = (x - xc)/distx;
  const Real dy = (y - yc)/disty;

#ifndef NDEBUG
  // totalorder is only used in the assertion below, so
  // we avoid declaring it when asserts are not active.
  const unsigned int totalorder = order + elem->p_level();
#endif
  libmesh_assert_less (i, (totalorder+1)*(totalorder+2)/2);

  // monomials. since they are hierarchic we only need one case block.

  switch (j)
    {
      // d()/dx
    case 0:
      {
        switch (i)
          {
            // constants
          case 0:
            return 0.;

            // linears
          case 1:
            return 1./distx;

          case 2:
            return 0.;

            // quadratics
          case 3:
            return 2.*dx/distx;

          case 4:
            return dy/distx;

          case 5:
            return 0.;

            // cubics
          case 6:
            return 3.*dx*dx/distx;

          case 7:
            return 2.*dx*dy/distx;

          case 8:
            return dy*dy/distx;

          case 9:
            return 0.;

            // quartics
          case 10:
            return 4.*dx*dx*dx/distx;

          case 11:
            return 3.*dx*dx*dy/distx;

          case 12:
            return 2.*dx*dy*dy/distx;

          case 13:
            return dy*dy*dy/distx;

          case 14:
            return 0.;

          default:
            unsigned int o = 0;
            for (; i >= (o+1)*(o+2)/2; o++) { }
            unsigned int i2 = i - (o*(o+1)/2);
            Real val = o - i2;
            for (unsigned int index=i2+1; index < o; index++)
              val *= dx;
            for (unsigned int index=0; index != i2; index++)
              val *= dy;
            return val/distx;
          }
      }


      // d()/dy
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
            return 1./disty;

            // quadratics
          case 3:
            return 0.;

          case 4:
            return dx/disty;

          case 5:
            return 2.*dy/disty;

            // cubics
          case 6:
            return 0.;

          case 7:
            return dx*dx/disty;

          case 8:
            return 2.*dx*dy/disty;

          case 9:
            return 3.*dy*dy/disty;

            // quartics
          case 10:
            return 0.;

          case 11:
            return dx*dx*dx/disty;

          case 12:
            return 2.*dx*dx*dy/disty;

          case 13:
            return 3.*dx*dy*dy/disty;

          case 14:
            return 4.*dy*dy*dy/disty;

          default:
            unsigned int o = 0;
            for (; i >= (o+1)*(o+2)/2; o++) { }
            unsigned int i2 = i - (o*(o+1)/2);
            Real val = i2;
            for (unsigned int index=i2; index != o; index++)
              val *= dx;
            for (unsigned int index=1; index <= i2; index++)
              val *= dy;
            return val/disty;
          }
      }


    default:
      libmesh_error_msg("Invalid j = " << j);
    }

  libmesh_error_msg("We'll never get here!");
  return 0.;

#endif
}



template <>
Real FE<2,XYZ>::shape_second_deriv(const ElemType,
                                   const Order,
                                   const unsigned int,
                                   const unsigned int,
                                   const Point &)
{
  libmesh_error_msg("XYZ polynomials require the element \nbecause the centroid is needed.");
  return 0.;
}



template <>
Real FE<2,XYZ>::shape_second_deriv(const Elem * elem,
                                   const Order libmesh_dbg_var(order),
                                   const unsigned int i,
                                   const unsigned int j,
                                   const Point & point_in)
{
#if LIBMESH_DIM > 1

  libmesh_assert_less_equal (j, 2);
  libmesh_assert(elem);

  Point centroid = elem->centroid();
  Point max_distance = Point(0.,0.,0.);
  for (unsigned int p = 0; p < elem->n_nodes(); p++)
    for (unsigned int d = 0; d < 2; d++)
      {
        const Real distance = std::abs(centroid(d) - elem->point(p)(d));
        max_distance(d) = std::max(distance, max_distance(d));
      }

  const Real x  = point_in(0);
  const Real y  = point_in(1);
  const Real xc = centroid(0);
  const Real yc = centroid(1);
  const Real distx = max_distance(0);
  const Real disty = max_distance(1);
  const Real dx = (x - xc)/distx;
  const Real dy = (y - yc)/disty;
  const Real dist2x = pow(distx,2.);
  const Real dist2y = pow(disty,2.);
  const Real distxy = distx * disty;

#ifndef NDEBUG
  // totalorder is only used in the assertion below, so
  // we avoid declaring it when asserts are not active.
  const unsigned int totalorder = order + elem->p_level();
#endif
  libmesh_assert_less (i, (totalorder+1)*(totalorder+2)/2);

  // monomials. since they are hierarchic we only need one case block.

  switch (j)
    {
      // d^2()/dx^2
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
            return 2./dist2x;

          case 4:
          case 5:
            return 0.;

            // cubics
          case 6:
            return 6.*dx/dist2x;

          case 7:
            return 2.*dy/dist2x;

          case 8:
          case 9:
            return 0.;

            // quartics
          case 10:
            return 12.*dx*dx/dist2x;

          case 11:
            return 6.*dx*dy/dist2x;

          case 12:
            return 2.*dy*dy/dist2x;

          case 13:
          case 14:
            return 0.;

          default:
            unsigned int o = 0;
            for (; i >= (o+1)*(o+2)/2; o++) { }
            unsigned int i2 = i - (o*(o+1)/2);
            Real val = (o - i2) * (o - i2 - 1);
            for (unsigned int index=i2+2; index < o; index++)
              val *= dx;
            for (unsigned int index=0; index != i2; index++)
              val *= dy;
            return val/dist2x;
          }
      }

      // d^2()/dxdy
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
            return 1./distxy;

          case 5:
            return 0.;

            // cubics
          case 6:
            return 0.;
          case 7:
            return 2.*dx/distxy;

          case 8:
            return 2.*dy/distxy;

          case 9:
            return 0.;

            // quartics
          case 10:
            return 0.;

          case 11:
            return 3.*dx*dx/distxy;

          case 12:
            return 4.*dx*dy/distxy;

          case 13:
            return 3.*dy*dy/distxy;

          case 14:
            return 0.;

          default:
            unsigned int o = 0;
            for (; i >= (o+1)*(o+2)/2; o++) { }
            unsigned int i2 = i - (o*(o+1)/2);
            Real val = (o - i2) * i2;
            for (unsigned int index=i2+1; index < o; index++)
              val *= dx;
            for (unsigned int index=1; index < i2; index++)
              val *= dy;
            return val/distxy;
          }
      }

      // d^2()/dy^2
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
            return 2./dist2y;

            // cubics
          case 6:
            return 0.;

          case 7:
            return 0.;

          case 8:
            return 2.*dx/dist2y;

          case 9:
            return 6.*dy/dist2y;

            // quartics
          case 10:
          case 11:
            return 0.;

          case 12:
            return 2.*dx*dx/dist2y;

          case 13:
            return 6.*dx*dy/dist2y;

          case 14:
            return 12.*dy*dy/dist2y;

          default:
            unsigned int o = 0;
            for (; i >= (o+1)*(o+2)/2; o++) { }
            unsigned int i2 = i - (o*(o+1)/2);
            Real val = i2 * (i2 - 1);
            for (unsigned int index=i2; index != o; index++)
              val *= dx;
            for (unsigned int index=2; index < i2; index++)
              val *= dy;
            return val/dist2y;
          }
      }

    default:
      libmesh_error_msg("Invalid shape function derivative j = " << j);
    }

  libmesh_error_msg("We'll never get here!");
  return 0.;

#endif
}

} // namespace libMesh
