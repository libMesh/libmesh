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
#include "libmesh/number_lookups.h"
#include "libmesh/utility.h"


namespace libMesh
{

template <>
Real FE<2,HIERARCHIC>::shape(const ElemType,
                             const Order,
                             const unsigned int,
                             const Point &)
{
  libmesh_error_msg("Hierarchic polynomials require the element type \nbecause edge orientation is needed.");
  return 0.;
}



template <>
Real FE<2,HIERARCHIC>::shape(const Elem * elem,
                             const Order order,
                             const unsigned int i,
                             const Point & p)
{
  libmesh_assert(elem);

  const Order totalorder = static_cast<Order>(order+elem->p_level());
  libmesh_assert_greater (totalorder, 0);

  switch (elem->type())
    {
    case TRI3:
    case TRISHELL3:
    case TRI6:
      {
        const Real zeta1 = p(0);
        const Real zeta2 = p(1);
        const Real zeta0 = 1. - zeta1 - zeta2;

        libmesh_assert_less (i, (totalorder+1u)*(totalorder+2u)/2);
        libmesh_assert (elem->type() == TRI6 || totalorder < 2);

        // Vertex DoFs
        if (i == 0)
          return zeta0;
        else if (i == 1)
          return zeta1;
        else if (i == 2)
          return zeta2;
        // Edge DoFs
        else if (i < totalorder + 2u)
          {
            // Avoid returning NaN on vertices!
            if (zeta0 + zeta1 == 0.)
              return 0.;

            const unsigned int basisorder = i - 1;
            // Get factors to account for edge-flipping
            Real f0 = 1;
            if (basisorder%2 && (elem->point(0) > elem->point(1)))
              f0 = -1.;

            Real edgeval = (zeta1 - zeta0) / (zeta1 + zeta0);
            Real crossfunc = zeta0 + zeta1;
            for (unsigned int n=1; n != basisorder; ++n)
              crossfunc *= (zeta0 + zeta1);

            return f0 * crossfunc *
              FE<1,HIERARCHIC>::shape(EDGE3, totalorder,
                                      basisorder, edgeval);
          }
        else if (i < 2u*totalorder + 1)
          {
            // Avoid returning NaN on vertices!
            if (zeta1 + zeta2 == 0.)
              return 0.;

            const unsigned int basisorder = i - totalorder;
            // Get factors to account for edge-flipping
            Real f1 = 1;
            if (basisorder%2 && (elem->point(1) > elem->point(2)))
              f1 = -1.;

            Real edgeval = (zeta2 - zeta1) / (zeta2 + zeta1);
            Real crossfunc = zeta2 + zeta1;
            for (unsigned int n=1; n != basisorder; ++n)
              crossfunc *= (zeta2 + zeta1);

            return f1 * crossfunc *
              FE<1,HIERARCHIC>::shape(EDGE3, totalorder,
                                      basisorder, edgeval);
          }
        else if (i < 3u*totalorder)
          {
            // Avoid returning NaN on vertices!
            if (zeta0 + zeta2 == 0.)
              return 0.;

            const unsigned int basisorder = i - (2u*totalorder) + 1;
            // Get factors to account for edge-flipping
            Real f2 = 1;
            if (basisorder%2 && (elem->point(2) > elem->point(0)))
              f2 = -1.;

            Real edgeval = (zeta0 - zeta2) / (zeta0 + zeta2);
            Real crossfunc = zeta0 + zeta2;
            for (unsigned int n=1; n != basisorder; ++n)
              crossfunc *= (zeta0 + zeta2);

            return f2 * crossfunc *
              FE<1,HIERARCHIC>::shape(EDGE3, totalorder,
                                      basisorder, edgeval);
          }
        // Interior DoFs
        else
          {
            const unsigned int basisnum = i - (3u*totalorder);
            unsigned int exp0 = triangular_number_column[basisnum] + 1;
            unsigned int exp1 = triangular_number_row[basisnum] + 1 -
              triangular_number_column[basisnum];

            Real returnval = 1;
            for (unsigned int n = 0; n != exp0; ++n)
              returnval *= zeta0;
            for (unsigned int n = 0; n != exp1; ++n)
              returnval *= zeta1;
            returnval *= zeta2;
            return returnval;
          }
      }

      // Hierarchic shape functions on the quadrilateral.
    case QUAD4:
    case QUADSHELL4:
      libmesh_assert_less (totalorder, 2);
      libmesh_fallthrough();
    case QUAD8:
    case QUADSHELL8:
    case QUAD9:
      {
        // Compute quad shape functions as a tensor-product
        const Real xi  = p(0);
        const Real eta = p(1);

        libmesh_assert_less (i, (totalorder+1u)*(totalorder+1u));

        // Example i, i0, i1 values for totalorder = 5:
        //                                    0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35
        //  static const unsigned int i0[] = {0, 1, 1, 0, 2, 3, 4, 5, 1, 1, 1, 1, 2, 3, 4, 5, 0, 0, 0, 0, 2, 3, 3, 2, 4, 4, 4, 3, 2, 5, 5, 5, 5, 4, 3, 2};
        //  static const unsigned int i1[] = {0, 0, 1, 1, 0, 0, 0, 0, 2, 3, 4, 5, 1, 1, 1, 1, 2, 3, 4, 5, 2, 2, 3, 3, 2, 3, 4, 4, 4, 2, 3, 4, 5, 5, 5, 5};

        unsigned int i0, i1;

        // Vertex DoFs
        if (i == 0)
          { i0 = 0; i1 = 0; }
        else if (i == 1)
          { i0 = 1; i1 = 0; }
        else if (i == 2)
          { i0 = 1; i1 = 1; }
        else if (i == 3)
          { i0 = 0; i1 = 1; }
        // Edge DoFs
        else if (i < totalorder + 3u)
          { i0 = i - 2; i1 = 0; }
        else if (i < 2u*totalorder + 2)
          { i0 = 1; i1 = i - totalorder - 1; }
        else if (i < 3u*totalorder + 1)
          { i0 = i - 2u*totalorder; i1 = 1; }
        else if (i < 4u*totalorder)
          { i0 = 0; i1 = i - 3u*totalorder + 1; }
        // Interior DoFs
        else
          {
            unsigned int basisnum = i - 4*totalorder;
            i0 = square_number_column[basisnum] + 2;
            i1 = square_number_row[basisnum] + 2;
          }

        // Flip odd degree of freedom values if necessary
        // to keep continuity on sides
        Real f = 1.;

        if ((i0%2) && (i0 > 2) && (i1 == 0))
          f = (elem->point(0) > elem->point(1))?-1.:1.;
        else if ((i0%2) && (i0>2) && (i1 == 1))
          f = (elem->point(3) > elem->point(2))?-1.:1.;
        else if ((i0 == 0) && (i1%2) && (i1>2))
          f = (elem->point(0) > elem->point(3))?-1.:1.;
        else if ((i0 == 1) && (i1%2) && (i1>2))
          f = (elem->point(1) > elem->point(2))?-1.:1.;

        return f*(FE<1,HIERARCHIC>::shape(EDGE3, totalorder, i0, xi)*
                  FE<1,HIERARCHIC>::shape(EDGE3, totalorder, i1, eta));
      }

    default:
      libmesh_error_msg("ERROR: Unsupported element type = " << elem->type());
    }

  return 0.;
}



template <>
Real FE<2,HIERARCHIC>::shape_deriv(const ElemType,
                                   const Order,
                                   const unsigned int,
                                   const unsigned int,
                                   const Point &)
{
  libmesh_error_msg("Hierarchic polynomials require the element type \nbecause edge orientation is needed.");
  return 0.;
}



template <>
Real FE<2,HIERARCHIC>::shape_deriv(const Elem * elem,
                                   const Order order,
                                   const unsigned int i,
                                   const unsigned int j,
                                   const Point & p)
{
  libmesh_assert(elem);

  const ElemType type = elem->type();

  const Order totalorder = static_cast<Order>(order+elem->p_level());

  libmesh_assert_greater (totalorder, 0);

  switch (type)
    {
      // 1st & 2nd-order Hierarchics.
    case TRI3:
    case TRISHELL3:
    case TRI6:
      {
        const Real eps = 1.e-6;

        libmesh_assert_less (j, 2);

        switch (j)
          {
            //  d()/dxi
          case 0:
            {
              const Point pp(p(0)+eps, p(1));
              const Point pm(p(0)-eps, p(1));

              return (FE<2,HIERARCHIC>::shape(elem, order, i, pp) -
                      FE<2,HIERARCHIC>::shape(elem, order, i, pm))/2./eps;
            }

            // d()/deta
          case 1:
            {
              const Point pp(p(0), p(1)+eps);
              const Point pm(p(0), p(1)-eps);

              return (FE<2,HIERARCHIC>::shape(elem, order, i, pp) -
                      FE<2,HIERARCHIC>::shape(elem, order, i, pm))/2./eps;
            }

          default:
            libmesh_error_msg("Invalid derivative index j = " << j);
          }
      }

    case QUAD4:
    case QUADSHELL4:
      libmesh_assert_less (totalorder, 2);
      libmesh_fallthrough();
    case QUAD8:
    case QUADSHELL8:
    case QUAD9:
      {
        // Compute quad shape functions as a tensor-product
        const Real xi  = p(0);
        const Real eta = p(1);

        libmesh_assert_less (i, (totalorder+1u)*(totalorder+1u));

        // Example i, i0, i1 values for totalorder = 5:
        //                                    0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35
        //  static const unsigned int i0[] = {0, 1, 1, 0, 2, 3, 4, 5, 1, 1, 1, 1, 2, 3, 4, 5, 0, 0, 0, 0, 2, 3, 4, 5, 2, 3, 4, 5, 2, 3, 4, 5, 2, 3, 4, 5};
        //  static const unsigned int i1[] = {0, 0, 1, 1, 0, 0, 0, 0, 2, 3, 4, 5, 1, 1, 1, 1, 2, 3, 4, 5, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5};

        unsigned int i0, i1;

        // Vertex DoFs
        if (i == 0)
          { i0 = 0; i1 = 0; }
        else if (i == 1)
          { i0 = 1; i1 = 0; }
        else if (i == 2)
          { i0 = 1; i1 = 1; }
        else if (i == 3)
          { i0 = 0; i1 = 1; }
        // Edge DoFs
        else if (i < totalorder + 3u)
          { i0 = i - 2; i1 = 0; }
        else if (i < 2u*totalorder + 2)
          { i0 = 1; i1 = i - totalorder - 1; }
        else if (i < 3u*totalorder + 1u)
          { i0 = i - 2u*totalorder; i1 = 1; }
        else if (i < 4u*totalorder)
          { i0 = 0; i1 = i - 3u*totalorder + 1; }
        // Interior DoFs
        else
          {
            unsigned int basisnum = i - 4*totalorder;
            i0 = square_number_column[basisnum] + 2;
            i1 = square_number_row[basisnum] + 2;
          }

        // Flip odd degree of freedom values if necessary
        // to keep continuity on sides
        Real f = 1.;

        if ((i0%2) && (i0 > 2) && (i1 == 0))
          f = (elem->point(0) > elem->point(1))?-1.:1.;
        else if ((i0%2) && (i0>2) && (i1 == 1))
          f = (elem->point(3) > elem->point(2))?-1.:1.;
        else if ((i0 == 0) && (i1%2) && (i1>2))
          f = (elem->point(0) > elem->point(3))?-1.:1.;
        else if ((i0 == 1) && (i1%2) && (i1>2))
          f = (elem->point(1) > elem->point(2))?-1.:1.;

        switch (j)
          {
            // d()/dxi
          case 0:
            return f*(FE<1,HIERARCHIC>::shape_deriv(EDGE3, totalorder, i0, 0, xi)*
                      FE<1,HIERARCHIC>::shape      (EDGE3, totalorder, i1,    eta));

            // d()/deta
          case 1:
            return f*(FE<1,HIERARCHIC>::shape      (EDGE3, totalorder, i0,    xi)*
                      FE<1,HIERARCHIC>::shape_deriv(EDGE3, totalorder, i1, 0, eta));

          default:
            libmesh_error_msg("Invalid derivative index j = " << j);
          }
      }

    default:
      libmesh_error_msg("ERROR: Unsupported element type = " << type);
    }

  return 0.;
}



template <>
Real FE<2,HIERARCHIC>::shape_second_deriv(const ElemType,
                                          const Order,
                                          const unsigned int,
                                          const unsigned int,
                                          const Point &)
{
  libmesh_error_msg("Hierarchic polynomials require the element type \nbecause edge orientation is needed.");
  return 0.;
}



template <>
Real FE<2,HIERARCHIC>::shape_second_deriv(const Elem * elem,
                                          const Order order,
                                          const unsigned int i,
                                          const unsigned int j,
                                          const Point & p)
{
  libmesh_assert(elem);

  // I have been lazy here and am using finite differences
  // to compute the derivatives!
  const Real eps = 1.e-6;
  Point pp, pm;
  unsigned int prevj = libMesh::invalid_uint;

  switch (j)
    {
      //  d^2()/dxi^2
    case 0:
      {
        pp = Point(p(0)+eps, p(1));
        pm = Point(p(0)-eps, p(1));
        prevj = 0;
        break;
      }

      // d^2()/dxideta
    case 1:
      {
        pp = Point(p(0), p(1)+eps);
        pm = Point(p(0), p(1)-eps);
        prevj = 0;
        break;
      }

      // d^2()/deta^2
    case 2:
      {
        pp = Point(p(0), p(1)+eps);
        pm = Point(p(0), p(1)-eps);
        prevj = 1;
        break;
      }
    default:
      libmesh_error_msg("Invalid derivative index j = " << j);
    }

  return (FE<2,HIERARCHIC>::shape_deriv(elem, order, i, prevj, pp) -
          FE<2,HIERARCHIC>::shape_deriv(elem, order, i, prevj, pm)
          )/2./eps;
}

} // namespace libMesh
