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

#ifndef LIBMESH_FE_BERNSTEIN_SHAPE_2D_IMPL_H
#define LIBMESH_FE_BERNSTEIN_SHAPE_2D_IMPL_H

// Local includes
#include "libmesh/libmesh_config.h"
#ifdef LIBMESH_ENABLE_HIGHER_ORDER_SHAPES

#include "libmesh/fe.h"
#include "libmesh/elem.h"
#include "libmesh/number_lookups.h"
#include "libmesh/utility.h"


namespace libMesh
{


template <typename RealType>
typename FEShim<2,BERNSTEIN,RealType>::OutputShape FEShim<2,BERNSTEIN,RealType>::shape(const ElemType,
                            const Order,
                            const unsigned int,
                            const Point &)
{
  libmesh_error_msg("Bernstein polynomials require the element type \nbecause edge orientation is needed.");
  return 0.;
}



template <typename RealType>
typename FEShim<2,BERNSTEIN,RealType>::OutputShape FEShim<2,BERNSTEIN,RealType>::shape(const Elem * elem,
                            const Order order,
                            const unsigned int i,
                            const Point & p,
                            const bool add_p_level)
{
  libmesh_assert(elem);

  const ElemType type = elem->type();

  const Order totalorder =
    static_cast<Order>(order + add_p_level * elem->p_level());

  // Declare that we are using our own special power function
  // from the Utility namespace.  This saves typing later.
  using Utility::pow;

  switch (type)
    {
      // Hierarchic shape functions on the quadrilateral.
    case QUAD4:
    case QUADSHELL4:
    case QUAD9:
      {
        // Compute quad shape functions as a tensor-product
        const auto xi  = p(0);
        const auto eta = p(1);

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
        // Interior DoFs. Use Roy's number look up
        else
          {
            unsigned int basisnum = i - 4*totalorder;
            i0 = square_number_column[basisnum] + 2;
            i1 = square_number_row[basisnum] + 2;
          }


        // Flip odd degree of freedom values if necessary
        // to keep continuity on sides.
        if     ((i>= 4                 && i<= 4+  totalorder-2u) && elem->point(0) > elem->point(1)) i0=totalorder+2-i0;//
        else if ((i>= 4+  totalorder-1u && i<= 4+2*totalorder-3u) && elem->point(1) > elem->point(2)) i1=totalorder+2-i1;
        else if ((i>= 4+2*totalorder-2u && i<= 4+3*totalorder-4u) && elem->point(3) > elem->point(2)) i0=totalorder+2-i0;
        else if ((i>= 4+3*totalorder-3u && i<= 4+4*totalorder-5u) && elem->point(0) > elem->point(3)) i1=totalorder+2-i1;


        return (FEShim<1,BERNSTEIN,RealType>::shape(EDGE3, totalorder, i0, xi)*
                FEShim<1,BERNSTEIN,RealType>::shape(EDGE3, totalorder, i1, eta));
      }
      // handle serendipity QUAD8 element separately
    case QUAD8:
    case QUADSHELL8:
      {
        libmesh_assert_less (totalorder, 3);

        const auto xi  = p(0);
        const auto eta = p(1);

        libmesh_assert_less (i, 8);

        //                                0  1  2  3  4  5  6  7  8
        static const unsigned int i0[] = {0, 1, 1, 0, 2, 1, 2, 0, 2};
        static const unsigned int i1[] = {0, 0, 1, 1, 0, 2, 1, 2, 2};
        static const Real scal[] = {-0.25, -0.25, -0.25, -0.25, 0.5, 0.5, 0.5, 0.5};

        //B_t,i0(i)|xi * B_s,i1(i)|eta + scal(i) * B_t,2|xi * B_t,2|eta
        return (FEShim<1,BERNSTEIN,RealType>::shape(EDGE3, totalorder, i0[i], xi)*
                FEShim<1,BERNSTEIN,RealType>::shape(EDGE3, totalorder, i1[i], eta)
                +scal[i]*
                FEShim<1,BERNSTEIN,RealType>::shape(EDGE3, totalorder, i0[8], xi)*
                FEShim<1,BERNSTEIN,RealType>::shape(EDGE3, totalorder, i1[8], eta));

      }

    case TRI3:
    case TRISHELL3:
      libmesh_assert_less (totalorder, 2);
      libmesh_fallthrough();
    case TRI6:
      switch (totalorder)
        {
        case FIRST:
          {
            const auto x=p(0);
            const auto y=p(1);
            const auto r=1.-x-y;

            libmesh_assert_less (i, 3);

            switch(i)
              {
              case 0: return r;  //f0,0,1
              case 1: return x;  //f0,1,1
              case 2: return y;  //f1,0,1

              default:
                libmesh_error_msg("Invalid shape function index i = " << i);
              }
          }
        case SECOND:
          {
            const auto x=p(0);
            const auto y=p(1);
            const auto r=1.-x-y;

            libmesh_assert_less (i, 6);

            switch(i)
              {
              case 0: return r*r;
              case 1: return x*x;
              case 2: return y*y;

              case 3: return 2.*x*r;
              case 4: return 2.*x*y;
              case 5: return 2.*r*y;

              default:
                libmesh_error_msg("Invalid shape function index i = " << i);
              }
          }
        case THIRD:
          {
            const auto x=p(0);
            const auto y=p(1);
            const auto r=1.-x-y;
            libmesh_assert_less (i, 10);

            unsigned int shape=i;


            if ((i==3||i==4) && elem->point(0) > elem->point(1)) shape=7-i;
            if ((i==5||i==6) && elem->point(1) > elem->point(2)) shape=11-i;
            if ((i==7||i==8) && elem->point(0) > elem->point(2)) shape=15-i;

            switch(shape)
              {
              case 0: return r*r*r;
              case 1: return x*x*x;
              case 2: return y*y*y;

              case 3: return 3.*x*r*r;
              case 4: return 3.*x*x*r;

              case 5: return 3.*y*x*x;
              case 6: return 3.*y*y*x;

              case 7: return 3.*y*r*r;
              case 8: return 3.*y*y*r;

              case 9: return 6.*x*y*r;

              default:
                libmesh_error_msg("Invalid shape function index shape = " << shape);
              }
          }
        case FOURTH:
          {
            const auto x=p(0);
            const auto y=p(1);
            const auto r=1-x-y;
            unsigned int shape=i;

            libmesh_assert_less (i, 15);

            if ((i==3||i== 5) && elem->point(0) > elem->point(1)) shape=8-i;
            if ((i==6||i== 8) && elem->point(1) > elem->point(2)) shape=14-i;
            if ((i==9||i==11) && elem->point(0) > elem->point(2)) shape=20-i;


            switch(shape)
              {
                // point functions
              case  0: return r*r*r*r;
              case  1: return x*x*x*x;
              case  2: return y*y*y*y;

                // edge functions
              case  3: return 4.*x*r*r*r;
              case  4: return 6.*x*x*r*r;
              case  5: return 4.*x*x*x*r;

              case  6: return 4.*y*x*x*x;
              case  7: return 6.*y*y*x*x;
              case  8: return 4.*y*y*y*x;

              case  9: return 4.*y*r*r*r;
              case 10: return 6.*y*y*r*r;
              case 11: return 4.*y*y*y*r;

                // inner functions
              case 12: return 12.*x*y*r*r;
              case 13: return 12.*x*x*y*r;
              case 14: return 12.*x*y*y*r;

              default:
                libmesh_error_msg("Invalid shape function index shape = " << shape);
              }
          }
        case FIFTH:
          {
            const auto x=p(0);
            const auto y=p(1);
            const auto r=1-x-y;
            unsigned int shape=i;

            libmesh_assert_less (i, 21);

            if ((i>= 3&&i<= 6) && elem->point(0) > elem->point(1)) shape=9-i;
            if ((i>= 7&&i<=10) && elem->point(1) > elem->point(2)) shape=17-i;
            if ((i>=11&&i<=14) && elem->point(0) > elem->point(2)) shape=25-i;

            switch(shape)
              {
                //point functions
              case  0: return pow<5>(r);
              case  1: return pow<5>(x);
              case  2: return pow<5>(y);

                //edge functions
              case  3: return  5.*x        *pow<4>(r);
              case  4: return 10.*pow<2>(x)*pow<3>(r);
              case  5: return 10.*pow<3>(x)*pow<2>(r);
              case  6: return  5.*pow<4>(x)*r;

              case  7: return  5.*y   *pow<4>(x);
              case  8: return 10.*pow<2>(y)*pow<3>(x);
              case  9: return 10.*pow<3>(y)*pow<2>(x);
              case 10: return  5.*pow<4>(y)*x;

              case 11: return  5.*y   *pow<4>(r);
              case 12: return 10.*pow<2>(y)*pow<3>(r);
              case 13: return 10.*pow<3>(y)*pow<2>(r);
              case 14: return  5.*pow<4>(y)*r;

                //inner functions
              case 15: return 20.*x*y*pow<3>(r);
              case 16: return 30.*x*pow<2>(y)*pow<2>(r);
              case 17: return 30.*pow<2>(x)*y*pow<2>(r);
              case 18: return 20.*x*pow<3>(y)*r;
              case 19: return 20.*pow<3>(x)*y*r;
              case 20: return 30.*pow<2>(x)*pow<2>(y)*r;

              default:
                libmesh_error_msg("Invalid shape function index shape = " << shape);
              }
          }
        case SIXTH:
          {
            const auto x=p(0);
            const auto y=p(1);
            const auto r=1-x-y;
            unsigned int shape=i;

            libmesh_assert_less (i, 28);

            if ((i>= 3&&i<= 7) && elem->point(0) > elem->point(1)) shape=10-i;
            if ((i>= 8&&i<=12) && elem->point(1) > elem->point(2)) shape=20-i;
            if ((i>=13&&i<=17) && elem->point(0) > elem->point(2)) shape=30-i;

            switch(shape)
              {
                //point functions
              case  0: return pow<6>(r);
              case  1: return pow<6>(x);
              case  2: return pow<6>(y);

                //edge functions
              case  3: return  6.*x        *pow<5>(r);
              case  4: return 15.*pow<2>(x)*pow<4>(r);
              case  5: return 20.*pow<3>(x)*pow<3>(r);
              case  6: return 15.*pow<4>(x)*pow<2>(r);
              case  7: return  6.*pow<5>(x)*r;

              case  8: return  6.*y        *pow<5>(x);
              case  9: return 15.*pow<2>(y)*pow<4>(x);
              case 10: return 20.*pow<3>(y)*pow<3>(x);
              case 11: return 15.*pow<4>(y)*pow<2>(x);
              case 12: return  6.*pow<5>(y)*x;

              case 13: return  6.*y        *pow<5>(r);
              case 14: return 15.*pow<2>(y)*pow<4>(r);
              case 15: return 20.*pow<3>(y)*pow<3>(r);
              case 16: return 15.*pow<4>(y)*pow<2>(r);
              case 17: return  6.*pow<5>(y)*r;

                //inner functions
              case 18: return 30.*x*y*pow<4>(r);
              case 19: return 60.*x*pow<2>(y)*pow<3>(r);
              case 20: return 60.*  pow<2>(x)*y*pow<3>(r);
              case 21: return 60.*x*pow<3>(y)*pow<2>(r);
              case 22: return 60.*pow<3>(x)*y*pow<2>(r);
              case 23: return 90.*pow<2>(x)*pow<2>(y)*pow<2>(r);
              case 24: return 30.*x*pow<4>(y)*r;
              case 25: return 60.*pow<2>(x)*pow<3>(y)*r;
              case 26: return 60.*pow<3>(x)*pow<2>(y)*r;
              case 27: return 30.*pow<4>(x)*y*r;

              default:
                libmesh_error_msg("Invalid shape function index shape = " << shape);
              } // switch shape
          } // case TRI6
        default:
          libmesh_error_msg("Invalid totalorder = " << totalorder);
        } // switch order

    default:
      libmesh_error_msg("ERROR: Unsupported element type = " << type);
    } // switch type
}



template <typename RealType>
typename FEShim<2,BERNSTEIN,RealType>::OutputShape FEShim<2,BERNSTEIN,RealType>::shape_deriv(const ElemType,
                                  const Order,
                                  const unsigned int,
                                  const unsigned int,
                                  const Point &)
{
  libmesh_error_msg("Bernstein polynomials require the element type \nbecause edge orientation is needed.");
  return 0.;
}



template <typename RealType>
typename FEShim<2,BERNSTEIN,RealType>::OutputShape FEShim<2,BERNSTEIN,RealType>::shape_deriv(const Elem * elem,
                                  const Order order,
                                  const unsigned int i,
                                  const unsigned int j,
                                  const Point & p,
                                  const bool add_p_level)
{
  libmesh_assert(elem);

  const ElemType type = elem->type();

  const Order totalorder =
    static_cast<Order>(order + add_p_level * elem->p_level());

  switch (type)
    {
      // Hierarchic shape functions on the quadrilateral.
    case QUAD4:
    case QUAD9:
      {
        // Compute quad shape functions as a tensor-product
        const auto xi  = p(0);
        const auto eta = p(1);

        libmesh_assert_less (i, (totalorder+1u)*(totalorder+1u));

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
        if      ((i>= 4                 && i<= 4+  totalorder-2u) && elem->point(0) > elem->point(1)) i0=totalorder+2-i0;
        else if ((i>= 4+  totalorder-1u && i<= 4+2*totalorder-3u) && elem->point(1) > elem->point(2)) i1=totalorder+2-i1;
        else if ((i>= 4+2*totalorder-2u && i<= 4+3*totalorder-4u) && elem->point(3) > elem->point(2)) i0=totalorder+2-i0;
        else if ((i>= 4+3*totalorder-3u && i<= 4+4*totalorder-5u) && elem->point(0) > elem->point(3)) i1=totalorder+2-i1;

        switch (j)
          {
            // d()/dxi
          case 0:
            return (FEShim<1,BERNSTEIN,RealType>::shape_deriv(EDGE3, totalorder, i0, 0, xi)*
                    FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i1,    eta));

            // d()/deta
          case 1:
            return (FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i0,    xi)*
                    FEShim<1,BERNSTEIN,RealType>::shape_deriv(EDGE3, totalorder, i1, 0, eta));

          default:
            libmesh_error_msg("Invalid shape function derivative j = " << j);
          }
      }

      // Bernstein shape functions on the 8-noded quadrilateral
      // is handled separately.
    case QUAD8:
    case QUADSHELL8:
      {
        libmesh_assert_less (totalorder, 3);

        const auto xi  = p(0);
        const auto eta = p(1);

        libmesh_assert_less (i, 8);

        //                                0  1  2  3  4  5  6  7  8
        static const unsigned int i0[] = {0, 1, 1, 0, 2, 1, 2, 0, 2};
        static const unsigned int i1[] = {0, 0, 1, 1, 0, 2, 1, 2, 2};
        static const Real scal[] = {-0.25, -0.25, -0.25, -0.25, 0.5, 0.5, 0.5, 0.5};
        switch (j)
          {
            // d()/dxi
          case 0:
            return (FEShim<1,BERNSTEIN,RealType>::shape_deriv(EDGE3, totalorder, i0[i], 0, xi)*
                    FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i1[i],    eta)
                    +scal[i]*
                    FEShim<1,BERNSTEIN,RealType>::shape_deriv(EDGE3, totalorder, i0[8], 0, xi)*
                    FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i1[8],    eta));

            // d()/deta
          case 1:
            return (FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i0[i],    xi)*
                    FEShim<1,BERNSTEIN,RealType>::shape_deriv(EDGE3, totalorder, i1[i], 0, eta)
                    +scal[i]*
                    FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i0[8],    xi)*
                    FEShim<1,BERNSTEIN,RealType>::shape_deriv(EDGE3, totalorder, i1[8], 0, eta));

          default:
            libmesh_error_msg("Invalid shape function derivative j = " << j);
          }
      }

    case TRI3:
    case TRISHELL3:
      libmesh_assert_less (totalorder, 2);
      libmesh_fallthrough();
    case TRI6:
      {
        // I have been lazy here and am using finite differences
        // to compute the derivatives!
        const Real eps = 1.e-4;

        switch (j)
          {
            //  d()/dxi
          case 0:
            {
              const Point pp(p(0)+eps, p(1));
              const Point pm(p(0)-eps, p(1));

              return (FEShim<2,BERNSTEIN,RealType>::shape(elem, totalorder, i, pp) -
                      FEShim<2,BERNSTEIN,RealType>::shape(elem, totalorder, i, pm))/2./eps;
            }

            // d()/deta
          case 1:
            {
              const Point pp(p(0), p(1)+eps);
              const Point pm(p(0), p(1)-eps);

              return (FEShim<2,BERNSTEIN,RealType>::shape(elem, totalorder, i, pp) -
                      FEShim<2,BERNSTEIN,RealType>::shape(elem, totalorder, i, pm))/2./eps;
            }


          default:
            libmesh_error_msg("Invalid shape function derivative j = " << j);
          }
      }

    default:
      libmesh_error_msg("ERROR: Unsupported element type = " << type);
    }
}



#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
template <typename RealType>
typename FEShim<2,BERNSTEIN,RealType>::OutputShape FEShim<2,BERNSTEIN,RealType>::shape_second_deriv(const ElemType,
                                         const Order,
                                         const unsigned int,
                                         const unsigned int,
                                         const Point &)
{
  libmesh_error_msg("Bernstein polynomials require the element type \nbecause edge orientation is needed.");
  return 0.;
}



template <typename RealType>
typename FEShim<2,BERNSTEIN,RealType>::OutputShape FEShim<2,BERNSTEIN,RealType>::shape_second_deriv(const Elem * elem,
                                         const Order order,
                                         const unsigned int i,
                                         const unsigned int j,
                                         const Point & p,
                                         const bool add_p_level)
{
  libmesh_assert(elem);

  const ElemType type = elem->type();

  const Order totalorder =
    static_cast<Order>(order + add_p_level * elem->p_level());

  switch (type)
    {
      // Hierarchic shape functions on the quadrilateral.
    case QUAD4:
    case QUAD9:
      {
        // Compute quad shape functions as a tensor-product
        const auto xi  = p(0);
        const auto eta = p(1);

        libmesh_assert_less (i, (totalorder+1u)*(totalorder+1u));

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
        if      ((i>= 4                 && i<= 4+  totalorder-2u) && elem->point(0) > elem->point(1)) i0=totalorder+2-i0;
        else if ((i>= 4+  totalorder-1u && i<= 4+2*totalorder-3u) && elem->point(1) > elem->point(2)) i1=totalorder+2-i1;
        else if ((i>= 4+2*totalorder-2u && i<= 4+3*totalorder-4u) && elem->point(3) > elem->point(2)) i0=totalorder+2-i0;
        else if ((i>= 4+3*totalorder-3u && i<= 4+4*totalorder-5u) && elem->point(0) > elem->point(3)) i1=totalorder+2-i1;

        switch (j)
          {
            // d^2() / dxi^2
          case 0:
            return (FEShim<1,BERNSTEIN,RealType>::shape_second_deriv(EDGE3, totalorder, i0, 0, xi)*
                    FEShim<1,BERNSTEIN,RealType>::shape             (EDGE3, totalorder, i1,    eta));

            // d^2() / dxi deta
          case 1:
            return (FEShim<1,BERNSTEIN,RealType>::shape_deriv(EDGE3, totalorder, i0, 0, xi)*
                    FEShim<1,BERNSTEIN,RealType>::shape_deriv(EDGE3, totalorder, i1, 0, eta));

            // d^2() / deta^2
          case 2:
            return (FEShim<1,BERNSTEIN,RealType>::shape             (EDGE3, totalorder, i0,    xi)*
                    FEShim<1,BERNSTEIN,RealType>::shape_second_deriv(EDGE3, totalorder, i1, 0, eta));

          default:
            libmesh_error_msg("Invalid shape function derivative j = " << j);
          }
      }

      // Going to be lazy again about the hard cases.
    case TRI3:
    case TRISHELL3:
      libmesh_assert_less (totalorder, 2);
      libmesh_fallthrough();
    case QUAD8:
    case QUADSHELL8:
    case TRI6:
      {
        // I have been lazy here and am using finite differences
        // to compute the derivatives!
        const Real eps = 1.e-4;

        switch (j)
          {
            //  d^2() / dxi^2
          case 0:
            {
              const Point pp(p(0)+eps, p(1));
              const Point pm(p(0)-eps, p(1));

              return (FEShim<2,BERNSTEIN,RealType>::shape_deriv(elem, totalorder, i, 0, pp) -
                      FEShim<2,BERNSTEIN,RealType>::shape_deriv(elem, totalorder, i, 0, pm))/2./eps;
            }

            // d^2() / dxi deta
          case 1:
            {
              const Point pp(p(0), p(1)+eps);
              const Point pm(p(0), p(1)-eps);

              return (FEShim<2,BERNSTEIN,RealType>::shape_deriv(elem, totalorder, i, 0, pp) -
                      FEShim<2,BERNSTEIN,RealType>::shape_deriv(elem, totalorder, i, 0, pm))/2./eps;
            }

            // d^2() / deta^2
          case 2:
            {
              const Point pp(p(0), p(1)+eps);
              const Point pm(p(0), p(1)-eps);

              return (FEShim<2,BERNSTEIN,RealType>::shape_deriv(elem, totalorder, i, 1, pp) -
                      FEShim<2,BERNSTEIN,RealType>::shape_deriv(elem, totalorder, i, 1, pm))/2./eps;
            }

          default:
            libmesh_error_msg("Invalid shape function derivative j = " << j);
          }
      }

    default:
      libmesh_error_msg("ERROR: Unsupported element type = " << type);
    }
}

#endif

} // namespace libMesh


#endif// LIBMESH_ENABLE_HIGHER_ORDER_SHAPES

#endif // LIBMESH_FE_BERNSTEIN_SHAPE_2D_IMPL_H
