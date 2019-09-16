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



// C++ includes
#include <cstdlib> // *must* precede <cmath> for proper std:abs() on PGI, Sun Studio CC
#include <cmath> // for std::sqrt


// Local includes
#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_ENABLE_HIGHER_ORDER_SHAPES

#include "libmesh/fe.h"
#include "libmesh/elem.h"
#include "libmesh/utility.h"


// Anonymous namespace to hold static std::sqrt values
namespace
{
using libMesh::Real;

static const Real sqrt2  = std::sqrt(2.);
static const Real sqrt6  = std::sqrt(6.);
static const Real sqrt10 = std::sqrt(10.);
static const Real sqrt14 = std::sqrt(14.);
static const Real sqrt22 = std::sqrt(22.);
static const Real sqrt26 = std::sqrt(26.);
}


namespace libMesh
{

template <>
Real FE<2,SZABAB>::shape(const ElemType,
                         const Order,
                         const unsigned int,
                         const Point &)
{
  libmesh_error_msg("Szabo-Babuska polynomials require the element type \nbecause edge orientation is needed.");
  return 0.;
}



template <>
Real FE<2,SZABAB>::shape(const Elem * elem,
                         const Order order,
                         const unsigned int i,
                         const Point & p,
                         const bool add_p_level)
{
  libmesh_assert(elem);

  const ElemType type = elem->type();

  const Order totalorder = static_cast<Order>(order + add_p_level * elem->p_level());

  // Declare that we are using our own special power function
  // from the Utility namespace.  This saves typing later.
  using Utility::pow;

  switch (totalorder)
    {
      // 1st & 2nd-order Szabo-Babuska.
    case FIRST:
    case SECOND:
      {
        switch (type)
          {

            // Szabo-Babuska shape functions on the triangle.
          case TRI3:
          case TRI6:
            {
              const Real l1 = 1-p(0)-p(1);
              const Real l2 = p(0);
              const Real l3 = p(1);

              libmesh_assert_less (i, 6);

              switch (i)
                {
                case 0: return l1;
                case 1: return l2;
                case 2: return l3;

                case 3: return l1*l2*(-4.*sqrt6);
                case 4: return l2*l3*(-4.*sqrt6);
                case 5: return l3*l1*(-4.*sqrt6);

                default:
                  libmesh_error_msg("Invalid i = " << i);
                }
            }


            // Szabo-Babuska shape functions on the quadrilateral.
          case QUAD4:
          case QUAD8:
          case QUAD9:
            {
              // Compute quad shape functions as a tensor-product
              const Real xi  = p(0);
              const Real eta = p(1);

              libmesh_assert_less (i, 9);

              //                                0  1  2  3  4  5  6  7  8
              static const unsigned int i0[] = {0, 1, 1, 0, 2, 1, 2, 0, 2};
              static const unsigned int i1[] = {0, 0, 1, 1, 0, 2, 1, 2, 2};

              return (FE<1,SZABAB>::shape(EDGE3, totalorder, i0[i], xi)*
                      FE<1,SZABAB>::shape(EDGE3, totalorder, i1[i], eta));

            }

          default:
            libmesh_error_msg("Invalid element type = " << type);
          }
      }


      // 3rd-order Szabo-Babuska.
    case THIRD:
      {
        switch (type)
          {

            // Szabo-Babuska shape functions on the triangle.
          case TRI6:
            {
              Real l1 = 1-p(0)-p(1);
              Real l2 = p(0);
              Real l3 = p(1);

              Real f=1;

              libmesh_assert_less (i, 10);


              if (i==4 && (elem->point(0) > elem->point(1)))f=-1;
              if (i==6 && (elem->point(1) > elem->point(2)))f=-1;
              if (i==8 && (elem->point(2) > elem->point(0)))f=-1;


              switch (i)
                {
                  //nodal modes
                case 0: return l1;
                case 1: return l2;
                case 2: return l3;

                  //side modes
                case 3: return   l1*l2*(-4.*sqrt6);
                case 4: return f*l1*l2*(-4.*sqrt10)*(l2-l1);

                case 5: return   l2*l3*(-4.*sqrt6);
                case 6: return f*l2*l3*(-4.*sqrt10)*(l3-l2);

                case 7: return   l3*l1*(-4.*sqrt6);
                case 8: return f*l3*l1*(-4.*sqrt10)*(l1-l3);

                  //internal modes
                case 9: return l1*l2*l3;

                default:
                  libmesh_error_msg("Invalid i = " << i);
                }
            }


            // Szabo-Babuska shape functions on the quadrilateral.
          case QUAD8:
          case QUAD9:
            {
              // Compute quad shape functions as a tensor-product
              const Real xi  = p(0);
              const Real eta = p(1);

              libmesh_assert_less (i, 16);

              //                                0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15
              static const unsigned int i0[] = {0,  1,  1,  0,  2,  3,  1,  1,  2,  3,  0,  0,  2,  3,  2,  3};
              static const unsigned int i1[] = {0,  0,  1,  1,  0,  0,  2,  3,  1,  1,  2,  3,  2,  2,  3,  3};

              Real f=1.;

              // take care of edge orientation, this is needed at
              // edge shapes with (y=0)-asymmetric 1D shapes, these have
              // one 1D shape index being 0 or 1, the other one being odd and >=3

              switch(i)
                {
                case  5: // edge 0 points
                  if (elem->point(0) > elem->point(1))f = -1.;
                  break;
                case  7: // edge 1 points
                  if (elem->point(1) > elem->point(2))f = -1.;
                  break;
                case  9: // edge 2 points
                  if (elem->point(3) > elem->point(2))f = -1.;
                  break;
                case 11: // edge 3 points
                  if (elem->point(0) > elem->point(3))f = -1.;
                  break;

                default:
                  // Everything else keeps f=1
                  break;
                }

              return f*(FE<1,SZABAB>::shape(EDGE3, totalorder, i0[i], xi)*
                        FE<1,SZABAB>::shape(EDGE3, totalorder, i1[i], eta));
            }

          default:
            libmesh_error_msg("Invalid element type = " << type);
          }
      }




      // 4th-order Szabo-Babuska.
    case FOURTH:
      {
        switch (type)
          {
            // Szabo-Babuska shape functions on the triangle.
          case TRI6:
            {
              Real l1 = 1-p(0)-p(1);
              Real l2 = p(0);
              Real l3 = p(1);

              Real f=1;

              libmesh_assert_less (i, 15);


              if (i== 4 && (elem->point(0) > elem->point(1)))f=-1;
              if (i== 7 && (elem->point(1) > elem->point(2)))f=-1;
              if (i==10 && (elem->point(2) > elem->point(0)))f=-1;


              switch (i)
                {
                  //nodal modes
                case  0: return l1;
                case  1: return l2;
                case  2: return l3;

                  //side modes
                case  3: return   l1*l2*(-4.*sqrt6);
                case  4: return f*l1*l2*(-4.*sqrt10)*(l2-l1);
                case  5: return   l1*l2*(-sqrt14)*(5.*pow<2>(l2-l1)-1);

                case  6: return   l2*l3*(-4.*sqrt6);
                case  7: return f*l2*l3*(-4.*sqrt10)*(l3-l2);
                case  8: return   l2*l3*(-sqrt14)*(5.*pow<2>(l3-l2)-1);

                case  9: return   l3*l1*(-4.*sqrt6);
                case 10: return f*l3*l1*(-4.*sqrt10)*(l1-l3);
                case 11: return   l3*l1*(-sqrt14)*(5.*pow<2>(l1-l3)-1);

                  //internal modes
                case 12: return l1*l2*l3;

                case 13: return l1*l2*l3*(l2-l1);
                case 14: return l1*l2*l3*(2*l3-1);

                default:
                  libmesh_error_msg("Invalid i = " << i);
                }
            }


            // Szabo-Babuska shape functions on the quadrilateral.
          case QUAD8:
          case QUAD9:
            {
              // Compute quad shape functions as a tensor-product
              const Real xi  = p(0);
              const Real eta = p(1);

              libmesh_assert_less (i, 25);

              //                                0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24
              static const unsigned int i0[] = {0, 1, 1, 0, 2, 3, 4, 1, 1, 1, 2, 3, 4, 0, 0, 0, 2, 3, 4, 2, 3, 4, 2, 3, 4};
              static const unsigned int i1[] = {0, 0, 1, 1, 0, 0, 0, 2, 3, 4, 1, 1, 1, 2, 3, 4, 2, 2, 2, 3, 3, 3, 4, 4, 4};

              Real f=1.;

              switch(i)
                {
                case  5: // edge 0 points
                  if (elem->point(0) > elem->point(1))f = -1.;
                  break;
                case  8: // edge 1 points
                  if (elem->point(1) > elem->point(2))f = -1.;
                  break;
                case 11: // edge 2 points
                  if (elem->point(3) > elem->point(2))f = -1.;
                  break;
                case 14: // edge 3 points
                  if (elem->point(0) > elem->point(3))f = -1.;
                  break;

                default:
                  // Everything else keeps f=1
                  break;
                }

              return f*(FE<1,SZABAB>::shape(EDGE3, totalorder, i0[i], xi)*
                        FE<1,SZABAB>::shape(EDGE3, totalorder, i1[i], eta));
            }

          default:
            libmesh_error_msg("Invalid element type = " << type);
          }
      }




      // 5th-order Szabo-Babuska.
    case FIFTH:
      {
        switch (type)
          {
            // Szabo-Babuska shape functions on the triangle.
          case TRI6:
            {
              Real l1 = 1-p(0)-p(1);
              Real l2 = p(0);
              Real l3 = p(1);

              const Real x=l2-l1;
              const Real y=2.*l3-1;

              Real f=1;

              libmesh_assert_less (i, 21);


              if ((i== 4||i== 6) && (elem->point(0) > elem->point(1)))f=-1;
              if ((i== 8||i==10) && (elem->point(1) > elem->point(2)))f=-1;
              if ((i==12||i==14) && (elem->point(2) > elem->point(0)))f=-1;


              switch (i)
                {
                  //nodal modes
                case  0: return l1;
                case  1: return l2;
                case  2: return l3;

                  //side modes
                case  3: return   l1*l2*(-4.*sqrt6);
                case  4: return f*l1*l2*(-4.*sqrt10)*(l2-l1);
                case  5: return   -sqrt14*l1*l2*(5.0*l1*l1-1.0+(-10.0*l1+5.0*l2)*l2);
                case  6: return f*(-sqrt2)*l1*l2*((9.-21.*l1*l1)*l1+(-9.+63.*l1*l1+(-63.*l1+21.*l2)*l2)*l2);

                case  7: return   l2*l3*(-4.*sqrt6);
                case  8: return f*l2*l3*(-4.*sqrt10)*(l3-l2);
                case  9: return   -sqrt14*l2*l3*(5.0*l3*l3-1.0+(-10.0*l3+5.0*l2)*l2);
                case 10: return -f*sqrt2*l2*l3*((-9.0+21.0*l3*l3)*l3+(-63.0*l3*l3+9.0+(63.0*l3-21.0*l2)*l2)*l2);

                case 11: return   l3*l1*(-4.*sqrt6);
                case 12: return f*l3*l1*(-4.*sqrt10)*(l1-l3);
                case 13: return -sqrt14*l3*l1*(5.0*l3*l3-1.0+(-10.0*l3+5.0*l1)*l1);
                case 14: return f*(-sqrt2)*l3*l1*((9.0-21.0*l3*l3)*l3+(-9.0+63.0*l3*l3+(-63.0*l3+21.0*l1)*l1)*l1);

                  //internal modes
                case 15: return l1*l2*l3;

                case 16: return l1*l2*l3*x;
                case 17: return l1*l2*l3*y;

                case 18: return l1*l2*l3*(1.5*l1*l1-.5+(-3.0*l1+1.5*l2)*l2);
                case 19: return l1*l2*l3*x*y;
                case 20: return l1*l2*l3*(1.0+(-6.0+6.0*l3)*l3);

                default:
                  libmesh_error_msg("Invalid i = " << i);
                }
            } // case TRI6

            // Szabo-Babuska shape functions on the quadrilateral.
          case QUAD8:
          case QUAD9:
            {
              // Compute quad shape functions as a tensor-product
              const Real xi  = p(0);
              const Real eta = p(1);

              libmesh_assert_less (i, 36);

              //                                0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35
              static const unsigned int i0[] = {0, 1, 1, 0, 2, 3, 4, 5, 1, 1, 1, 1, 2, 3, 4, 5, 0, 0, 0, 0, 2, 3, 4, 5, 2, 3, 4, 5, 2, 3, 4, 5, 2, 3, 4, 5};
              static const unsigned int i1[] = {0, 0, 1, 1, 0, 0, 0, 0, 2, 3, 4, 5, 1, 1, 1, 1, 2, 3, 4, 5, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5};

              Real f=1.;

              switch(i)
                {
                case  5: // edge 0 points
                case  7:
                  if (elem->point(0) > elem->point(1))f = -1.;
                  break;
                case  9: // edge 1 points
                case 11:
                  if (elem->point(1) > elem->point(2))f = -1.;
                  break;
                case 13: // edge 2 points
                case 15:
                  if (elem->point(3) > elem->point(2))f = -1.;
                  break;
                case 14: // edge 3 points
                case 19:
                  if (elem->point(0) > elem->point(3))f = -1.;
                  break;

                default:
                  // Everything else keeps f=1
                  break;
                }

              return f*(FE<1,SZABAB>::shape(EDGE3, totalorder, i0[i], xi)*
                        FE<1,SZABAB>::shape(EDGE3, totalorder, i1[i], eta));

            } // case QUAD8/QUAD9

          default:
            libmesh_error_msg("Invalid element type = " << type);

          } // switch type

      } // case FIFTH

      // 6th-order Szabo-Babuska.
    case SIXTH:
      {
        switch (type)
          {
            // Szabo-Babuska shape functions on the triangle.
          case TRI6:
            {
              Real l1 = 1-p(0)-p(1);
              Real l2 = p(0);
              Real l3 = p(1);

              const Real x=l2-l1;
              const Real y=2.*l3-1;

              Real f=1;

              libmesh_assert_less (i, 28);


              if ((i== 4||i== 6) && (elem->point(0) > elem->point(1)))f=-1;
              if ((i== 9||i==11) && (elem->point(1) > elem->point(2)))f=-1;
              if ((i==14||i==16) && (elem->point(2) > elem->point(0)))f=-1;


              switch (i)
                {
                  //nodal modes
                case  0: return l1;
                case  1: return l2;
                case  2: return l3;

                  //side modes
                case  3: return   l1*l2*(-4.*sqrt6);
                case  4: return f*l1*l2*(-4.*sqrt10)*(l2-l1);
                case  5: return   -sqrt14*l1*l2*(5.0*l1*l1-1.0+(-10.0*l1+5.0*l2)*l2);
                case  6: return f*(-sqrt2)*l1*l2*((9.0-21.0*l1*l1)*l1+(-9.0+63.0*l1*l1+(-63.0*l1+21.0*l2)*l2)*l2);
                case  7: return   -sqrt22*l1*l2*(0.5+(-7.0+0.105E2*l1*l1)*l1*l1+((14.0-0.42E2*l1*l1)*l1+(-7.0+0.63E2*l1*l1+(-0.42E2*l1+0.105E2*l2)*l2)*l2)*l2);

                case  8: return   l2*l3*(-4.*sqrt6);
                case  9: return f*l2*l3*(-4.*sqrt10)*(l3-l2);
                case 10: return   -sqrt14*l2*l3*(5.0*l3*l3-1.0+(-10.0*l3+5.0*l2)*l2);
                case 11: return f*(-sqrt2)*l2*l3*((-9.0+21.0*l3*l3)*l3+(-63.0*l3*l3+9.0+(63.0*l3-21.0*l2)*l2)*l2);
                case 12: return   -sqrt22*l2*l3*(0.5+(-7.0+0.105E2*l3*l3)*l3*l3+((14.0-0.42E2*l3*l3)*l3+(-7.0+0.63E2*l3*l3+(-0.42E2*l3+0.105E2*l2)*l2)*l2)*l2);

                case 13: return   l3*l1*(-4.*sqrt6);
                case 14: return f*l3*l1*(-4.*sqrt10)*(l1-l3);
                case 15: return   -sqrt14*l3*l1*(5.0*l3*l3-1.0+(-10.0*l3+5.0*l1)*l1);
                case 16: return f*(-sqrt2)*l3*l1*((9.0-21.0*l3*l3)*l3+(-9.0+63.0*l3*l3+(-63.0*l3+21.0*l1)*l1)*l1);
                case 17: return   -sqrt22*l3*l1*(0.5+(-7.0+0.105E2*l3*l3)*l3*l3+((14.0-0.42E2*l3*l3)*l3+(-7.0+0.63E2*l3*l3+(-0.42E2*l3+0.105E2*l1)*l1)*l1)*l1);



                  //internal modes
                case 18: return l1*l2*l3;

                case 19: return l1*l2*l3*x;
                case 20: return l1*l2*l3*y;

                case 21: return 0.5*l1*l2*l3*(3.0*l1*l1-1.0+(-6.0*l1+3.0*l2)*l2);
                case 22: return l1*l2*l3*(l2-l1)*(2.0*l3-1.0);
                case 23: return 0.5*l1*l2*l3*(2.0+(-12.0+12.0*l3)*l3);
                case 24: return 0.5*l1*l2*l3*((3.0-5.0*l1*l1)*l1+(-3.0+15.0*l1*l1+(-15.0*l1+5.0*l2)*l2)*l2);
                case 25: return 0.5*l1*l2*l3*(3.0*l1*l1-1.0+(-6.0*l1+3.0*l2)*l2)*(2.0*l3-1.0);
                case 26: return 0.5*l1*l2*l3*(2.0+(-12.0+12.0*l3)*l3)*(l2-l1);
                case 27: return 0.5*l1*l2*l3*(-2.0+(24.0+(-60.0+40.0*l3)*l3)*l3);


                default:
                  libmesh_error_msg("Invalid i = " << i);
                }
            } // case TRI6

            // Szabo-Babuska shape functions on the quadrilateral.
          case QUAD8:
          case QUAD9:
            {
              // Compute quad shape functions as a tensor-product
              const Real xi  = p(0);
              const Real eta = p(1);

              libmesh_assert_less (i, 49);

              //                                0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48
              static const unsigned int i0[] = {0, 1, 1, 0, 2, 3, 4, 5, 6, 1, 1, 1, 1, 1, 2, 3, 4, 5, 6, 0, 0, 0, 0, 0, 2, 3, 4, 5, 6, 2, 3, 4, 5, 6, 2, 3, 4, 5, 6, 2, 3, 4, 5, 6, 2, 3, 4, 5, 6};
              static const unsigned int i1[] = {0, 0, 1, 1, 0, 0, 0, 0, 0, 2, 3, 4, 5, 6, 1, 1, 1, 1, 1, 2, 3, 4, 5, 6, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6};

              Real f=1.;

              switch(i)
                {
                case  5: // edge 0 points
                case  7:
                  if (elem->point(0) > elem->point(1))f = -1.;
                  break;
                case 10: // edge 1 points
                case 12:
                  if (elem->point(1) > elem->point(2))f = -1.;
                  break;
                case 15: // edge 2 points
                case 17:
                  if (elem->point(3) > elem->point(2))f = -1.;
                  break;
                case 20: // edge 3 points
                case 22:
                  if (elem->point(0) > elem->point(3))f = -1.;
                  break;

                default:
                  // Everything else keeps f=1
                  break;
                }

              return f*(FE<1,SZABAB>::shape(EDGE3, totalorder, i0[i], xi)*
                        FE<1,SZABAB>::shape(EDGE3, totalorder, i1[i], eta));

            } // case QUAD8/QUAD9

          default:
            libmesh_error_msg("Invalid element type = " << type);

          } // switch type

      } // case SIXTH


      // 7th-order Szabo-Babuska.
    case SEVENTH:
      {
        switch (type)
          {
            // Szabo-Babuska shape functions on the triangle.
          case TRI6:
            {

              Real l1 = 1-p(0)-p(1);
              Real l2 = p(0);
              Real l3 = p(1);

              const Real x=l2-l1;
              const Real y=2.*l3-1.;

              Real f=1;

              libmesh_assert_less (i, 36);


              if ((i>= 4&&i<= 8) && (elem->point(0) > elem->point(1)))f=-1;
              if ((i>=10&&i<=14) && (elem->point(1) > elem->point(2)))f=-1;
              if ((i>=16&&i<=20) && (elem->point(2) > elem->point(0)))f=-1;


              switch (i)
                {
                  //nodal modes
                case  0: return l1;
                case  1: return l2;
                case  2: return l3;

                  //side modes
                case  3: return   l1*l2*(-4.*sqrt6);
                case  4: return f*l1*l2*(-4.*sqrt10)*(l2-l1);

                case  5: return   -sqrt14*l1*l2*(5.0*l1*l1-1.0+(-10.0*l1+5.0*l2)*l2);
                case  6: return f*-sqrt2*l1*l2*((9.0-21.0*l1*l1)*l1+(-9.0+63.0*l1*l1+(-63.0*l1+21.0*l2)*l2)*l2);
                case  7: return   -sqrt22*l1*l2*(0.5+(-7.0+0.105E2*l1*l1)*l1*l1+((14.0-0.42E2*l1*l1)*l1+(-7.0+0.63E2*l1*l1+(-0.42E2*l1+0.105E2*l2)*l2)*l2)*l2);
                case  8: return f*-sqrt26*l1*l2*((-0.25E1+(15.0-0.165E2*l1*l1)*l1*l1)*l1+(0.25E1+(-45.0+0.825E2*l1*l1)*l1*l1+((45.0-0.165E3*l1*l1)*l1+(-15.0+0.165E3*l1*l1+(-0.825E2*l1+0.165E2*l2)*l2)*l2)*l2)*l2);

                case  9: return   l2*l3*(-4.*sqrt6);
                case 10: return f*l2*l3*(-4.*sqrt10)*(l3-l2);

                case 11: return   -sqrt14*l2*l3*(5.0*l3*l3-1.0+(-10.0*l3+5.0*l2)*l2);
                case 12: return f*-sqrt2*l2*l3*((-9.0+21.0*l3*l3)*l3+(-63.0*l3*l3+9.0+(63.0*l3-21.0*l2)*l2)*l2);
                case 13: return   -sqrt22*l2*l3*(0.5+(-7.0+0.105E2*l3*l3)*l3*l3+((14.0-0.42E2*l3*l3)*l3+(-7.0+0.63E2*l3*l3+(-0.42E2*l3+0.105E2*l2)*l2)*l2)*l2);
                case 14: return f*-sqrt26*l2*l3*((0.25E1+(-15.0+0.165E2*l3*l3)*l3*l3)*l3+(-0.25E1+(45.0-0.825E2*l3*l3)*l3*l3+((-45.0+0.165E3*l3*l3)*l3+(15.0-0.165E3*l3*l3+(0.825E2*l3-0.165E2*l2)*l2)*l2)*l2)*l2);

                case 15: return   l3*l1*(-4.*sqrt6);
                case 16: return f*l3*l1*(-4.*sqrt10)*(l1-l3);

                case 17: return   -sqrt14*l3*l1*(5.0*l3*l3-1.0+(-10.0*l3+5.0*l1)*l1);
                case 18: return -f*sqrt2*l3*l1*((9.-21.*l3*l3)*l3+(-9.+63.*l3*l3+(-63.*l3+21.*l1)*l1)*l1);
                case 19: return   -sqrt22*l3*l1*(0.5+(-7.0+0.105E2*l3*l3)*l3*l3+((14.0-0.42E2*l3*l3)*l3+(-7.0+0.63E2*l3*l3+(-0.42E2*l3+0.105E2*l1)*l1)*l1)*l1);
                case 20: return f*-sqrt26*l3*l1*((-0.25E1+(15.0-0.165E2*l3*l3)*l3*l3)*l3+(0.25E1+(-45.0+0.825E2*l3*l3)*l3*l3+((45.0-0.165E3*l3*l3)*l3+(-15.0+0.165E3*l3*l3+(-0.825E2*l3+0.165E2*l1)*l1)*l1)*l1)*l1);


                  //internal modes
                case 21: return l1*l2*l3;

                case 22: return l1*l2*l3*x;
                case 23: return l1*l2*l3*y;

                case 24: return l1*l2*l3*0.5*(3.*pow<2>(x)-1.);
                case 25: return l1*l2*l3*x*y;
                case 26: return l1*l2*l3*0.5*(3.*pow<2>(y)-1.);

                case 27: return 0.5*l1*l2*l3*((3.0-5.0*l1*l1)*l1+(-3.0+15.0*l1*l1+(-15.0*l1+5.0*l2)*l2)*l2);
                case 28: return 0.5*l1*l2*l3*(3.0*l1*l1-1.0+(-6.0*l1+3.0*l2)*l2)*(2.0*l3-1.0);
                case 29: return 0.5*l1*l2*l3*(2.0+(-12.0+12.0*l3)*l3)*(l2-l1);
                case 30: return 0.5*l1*l2*l3*(-2.0+(24.0+(-60.0+40.0*l3)*l3)*l3);
                case 31: return 0.125*l1*l2*l3*((-15.0+(70.0-63.0*l1*l1)*l1*l1)*l1+(15.0+(-210.0+315.0*l1*l1)*l1*l1+((210.0-630.0*l1*l1)*l1+(-70.0+630.0*l1*l1+(-315.0*l1+63.0*l2)*l2)*l2)*l2)*l2);
                case 32: return 0.5*l1*l2*l3*((3.0-5.0*l1*l1)*l1+(-3.0+15.0*l1*l1+(-15.0*l1+5.0*l2)*l2)*l2)*(2.0*l3-1.0);
                case 33: return 0.25*l1*l2*l3*(3.0*l1*l1-1.0+(-6.0*l1+3.0*l2)*l2)*(2.0+(-12.0+12.0*l3)*l3);
                case 34: return 0.5*l1*l2*l3*(-2.0+(24.0+(-60.0+40.0*l3)*l3)*l3)*(l2-l1);
                case 35: return 0.125*l1*l2*l3*(-8.0+(240.0+(-1680.0+(4480.0+(-5040.0+2016.0*l3)*l3)*l3)*l3)*l3);

                default:
                  libmesh_error_msg("Invalid i = " << i);
                }
            } // case TRI6

            // Szabo-Babuska shape functions on the quadrilateral.
          case QUAD8:
          case QUAD9:
            {
              // Compute quad shape functions as a tensor-product
              const Real xi  = p(0);
              const Real eta = p(1);

              libmesh_assert_less (i, 64);

              //                                0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63
              static const unsigned int i0[] = {0, 1, 1, 0, 2, 3, 4, 5, 6, 7, 1, 1, 1, 1, 1, 1, 2, 3, 4, 5, 6, 7, 0, 0, 0, 0, 0, 0, 2, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7};
              static const unsigned int i1[] = {0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 2, 3, 4, 5, 6, 7, 1, 1, 1, 1, 1, 1, 2, 3, 4, 5, 6, 7, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7};

              Real f=1.;

              switch(i)
                {
                case  5: // edge 0 points
                case  7:
                case  9:
                  if (elem->point(0) > elem->point(1))f = -1.;
                  break;
                case 11: // edge 1 points
                case 13:
                case 15:
                  if (elem->point(1) > elem->point(2))f = -1.;
                  break;
                case 17: // edge 2 points
                case 19:
                case 21:
                  if (elem->point(3) > elem->point(2))f = -1.;
                  break;
                case 23: // edge 3 points
                case 25:
                case 27:
                  if (elem->point(0) > elem->point(3))f = -1.;
                  break;

                default:
                  // Everything else keeps f=1
                  break;
                }

              return f*(FE<1,SZABAB>::shape(EDGE3, totalorder, i0[i], xi)*
                        FE<1,SZABAB>::shape(EDGE3, totalorder, i1[i], eta));

            } // case QUAD8/QUAD9

          default:
            libmesh_error_msg("Invalid element type = " << type);

          } // switch type

      } // case SEVENTH


      // by default throw an error
    default:
      libmesh_error_msg("ERROR: Unsupported polynomial order!");
    } // switch order
}





template <>
Real FE<2,SZABAB>::shape_deriv(const ElemType,
                               const Order,
                               const unsigned int,
                               const unsigned int,
                               const Point &)
{
  libmesh_error_msg("Szabo-Babuska polynomials require the element type \nbecause edge orientation is needed.");
  return 0.;
}



template <>
Real FE<2,SZABAB>::shape_deriv(const Elem * elem,
                               const Order order,
                               const unsigned int i,
                               const unsigned int j,
                               const Point & p,
                               const bool add_p_level)
{
  libmesh_assert(elem);

  const ElemType type = elem->type();

  const Order totalorder = static_cast<Order>(order + add_p_level * elem->p_level());

  switch (totalorder)
    {

      // 1st & 2nd-order Szabo-Babuska.
    case FIRST:
    case SECOND:
      {
        switch (type)
          {

            // Szabo-Babuska shape functions on the triangle.
          case TRI3:
          case TRI6:
            {
              // Here we use finite differences to compute the derivatives!
              const Real eps = 1.e-6;

              libmesh_assert_less (i, 6);
              libmesh_assert_less (j, 2);

              switch (j)
                {
                  //  d()/dxi
                case 0:
                  {
                    const Point pp(p(0)+eps, p(1));
                    const Point pm(p(0)-eps, p(1));

                    return (FE<2,SZABAB>::shape(elem, order, i, pp) -
                            FE<2,SZABAB>::shape(elem, order, i, pm))/2./eps;
                  }

                  // d()/deta
                case 1:
                  {
                    const Point pp(p(0), p(1)+eps);
                    const Point pm(p(0), p(1)-eps);

                    return (FE<2,SZABAB>::shape(elem, order, i, pp) -
                            FE<2,SZABAB>::shape(elem, order, i, pm))/2./eps;
                  }

                default:
                  libmesh_error_msg("Invalid j = " << j);
                }
            }



            // Szabo-Babuska shape functions on the quadrilateral.
          case QUAD4:
          case QUAD8:
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
                  return (FE<1,SZABAB>::shape_deriv(EDGE3, totalorder, i0[i], 0, xi)*
                          FE<1,SZABAB>::shape      (EDGE3, totalorder, i1[i],    eta));

                  // d()/deta
                case 1:
                  return (FE<1,SZABAB>::shape      (EDGE3, totalorder, i0[i],    xi)*
                          FE<1,SZABAB>::shape_deriv(EDGE3, totalorder, i1[i], 0, eta));

                default:
                  libmesh_error_msg("Invalid j = " << j);
                }
            }

          default:
            libmesh_error_msg("Invalid element type = " << type);
          }
      }



      // 3rd-order Szabo-Babuska.
    case THIRD:
      {
        switch (type)
          {
            // Szabo-Babuska shape functions on the triangle.
          case TRI6:
            {
              // Here we use finite differences to compute the derivatives!
              const Real eps = 1.e-6;

              libmesh_assert_less (i, 10);
              libmesh_assert_less (j, 2);

              switch (j)
                {
                  //  d()/dxi
                case 0:
                  {
                    const Point pp(p(0)+eps, p(1));
                    const Point pm(p(0)-eps, p(1));

                    return (FE<2,SZABAB>::shape(elem, order, i, pp) -
                            FE<2,SZABAB>::shape(elem, order, i, pm))/2./eps;
                  }

                  // d()/deta
                case 1:
                  {
                    const Point pp(p(0), p(1)+eps);
                    const Point pm(p(0), p(1)-eps);

                    return (FE<2,SZABAB>::shape(elem, order, i, pp) -
                            FE<2,SZABAB>::shape(elem, order, i, pm))/2./eps;
                  }


                default:
                  libmesh_error_msg("Invalid j = " << j);
                }
            }


            // Szabo-Babuska shape functions on the quadrilateral.
          case QUAD8:
          case QUAD9:
            {
              // Compute quad shape functions as a tensor-product
              const Real xi  = p(0);
              const Real eta = p(1);

              libmesh_assert_less (i, 16);

              //                                0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15
              static const unsigned int i0[] = {0,  1,  1,  0,  2,  3,  1,  1,  2,  3,  0,  0,  2,  3,  2,  3};
              static const unsigned int i1[] = {0,  0,  1,  1,  0,  0,  2,  3,  1,  1,  2,  3,  2,  2,  3,  3};

              Real f=1.;

              switch(i)
                {
                case  5: // edge 0 points
                  if (elem->point(0) > elem->point(1))f = -1.;
                  break;
                case  7: // edge 1 points
                  if (elem->point(1) > elem->point(2))f = -1.;
                  break;
                case  9: // edge 2 points
                  if (elem->point(3) > elem->point(2))f = -1.;
                  break;
                case 11: // edge 3 points
                  if (elem->point(0) > elem->point(3))f = -1.;
                  break;

                default:
                  // Everything else keeps f=1
                  break;
                }


              switch (j)
                {
                  // d()/dxi
                case 0:
                  return f*(FE<1,SZABAB>::shape_deriv(EDGE3, totalorder, i0[i], 0, xi)*
                            FE<1,SZABAB>::shape      (EDGE3, totalorder, i1[i],    eta));

                  // d()/deta
                case 1:
                  return f*(FE<1,SZABAB>::shape      (EDGE3, totalorder, i0[i],    xi)*
                            FE<1,SZABAB>::shape_deriv(EDGE3, totalorder, i1[i], 0, eta));

                default:
                  libmesh_error_msg("Invalid j = " << j);
                }
            }

          default:
            libmesh_error_msg("Invalid element type = " << type);
          }
      }




      // 4th-order Szabo-Babuska.
    case FOURTH:
      {
        switch (type)
          {

            // Szabo-Babuska shape functions on the triangle.
          case TRI6:
            {
              // Here we use finite differences to compute the derivatives!
              const Real eps = 1.e-6;

              libmesh_assert_less (i, 15);
              libmesh_assert_less (j, 2);

              switch (j)
                {
                  //  d()/dxi
                case 0:
                  {
                    const Point pp(p(0)+eps, p(1));
                    const Point pm(p(0)-eps, p(1));

                    return (FE<2,SZABAB>::shape(elem, order, i, pp) -
                            FE<2,SZABAB>::shape(elem, order, i, pm))/2./eps;
                  }

                  // d()/deta
                case 1:
                  {
                    const Point pp(p(0), p(1)+eps);
                    const Point pm(p(0), p(1)-eps);

                    return (FE<2,SZABAB>::shape(elem, order, i, pp) -
                            FE<2,SZABAB>::shape(elem, order, i, pm))/2./eps;
                  }


                default:
                  libmesh_error_msg("Invalid j = " << j);
                }
            }



            // Szabo-Babuska shape functions on the quadrilateral.
          case QUAD8:
          case QUAD9:
            {
              // Compute quad shape functions as a tensor-product
              const Real xi  = p(0);
              const Real eta = p(1);

              libmesh_assert_less (i, 25);

              //                                0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24
              static const unsigned int i0[] = {0, 1, 1, 0, 2, 3, 4, 1, 1, 1, 2, 3, 4, 0, 0, 0, 2, 3, 4, 2, 3, 4, 2, 3, 4};
              static const unsigned int i1[] = {0, 0, 1, 1, 0, 0, 0, 2, 3, 4, 1, 1, 1, 2, 3, 4, 2, 2, 2, 3, 3, 3, 4, 4, 4};

              Real f=1.;

              switch(i)
                {
                case  5: // edge 0 points
                  if (elem->point(0) > elem->point(1))f = -1.;
                  break;
                case  8: // edge 1 points
                  if (elem->point(1) > elem->point(2))f = -1.;
                  break;
                case 11: // edge 2 points
                  if (elem->point(3) > elem->point(2))f = -1.;
                  break;
                case 14: // edge 3 points
                  if (elem->point(0) > elem->point(3))f = -1.;
                  break;

                default:
                  // Everything else keeps f=1
                  break;
                }


              switch (j)
                {
                  // d()/dxi
                case 0:
                  return f*(FE<1,SZABAB>::shape_deriv(EDGE3, totalorder, i0[i], 0, xi)*
                            FE<1,SZABAB>::shape      (EDGE3, totalorder, i1[i],    eta));

                  // d()/deta
                case 1:
                  return f*(FE<1,SZABAB>::shape      (EDGE3, totalorder, i0[i],    xi)*
                            FE<1,SZABAB>::shape_deriv(EDGE3, totalorder, i1[i], 0, eta));

                default:
                  libmesh_error_msg("Invalid j = " << j);
                }
            }

          default:
            libmesh_error_msg("Invalid element type = " << type);
          }
      }




      // 5th-order Szabo-Babuska.
    case FIFTH:
      {
        // Szabo-Babuska shape functions on the quadrilateral.
        switch (type)
          {

            // Szabo-Babuska shape functions on the triangle.
          case TRI6:
            {
              // Here we use finite differences to compute the derivatives!
              const Real eps = 1.e-6;

              libmesh_assert_less (i, 21);
              libmesh_assert_less (j, 2);

              switch (j)
                {
                  //  d()/dxi
                case 0:
                  {
                    const Point pp(p(0)+eps, p(1));
                    const Point pm(p(0)-eps, p(1));

                    return (FE<2,SZABAB>::shape(elem, order, i, pp) -
                            FE<2,SZABAB>::shape(elem, order, i, pm))/2./eps;
                  }

                  // d()/deta
                case 1:
                  {
                    const Point pp(p(0), p(1)+eps);
                    const Point pm(p(0), p(1)-eps);

                    return (FE<2,SZABAB>::shape(elem, order, i, pp) -
                            FE<2,SZABAB>::shape(elem, order, i, pm))/2./eps;
                  }

                default:
                  libmesh_error_msg("Invalid j = " << j);
                }
            }



          case QUAD8:
          case QUAD9:
            {
              // Compute quad shape functions as a tensor-product
              const Real xi  = p(0);
              const Real eta = p(1);

              libmesh_assert_less (i, 36);

              //                                0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35
              static const unsigned int i0[] = {0, 1, 1, 0, 2, 3, 4, 5, 1, 1, 1, 1, 2, 3, 4, 5, 0, 0, 0, 0, 2, 3, 4, 5, 2, 3, 4, 5, 2, 3, 4, 5, 2, 3, 4, 5};
              static const unsigned int i1[] = {0, 0, 1, 1, 0, 0, 0, 0, 2, 3, 4, 5, 1, 1, 1, 1, 2, 3, 4, 5, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5};

              Real f=1.;

              switch(i)
                {
                case  5: // edge 0 points
                case  7:
                  if (elem->point(0) > elem->point(1))f = -1.;
                  break;
                case  9: // edge 1 points
                case 11:
                  if (elem->point(1) > elem->point(2))f = -1.;
                  break;
                case 13: // edge 2 points
                case 15:
                  if (elem->point(3) > elem->point(2))f = -1.;
                  break;
                case 14: // edge 3 points
                case 19:
                  if (elem->point(0) > elem->point(3))f = -1.;
                  break;

                default:
                  // Everything else keeps f=1
                  break;
                }


              switch (j)
                {
                  // d()/dxi
                case 0:
                  return f*(FE<1,SZABAB>::shape_deriv(EDGE3, totalorder, i0[i], 0, xi)*
                            FE<1,SZABAB>::shape      (EDGE3, totalorder, i1[i],    eta));

                  // d()/deta
                case 1:
                  return f*(FE<1,SZABAB>::shape      (EDGE3, totalorder, i0[i],    xi)*
                            FE<1,SZABAB>::shape_deriv(EDGE3, totalorder, i1[i], 0, eta));

                default:
                  libmesh_error_msg("Invalid j = " << j);
                }
            }

          default:
            libmesh_error_msg("Invalid element type = " << type);
          }
      }


      // 6th-order Szabo-Babuska.
    case SIXTH:
      {
        // Szabo-Babuska shape functions on the quadrilateral.
        switch (type)
          {

            // Szabo-Babuska shape functions on the triangle.
          case TRI6:
            {
              // Here we use finite differences to compute the derivatives!
              const Real eps = 1.e-6;

              libmesh_assert_less (i, 28);
              libmesh_assert_less (j, 2);

              switch (j)
                {
                  //  d()/dxi
                case 0:
                  {
                    const Point pp(p(0)+eps, p(1));
                    const Point pm(p(0)-eps, p(1));

                    return (FE<2,SZABAB>::shape(elem, order, i, pp) -
                            FE<2,SZABAB>::shape(elem, order, i, pm))/2./eps;
                  }

                  // d()/deta
                case 1:
                  {
                    const Point pp(p(0), p(1)+eps);
                    const Point pm(p(0), p(1)-eps);

                    return (FE<2,SZABAB>::shape(elem, order, i, pp) -
                            FE<2,SZABAB>::shape(elem, order, i, pm))/2./eps;
                  }

                default:
                  libmesh_error_msg("Invalid j = " << j);
                }
            }



          case QUAD8:
          case QUAD9:
            {
              // Compute quad shape functions as a tensor-product
              const Real xi  = p(0);
              const Real eta = p(1);

              libmesh_assert_less (i, 49);

              //                                0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48
              static const unsigned int i0[] = {0, 1, 1, 0, 2, 3, 4, 5, 6, 1, 1, 1, 1, 1, 2, 3, 4, 5, 6, 0, 0, 0, 0, 0, 2, 3, 4, 5, 6, 2, 3, 4, 5, 6, 2, 3, 4, 5, 6, 2, 3, 4, 5, 6, 2, 3, 4, 5, 6};
              static const unsigned int i1[] = {0, 0, 1, 1, 0, 0, 0, 0, 0, 2, 3, 4, 5, 6, 1, 1, 1, 1, 1, 2, 3, 4, 5, 6, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6};

              Real f=1.;

              switch(i)
                {
                case  5: // edge 0 points
                case  7:
                  if (elem->point(0) > elem->point(1))f = -1.;
                  break;
                case 10: // edge 1 points
                case 12:
                  if (elem->point(1) > elem->point(2))f = -1.;
                  break;
                case 15: // edge 2 points
                case 17:
                  if (elem->point(3) > elem->point(2))f = -1.;
                  break;
                case 20: // edge 3 points
                case 22:
                  if (elem->point(0) > elem->point(3))f = -1.;
                  break;

                default:
                  // Everything else keeps f=1
                  break;
                }


              switch (j)
                {
                  // d()/dxi
                case 0:
                  return f*(FE<1,SZABAB>::shape_deriv(EDGE3, totalorder, i0[i], 0, xi)*
                            FE<1,SZABAB>::shape      (EDGE3, totalorder, i1[i],    eta));

                  // d()/deta
                case 1:
                  return f*(FE<1,SZABAB>::shape      (EDGE3, totalorder, i0[i],    xi)*
                            FE<1,SZABAB>::shape_deriv(EDGE3, totalorder, i1[i], 0, eta));

                default:
                  libmesh_error_msg("Invalid j = " << j);
                }
            }

          default:
            libmesh_error_msg("Invalid element type = " << type);
          }
      }


      // 7th-order Szabo-Babuska.
    case SEVENTH:
      {
        // Szabo-Babuska shape functions on the quadrilateral.
        switch (type)
          {

            // Szabo-Babuska shape functions on the triangle.
          case TRI6:
            {
              // Here we use finite differences to compute the derivatives!
              const Real eps = 1.e-6;

              libmesh_assert_less (i, 36);
              libmesh_assert_less (j, 2);

              switch (j)
                {
                  //  d()/dxi
                case 0:
                  {
                    const Point pp(p(0)+eps, p(1));
                    const Point pm(p(0)-eps, p(1));

                    return (FE<2,SZABAB>::shape(elem, order, i, pp) -
                            FE<2,SZABAB>::shape(elem, order, i, pm))/2./eps;
                  }

                  // d()/deta
                case 1:
                  {
                    const Point pp(p(0), p(1)+eps);
                    const Point pm(p(0), p(1)-eps);

                    return (FE<2,SZABAB>::shape(elem, order, i, pp) -
                            FE<2,SZABAB>::shape(elem, order, i, pm))/2./eps;
                  }

                default:
                  libmesh_error_msg("Invalid j = " << j);
                }
            }



          case QUAD8:
          case QUAD9:
            {
              // Compute quad shape functions as a tensor-product
              const Real xi  = p(0);
              const Real eta = p(1);

              libmesh_assert_less (i, 64);

              //                                0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63
              static const unsigned int i0[] = {0, 1, 1, 0, 2, 3, 4, 5, 6, 7, 1, 1, 1, 1, 1, 1, 2, 3, 4, 5, 6, 7, 0, 0, 0, 0, 0, 0, 2, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7};
              static const unsigned int i1[] = {0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 2, 3, 4, 5, 6, 7, 1, 1, 1, 1, 1, 1, 2, 3, 4, 5, 6, 7, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7};

              Real f=1.;

              switch(i)
                {
                case  5: // edge 0 points
                case  7:
                case  9:
                  if (elem->point(0) > elem->point(1))f = -1.;
                  break;
                case 11: // edge 1 points
                case 13:
                case 15:
                  if (elem->point(1) > elem->point(2))f = -1.;
                  break;
                case 17: // edge 2 points
                case 19:
                case 21:
                  if (elem->point(3) > elem->point(2))f = -1.;
                  break;
                case 23: // edge 3 points
                case 25:
                case 27:
                  if (elem->point(0) > elem->point(3))f = -1.;
                  break;

                default:
                  // Everything else keeps f=1
                  break;
                }


              switch (j)
                {
                  // d()/dxi
                case 0:
                  return f*(FE<1,SZABAB>::shape_deriv(EDGE3, totalorder, i0[i], 0, xi)*
                            FE<1,SZABAB>::shape      (EDGE3, totalorder, i1[i],    eta));

                  // d()/deta
                case 1:
                  return f*(FE<1,SZABAB>::shape      (EDGE3, totalorder, i0[i],    xi)*
                            FE<1,SZABAB>::shape_deriv(EDGE3, totalorder, i1[i], 0, eta));

                default:
                  libmesh_error_msg("Invalid j = " << j);
                }
            }

          default:
            libmesh_error_msg("Invalid element type = " << type);
          }
      }



      // by default throw an error;call the orientation-independent shape functions
    default:
      libmesh_error_msg("ERROR: Unsupported polynomial order!");
    }
}


#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

template <>
Real FE<2,SZABAB>::shape_second_deriv(const ElemType,
                                      const Order,
                                      const unsigned int,
                                      const unsigned int,
                                      const Point &)
{
  static bool warning_given = false;

  if (!warning_given)
    libMesh::err << "Second derivatives for Szabab elements "
                 << " are not yet implemented!"
                 << std::endl;

  warning_given = true;
  return 0.;
}



template <>
Real FE<2,SZABAB>::shape_second_deriv(const Elem *,
                                      const Order,
                                      const unsigned int,
                                      const unsigned int,
                                      const Point &,
                                      const bool)
{
  static bool warning_given = false;

  if (!warning_given)
    libMesh::err << "Second derivatives for Szabab elements "
                 << " are not yet implemented!"
                 << std::endl;

  warning_given = true;
  return 0.;
}

#endif

} // namespace libMesh

#endif// LIBMESH_ENABLE_HIGHER_ORDER_SHAPES
