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

#ifndef LIBMESH_FE_BERNSTEIN_SHAPE_3D_IMPL_H
#define LIBMESH_FE_BERNSTEIN_SHAPE_3D_IMPL_H

// Local includes
#include "libmesh/libmesh_config.h"
#ifdef LIBMESH_ENABLE_HIGHER_ORDER_SHAPES

#include "libmesh/fe.h"
#include "libmesh/elem.h"


namespace libMesh
{



template <typename RealType>
typename FEShim<3,BERNSTEIN,RealType>::OutputShape FEShim<3,BERNSTEIN,RealType>::shape(const ElemType,
                            const Order,
                            const unsigned int,
                            const Point &)
{
  libmesh_error_msg("Bernstein polynomials require the element type \nbecause edge and face orientation is needed.");
  return 0.;
}



template <typename RealType>
typename FEShim<3,BERNSTEIN,RealType>::OutputShape FEShim<3,BERNSTEIN,RealType>::shape(const Elem * elem,
                            const Order order,
                            const unsigned int i,
                            const Point & p,
                            const bool add_p_level)
{

#if LIBMESH_DIM == 3

  libmesh_assert(elem);
  const ElemType type = elem->type();

  const Order totalorder =
    static_cast<Order>(order + add_p_level * elem->p_level());

  switch (totalorder)
    {
      // 1st order Bernstein.
    case FIRST:
      {
        switch (type)
          {

            // Bernstein shape functions on the tetrahedron.
          case TET4:
          case TET10:
            {
              libmesh_assert_less (i, 4);

              // Area coordinates
              const auto zeta1 = p(0);
              const auto zeta2 = p(1);
              const auto zeta3 = p(2);
              const auto zeta0 = 1. - zeta1 - zeta2 - zeta3;

              switch(i)
                {
                case  0:  return zeta0;
                case  1:  return zeta1;
                case  2:  return zeta2;
                case  3:  return zeta3;

                default:
                  libmesh_error_msg("Invalid shape function index i = " << i);
                }
            }

            // Bernstein shape functions on the hexahedral.
          case HEX8:
          case HEX20:
          case HEX27:
            {
              libmesh_assert_less (i, 8);

              // Compute hex shape functions as a tensor-product
              const auto xi   = p(0);
              const auto eta  = p(1);
              const auto zeta = p(2);

              // The only way to make any sense of this
              // is to look at the mgflo/mg2/mgf documentation
              // and make the cut-out cube!
              //                                0  1  2  3  4  5  6  7
              static const unsigned int i0[] = {0, 1, 1, 0, 0, 1, 1, 0};
              static const unsigned int i1[] = {0, 0, 1, 1, 0, 0, 1, 1};
              static const unsigned int i2[] = {0, 0, 0, 0, 1, 1, 1, 1};

              return (FEShim<1,BERNSTEIN,RealType>::shape(EDGE3, totalorder, i0[i], xi)*
                      FEShim<1,BERNSTEIN,RealType>::shape(EDGE3, totalorder, i1[i], eta)*
                      FEShim<1,BERNSTEIN,RealType>::shape(EDGE3, totalorder, i2[i], zeta));
            }


          default:
            libmesh_error_msg("Invalid element type = " << type);
          }
      }




    case SECOND:
      {
        switch (type)
          {

            // Bernstein shape functions on the tetrahedron.
          case TET10:
            {
              libmesh_assert_less (i, 10);

              // Area coordinates
              const auto zeta1 = p(0);
              const auto zeta2 = p(1);
              const auto zeta3 = p(2);
              const auto zeta0 = 1. - zeta1 - zeta2 - zeta3;

              switch(i)
                {
                case  0:  return zeta0*zeta0;
                case  1:  return zeta1*zeta1;
                case  2:  return zeta2*zeta2;
                case  3:  return zeta3*zeta3;
                case  4:  return 2.*zeta0*zeta1;
                case  5:  return 2.*zeta1*zeta2;
                case  6:  return 2.*zeta0*zeta2;
                case  7:  return 2.*zeta3*zeta0;
                case  8:  return 2.*zeta1*zeta3;
                case  9:  return 2.*zeta2*zeta3;

                default:
                  libmesh_error_msg("Invalid shape function index i = " << i);
                }
            }

            // Bernstein shape functions on the 20-noded hexahedral.
          case HEX20:
            {
              libmesh_assert_less (i, 20);

              // Compute hex shape functions as a tensor-product
              const auto xi   = p(0);
              const auto eta  = p(1);
              const auto zeta = p(2);

              // The only way to make any sense of this
              // is to look at the mgflo/mg2/mgf documentation
              // and make the cut-out cube!
              //                                0      1      2      3      4      5      6      7      8      9      10     11     12     13     14     15     16     17     18     19  20 21 22 23 24 25 26
              static const unsigned int i0[] = {0,     1,     1,     0,     0,     1,     1,     0,     2,     1,     2,     0,     0,     1,     1,     0,     2,     1,     2,     0,  2, 2, 1, 2, 0, 2, 2};
              static const unsigned int i1[] = {0,     0,     1,     1,     0,     0,     1,     1,     0,     2,     1,     2,     0,     0,     1,     1,     0,     2,     1,     2,  2, 0, 2, 1, 2, 2, 2};
              static const unsigned int i2[] = {0,     0,     0,     0,     1,     1,     1,     1,     0,     0,     0,     0,     2,     2,     2,     2,     1,     1,     1,     1,  0, 2, 2, 2, 2, 1, 2};
              //To compute the hex20 shape functions the original shape functions for hex27 are used.
              //scalx[i] tells how often the original x-th shape function has to be added to the original i-th shape function
              //to compute the new i-th shape function for hex20
              //example: B_0^HEX20 = B_0^HEX27 - 0.25*B_20^HEX27 - 0.25*B_21^HEX27 + 0*B_22^HEX27 + 0*B_23^HEX27 - 0.25*B_24^HEX27 + 0*B_25^HEX27 - 0.25*B_26^HEX27
              //         B_0^HEX20 = B_0^HEX27 + scal20[0]*B_20^HEX27 + scal21[0]*B_21^HEX27 + ...
              static const Real scal20[] =     {-0.25, -0.25, -0.25, -0.25, 0,     0,     0,     0,     0.5,   0.5,   0.5,   0.5,   0,     0,     0,     0,     0,     0,     0,     0};
              static const Real scal21[] =     {-0.25, -0.25, 0,     0,     -0.25, -0.25, 0,     0,     0.5,   0,     0,     0,     0.5,   0.5,   0,     0,     0.5,   0,     0,     0};
              static const Real scal22[] =     {0,     -0.25, -0.25, 0,     0,     -0.25, -0.25, 0,     0,     0.5,   0,     0,     0,     0.5,   0.5,   0,     0,     0.5,   0,     0};
              static const Real scal23[] =     {0,     0,     -0.25, -0.25, 0,     0,     -0.25, -0.25, 0,     0,     0.5,   0,     0,     0,     0.5,   0.5,   0,     0,     0.5,   0};
              static const Real scal24[] =     {-0.25, 0,     0,     -0.25, -0.25, 0,     0,     -0.25, 0,     0,     0,     0.5,   0.5,   0,     0,     0.5,   0,     0,     0,     0.5};
              static const Real scal25[] =     {0,     0,     0,     0,     -0.25, -0.25, -0.25, -0.25, 0,     0,     0,     0,     0,     0,     0,     0,     0.5,   0.5,   0.5,   0.5};
              static const Real scal26[] =     {-0.25, -0.25, -0.25, -0.25, -0.25, -0.25, -0.25, -0.25, 0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25};

              return (FEShim<1,BERNSTEIN,RealType>::shape(EDGE3, totalorder, i0[i], xi)*
                      FEShim<1,BERNSTEIN,RealType>::shape(EDGE3, totalorder, i1[i], eta)*
                      FEShim<1,BERNSTEIN,RealType>::shape(EDGE3, totalorder, i2[i], zeta)
                      +scal20[i]*
                      FEShim<1,BERNSTEIN,RealType>::shape(EDGE3, totalorder, i0[20], xi)*
                      FEShim<1,BERNSTEIN,RealType>::shape(EDGE3, totalorder, i1[20], eta)*
                      FEShim<1,BERNSTEIN,RealType>::shape(EDGE3, totalorder, i2[20], zeta)
                      +scal21[i]*
                      FEShim<1,BERNSTEIN,RealType>::shape(EDGE3, totalorder, i0[21], xi)*
                      FEShim<1,BERNSTEIN,RealType>::shape(EDGE3, totalorder, i1[21], eta)*
                      FEShim<1,BERNSTEIN,RealType>::shape(EDGE3, totalorder, i2[21], zeta)
                      +scal22[i]*
                      FEShim<1,BERNSTEIN,RealType>::shape(EDGE3, totalorder, i0[22], xi)*
                      FEShim<1,BERNSTEIN,RealType>::shape(EDGE3, totalorder, i1[22], eta)*
                      FEShim<1,BERNSTEIN,RealType>::shape(EDGE3, totalorder, i2[22], zeta)
                      +scal23[i]*
                      FEShim<1,BERNSTEIN,RealType>::shape(EDGE3, totalorder, i0[23], xi)*
                      FEShim<1,BERNSTEIN,RealType>::shape(EDGE3, totalorder, i1[23], eta)*
                      FEShim<1,BERNSTEIN,RealType>::shape(EDGE3, totalorder, i2[23], zeta)
                      +scal24[i]*
                      FEShim<1,BERNSTEIN,RealType>::shape(EDGE3, totalorder, i0[24], xi)*
                      FEShim<1,BERNSTEIN,RealType>::shape(EDGE3, totalorder, i1[24], eta)*
                      FEShim<1,BERNSTEIN,RealType>::shape(EDGE3, totalorder, i2[24], zeta)
                      +scal25[i]*
                      FEShim<1,BERNSTEIN,RealType>::shape(EDGE3, totalorder, i0[25], xi)*
                      FEShim<1,BERNSTEIN,RealType>::shape(EDGE3, totalorder, i1[25], eta)*
                      FEShim<1,BERNSTEIN,RealType>::shape(EDGE3, totalorder, i2[25], zeta)
                      +scal26[i]*
                      FEShim<1,BERNSTEIN,RealType>::shape(EDGE3, totalorder, i0[26], xi)*
                      FEShim<1,BERNSTEIN,RealType>::shape(EDGE3, totalorder, i1[26], eta)*
                      FEShim<1,BERNSTEIN,RealType>::shape(EDGE3, totalorder, i2[26], zeta));
            }

            // Bernstein shape functions on the hexahedral.
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

              return (FEShim<1,BERNSTEIN,RealType>::shape(EDGE3, totalorder, i0[i], xi)*
                      FEShim<1,BERNSTEIN,RealType>::shape(EDGE3, totalorder, i1[i], eta)*
                      FEShim<1,BERNSTEIN,RealType>::shape(EDGE3, totalorder, i2[i], zeta));
            }


          default:
            libmesh_error_msg("Invalid element type = " << type);
          }

      }



      // 3rd-order Bernstein.
    case THIRD:
      {
        switch (type)
          {

            //     // Bernstein shape functions on the tetrahedron.
            //   case TET10:
            //     {
            //       libmesh_assert_less (i, 20);

            //       // Area coordinates
            //       const auto zeta1 = p(0);
            //       const auto zeta2 = p(1);
            //       const auto zeta3 = p(2);
            //       const auto zeta0 = 1. - zeta1 - zeta2 - zeta3;


            //       unsigned int shape=i;

            //       // handle the edge orientation

            //       if ((i== 4||i== 5) && elem->node_id(0) > elem->node_id(1))shape= 9-i;   //Edge 0
            //       if ((i== 6||i== 7) && elem->node_id(1) > elem->node_id(2))shape=13-i;   //Edge 1
            //       if ((i== 8||i== 9) && elem->node_id(0) > elem->node_id(2))shape=17-i;   //Edge 2
            //       if ((i==10||i==11) && elem->node_id(0) > elem->node_id(3))shape=21-i;   //Edge 3
            //       if ((i==12||i==13) && elem->node_id(1) > elem->node_id(3))shape=25-i;   //Edge 4
            //       if ((i==14||i==15) && elem->node_id(2) > elem->node_id(3))shape=29-i;   //Edge 5

            //       // No need to handle face orientation in 3rd order.


            //       switch(shape)
            // {
            //   //point function
            // case  0:  return zeta0*zeta0*zeta0;
            // case  1:  return zeta1*zeta1*zeta1;
            // case  2:  return zeta2*zeta2*zeta2;
            // case  3:  return zeta3*zeta3*zeta3;

            //   //edge functions
            // case  4:  return 3.*zeta0*zeta0*zeta1;
            // case  5:  return 3.*zeta1*zeta1*zeta0;

            // case  6:  return 3.*zeta1*zeta1*zeta2;
            // case  7:  return 3.*zeta2*zeta2*zeta1;

            // case  8:  return 3.*zeta0*zeta0*zeta2;
            // case  9:  return 3.*zeta2*zeta2*zeta0;

            // case 10:  return 3.*zeta0*zeta0*zeta3;
            // case 11:  return 3.*zeta3*zeta3*zeta0;

            // case 12:  return 3.*zeta1*zeta1*zeta3;
            // case 13:  return 3.*zeta3*zeta3*zeta1;

            // case 14:  return 3.*zeta2*zeta2*zeta3;
            // case 15:  return 3.*zeta3*zeta3*zeta2;

            //   //face functions
            // case 16:  return 6.*zeta0*zeta1*zeta2;
            // case 17:  return 6.*zeta0*zeta1*zeta3;
            // case 18:  return 6.*zeta1*zeta2*zeta3;
            // case 19:  return 6.*zeta2*zeta0*zeta3;

            // default:
            // libmesh_error_msg("Invalid shape function index i = " << i);
            // }
            //     }


            // Bernstein shape functions on the hexahedral.
          case HEX27:
            {
              libmesh_assert_less (i, 64);

              // Compute hex shape functions as a tensor-product
              const auto xi    = p(0);
              const auto eta   = p(1);
              const auto zeta  = p(2);
              auto xi_mapped   = p(0);
              auto eta_mapped  = p(1);
              auto zeta_mapped = p(2);

              // The only way to make any sense of this
              // is to look at the mgflo/mg2/mgf documentation
              // and make the cut-out cube!
              //  Nodes                         0  1  2  3  4  5  6  7  8  8  9  9 10 10 11 11 12 12 13 13 14 14 15 15 16 16 17 17 18 18 19 19 20 20 20 20 21 21 21 21 22 22 22 22 23 23 23 23 24 24 24 24 25 25 25 25 26 26 26 26 26 26 26 26
              //  DOFS                          0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 18 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 60 62 63
              static const unsigned int i0[] = {0, 1, 1, 0, 0, 1, 1, 0, 2, 3, 1, 1, 2, 3, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 2, 3, 1, 1, 2, 3, 0, 0, 2, 3, 2, 3, 2, 3, 2, 3, 1, 1, 1, 1, 2, 3, 2, 3, 0, 0, 0, 0, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3};
              static const unsigned int i1[] = {0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 2, 3, 1, 1, 2, 3, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 2, 3, 1, 1, 2, 3, 2, 2, 3, 3, 0, 0, 0, 0, 2, 3, 2, 3, 1, 1, 1, 1, 2, 3, 2, 3, 2, 2, 3, 3, 2, 2, 3, 3, 2, 2, 3, 3};
              static const unsigned int i2[] = {0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 2, 3, 2, 3, 2, 3, 2, 3, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 2, 2, 3, 3, 2, 2, 3, 3, 2, 2, 3, 3, 2, 2, 3, 3, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3};



              // handle the edge orientation
              {
                // Edge 0
                if ((i0[i] >= 2) && (i1[i] == 0) && (i2[i] == 0))
                  {
                    if (elem->point(0) != std::min(elem->point(0), elem->point(1)))
                      xi_mapped = -xi;
                  }
                // Edge 1
                else if ((i0[i] == 1) && (i1[i] >= 2) && (i2[i] == 0))
                  {
                    if (elem->point(1) != std::min(elem->point(1), elem->point(2)))
                      eta_mapped = -eta;
                  }
                // Edge 2
                else if ((i0[i] >= 2) && (i1[i] == 1) && (i2[i] == 0))
                  {
                    if (elem->point(3) != std::min(elem->point(3), elem->point(2)))
                      xi_mapped = -xi;
                  }
                // Edge 3
                else if ((i0[i] == 0) && (i1[i] >= 2) && (i2[i] == 0))
                  {
                    if (elem->point(0) != std::min(elem->point(0), elem->point(3)))
                      eta_mapped = -eta;
                  }
                // Edge 4
                else if ((i0[i] == 0) && (i1[i] == 0) && (i2[i] >=2 ))
                  {
                    if (elem->point(0) != std::min(elem->point(0), elem->point(4)))
                      zeta_mapped = -zeta;
                  }
                // Edge 5
                else if ((i0[i] == 1) && (i1[i] == 0) && (i2[i] >=2 ))
                  {
                    if (elem->point(1) != std::min(elem->point(1), elem->point(5)))
                      zeta_mapped = -zeta;
                  }
                // Edge 6
                else if ((i0[i] == 1) && (i1[i] == 1) && (i2[i] >=2 ))
                  {
                    if (elem->point(2) != std::min(elem->point(2), elem->point(6)))
                      zeta_mapped = -zeta;
                  }
                // Edge 7
                else if ((i0[i] == 0) && (i1[i] == 1) && (i2[i] >=2 ))
                  {
                    if (elem->point(3) != std::min(elem->point(3), elem->point(7)))
                      zeta_mapped = -zeta;
                  }
                // Edge 8
                else if ((i0[i] >=2 ) && (i1[i] == 0) && (i2[i] == 1))
                  {
                    if (elem->point(4) != std::min(elem->point(4), elem->point(5)))
                      xi_mapped = -xi;
                  }
                // Edge 9
                else if ((i0[i] == 1) && (i1[i] >=2 ) && (i2[i] == 1))
                  {
                    if (elem->point(5) != std::min(elem->point(5), elem->point(6)))
                      eta_mapped = -eta;
                  }
                // Edge 10
                else if ((i0[i] >=2 ) && (i1[i] == 1) && (i2[i] == 1))
                  {
                    if (elem->point(7) != std::min(elem->point(7), elem->point(6)))
                      xi_mapped = -xi;
                  }
                // Edge 11
                else if ((i0[i] == 0) && (i1[i] >=2 ) && (i2[i] == 1))
                  {
                    if (elem->point(4) != std::min(elem->point(4), elem->point(7)))
                      eta_mapped = -eta;
                  }
              }


              // handle the face orientation
              {
                // Face 0
                if ((i2[i] == 0) && (i0[i] >= 2) && (i1[i] >= 2))
                  {
                    const Point min_point = std::min(elem->point(1),
                                                     std::min(elem->point(2),
                                                              std::min(elem->point(0),
                                                                       elem->point(3))));
                    if (elem->point(0) == min_point)
                      if (elem->point(1) == std::min(elem->point(1), elem->point(3)))
                        {
                          // Case 1
                          xi_mapped  = xi;
                          eta_mapped = eta;
                        }
                      else
                        {
                          // Case 2
                          xi_mapped  = eta;
                          eta_mapped = xi;
                        }

                    else if (elem->point(3) == min_point)
                      if (elem->point(0) == std::min(elem->point(0), elem->point(2)))
                        {
                          // Case 3
                          xi_mapped  = -eta;
                          eta_mapped = xi;
                        }
                      else
                        {
                          // Case 4
                          xi_mapped  = xi;
                          eta_mapped = -eta;
                        }

                    else if (elem->point(2) == min_point)
                      if (elem->point(3) == std::min(elem->point(3), elem->point(1)))
                        {
                          // Case 5
                          xi_mapped  = -xi;
                          eta_mapped = -eta;
                        }
                      else
                        {
                          // Case 6
                          xi_mapped  = -eta;
                          eta_mapped = -xi;
                        }

                    else if (elem->point(1) == min_point)
                      {
                        if (elem->point(2) == std::min(elem->point(2), elem->point(0)))
                          {
                            // Case 7
                            xi_mapped  = eta;
                            eta_mapped = -xi;
                          }
                        else
                          {
                            // Case 8
                            xi_mapped  = -xi;
                            eta_mapped = eta;
                          }
                      }
                  }


                // Face 1
                else if ((i1[i] == 0) && (i0[i] >= 2) && (i2[i] >= 2))
                  {
                    const Point min_point = std::min(elem->point(0),
                                                     std::min(elem->point(1),
                                                              std::min(elem->point(5),
                                                                       elem->point(4))));
                    if (elem->point(0) == min_point)
                      if (elem->point(1) == std::min(elem->point(1), elem->point(4)))
                        {
                          // Case 1
                          xi_mapped   = xi;
                          zeta_mapped = zeta;
                        }
                      else
                        {
                          // Case 2
                          xi_mapped   = zeta;
                          zeta_mapped = xi;
                        }

                    else if (elem->point(1) == min_point)
                      if (elem->point(5) == std::min(elem->point(5), elem->point(0)))
                        {
                          // Case 3
                          xi_mapped   = zeta;
                          zeta_mapped = -xi;
                        }
                      else
                        {
                          // Case 4
                          xi_mapped   = -xi;
                          zeta_mapped = zeta;
                        }

                    else if (elem->point(5) == min_point)
                      if (elem->point(4) == std::min(elem->point(4), elem->point(1)))
                        {
                          // Case 5
                          xi_mapped   = -xi;
                          zeta_mapped = -zeta;
                        }
                      else
                        {
                          // Case 6
                          xi_mapped   = -zeta;
                          zeta_mapped = -xi;
                        }

                    else if (elem->point(4) == min_point)
                      {
                        if (elem->point(0) == std::min(elem->point(0), elem->point(5)))
                          {
                            // Case 7
                            xi_mapped   = -xi;
                            zeta_mapped = zeta;
                          }
                        else
                          {
                            // Case 8
                            xi_mapped   = xi;
                            zeta_mapped = -zeta;
                          }
                      }
                  }


                // Face 2
                else if ((i0[i] == 1) && (i1[i] >= 2) && (i2[i] >= 2))
                  {
                    const Point min_point = std::min(elem->point(1),
                                                     std::min(elem->point(2),
                                                              std::min(elem->point(6),
                                                                       elem->point(5))));
                    if (elem->point(1) == min_point)
                      if (elem->point(2) == std::min(elem->point(2), elem->point(5)))
                        {
                          // Case 1
                          eta_mapped  = eta;
                          zeta_mapped = zeta;
                        }
                      else
                        {
                          // Case 2
                          eta_mapped  = zeta;
                          zeta_mapped = eta;
                        }

                    else if (elem->point(2) == min_point)
                      if (elem->point(6) == std::min(elem->point(6), elem->point(1)))
                        {
                          // Case 3
                          eta_mapped  = zeta;
                          zeta_mapped = -eta;
                        }
                      else
                        {
                          // Case 4
                          eta_mapped  = -eta;
                          zeta_mapped = zeta;
                        }

                    else if (elem->point(6) == min_point)
                      if (elem->point(5) == std::min(elem->point(5), elem->point(2)))
                        {
                          // Case 5
                          eta_mapped  = -eta;
                          zeta_mapped = -zeta;
                        }
                      else
                        {
                          // Case 6
                          eta_mapped  = -zeta;
                          zeta_mapped = -eta;
                        }

                    else if (elem->point(5) == min_point)
                      {
                        if (elem->point(1) == std::min(elem->point(1), elem->point(6)))
                          {
                            // Case 7
                            eta_mapped  = -zeta;
                            zeta_mapped = eta;
                          }
                        else
                          {
                            // Case 8
                            eta_mapped   = eta;
                            zeta_mapped = -zeta;
                          }
                      }
                  }


                // Face 3
                else if ((i1[i] == 1) && (i0[i] >= 2) && (i2[i] >= 2))
                  {
                    const Point min_point = std::min(elem->point(2),
                                                     std::min(elem->point(3),
                                                              std::min(elem->point(7),
                                                                       elem->point(6))));
                    if (elem->point(3) == min_point)
                      if (elem->point(2) == std::min(elem->point(2), elem->point(7)))
                        {
                          // Case 1
                          xi_mapped   = xi;
                          zeta_mapped = zeta;
                        }
                      else
                        {
                          // Case 2
                          xi_mapped   = zeta;
                          zeta_mapped = xi;
                        }

                    else if (elem->point(7) == min_point)
                      if (elem->point(3) == std::min(elem->point(3), elem->point(6)))
                        {
                          // Case 3
                          xi_mapped   = -zeta;
                          zeta_mapped = xi;
                        }
                      else
                        {
                          // Case 4
                          xi_mapped   = xi;
                          zeta_mapped = -zeta;
                        }

                    else if (elem->point(6) == min_point)
                      if (elem->point(7) == std::min(elem->point(7), elem->point(2)))
                        {
                          // Case 5
                          xi_mapped   = -xi;
                          zeta_mapped = -zeta;
                        }
                      else
                        {
                          // Case 6
                          xi_mapped   = -zeta;
                          zeta_mapped = -xi;
                        }

                    else if (elem->point(2) == min_point)
                      {
                        if (elem->point(6) == std::min(elem->point(3), elem->point(6)))
                          {
                            // Case 7
                            xi_mapped   = zeta;
                            zeta_mapped = -xi;
                          }
                        else
                          {
                            // Case 8
                            xi_mapped   = -xi;
                            zeta_mapped = zeta;
                          }
                      }
                  }


                // Face 4
                else if ((i0[i] == 0) && (i1[i] >= 2) && (i2[i] >= 2))
                  {
                    const Point min_point = std::min(elem->point(3),
                                                     std::min(elem->point(0),
                                                              std::min(elem->point(4),
                                                                       elem->point(7))));
                    if (elem->point(0) == min_point)
                      if (elem->point(3) == std::min(elem->point(3), elem->point(4)))
                        {
                          // Case 1
                          eta_mapped  = eta;
                          zeta_mapped = zeta;
                        }
                      else
                        {
                          // Case 2
                          eta_mapped  = zeta;
                          zeta_mapped = eta;
                        }

                    else if (elem->point(4) == min_point)
                      if (elem->point(0) == std::min(elem->point(0), elem->point(7)))
                        {
                          // Case 3
                          eta_mapped  = -zeta;
                          zeta_mapped = eta;
                        }
                      else
                        {
                          // Case 4
                          eta_mapped  = eta;
                          zeta_mapped = -zeta;
                        }

                    else if (elem->point(7) == min_point)
                      if (elem->point(4) == std::min(elem->point(4), elem->point(3)))
                        {
                          // Case 5
                          eta_mapped  = -eta;
                          zeta_mapped = -zeta;
                        }
                      else
                        {
                          // Case 6
                          eta_mapped  = -zeta;
                          zeta_mapped = -eta;
                        }

                    else if (elem->point(3) == min_point)
                      {
                        if (elem->point(7) == std::min(elem->point(7), elem->point(0)))
                          {
                            // Case 7
                            eta_mapped   = zeta;
                            zeta_mapped = -eta;
                          }
                        else
                          {
                            // Case 8
                            eta_mapped  = -eta;
                            zeta_mapped = zeta;
                          }
                      }
                  }


                // Face 5
                else if ((i2[i] == 1) && (i0[i] >= 2) && (i1[i] >= 2))
                  {
                    const Point min_point = std::min(elem->point(4),
                                                     std::min(elem->point(5),
                                                              std::min(elem->point(6),
                                                                       elem->point(7))));
                    if (elem->point(4) == min_point)
                      if (elem->point(5) == std::min(elem->point(5), elem->point(7)))
                        {
                          // Case 1
                          xi_mapped  = xi;
                          eta_mapped = eta;
                        }
                      else
                        {
                          // Case 2
                          xi_mapped  = eta;
                          eta_mapped = xi;
                        }

                    else if (elem->point(5) == min_point)
                      if (elem->point(6) == std::min(elem->point(6), elem->point(4)))
                        {
                          // Case 3
                          xi_mapped  = eta;
                          eta_mapped = -xi;
                        }
                      else
                        {
                          // Case 4
                          xi_mapped  = -xi;
                          eta_mapped = eta;
                        }

                    else if (elem->point(6) == min_point)
                      if (elem->point(7) == std::min(elem->point(7), elem->point(5)))
                        {
                          // Case 5
                          xi_mapped  = -xi;
                          eta_mapped = -eta;
                        }
                      else
                        {
                          // Case 6
                          xi_mapped  = -eta;
                          eta_mapped = -xi;
                        }

                    else if (elem->point(7) == min_point)
                      {
                        if (elem->point(4) == std::min(elem->point(4), elem->point(6)))
                          {
                            // Case 7
                            xi_mapped  = -eta;
                            eta_mapped = xi;
                          }
                        else
                          {
                            // Case 8
                            xi_mapped  = xi;
                            eta_mapped = eta;
                          }
                      }
                  }
              }

              return (FEShim<1,BERNSTEIN,RealType>::shape(EDGE3, totalorder, i0[i], xi_mapped)*
                      FEShim<1,BERNSTEIN,RealType>::shape(EDGE3, totalorder, i1[i], eta_mapped)*
                      FEShim<1,BERNSTEIN,RealType>::shape(EDGE3, totalorder, i2[i], zeta_mapped));
            }


          default:
            libmesh_error_msg("Invalid element type = " << type);
          } //case HEX27

      }//case THIRD


      // 4th-order Bernstein.
    case FOURTH:
      {
        switch (type)
          {

            // Bernstein shape functions on the hexahedral.
          case HEX27:
            {
              libmesh_assert_less (i, 125);

              // Compute hex shape functions as a tensor-product
              const auto xi    = p(0);
              const auto eta   = p(1);
              const auto zeta  = p(2);
              auto xi_mapped   = p(0);
              auto eta_mapped  = p(1);
              auto zeta_mapped = p(2);

              // The only way to make any sense of this
              // is to look at the mgflo/mg2/mgf documentation
              // and make the cut-out cube!
              //  Nodes                         0  1  2  3  4  5  6  7  8  8  9  9 10 10 11 11 12 12 13 13 14 14 15 15 16 16 17 17 18 18 19 19 20 20 20 20 21 21 21 21 22 22 22 22 23 23 23 23 24 24 24 24 25 25 25 25 26 26 26 26 26 26 26 26
              //  DOFS                          0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 18 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 60 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24
              static const unsigned int i0[] = {0, 1, 1, 0, 0, 1, 1, 0, 2, 3, 4, 1, 1, 1, 2, 3, 4, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 2, 3, 4, 1, 1, 1, 2, 3, 4, 0, 0, 0, 2, 3, 4, 2, 3, 4, 2, 3, 4, 2, 3, 4, 2, 3, 4, 2, 3, 4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 4, 2, 3, 4, 2, 3, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 3, 4, 2, 3, 4, 2, 3, 4, 2, 3, 4, 2, 3, 4, 2, 3, 4, 2, 3, 4, 2, 3, 4, 2, 3, 4, 2, 3, 4, 2, 3, 4, 2, 3, 4};
              static const unsigned int i1[] = {0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 2, 3, 4, 1, 1, 1, 2, 3, 4, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 2, 3, 4, 1, 1, 1, 2, 3, 4, 2, 2, 2, 3, 3, 3, 4, 4, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 3, 4, 2, 3, 4, 2, 3, 4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 4, 2, 3, 4, 2, 3, 4, 2, 2, 2, 3, 3, 3, 4, 4, 4, 2, 2, 2, 3, 3, 3, 4, 4, 4, 2, 2, 2, 3, 3, 3, 4, 4, 4, 2, 2, 2, 3, 3, 3, 4, 4, 4};
              static const unsigned int i2[] = {0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 3, 4, 2, 3, 4, 2, 3, 4, 2, 3, 4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 3, 3, 3, 4, 4, 4, 2, 2, 2, 3, 3, 3, 4, 4, 4, 2, 2, 2, 3, 3, 3, 4, 4, 4, 2, 2, 2, 3, 3, 3, 4, 4, 4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4};



              // handle the edge orientation
              {
                // Edge 0
                if ((i0[i] >= 2) && (i1[i] == 0) && (i2[i] == 0))
                  {
                    if (elem->point(0) != std::min(elem->point(0), elem->point(1)))
                      xi_mapped = -xi;
                  }
                // Edge 1
                else if ((i0[i] == 1) && (i1[i] >= 2) && (i2[i] == 0))
                  {
                    if (elem->point(1) != std::min(elem->point(1), elem->point(2)))
                      eta_mapped = -eta;
                  }
                // Edge 2
                else if ((i0[i] >= 2) && (i1[i] == 1) && (i2[i] == 0))
                  {
                    if (elem->point(3) != std::min(elem->point(3), elem->point(2)))
                      xi_mapped = -xi;
                  }
                // Edge 3
                else if ((i0[i] == 0) && (i1[i] >= 2) && (i2[i] == 0))
                  {
                    if (elem->point(0) != std::min(elem->point(0), elem->point(3)))
                      eta_mapped = -eta;
                  }
                // Edge 4
                else if ((i0[i] == 0) && (i1[i] == 0) && (i2[i] >=2 ))
                  {
                    if (elem->point(0) != std::min(elem->point(0), elem->point(4)))
                      zeta_mapped = -zeta;
                  }
                // Edge 5
                else if ((i0[i] == 1) && (i1[i] == 0) && (i2[i] >=2 ))
                  {
                    if (elem->point(1) != std::min(elem->point(1), elem->point(5)))
                      zeta_mapped = -zeta;
                  }
                // Edge 6
                else if ((i0[i] == 1) && (i1[i] == 1) && (i2[i] >=2 ))
                  {
                    if (elem->point(2) != std::min(elem->point(2), elem->point(6)))
                      zeta_mapped = -zeta;
                  }
                // Edge 7
                else if ((i0[i] == 0) && (i1[i] == 1) && (i2[i] >=2 ))
                  {
                    if (elem->point(3) != std::min(elem->point(3), elem->point(7)))
                      zeta_mapped = -zeta;
                  }
                // Edge 8
                else if ((i0[i] >=2 ) && (i1[i] == 0) && (i2[i] == 1))
                  {
                    if (elem->point(4) != std::min(elem->point(4), elem->point(5)))
                      xi_mapped = -xi;
                  }
                // Edge 9
                else if ((i0[i] == 1) && (i1[i] >=2 ) && (i2[i] == 1))
                  {
                    if (elem->point(5) != std::min(elem->point(5), elem->point(6)))
                      eta_mapped = -eta;
                  }
                // Edge 10
                else if ((i0[i] >=2 ) && (i1[i] == 1) && (i2[i] == 1))
                  {
                    if (elem->point(7) != std::min(elem->point(7), elem->point(6)))
                      xi_mapped = -xi;
                  }
                // Edge 11
                else if ((i0[i] == 0) && (i1[i] >=2 ) && (i2[i] == 1))
                  {
                    if (elem->point(4) != std::min(elem->point(4), elem->point(7)))
                      eta_mapped = -eta;
                  }
              }


              // handle the face orientation
              {
                // Face 0
                if ((i2[i] == 0) && (i0[i] >= 2) && (i1[i] >= 2))
                  {
                    const Point min_point = std::min(elem->point(1),
                                                     std::min(elem->point(2),
                                                              std::min(elem->point(0),
                                                                       elem->point(3))));
                    if (elem->point(0) == min_point)
                      if (elem->point(1) == std::min(elem->point(1), elem->point(3)))
                        {
                          // Case 1
                          xi_mapped  = xi;
                          eta_mapped = eta;
                        }
                      else
                        {
                          // Case 2
                          xi_mapped  = eta;
                          eta_mapped = xi;
                        }

                    else if (elem->point(3) == min_point)
                      if (elem->point(0) == std::min(elem->point(0), elem->point(2)))
                        {
                          // Case 3
                          xi_mapped  = -eta;
                          eta_mapped = xi;
                        }
                      else
                        {
                          // Case 4
                          xi_mapped  = xi;
                          eta_mapped = -eta;
                        }

                    else if (elem->point(2) == min_point)
                      if (elem->point(3) == std::min(elem->point(3), elem->point(1)))
                        {
                          // Case 5
                          xi_mapped  = -xi;
                          eta_mapped = -eta;
                        }
                      else
                        {
                          // Case 6
                          xi_mapped  = -eta;
                          eta_mapped = -xi;
                        }

                    else if (elem->point(1) == min_point)
                      {
                        if (elem->point(2) == std::min(elem->point(2), elem->point(0)))
                          {
                            // Case 7
                            xi_mapped  = eta;
                            eta_mapped = -xi;
                          }
                        else
                          {
                            // Case 8
                            xi_mapped  = -xi;
                            eta_mapped = eta;
                          }
                      }
                  }


                // Face 1
                else if ((i1[i] == 0) && (i0[i] >= 2) && (i2[i] >= 2))
                  {
                    const Point min_point = std::min(elem->point(0),
                                                     std::min(elem->point(1),
                                                              std::min(elem->point(5),
                                                                       elem->point(4))));
                    if (elem->point(0) == min_point)
                      if (elem->point(1) == std::min(elem->point(1), elem->point(4)))
                        {
                          // Case 1
                          xi_mapped   = xi;
                          zeta_mapped = zeta;
                        }
                      else
                        {
                          // Case 2
                          xi_mapped   = zeta;
                          zeta_mapped = xi;
                        }

                    else if (elem->point(1) == min_point)
                      if (elem->point(5) == std::min(elem->point(5), elem->point(0)))
                        {
                          // Case 3
                          xi_mapped   = zeta;
                          zeta_mapped = -xi;
                        }
                      else
                        {
                          // Case 4
                          xi_mapped   = -xi;
                          zeta_mapped = zeta;
                        }

                    else if (elem->point(5) == min_point)
                      if (elem->point(4) == std::min(elem->point(4), elem->point(1)))
                        {
                          // Case 5
                          xi_mapped   = -xi;
                          zeta_mapped = -zeta;
                        }
                      else
                        {
                          // Case 6
                          xi_mapped   = -zeta;
                          zeta_mapped = -xi;
                        }

                    else if (elem->point(4) == min_point)
                      {
                        if (elem->point(0) == std::min(elem->point(0), elem->point(5)))
                          {
                            // Case 7
                            xi_mapped   = -xi;
                            zeta_mapped = zeta;
                          }
                        else
                          {
                            // Case 8
                            xi_mapped   = xi;
                            zeta_mapped = -zeta;
                          }
                      }
                  }


                // Face 2
                else if ((i0[i] == 1) && (i1[i] >= 2) && (i2[i] >= 2))
                  {
                    const Point min_point = std::min(elem->point(1),
                                                     std::min(elem->point(2),
                                                              std::min(elem->point(6),
                                                                       elem->point(5))));
                    if (elem->point(1) == min_point)
                      if (elem->point(2) == std::min(elem->point(2), elem->point(5)))
                        {
                          // Case 1
                          eta_mapped  = eta;
                          zeta_mapped = zeta;
                        }
                      else
                        {
                          // Case 2
                          eta_mapped  = zeta;
                          zeta_mapped = eta;
                        }

                    else if (elem->point(2) == min_point)
                      if (elem->point(6) == std::min(elem->point(6), elem->point(1)))
                        {
                          // Case 3
                          eta_mapped  = zeta;
                          zeta_mapped = -eta;
                        }
                      else
                        {
                          // Case 4
                          eta_mapped  = -eta;
                          zeta_mapped = zeta;
                        }

                    else if (elem->point(6) == min_point)
                      if (elem->point(5) == std::min(elem->point(5), elem->point(2)))
                        {
                          // Case 5
                          eta_mapped  = -eta;
                          zeta_mapped = -zeta;
                        }
                      else
                        {
                          // Case 6
                          eta_mapped  = -zeta;
                          zeta_mapped = -eta;
                        }

                    else if (elem->point(5) == min_point)
                      {
                        if (elem->point(1) == std::min(elem->point(1), elem->point(6)))
                          {
                            // Case 7
                            eta_mapped  = -zeta;
                            zeta_mapped = eta;
                          }
                        else
                          {
                            // Case 8
                            eta_mapped   = eta;
                            zeta_mapped = -zeta;
                          }
                      }
                  }


                // Face 3
                else if ((i1[i] == 1) && (i0[i] >= 2) && (i2[i] >= 2))
                  {
                    const Point min_point = std::min(elem->point(2),
                                                     std::min(elem->point(3),
                                                              std::min(elem->point(7),
                                                                       elem->point(6))));
                    if (elem->point(3) == min_point)
                      if (elem->point(2) == std::min(elem->point(2), elem->point(7)))
                        {
                          // Case 1
                          xi_mapped   = xi;
                          zeta_mapped = zeta;
                        }
                      else
                        {
                          // Case 2
                          xi_mapped   = zeta;
                          zeta_mapped = xi;
                        }

                    else if (elem->point(7) == min_point)
                      if (elem->point(3) == std::min(elem->point(3), elem->point(6)))
                        {
                          // Case 3
                          xi_mapped   = -zeta;
                          zeta_mapped = xi;
                        }
                      else
                        {
                          // Case 4
                          xi_mapped   = xi;
                          zeta_mapped = -zeta;
                        }

                    else if (elem->point(6) == min_point)
                      if (elem->point(7) == std::min(elem->point(7), elem->point(2)))
                        {
                          // Case 5
                          xi_mapped   = -xi;
                          zeta_mapped = -zeta;
                        }
                      else
                        {
                          // Case 6
                          xi_mapped   = -zeta;
                          zeta_mapped = -xi;
                        }

                    else if (elem->point(2) == min_point)
                      {
                        if (elem->point(6) == std::min(elem->point(3), elem->point(6)))
                          {
                            // Case 7
                            xi_mapped   = zeta;
                            zeta_mapped = -xi;
                          }
                        else
                          {
                            // Case 8
                            xi_mapped   = -xi;
                            zeta_mapped = zeta;
                          }
                      }
                  }


                // Face 4
                else if ((i0[i] == 0) && (i1[i] >= 2) && (i2[i] >= 2))
                  {
                    const Point min_point = std::min(elem->point(3),
                                                     std::min(elem->point(0),
                                                              std::min(elem->point(4),
                                                                       elem->point(7))));
                    if (elem->point(0) == min_point)
                      if (elem->point(3) == std::min(elem->point(3), elem->point(4)))
                        {
                          // Case 1
                          eta_mapped  = eta;
                          zeta_mapped = zeta;
                        }
                      else
                        {
                          // Case 2
                          eta_mapped  = zeta;
                          zeta_mapped = eta;
                        }

                    else if (elem->point(4) == min_point)
                      if (elem->point(0) == std::min(elem->point(0), elem->point(7)))
                        {
                          // Case 3
                          eta_mapped  = -zeta;
                          zeta_mapped = eta;
                        }
                      else
                        {
                          // Case 4
                          eta_mapped  = eta;
                          zeta_mapped = -zeta;
                        }

                    else if (elem->point(7) == min_point)
                      if (elem->point(4) == std::min(elem->point(4), elem->point(3)))
                        {
                          // Case 5
                          eta_mapped  = -eta;
                          zeta_mapped = -zeta;
                        }
                      else
                        {
                          // Case 6
                          eta_mapped  = -zeta;
                          zeta_mapped = -eta;
                        }

                    else if (elem->point(3) == min_point)
                      {
                        if (elem->point(7) == std::min(elem->point(7), elem->point(0)))
                          {
                            // Case 7
                            eta_mapped   = zeta;
                            zeta_mapped = -eta;
                          }
                        else
                          {
                            // Case 8
                            eta_mapped  = -eta;
                            zeta_mapped = zeta;
                          }
                      }
                  }


                // Face 5
                else if ((i2[i] == 1) && (i0[i] >= 2) && (i1[i] >= 2))
                  {
                    const Point min_point = std::min(elem->point(4),
                                                     std::min(elem->point(5),
                                                              std::min(elem->point(6),
                                                                       elem->point(7))));
                    if (elem->point(4) == min_point)
                      if (elem->point(5) == std::min(elem->point(5), elem->point(7)))
                        {
                          // Case 1
                          xi_mapped  = xi;
                          eta_mapped = eta;
                        }
                      else
                        {
                          // Case 2
                          xi_mapped  = eta;
                          eta_mapped = xi;
                        }

                    else if (elem->point(5) == min_point)
                      if (elem->point(6) == std::min(elem->point(6), elem->point(4)))
                        {
                          // Case 3
                          xi_mapped  = eta;
                          eta_mapped = -xi;
                        }
                      else
                        {
                          // Case 4
                          xi_mapped  = -xi;
                          eta_mapped = eta;
                        }

                    else if (elem->point(6) == min_point)
                      if (elem->point(7) == std::min(elem->point(7), elem->point(5)))
                        {
                          // Case 5
                          xi_mapped  = -xi;
                          eta_mapped = -eta;
                        }
                      else
                        {
                          // Case 6
                          xi_mapped  = -eta;
                          eta_mapped = -xi;
                        }

                    else if (elem->point(7) == min_point)
                      {
                        if (elem->point(4) == std::min(elem->point(4), elem->point(6)))
                          {
                            // Case 7
                            xi_mapped  = -eta;
                            eta_mapped = xi;
                          }
                        else
                          {
                            // Case 8
                            xi_mapped  = xi;
                            eta_mapped = eta;
                          }
                      }
                  }


              }


              return (FEShim<1,BERNSTEIN,RealType>::shape(EDGE3, totalorder, i0[i], xi_mapped)*
                      FEShim<1,BERNSTEIN,RealType>::shape(EDGE3, totalorder, i1[i], eta_mapped)*
                      FEShim<1,BERNSTEIN,RealType>::shape(EDGE3, totalorder, i2[i], zeta_mapped));
            }


          default:
            libmesh_error_msg("Invalid element type = " << type);
          }
      }


    default:
      libmesh_error_msg("Invalid totalorder = " << totalorder);
    }
#else // LIBMESH_DIM != 3
  libmesh_ignore(elem, order, i, p, add_p_level);
  libmesh_not_implemented();
#endif
}




template <typename RealType>
typename FEShim<3,BERNSTEIN,RealType>::OutputShape FEShim<3,BERNSTEIN,RealType>::shape_deriv(const ElemType,
                                  const Order,
                                  const unsigned int,
                                  const unsigned int,
                                  const Point & )
{
  libmesh_error_msg("Bernstein polynomials require the element type \nbecause edge and face orientation is needed.");
  return 0.;
}



template <typename RealType>
typename FEShim<3,BERNSTEIN,RealType>::OutputShape FEShim<3,BERNSTEIN,RealType>::shape_deriv(const Elem * elem,
                                  const Order order,
                                  const unsigned int i,
                                  const unsigned int j,
                                  const Point & p,
                                  const bool add_p_level)
{

#if LIBMESH_DIM == 3
  libmesh_assert(elem);
  const ElemType type = elem->type();

  const Order totalorder =
    static_cast<Order>(order + add_p_level * elem->p_level());

  libmesh_assert_less (j, 3);

  switch (totalorder)
    {
      // 1st order Bernstein.
    case FIRST:
      {
        switch (type)
          {
            // Bernstein shape functions on the tetrahedron.
          case TET4:
          case TET10:
            {
              // I have been lazy here and am using finite differences
              // to compute the derivatives!
              const Real eps = 1.e-4;

              libmesh_assert_less (i, 4);


              switch (j)
                {
                  //  d()/dxi
                case 0:
                  {
                    const Point pp(p(0)+eps, p(1), p(2));
                    const Point pm(p(0)-eps, p(1), p(2));

                    return (FEShim<3,BERNSTEIN,RealType>::shape(elem, order, i, pp) -
                            FEShim<3,BERNSTEIN,RealType>::shape(elem, order, i, pm))/2./eps;
                  }

                  // d()/deta
                case 1:
                  {
                    const Point pp(p(0), p(1)+eps, p(2));
                    const Point pm(p(0), p(1)-eps, p(2));

                    return (FEShim<3,BERNSTEIN,RealType>::shape(elem, order, i, pp) -
                            FEShim<3,BERNSTEIN,RealType>::shape(elem, order, i, pm))/2./eps;
                  }
                  // d()/dzeta
                case 2:
                  {
                    const Point pp(p(0), p(1), p(2)+eps);
                    const Point pm(p(0), p(1), p(2)-eps);

                    return (FEShim<3,BERNSTEIN,RealType>::shape(elem, order, i, pp) -
                            FEShim<3,BERNSTEIN,RealType>::shape(elem, order, i, pm))/2./eps;
                  }
                default:
                  libmesh_error_msg("Invalid derivative index j = " << j);
                }
            }




            // Bernstein shape functions on the hexahedral.
          case HEX8:
          case HEX20:
          case HEX27:
            {
              libmesh_assert_less (i, 8);

              // Compute hex shape functions as a tensor-product
              const auto xi   = p(0);
              const auto eta  = p(1);
              const auto zeta = p(2);

              // The only way to make any sense of this
              // is to look at the mgflo/mg2/mgf documentation
              // and make the cut-out cube!
              //                                0  1  2  3  4  5  6  7
              static const unsigned int i0[] = {0, 1, 1, 0, 0, 1, 1, 0};
              static const unsigned int i1[] = {0, 0, 1, 1, 0, 0, 1, 1};
              static const unsigned int i2[] = {0, 0, 0, 0, 1, 1, 1, 1};

              switch (j)
                {
                  // d()/dxi
                case 0:
                  return (FEShim<1,BERNSTEIN,RealType>::shape_deriv(EDGE3, totalorder, i0[i], 0, xi)*
                          FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i1[i],    eta)*
                          FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i2[i],    zeta));

                  // d()/deta
                case 1:
                  return (FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i0[i],     xi)*
                          FEShim<1,BERNSTEIN,RealType>::shape_deriv(EDGE3, totalorder, i1[i], 0, eta)*
                          FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i2[i],    zeta));

                  // d()/dzeta
                case 2:
                  return (FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i0[i],    xi)*
                          FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i1[i],    eta)*
                          FEShim<1,BERNSTEIN,RealType>::shape_deriv(EDGE3, totalorder, i2[i], 0, zeta));

                default:
                  libmesh_error_msg("Invalid derivative index j = " << j);
                }
            }

          default:
            libmesh_error_msg("Invalid element type = " << type);
          }
      }




    case SECOND:
      {
        switch (type)
          {
            // Bernstein shape functions on the tetrahedron.
          case TET10:
            {
              // I have been lazy here and am using finite differences
              // to compute the derivatives!
              const Real eps = 1.e-4;

              libmesh_assert_less (i, 10);


              switch (j)
                {
                  //  d()/dxi
                case 0:
                  {
                    const Point pp(p(0)+eps, p(1), p(2));
                    const Point pm(p(0)-eps, p(1), p(2));

                    return (FEShim<3,BERNSTEIN,RealType>::shape(elem, order, i, pp) -
                            FEShim<3,BERNSTEIN,RealType>::shape(elem, order, i, pm))/2./eps;
                  }

                  // d()/deta
                case 1:
                  {
                    const Point pp(p(0), p(1)+eps, p(2));
                    const Point pm(p(0), p(1)-eps, p(2));

                    return (FEShim<3,BERNSTEIN,RealType>::shape(elem, order, i, pp) -
                            FEShim<3,BERNSTEIN,RealType>::shape(elem, order, i, pm))/2./eps;
                  }
                  // d()/dzeta
                case 2:
                  {
                    const Point pp(p(0), p(1), p(2)+eps);
                    const Point pm(p(0), p(1), p(2)-eps);

                    return (FEShim<3,BERNSTEIN,RealType>::shape(elem, order, i, pp) -
                            FEShim<3,BERNSTEIN,RealType>::shape(elem, order, i, pm))/2./eps;
                  }
                default:
                  libmesh_error_msg("Invalid derivative index j = " << j);
                }
            }

            // Bernstein shape functions on the hexahedral.
          case HEX20:
            {
              libmesh_assert_less (i, 20);

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
              static const Real scal20[] =     {-0.25, -0.25, -0.25, -0.25, 0,     0,     0,     0,     0.5,   0.5,   0.5,   0.5,   0,     0,     0,     0,     0,     0,     0,     0};
              static const Real scal21[] =     {-0.25, -0.25, 0,     0,     -0.25, -0.25, 0,     0,     0.5,   0,     0,     0,     0.5,   0.5,   0,     0,     0.5,   0,     0,     0};
              static const Real scal22[] =     {0,     -0.25, -0.25, 0,     0,     -0.25, -0.25, 0,     0,     0.5,   0,     0,     0,     0.5,   0.5,   0,     0,     0.5,   0,     0};
              static const Real scal23[] =     {0,     0,     -0.25, -0.25, 0,     0,     -0.25, -0.25, 0,     0,     0.5,   0,     0,     0,     0.5,   0.5,   0,     0,     0.5,   0};
              static const Real scal24[] =     {-0.25, 0,     0,     -0.25, -0.25, 0,     0,     -0.25, 0,     0,     0,     0.5,   0.5,   0,     0,     0.5,   0,     0,     0,     0.5};
              static const Real scal25[] =     {0,     0,     0,     0,     -0.25, -0.25, -0.25, -0.25, 0,     0,     0,     0,     0,     0,     0,     0,     0.5,   0.5,   0.5,   0.5};
              static const Real scal26[] =     {-0.25, -0.25, -0.25, -0.25, -0.25, -0.25, -0.25, -0.25, 0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25};

              switch (j)
                {
                  // d()/dxi
                case 0:
                  return (FEShim<1,BERNSTEIN,RealType>::shape_deriv(EDGE3, totalorder, i0[i], 0, xi)*
                          FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i1[i],    eta)*
                          FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i2[i],    zeta)
                          +scal20[i]*
                          FEShim<1,BERNSTEIN,RealType>::shape_deriv(EDGE3, totalorder, i0[20], 0, xi)*
                          FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i1[20],    eta)*
                          FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i2[20],    zeta)
                          +scal21[i]*
                          FEShim<1,BERNSTEIN,RealType>::shape_deriv(EDGE3, totalorder, i0[21], 0, xi)*
                          FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i1[21],    eta)*
                          FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i2[21],    zeta)
                          +scal22[i]*
                          FEShim<1,BERNSTEIN,RealType>::shape_deriv(EDGE3, totalorder, i0[22], 0, xi)*
                          FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i1[22],    eta)*
                          FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i2[22],    zeta)
                          +scal23[i]*
                          FEShim<1,BERNSTEIN,RealType>::shape_deriv(EDGE3, totalorder, i0[23], 0, xi)*
                          FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i1[23],    eta)*
                          FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i2[23],    zeta)
                          +scal24[i]*
                          FEShim<1,BERNSTEIN,RealType>::shape_deriv(EDGE3, totalorder, i0[24], 0, xi)*
                          FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i1[24],    eta)*
                          FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i2[24],    zeta)
                          +scal25[i]*
                          FEShim<1,BERNSTEIN,RealType>::shape_deriv(EDGE3, totalorder, i0[25], 0, xi)*
                          FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i1[25],    eta)*
                          FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i2[25],    zeta)
                          +scal26[i]*
                          FEShim<1,BERNSTEIN,RealType>::shape_deriv(EDGE3, totalorder, i0[26], 0, xi)*
                          FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i1[26],    eta)*
                          FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i2[26],    zeta));

                  // d()/deta
                case 1:
                  return (FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i0[i],     xi)*
                          FEShim<1,BERNSTEIN,RealType>::shape_deriv(EDGE3, totalorder, i1[i], 0, eta)*
                          FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i2[i],    zeta)
                          +scal20[i]*
                          FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i0[20],     xi)*
                          FEShim<1,BERNSTEIN,RealType>::shape_deriv(EDGE3, totalorder, i1[20], 0, eta)*
                          FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i2[20],    zeta)
                          +scal21[i]*
                          FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i0[21],     xi)*
                          FEShim<1,BERNSTEIN,RealType>::shape_deriv(EDGE3, totalorder, i1[21], 0, eta)*
                          FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i2[21],    zeta)
                          +scal22[i]*
                          FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i0[22],     xi)*
                          FEShim<1,BERNSTEIN,RealType>::shape_deriv(EDGE3, totalorder, i1[22], 0, eta)*
                          FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i2[22],    zeta)
                          +scal23[i]*
                          FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i0[23],     xi)*
                          FEShim<1,BERNSTEIN,RealType>::shape_deriv(EDGE3, totalorder, i1[23], 0, eta)*
                          FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i2[23],    zeta)
                          +scal24[i]*
                          FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i0[24],     xi)*
                          FEShim<1,BERNSTEIN,RealType>::shape_deriv(EDGE3, totalorder, i1[24], 0, eta)*
                          FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i2[24],    zeta)
                          +scal25[i]*
                          FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i0[25],     xi)*
                          FEShim<1,BERNSTEIN,RealType>::shape_deriv(EDGE3, totalorder, i1[25], 0, eta)*
                          FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i2[25],    zeta)
                          +scal26[i]*
                          FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i0[26],     xi)*
                          FEShim<1,BERNSTEIN,RealType>::shape_deriv(EDGE3, totalorder, i1[26], 0, eta)*
                          FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i2[26],    zeta));

                  // d()/dzeta
                case 2:
                  return (FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i0[i],    xi)*
                          FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i1[i],    eta)*
                          FEShim<1,BERNSTEIN,RealType>::shape_deriv(EDGE3, totalorder, i2[i], 0, zeta)
                          +scal20[i]*
                          FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i0[20],    xi)*
                          FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i1[20],    eta)*
                          FEShim<1,BERNSTEIN,RealType>::shape_deriv(EDGE3, totalorder, i2[20], 0, zeta)
                          +scal21[i]*
                          FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i0[21],    xi)*
                          FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i1[21],    eta)*
                          FEShim<1,BERNSTEIN,RealType>::shape_deriv(EDGE3, totalorder, i2[21], 0, zeta)
                          +scal22[i]*
                          FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i0[22],    xi)*
                          FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i1[22],    eta)*
                          FEShim<1,BERNSTEIN,RealType>::shape_deriv(EDGE3, totalorder, i2[22], 0, zeta)
                          +scal23[i]*
                          FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i0[23],    xi)*
                          FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i1[23],    eta)*
                          FEShim<1,BERNSTEIN,RealType>::shape_deriv(EDGE3, totalorder, i2[23], 0, zeta)
                          +scal24[i]*
                          FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i0[24],    xi)*
                          FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i1[24],    eta)*
                          FEShim<1,BERNSTEIN,RealType>::shape_deriv(EDGE3, totalorder, i2[24], 0, zeta)
                          +scal25[i]*
                          FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i0[25],    xi)*
                          FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i1[25],    eta)*
                          FEShim<1,BERNSTEIN,RealType>::shape_deriv(EDGE3, totalorder, i2[25], 0, zeta)
                          +scal26[i]*
                          FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i0[26],    xi)*
                          FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i1[26],    eta)*
                          FEShim<1,BERNSTEIN,RealType>::shape_deriv(EDGE3, totalorder, i2[26], 0, zeta));

                default:
                  libmesh_error_msg("Invalid derivative index j = " << j);
                }
            }

            // Bernstein shape functions on the hexahedral.
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

              switch (j)
                {
                  // d()/dxi
                case 0:
                  return (FEShim<1,BERNSTEIN,RealType>::shape_deriv(EDGE3, totalorder, i0[i], 0, xi)*
                          FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i1[i],    eta)*
                          FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i2[i],    zeta));

                  // d()/deta
                case 1:
                  return (FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i0[i],     xi)*
                          FEShim<1,BERNSTEIN,RealType>::shape_deriv(EDGE3, totalorder, i1[i], 0, eta)*
                          FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i2[i],    zeta));

                  // d()/dzeta
                case 2:
                  return (FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i0[i],    xi)*
                          FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i1[i],    eta)*
                          FEShim<1,BERNSTEIN,RealType>::shape_deriv(EDGE3, totalorder, i2[i], 0, zeta));

                default:
                  libmesh_error_msg("Invalid derivative index j = " << j);
                }
            }


          default:
            libmesh_error_msg("Invalid element type = " << type);
          }
      }



      // 3rd-order Bernstein.
    case THIRD:
      {
        switch (type)
          {

            //     // Bernstein shape functions derivatives.
            //   case TET10:
            //     {
            //       // I have been lazy here and am using finite differences
            //       // to compute the derivatives!
            //       const Real eps = 1.e-4;

            //       libmesh_assert_less (i, 20);
            //       libmesh_assert_less (j, 3);

            // switch (j)
            // {
            //   //  d()/dxi
            // case 0:
            //   {
            //     const Point pp(p(0)+eps, p(1), p(2));
            //     const Point pm(p(0)-eps, p(1), p(2));

            //     return (FEShim<3,BERNSTEIN,RealType>::shape(elem, order, i, pp) -
            //     FEShim<3,BERNSTEIN,RealType>::shape(elem, order, i, pm))/2./eps;
            //   }

            //   // d()/deta
            // case 1:
            //   {
            //     const Point pp(p(0), p(1)+eps, p(2));
            //     const Point pm(p(0), p(1)-eps, p(2));

            //     return (FEShim<3,BERNSTEIN,RealType>::shape(elem, order, i, pp) -
            //     FEShim<3,BERNSTEIN,RealType>::shape(elem, order, i, pm))/2./eps;
            //   }
            //   // d()/dzeta
            // case 2:
            //   {
            //     const Point pp(p(0), p(1), p(2)+eps);
            //     const Point pm(p(0), p(1), p(2)-eps);

            //     return (FEShim<3,BERNSTEIN,RealType>::shape(elem, order, i, pp) -
            //     FEShim<3,BERNSTEIN,RealType>::shape(elem, order, i, pm))/2./eps;
            //   }
            // default:
            // libmesh_error_msg("Invalid derivative index j = " << j);
            // }


            //     }


            // Bernstein shape functions on the hexahedral.
          case HEX27:
            {
              // I have been lazy here and am using finite differences
              // to compute the derivatives!
              const Real eps = 1.e-4;

              libmesh_assert_less (i, 64);

              switch (j)
                {
                  //  d()/dxi
                case 0:
                  {
                    const Point pp(p(0)+eps, p(1), p(2));
                    const Point pm(p(0)-eps, p(1), p(2));

                    return (FEShim<3,BERNSTEIN,RealType>::shape(elem, order, i, pp) -
                            FEShim<3,BERNSTEIN,RealType>::shape(elem, order, i, pm))/2./eps;
                  }

                  // d()/deta
                case 1:
                  {
                    const Point pp(p(0), p(1)+eps, p(2));
                    const Point pm(p(0), p(1)-eps, p(2));

                    return (FEShim<3,BERNSTEIN,RealType>::shape(elem, order, i, pp) -
                            FEShim<3,BERNSTEIN,RealType>::shape(elem, order, i, pm))/2./eps;
                  }
                  // d()/dzeta
                case 2:
                  {
                    const Point pp(p(0), p(1), p(2)+eps);
                    const Point pm(p(0), p(1), p(2)-eps);

                    return (FEShim<3,BERNSTEIN,RealType>::shape(elem, order, i, pp) -
                            FEShim<3,BERNSTEIN,RealType>::shape(elem, order, i, pm))/2./eps;
                  }
                default:
                  libmesh_error_msg("Invalid derivative index j = " << j);
                }

            }

            //       // Compute hex shape functions as a tensor-product
            //       const auto xi    = p(0);
            //       const auto eta   = p(1);
            //       const auto zeta  = p(2);
            //       auto xi_mapped   = p(0);
            //       auto eta_mapped  = p(1);
            //       auto zeta_mapped = p(2);

            //       // The only way to make any sense of this
            //       // is to look at the mgflo/mg2/mgf documentation
            //       // and make the cut-out cube!
            //       //  Nodes                         0  1  2  3  4  5  6  7  8  8  9  9 10 10 11 11 12 12 13 13 14 14 15 15 16 16 17 17 18 18 19 19 20 20 20 20 21 21 21 21 22 22 22 22 23 23 23 23 24 24 24 24 25 25 25 25 26 26 26 26 26 26 26 26
            //       //  DOFS                          0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 18 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 60 62 63
            //       static const unsigned int i0[] = {0, 1, 1, 0, 0, 1, 1, 0, 2, 3, 1, 1, 2, 3, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 2, 3, 1, 1, 2, 3, 0, 0, 2, 3, 2, 3, 2, 3, 2, 3, 1, 1, 1, 1, 2, 3, 2, 3, 0, 0, 0, 0, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3};
            //       static const unsigned int i1[] = {0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 2, 3, 1, 1, 2, 3, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 2, 3, 1, 1, 2, 3, 2, 2, 3, 3, 0, 0, 0, 0, 2, 3, 2, 3, 1, 1, 1, 1, 2, 3, 2, 3, 2, 2, 3, 3, 2, 2, 3, 3, 2, 2, 3, 3};
            //       static const unsigned int i2[] = {0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 2, 3, 2, 3, 2, 3, 2, 3, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 2, 2, 3, 3, 2, 2, 3, 3, 2, 2, 3, 3, 2, 2, 3, 3, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3};



            //       // handle the edge orientation
            //       {
            // // Edge 0
            // if ((i1[i] == 0) && (i2[i] == 0))
            //   {
            //     if (elem->node_id(0) != std::min(elem->node_id(0), elem->node_id(1)))
            //       xi_mapped = -xi;
            //   }
            // // Edge 1
            // else if ((i0[i] == 1) && (i2[i] == 0))
            //   {
            //     if (elem->node_id(1) != std::min(elem->node_id(1), elem->node_id(2)))
            //       eta_mapped = -eta;
            //   }
            // // Edge 2
            // else if ((i1[i] == 1) && (i2[i] == 0))
            //   {
            //     if (elem->node_id(3) != std::min(elem->node_id(3), elem->node_id(2)))
            //       xi_mapped = -xi;
            //   }
            // // Edge 3
            // else if ((i0[i] == 0) && (i2[i] == 0))
            //   {
            //     if (elem->node_id(0) != std::min(elem->node_id(0), elem->node_id(3)))
            //       eta_mapped = -eta;
            //   }
            // // Edge 4
            // else if ((i0[i] == 0) && (i1[i] == 0))
            //   {
            //     if (elem->node_id(0) != std::min(elem->node_id(0), elem->node_id(4)))
            //       zeta_mapped = -zeta;
            //   }
            // // Edge 5
            // else if ((i0[i] == 1) && (i1[i] == 0))
            //   {
            //     if (elem->node_id(1) != std::min(elem->node_id(1), elem->node_id(5)))
            //       zeta_mapped = -zeta;
            //   }
            // // Edge 6
            // else if ((i0[i] == 1) && (i1[i] == 1))
            //   {
            //     if (elem->node_id(2) != std::min(elem->node_id(2), elem->node_id(6)))
            //       zeta_mapped = -zeta;
            //   }
            // // Edge 7
            // else if ((i0[i] == 0) && (i1[i] == 1))
            //   {
            //     if (elem->node_id(3) != std::min(elem->node_id(3), elem->node_id(7)))
            //       zeta_mapped = -zeta;
            //   }
            // // Edge 8
            // else if ((i1[i] == 0) && (i2[i] == 1))
            //   {
            //     if (elem->node_id(4) != std::min(elem->node_id(4), elem->node_id(5)))
            //       xi_mapped = -xi;
            //   }
            // // Edge 9
            // else if ((i0[i] == 1) && (i2[i] == 1))
            //   {
            //     if (elem->node_id(5) != std::min(elem->node_id(5), elem->node_id(6)))
            //       eta_mapped = -eta;
            //   }
            // // Edge 10
            // else if ((i1[i] == 1) && (i2[i] == 1))
            //   {
            //     if (elem->node_id(7) != std::min(elem->node_id(7), elem->node_id(6)))
            //       xi_mapped = -xi;
            //   }
            // // Edge 11
            // else if ((i0[i] == 0) && (i2[i] == 1))
            //   {
            //     if (elem->node_id(4) != std::min(elem->node_id(4), elem->node_id(7)))
            //       eta_mapped = -eta;
            //   }
            //       }


            //       // handle the face orientation
            //       {
            // // Face 0
            // if ((i2[i] == 0) && (i0[i] >= 2) && (i1[i] >= 2))
            //   {
            //     const unsigned int min_node = std::min(elem->node_id(1),
            //    std::min(elem->node_id(2),
            //     std::min(elem->node_id(0),
            //      elem->node_id(3))));
            //     if (elem->node_id(0) == min_node)
            //       if (elem->node_id(1) == std::min(elem->node_id(1), elem->node_id(3)))
            // {
            //   // Case 1
            //   xi_mapped  = xi;
            //   eta_mapped = eta;
            // }
            //       else
            // {
            //   // Case 2
            //   xi_mapped  = eta;
            //   eta_mapped = xi;
            // }

            //     else if (elem->node_id(3) == min_node)
            //       if (elem->node_id(0) == std::min(elem->node_id(0), elem->node_id(2)))
            // {
            //   // Case 3
            //   xi_mapped  = -eta;
            //   eta_mapped = xi;
            // }
            //       else
            // {
            //   // Case 4
            //   xi_mapped  = xi;
            //   eta_mapped = -eta;
            // }

            //     else if (elem->node_id(2) == min_node)
            //       if (elem->node_id(3) == std::min(elem->node_id(3), elem->node_id(1)))
            // {
            //   // Case 5
            //   xi_mapped  = -xi;
            //   eta_mapped = -eta;
            // }
            //       else
            // {
            //   // Case 6
            //   xi_mapped  = -eta;
            //   eta_mapped = -xi;
            // }

            //     else if (elem->node_id(1) == min_node)
            //       if (elem->node_id(2) == std::min(elem->node_id(2), elem->node_id(0)))
            // {
            //   // Case 7
            //   xi_mapped  = eta;
            //   eta_mapped = -xi;
            // }
            //       else
            // {
            //   // Case 8
            //   xi_mapped  = -xi;
            //   eta_mapped = eta;
            // }
            //   }


            // // Face 1
            // else if ((i1[i] == 0) && (i0[i] >= 2) && (i2[i] >= 2))
            //   {
            //     const unsigned int min_node = std::min(elem->node_id(0),
            //    std::min(elem->node_id(1),
            //     std::min(elem->node_id(5),
            //      elem->node_id(4))));
            //     if (elem->node_id(0) == min_node)
            //       if (elem->node_id(1) == std::min(elem->node_id(1), elem->node_id(4)))
            // {
            //   // Case 1
            //   xi_mapped   = xi;
            //   zeta_mapped = zeta;
            // }
            //       else
            // {
            //   // Case 2
            //   xi_mapped   = zeta;
            //   zeta_mapped = xi;
            // }

            //     else if (elem->node_id(1) == min_node)
            //       if (elem->node_id(5) == std::min(elem->node_id(5), elem->node_id(0)))
            // {
            //   // Case 3
            //   xi_mapped   = zeta;
            //   zeta_mapped = -xi;
            // }
            //       else
            // {
            //   // Case 4
            //   xi_mapped   = -xi;
            //   zeta_mapped = zeta;
            // }

            //     else if (elem->node_id(5) == min_node)
            //       if (elem->node_id(4) == std::min(elem->node_id(4), elem->node_id(1)))
            // {
            //   // Case 5
            //   xi_mapped   = -xi;
            //   zeta_mapped = -zeta;
            // }
            //       else
            // {
            //   // Case 6
            //   xi_mapped   = -zeta;
            //   zeta_mapped = -xi;
            // }

            //     else if (elem->node_id(4) == min_node)
            //       if (elem->node_id(0) == std::min(elem->node_id(0), elem->node_id(5)))
            // {
            //   // Case 7
            //   xi_mapped   = -xi;
            //   zeta_mapped = zeta;
            // }
            //       else
            // {
            //   // Case 8
            //   xi_mapped   = xi;
            //   zeta_mapped = -zeta;
            // }
            //   }


            // // Face 2
            // else if ((i0[i] == 1) && (i1[i] >= 2) && (i2[i] >= 2))
            //   {
            //     const unsigned int min_node = std::min(elem->node_id(1),
            //    std::min(elem->node_id(2),
            //     std::min(elem->node_id(6),
            //      elem->node_id(5))));
            //     if (elem->node_id(1) == min_node)
            //       if (elem->node_id(2) == std::min(elem->node_id(2), elem->node_id(5)))
            // {
            //   // Case 1
            //   eta_mapped  = eta;
            //   zeta_mapped = zeta;
            // }
            //       else
            // {
            //   // Case 2
            //   eta_mapped  = zeta;
            //   zeta_mapped = eta;
            // }

            //     else if (elem->node_id(2) == min_node)
            //       if (elem->node_id(6) == std::min(elem->node_id(6), elem->node_id(1)))
            // {
            //   // Case 3
            //   eta_mapped  = zeta;
            //   zeta_mapped = -eta;
            // }
            //       else
            // {
            //   // Case 4
            //   eta_mapped  = -eta;
            //   zeta_mapped = zeta;
            // }

            //     else if (elem->node_id(6) == min_node)
            //       if (elem->node_id(5) == std::min(elem->node_id(5), elem->node_id(2)))
            // {
            //   // Case 5
            //   eta_mapped  = -eta;
            //   zeta_mapped = -zeta;
            // }
            //       else
            // {
            //   // Case 6
            //   eta_mapped  = -zeta;
            //   zeta_mapped = -eta;
            // }

            //     else if (elem->node_id(5) == min_node)
            //       if (elem->node_id(1) == std::min(elem->node_id(1), elem->node_id(6)))
            // {
            //   // Case 7
            //   eta_mapped  = -zeta;
            //   zeta_mapped = eta;
            // }
            //       else
            // {
            //   // Case 8
            //   eta_mapped   = eta;
            //   zeta_mapped = -zeta;
            // }
            //   }


            // // Face 3
            // else if ((i1[i] == 1) && (i0[i] >= 2) && (i2[i] >= 2))
            //   {
            //     const unsigned int min_node = std::min(elem->node_id(2),
            //    std::min(elem->node_id(3),
            //     std::min(elem->node_id(7),
            //      elem->node_id(6))));
            //     if (elem->node_id(3) == min_node)
            //       if (elem->node_id(2) == std::min(elem->node_id(2), elem->node_id(7)))
            // {
            //   // Case 1
            //   xi_mapped   = xi;
            //   zeta_mapped = zeta;
            // }
            //       else
            // {
            //   // Case 2
            //   xi_mapped   = zeta;
            //   zeta_mapped = xi;
            // }

            //     else if (elem->node_id(7) == min_node)
            //       if (elem->node_id(3) == std::min(elem->node_id(3), elem->node_id(6)))
            // {
            //   // Case 3
            //   xi_mapped   = -zeta;
            //   zeta_mapped = xi;
            // }
            //       else
            // {
            //   // Case 4
            //   xi_mapped   = xi;
            //   zeta_mapped = -zeta;
            // }

            //     else if (elem->node_id(6) == min_node)
            //       if (elem->node_id(7) == std::min(elem->node_id(7), elem->node_id(2)))
            // {
            //   // Case 5
            //   xi_mapped   = -xi;
            //   zeta_mapped = -zeta;
            // }
            //       else
            // {
            //   // Case 6
            //   xi_mapped   = -zeta;
            //   zeta_mapped = -xi;
            // }

            //     else if (elem->node_id(2) == min_node)
            //       if (elem->node_id(6) == std::min(elem->node_id(3), elem->node_id(6)))
            // {
            //   // Case 7
            //   xi_mapped   = zeta;
            //   zeta_mapped = -xi;
            // }
            //       else
            // {
            //   // Case 8
            //   xi_mapped   = -xi;
            //   zeta_mapped = zeta;
            // }
            //   }


            // // Face 4
            // else if ((i0[i] == 0) && (i1[i] >= 2) && (i2[i] >= 2))
            //   {
            //     const unsigned int min_node = std::min(elem->node_id(3),
            //    std::min(elem->node_id(0),
            //     std::min(elem->node_id(4),
            //      elem->node_id(7))));
            //     if (elem->node_id(0) == min_node)
            //       if (elem->node_id(3) == std::min(elem->node_id(3), elem->node_id(4)))
            // {
            //   // Case 1
            //   eta_mapped  = eta;
            //   zeta_mapped = zeta;
            // }
            //       else
            // {
            //   // Case 2
            //   eta_mapped  = zeta;
            //   zeta_mapped = eta;
            // }

            //     else if (elem->node_id(4) == min_node)
            //       if (elem->node_id(0) == std::min(elem->node_id(0), elem->node_id(7)))
            // {
            //   // Case 3
            //   eta_mapped  = -zeta;
            //   zeta_mapped = eta;
            // }
            //       else
            // {
            //   // Case 4
            //   eta_mapped  = eta;
            //   zeta_mapped = -zeta;
            // }

            //     else if (elem->node_id(7) == min_node)
            //       if (elem->node_id(4) == std::min(elem->node_id(4), elem->node_id(3)))
            // {
            //   // Case 5
            //   eta_mapped  = -eta;
            //   zeta_mapped = -zeta;
            // }
            //       else
            // {
            //   // Case 6
            //   eta_mapped  = -zeta;
            //   zeta_mapped = -eta;
            // }

            //     else if (elem->node_id(3) == min_node)
            //       if (elem->node_id(7) == std::min(elem->node_id(7), elem->node_id(0)))
            // {
            //   // Case 7
            //   eta_mapped   = zeta;
            //   zeta_mapped = -eta;
            // }
            //       else
            // {
            //   // Case 8
            //   eta_mapped  = -eta;
            //   zeta_mapped = zeta;
            // }
            //   }


            // // Face 5
            // else if ((i2[i] == 1) && (i0[i] >= 2) && (i1[i] >= 2))
            //   {
            //     const unsigned int min_node = std::min(elem->node_id(4),
            //    std::min(elem->node_id(5),
            //     std::min(elem->node_id(6),
            //      elem->node_id(7))));
            //     if (elem->node_id(4) == min_node)
            //       if (elem->node_id(5) == std::min(elem->node_id(5), elem->node_id(7)))
            // {
            //   // Case 1
            //   xi_mapped  = xi;
            //   eta_mapped = eta;
            // }
            //       else
            // {
            //   // Case 2
            //   xi_mapped  = eta;
            //   eta_mapped = xi;
            // }

            //     else if (elem->node_id(5) == min_node)
            //       if (elem->node_id(6) == std::min(elem->node_id(6), elem->node_id(4)))
            // {
            //   // Case 3
            //   xi_mapped  = eta;
            //   eta_mapped = -xi;
            // }
            //       else
            // {
            //   // Case 4
            //   xi_mapped  = -xi;
            //   eta_mapped = eta;
            // }

            //     else if (elem->node_id(6) == min_node)
            //       if (elem->node_id(7) == std::min(elem->node_id(7), elem->node_id(5)))
            // {
            //   // Case 5
            //   xi_mapped  = -xi;
            //   eta_mapped = -eta;
            // }
            //       else
            // {
            //   // Case 6
            //   xi_mapped  = -eta;
            //   eta_mapped = -xi;
            // }

            //     else if (elem->node_id(7) == min_node)
            //       if (elem->node_id(4) == std::min(elem->node_id(4), elem->node_id(6)))
            // {
            //   // Case 7
            //   xi_mapped  = -eta;
            //   eta_mapped = xi;
            // }
            //       else
            // {
            //   // Case 8
            //   xi_mapped  = xi;
            //   eta_mapped = eta;
            // }
            //   }


            //       }



            //       libmesh_assert_less (j, 3);

            //       switch (j)
            // {
            //   // d()/dxi
            // case 0:
            //   return (FEShim<1,BERNSTEIN,RealType>::shape_deriv(EDGE3, totalorder, i0[i], 0, xi_mapped)*
            //   FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i1[i],    eta_mapped)*
            //   FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i2[i],    zeta_mapped));

            //   // d()/deta
            // case 1:
            //   return (FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i0[i],    xi_mapped)*
            //   FEShim<1,BERNSTEIN,RealType>::shape_deriv(EDGE3, totalorder, i1[i], 0, eta_mapped)*
            //   FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i2[i],    zeta_mapped));

            //   // d()/dzeta
            // case 2:
            //   return (FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i0[i],    xi_mapped)*
            //   FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i1[i],    eta_mapped)*
            //   FEShim<1,BERNSTEIN,RealType>::shape_deriv(EDGE3, totalorder, i2[i], 0, zeta_mapped));

            // default:
            // libmesh_error_msg("Invalid derivative index j = " << j);
            // }


          default:
            libmesh_error_msg("Invalid element type = " << type);
          }
      }

      // 4th-order Bernstein.
    case FOURTH:
      {
        switch (type)
          {

            // Bernstein shape functions derivatives on the hexahedral.
          case HEX27:
            {
              const Real eps = 1.e-4;

              libmesh_assert_less (i, 125);

              switch (j)
                {
                  //  d()/dxi
                case 0:
                  {
                    const Point pp(p(0)+eps, p(1), p(2));
                    const Point pm(p(0)-eps, p(1), p(2));

                    return (FEShim<3,BERNSTEIN,RealType>::shape(elem, order, i, pp) -
                            FEShim<3,BERNSTEIN,RealType>::shape(elem, order, i, pm))/2./eps;
                  }

                  // d()/deta
                case 1:
                  {
                    const Point pp(p(0), p(1)+eps, p(2));
                    const Point pm(p(0), p(1)-eps, p(2));

                    return (FEShim<3,BERNSTEIN,RealType>::shape(elem, order, i, pp) -
                            FEShim<3,BERNSTEIN,RealType>::shape(elem, order, i, pm))/2./eps;
                  }
                  // d()/dzeta
                case 2:
                  {
                    const Point pp(p(0), p(1), p(2)+eps);
                    const Point pm(p(0), p(1), p(2)-eps);

                    return (FEShim<3,BERNSTEIN,RealType>::shape(elem, order, i, pp) -
                            FEShim<3,BERNSTEIN,RealType>::shape(elem, order, i, pm))/2./eps;
                  }
                default:
                  libmesh_error_msg("Invalid derivative index j = " << j);
                }
            }

            //       // Compute hex shape functions as a tensor-product
            //       const auto xi    = p(0);
            //       const auto eta   = p(1);
            //       const auto zeta  = p(2);
            //       auto xi_mapped   = p(0);
            //       auto eta_mapped  = p(1);
            //       auto zeta_mapped = p(2);

            //       // The only way to make any sense of this
            //       // is to look at the mgflo/mg2/mgf documentation
            //       // and make the cut-out cube!
            //       //  Nodes                         0  1  2  3  4  5  6  7  8  8  9  9 10 10 11 11 12 12 13 13 14 14 15 15 16 16 17 17 18 18 19 19 20 20 20 20 21 21 21 21 22 22 22 22 23 23 23 23 24 24 24 24 25 25 25 25 26 26 26 26 26 26 26 26
            //       //  DOFS                          0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 18 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 60 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24
            //       static const unsigned int i0[] = {0, 1, 1, 0, 0, 1, 1, 0, 2, 3, 4, 1, 1, 1, 2, 3, 4, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 2, 3, 4, 1, 1, 1, 2, 3, 4, 0, 0, 0, 2, 3, 4, 2, 3, 4, 2, 3, 4, 2, 3, 4, 2, 3, 4, 2, 3, 4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 4, 2, 3, 4, 2, 3, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 3, 4, 2, 3, 4, 2, 3, 4, 2, 3, 4, 2, 3, 4, 2, 3, 4, 2, 3, 4, 2, 3, 4, 2, 3, 4, 2, 3, 4, 2, 3, 4, 2, 3, 4};
            //       static const unsigned int i1[] = {0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 2, 3, 4, 1, 1, 1, 2, 3, 4, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 2, 3, 4, 1, 1, 1, 2, 3, 4, 2, 2, 2, 3, 3, 3, 4, 4, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 3, 4, 2, 3, 4, 2, 3, 4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 4, 2, 3, 4, 2, 3, 4, 2, 2, 2, 3, 3, 3, 4, 4, 4, 2, 2, 2, 3, 3, 3, 4, 4, 4, 2, 2, 2, 3, 3, 3, 4, 4, 4, 2, 2, 2, 3, 3, 3, 4, 4, 4};
            //       static const unsigned int i2[] = {0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 3, 4, 2, 3, 4, 2, 3, 4, 2, 3, 4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 3, 3, 3, 4, 4, 4, 2, 2, 2, 3, 3, 3, 4, 4, 4, 2, 2, 2, 3, 3, 3, 4, 4, 4, 2, 2, 2, 3, 3, 3, 4, 4, 4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4};



            //       // handle the edge orientation
            //       {
            // // Edge 0
            // if ((i1[i] == 0) && (i2[i] == 0))
            //   {
            //     if (elem->node_id(0) != std::min(elem->node_id(0), elem->node_id(1)))
            //       xi_mapped = -xi;
            //   }
            // // Edge 1
            // else if ((i0[i] == 1) && (i2[i] == 0))
            //   {
            //     if (elem->node_id(1) != std::min(elem->node_id(1), elem->node_id(2)))
            //       eta_mapped = -eta;
            //   }
            // // Edge 2
            // else if ((i1[i] == 1) && (i2[i] == 0))
            //   {
            //     if (elem->node_id(3) != std::min(elem->node_id(3), elem->node_id(2)))
            //       xi_mapped = -xi;
            //   }
            // // Edge 3
            // else if ((i0[i] == 0) && (i2[i] == 0))
            //   {
            //     if (elem->node_id(0) != std::min(elem->node_id(0), elem->node_id(3)))
            //       eta_mapped = -eta;
            //   }
            // // Edge 4
            // else if ((i0[i] == 0) && (i1[i] == 0))
            //   {
            //     if (elem->node_id(0) != std::min(elem->node_id(0), elem->node_id(4)))
            //       zeta_mapped = -zeta;
            //   }
            // // Edge 5
            // else if ((i0[i] == 1) && (i1[i] == 0))
            //   {
            //     if (elem->node_id(1) != std::min(elem->node_id(1), elem->node_id(5)))
            //       zeta_mapped = -zeta;
            //   }
            // // Edge 6
            // else if ((i0[i] == 1) && (i1[i] == 1))
            //   {
            //     if (elem->node_id(2) != std::min(elem->node_id(2), elem->node_id(6)))
            //       zeta_mapped = -zeta;
            //   }
            // // Edge 7
            // else if ((i0[i] == 0) && (i1[i] == 1))
            //   {
            //     if (elem->node_id(3) != std::min(elem->node_id(3), elem->node_id(7)))
            //       zeta_mapped = -zeta;
            //   }
            // // Edge 8
            // else if ((i1[i] == 0) && (i2[i] == 1))
            //   {
            //     if (elem->node_id(4) != std::min(elem->node_id(4), elem->node_id(5)))
            //       xi_mapped = -xi;
            //   }
            // // Edge 9
            // else if ((i0[i] == 1) && (i2[i] == 1))
            //   {
            //     if (elem->node_id(5) != std::min(elem->node_id(5), elem->node_id(6)))
            //       eta_mapped = -eta;
            //   }
            // // Edge 10
            // else if ((i1[i] == 1) && (i2[i] == 1))
            //   {
            //     if (elem->node_id(7) != std::min(elem->node_id(7), elem->node_id(6)))
            //       xi_mapped = -xi;
            //   }
            // // Edge 11
            // else if ((i0[i] == 0) && (i2[i] == 1))
            //   {
            //     if (elem->node_id(4) != std::min(elem->node_id(4), elem->node_id(7)))
            //       eta_mapped = -eta;
            //   }
            //       }


            //       // handle the face orientation
            //       {
            // // Face 0
            // if ((i2[i] == 0) && (i0[i] >= 2) && (i1[i] >= 2))
            //   {
            //     const unsigned int min_node = std::min(elem->node_id(1),
            //    std::min(elem->node_id(2),
            //     std::min(elem->node_id(0),
            //      elem->node_id(3))));
            //     if (elem->node_id(0) == min_node)
            //       if (elem->node_id(1) == std::min(elem->node_id(1), elem->node_id(3)))
            // {
            //   // Case 1
            //   xi_mapped  = xi;
            //   eta_mapped = eta;
            // }
            //       else
            // {
            //   // Case 2
            //   xi_mapped  = eta;
            //   eta_mapped = xi;
            // }

            //     else if (elem->node_id(3) == min_node)
            //       if (elem->node_id(0) == std::min(elem->node_id(0), elem->node_id(2)))
            // {
            //   // Case 3
            //   xi_mapped  = -eta;
            //   eta_mapped = xi;
            // }
            //       else
            // {
            //   // Case 4
            //   xi_mapped  = xi;
            //   eta_mapped = -eta;
            // }

            //     else if (elem->node_id(2) == min_node)
            //       if (elem->node_id(3) == std::min(elem->node_id(3), elem->node_id(1)))
            // {
            //   // Case 5
            //   xi_mapped  = -xi;
            //   eta_mapped = -eta;
            // }
            //       else
            // {
            //   // Case 6
            //   xi_mapped  = -eta;
            //   eta_mapped = -xi;
            // }

            //     else if (elem->node_id(1) == min_node)
            //       if (elem->node_id(2) == std::min(elem->node_id(2), elem->node_id(0)))
            // {
            //   // Case 7
            //   xi_mapped  = eta;
            //   eta_mapped = -xi;
            // }
            //       else
            // {
            //   // Case 8
            //   xi_mapped  = -xi;
            //   eta_mapped = eta;
            // }
            //   }


            // // Face 1
            // else if ((i1[i] == 0) && (i0[i] >= 2) && (i2[i] >= 2))
            //   {
            //     const unsigned int min_node = std::min(elem->node_id(0),
            //    std::min(elem->node_id(1),
            //     std::min(elem->node_id(5),
            //      elem->node_id(4))));
            //     if (elem->node_id(0) == min_node)
            //       if (elem->node_id(1) == std::min(elem->node_id(1), elem->node_id(4)))
            // {
            //   // Case 1
            //   xi_mapped   = xi;
            //   zeta_mapped = zeta;
            // }
            //       else
            // {
            //   // Case 2
            //   xi_mapped   = zeta;
            //   zeta_mapped = xi;
            // }

            //     else if (elem->node_id(1) == min_node)
            //       if (elem->node_id(5) == std::min(elem->node_id(5), elem->node_id(0)))
            // {
            //   // Case 3
            //   xi_mapped   = zeta;
            //   zeta_mapped = -xi;
            // }
            //       else
            // {
            //   // Case 4
            //   xi_mapped   = -xi;
            //   zeta_mapped = zeta;
            // }

            //     else if (elem->node_id(5) == min_node)
            //       if (elem->node_id(4) == std::min(elem->node_id(4), elem->node_id(1)))
            // {
            //   // Case 5
            //   xi_mapped   = -xi;
            //   zeta_mapped = -zeta;
            // }
            //       else
            // {
            //   // Case 6
            //   xi_mapped   = -zeta;
            //   zeta_mapped = -xi;
            // }

            //     else if (elem->node_id(4) == min_node)
            //       if (elem->node_id(0) == std::min(elem->node_id(0), elem->node_id(5)))
            // {
            //   // Case 7
            //   xi_mapped   = -xi;
            //   zeta_mapped = zeta;
            // }
            //       else
            // {
            //   // Case 8
            //   xi_mapped   = xi;
            //   zeta_mapped = -zeta;
            // }
            //   }


            // // Face 2
            // else if ((i0[i] == 1) && (i1[i] >= 2) && (i2[i] >= 2))
            //   {
            //     const unsigned int min_node = std::min(elem->node_id(1),
            //    std::min(elem->node_id(2),
            //     std::min(elem->node_id(6),
            //      elem->node_id(5))));
            //     if (elem->node_id(1) == min_node)
            //       if (elem->node_id(2) == std::min(elem->node_id(2), elem->node_id(5)))
            // {
            //   // Case 1
            //   eta_mapped  = eta;
            //   zeta_mapped = zeta;
            // }
            //       else
            // {
            //   // Case 2
            //   eta_mapped  = zeta;
            //   zeta_mapped = eta;
            // }

            //     else if (elem->node_id(2) == min_node)
            //       if (elem->node_id(6) == std::min(elem->node_id(6), elem->node_id(1)))
            // {
            //   // Case 3
            //   eta_mapped  = zeta;
            //   zeta_mapped = -eta;
            // }
            //       else
            // {
            //   // Case 4
            //   eta_mapped  = -eta;
            //   zeta_mapped = zeta;
            // }

            //     else if (elem->node_id(6) == min_node)
            //       if (elem->node_id(5) == std::min(elem->node_id(5), elem->node_id(2)))
            // {
            //   // Case 5
            //   eta_mapped  = -eta;
            //   zeta_mapped = -zeta;
            // }
            //       else
            // {
            //   // Case 6
            //   eta_mapped  = -zeta;
            //   zeta_mapped = -eta;
            // }

            //     else if (elem->node_id(5) == min_node)
            //       if (elem->node_id(1) == std::min(elem->node_id(1), elem->node_id(6)))
            // {
            //   // Case 7
            //   eta_mapped  = -zeta;
            //   zeta_mapped = eta;
            // }
            //       else
            // {
            //   // Case 8
            //   eta_mapped   = eta;
            //   zeta_mapped = -zeta;
            // }
            //   }


            // // Face 3
            // else if ((i1[i] == 1) && (i0[i] >= 2) && (i2[i] >= 2))
            //   {
            //     const unsigned int min_node = std::min(elem->node_id(2),
            //    std::min(elem->node_id(3),
            //     std::min(elem->node_id(7),
            //      elem->node_id(6))));
            //     if (elem->node_id(3) == min_node)
            //       if (elem->node_id(2) == std::min(elem->node_id(2), elem->node_id(7)))
            // {
            //   // Case 1
            //   xi_mapped   = xi;
            //   zeta_mapped = zeta;
            // }
            //       else
            // {
            //   // Case 2
            //   xi_mapped   = zeta;
            //   zeta_mapped = xi;
            // }

            //     else if (elem->node_id(7) == min_node)
            //       if (elem->node_id(3) == std::min(elem->node_id(3), elem->node_id(6)))
            // {
            //   // Case 3
            //   xi_mapped   = -zeta;
            //   zeta_mapped = xi;
            // }
            //       else
            // {
            //   // Case 4
            //   xi_mapped   = xi;
            //   zeta_mapped = -zeta;
            // }

            //     else if (elem->node_id(6) == min_node)
            //       if (elem->node_id(7) == std::min(elem->node_id(7), elem->node_id(2)))
            // {
            //   // Case 5
            //   xi_mapped   = -xi;
            //   zeta_mapped = -zeta;
            // }
            //       else
            // {
            //   // Case 6
            //   xi_mapped   = -zeta;
            //   zeta_mapped = -xi;
            // }

            //     else if (elem->node_id(2) == min_node)
            //       if (elem->node_id(6) == std::min(elem->node_id(3), elem->node_id(6)))
            // {
            //   // Case 7
            //   xi_mapped   = zeta;
            //   zeta_mapped = -xi;
            // }
            //       else
            // {
            //   // Case 8
            //   xi_mapped   = -xi;
            //   zeta_mapped = zeta;
            // }
            //   }


            // // Face 4
            // else if ((i0[i] == 0) && (i1[i] >= 2) && (i2[i] >= 2))
            //   {
            //     const unsigned int min_node = std::min(elem->node_id(3),
            //    std::min(elem->node_id(0),
            //     std::min(elem->node_id(4),
            //      elem->node_id(7))));
            //     if (elem->node_id(0) == min_node)
            //       if (elem->node_id(3) == std::min(elem->node_id(3), elem->node_id(4)))
            // {
            //   // Case 1
            //   eta_mapped  = eta;
            //   zeta_mapped = zeta;
            // }
            //       else
            // {
            //   // Case 2
            //   eta_mapped  = zeta;
            //   zeta_mapped = eta;
            // }

            //     else if (elem->node_id(4) == min_node)
            //       if (elem->node_id(0) == std::min(elem->node_id(0), elem->node_id(7)))
            // {
            //   // Case 3
            //   eta_mapped  = -zeta;
            //   zeta_mapped = eta;
            // }
            //       else
            // {
            //   // Case 4
            //   eta_mapped  = eta;
            //   zeta_mapped = -zeta;
            // }

            //     else if (elem->node_id(7) == min_node)
            //       if (elem->node_id(4) == std::min(elem->node_id(4), elem->node_id(3)))
            // {
            //   // Case 5
            //   eta_mapped  = -eta;
            //   zeta_mapped = -zeta;
            // }
            //       else
            // {
            //   // Case 6
            //   eta_mapped  = -zeta;
            //   zeta_mapped = -eta;
            // }

            //     else if (elem->node_id(3) == min_node)
            //       if (elem->node_id(7) == std::min(elem->node_id(7), elem->node_id(0)))
            // {
            //   // Case 7
            //   eta_mapped   = zeta;
            //   zeta_mapped = -eta;
            // }
            //       else
            // {
            //   // Case 8
            //   eta_mapped  = -eta;
            //   zeta_mapped = zeta;
            // }
            //   }


            // // Face 5
            // else if ((i2[i] == 1) && (i0[i] >= 2) && (i1[i] >= 2))
            //   {
            //     const unsigned int min_node = std::min(elem->node_id(4),
            //    std::min(elem->node_id(5),
            //     std::min(elem->node_id(6),
            //      elem->node_id(7))));
            //     if (elem->node_id(4) == min_node)
            //       if (elem->node_id(5) == std::min(elem->node_id(5), elem->node_id(7)))
            // {
            //   // Case 1
            //   xi_mapped  = xi;
            //   eta_mapped = eta;
            // }
            //       else
            // {
            //   // Case 2
            //   xi_mapped  = eta;
            //   eta_mapped = xi;
            // }

            //     else if (elem->node_id(5) == min_node)
            //       if (elem->node_id(6) == std::min(elem->node_id(6), elem->node_id(4)))
            // {
            //   // Case 3
            //   xi_mapped  = eta;
            //   eta_mapped = -xi;
            // }
            //       else
            // {
            //   // Case 4
            //   xi_mapped  = -xi;
            //   eta_mapped = eta;
            // }

            //     else if (elem->node_id(6) == min_node)
            //       if (elem->node_id(7) == std::min(elem->node_id(7), elem->node_id(5)))
            // {
            //   // Case 5
            //   xi_mapped  = -xi;
            //   eta_mapped = -eta;
            // }
            //       else
            // {
            //   // Case 6
            //   xi_mapped  = -eta;
            //   eta_mapped = -xi;
            // }

            //     else if (elem->node_id(7) == min_node)
            //       if (elem->node_id(4) == std::min(elem->node_id(4), elem->node_id(6)))
            // {
            //   // Case 7
            //   xi_mapped  = -eta;
            //   eta_mapped = xi;
            // }
            //       else
            // {
            //   // Case 8
            //   xi_mapped  = xi;
            //   eta_mapped = eta;
            // }
            //   }


            //       }



            //       libmesh_assert_less (j, 3);

            //       switch (j)
            // {
            //   // d()/dxi
            // case 0:
            //   return (FEShim<1,BERNSTEIN,RealType>::shape_deriv(EDGE3, totalorder, i0[i], 0, xi_mapped)*
            //   FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i1[i],    eta_mapped)*
            //   FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i2[i],    zeta_mapped));

            //   // d()/deta
            // case 1:
            //   return (FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i0[i],    xi_mapped)*
            //   FEShim<1,BERNSTEIN,RealType>::shape_deriv(EDGE3, totalorder, i1[i], 0, eta_mapped)*
            //   FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i2[i],    zeta_mapped));

            //   // d()/dzeta
            // case 2:
            //   return (FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i0[i],    xi_mapped)*
            //   FEShim<1,BERNSTEIN,RealType>::shape      (EDGE3, totalorder, i1[i],    eta_mapped)*
            //   FEShim<1,BERNSTEIN,RealType>::shape_deriv(EDGE3, totalorder, i2[i], 0, zeta_mapped));

            // default:
            //   libmesh_error_msg("Invalid derivative index j = " << j);
            // }


          default:
            libmesh_error_msg("Invalid element type = " << type);
          }
      }


    default:
      libmesh_error_msg("Invalid totalorder = " << totalorder);
    }

#else // LIBMESH_DIM != 3
  libmesh_ignore(elem, order, i, j, p, add_p_level);
  libmesh_not_implemented();
#endif
}


#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

template <typename RealType>
typename FEShim<3,BERNSTEIN,RealType>::OutputShape FEShim<3,BERNSTEIN,RealType>::shape_second_deriv(const ElemType,
                                         const Order,
                                         const unsigned int,
                                         const unsigned int,
                                         const Point &)
{
  libmesh_error_msg("Bernstein polynomials require the element type \nbecause edge and face orientation is needed.");
  return 0.;
}



template <typename RealType>
typename FEShim<3,BERNSTEIN,RealType>::OutputShape FEShim<3,BERNSTEIN,RealType>::shape_second_deriv(const Elem * elem,
                                         const Order order,
                                         const unsigned int i,
                                         const unsigned int j,
                                         const Point & p,
                                         const bool add_p_level)
{

#if LIBMESH_DIM == 3
  libmesh_assert(elem);

  const Order totalorder =
    static_cast<Order>(order + add_p_level * elem->p_level());

  libmesh_assert_less (j, 6);

  {
    // I have been lazy here and am using finite differences
    // to compute the derivatives!
    const Real eps = 1.e-4;


    switch (j)
      {
        //  d^2()/dxi^2
      case 0:
        {
          const Point pp(p(0)+eps, p(1), p(2));
          const Point pm(p(0)-eps, p(1), p(2));

          return (FEShim<3,BERNSTEIN,RealType>::shape_deriv(elem, totalorder, i, 0, pp) -
                  FEShim<3,BERNSTEIN,RealType>::shape_deriv(elem, totalorder, i, 0, pm))/2./eps;
        }

        // d^2()/dxideta
      case 1:
        {
          const Point pp(p(0), p(1)+eps, p(2));
          const Point pm(p(0), p(1)-eps, p(2));

          return (FEShim<3,BERNSTEIN,RealType>::shape_deriv(elem, totalorder, i, 0, pp) -
                  FEShim<3,BERNSTEIN,RealType>::shape_deriv(elem, totalorder, i, 0, pm))/2./eps;
        }

        // d^2()/deta^2
      case 2:
        {
          const Point pp(p(0), p(1)+eps, p(2));
          const Point pm(p(0), p(1)-eps, p(2));

          return (FEShim<3,BERNSTEIN,RealType>::shape_deriv(elem, totalorder, i, 1, pp) -
                  FEShim<3,BERNSTEIN,RealType>::shape_deriv(elem, totalorder, i, 1, pm))/2./eps;
        }

        // d^2()/dxidzeta
      case 3:
        {
          const Point pp(p(0), p(1), p(2)+eps);
          const Point pm(p(0), p(1), p(2)-eps);

          return (FEShim<3,BERNSTEIN,RealType>::shape_deriv(elem, totalorder, i, 0, pp) -
                  FEShim<3,BERNSTEIN,RealType>::shape_deriv(elem, totalorder, i, 0, pm))/2./eps;
        }                  // d^2()/deta^2

        // d^2()/detadzeta
      case 4:
        {
          const Point pp(p(0), p(1), p(2)+eps);
          const Point pm(p(0), p(1), p(2)-eps);

          return (FEShim<3,BERNSTEIN,RealType>::shape_deriv(elem, totalorder, i, 1, pp) -
                  FEShim<3,BERNSTEIN,RealType>::shape_deriv(elem, totalorder, i, 1, pm))/2./eps;
        }

        // d^2()/dzeta^2
      case 5:
        {
          const Point pp(p(0), p(1), p(2)+eps);
          const Point pm(p(0), p(1), p(2)-eps);

          return (FEShim<3,BERNSTEIN,RealType>::shape_deriv(elem, totalorder, i, 2, pp) -
                  FEShim<3,BERNSTEIN,RealType>::shape_deriv(elem, totalorder, i, 2, pm))/2./eps;
        }

      default:
        libmesh_error_msg("Invalid derivative index j = " << j);
      }
  }

#else // LIBMESH_DIM != 3
  libmesh_ignore(elem, order, i, j, p, add_p_level);
  libmesh_not_implemented();
#endif
}

#endif

} // namespace libMesh



#endif //LIBMESH_ENABLE_HIGHER_ORDER_SHAPES

#endif // LIBMESH_FE_BERNSTEIN_SHAPE_3D_IMPL_H
