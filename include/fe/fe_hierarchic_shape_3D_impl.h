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

#ifndef LIBMESH_FE_HIERARCHIC_SHAPE_3D_IMPL_H
#define LIBMESH_FE_HIERARCHIC_SHAPE_3D_IMPL_H

// Local includes
#include "libmesh/fe.h"
#include "libmesh/elem.h"
#include "libmesh/number_lookups.h"

// Anonymous namespace for functions shared by HIERARCHIC and
// L2_HIERARCHIC implementations. Implementations appear at the bottom
// of this file.
namespace
{
using namespace libMesh;

template <typename RealType>
RealType fe_hierarchic_3D_shape(const ElemTempl<RealType> * elem,
                                const Order order,
                                const unsigned int i,
                                const PointTempl<RealType> & p,
                                const bool add_p_level);

template <typename RealType>
RealType fe_hierarchic_3D_shape_deriv(const ElemTempl<RealType> * elem,
                                      const Order order,
                                      const unsigned int i,
                                      const unsigned int j,
                                      const PointTempl<RealType> & p,
                                      const bool add_p_level);

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

template <typename RealType>
RealType fe_hierarchic_3D_shape_second_deriv(const ElemTempl<RealType> * elem,
                                             const Order order,
                                             const unsigned int i,
                                             const unsigned int j,
                                             const PointTempl<RealType> & p,
                                             const bool add_p_level);

#endif // LIBMESH_ENABLE_SECOND_DERIVATIVES

#if LIBMESH_DIM > 2
template <typename RealType>
PointTempl<RealType> get_min_point(const ElemTempl<RealType> * elem,
                                   unsigned int a,
                                   unsigned int b,
                                   unsigned int c,
                                   unsigned int d)
{
  return std::min(std::min(elem->point(a),elem->point(b)),
                  std::min(elem->point(c),elem->point(d)));
}

template <typename RealType>
void cube_indices(const ElemTempl<RealType> * elem,
                  const unsigned int totalorder,
                  const unsigned int i,
                  RealType & xi, RealType & eta, RealType & zeta,
                  unsigned int & i0,
                  unsigned int & i1,
                  unsigned int & i2)
{
  // The only way to make any sense of this
  // is to look at the mgflo/mg2/mgf documentation
  // and make the cut-out cube!
  // Example i0 and i1 values for totalorder = 3:
  // FIXME - these examples are incorrect now that we've got truly
  // hierarchic basis functions
  //     Nodes                         0  1  2  3  4  5  6  7  8  8  9  9 10 10 11 11 12 12 13 13 14 14 15 15 16 16 17 17 18 18 19 19 20 20 20 20 21 21 21 21 22 22 22 22 23 23 23 23 24 24 24 24 25 25 25 25 26 26 26 26 26 26 26 26
  //     DOFS                          0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 18 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 60 62 63
  // static const unsigned int i0[] = {0, 1, 1, 0, 0, 1, 1, 0, 2, 3, 1, 1, 2, 3, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 2, 3, 1, 1, 2, 3, 0, 0, 2, 3, 2, 3, 2, 3, 2, 3, 1, 1, 1, 1, 2, 3, 2, 3, 0, 0, 0, 0, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3};
  // static const unsigned int i1[] = {0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 2, 3, 1, 1, 2, 3, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 2, 3, 1, 1, 2, 3, 2, 2, 3, 3, 0, 0, 0, 0, 2, 3, 2, 3, 1, 1, 1, 1, 2, 3, 2, 3, 2, 2, 3, 3, 2, 2, 3, 3, 2, 2, 3, 3};
  // static const unsigned int i2[] = {0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 2, 3, 2, 3, 2, 3, 2, 3, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 2, 2, 3, 3, 2, 2, 3, 3, 2, 2, 3, 3, 2, 2, 3, 3, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3};

  // the number of DoFs per edge appears everywhere:
  const unsigned int e = totalorder - 1u;

  libmesh_assert_less (i, (totalorder+1u)*(totalorder+1u)*(totalorder+1u));

  RealType xi_saved = xi, eta_saved = eta, zeta_saved = zeta;

  // Vertices:
  if (i == 0)
    {
      i0 = 0;
      i1 = 0;
      i2 = 0;
    }
  else if (i == 1)
    {
      i0 = 1;
      i1 = 0;
      i2 = 0;
    }
  else if (i == 2)
    {
      i0 = 1;
      i1 = 1;
      i2 = 0;
    }
  else if (i == 3)
    {
      i0 = 0;
      i1 = 1;
      i2 = 0;
    }
  else if (i == 4)
    {
      i0 = 0;
      i1 = 0;
      i2 = 1;
    }
  else if (i == 5)
    {
      i0 = 1;
      i1 = 0;
      i2 = 1;
    }
  else if (i == 6)
    {
      i0 = 1;
      i1 = 1;
      i2 = 1;
    }
  else if (i == 7)
    {
      i0 = 0;
      i1 = 1;
      i2 = 1;
    }
  // Edge 0
  else if (i < 8 + e)
    {
      i0 = i - 6;
      i1 = 0;
      i2 = 0;
      if (elem->point(0) > elem->point(1))
        xi = -xi_saved;
    }
  // Edge 1
  else if (i < 8 + 2*e)
    {
      i0 = 1;
      i1 = i - e - 6;
      i2 = 0;
      if (elem->point(1) > elem->point(2))
        eta = -eta_saved;
    }
  // Edge 2
  else if (i < 8 + 3*e)
    {
      i0 = i - 2*e - 6;
      i1 = 1;
      i2 = 0;
      if (elem->point(3) > elem->point(2))
        xi = -xi_saved;
    }
  // Edge 3
  else if (i < 8 + 4*e)
    {
      i0 = 0;
      i1 = i - 3*e - 6;
      i2 = 0;
      if (elem->point(0) > elem->point(3))
        eta = -eta_saved;
    }
  // Edge 4
  else if (i < 8 + 5*e)
    {
      i0 = 0;
      i1 = 0;
      i2 = i - 4*e - 6;
      if (elem->point(0) > elem->point(4))
        zeta = -zeta_saved;
    }
  // Edge 5
  else if (i < 8 + 6*e)
    {
      i0 = 1;
      i1 = 0;
      i2 = i - 5*e - 6;
      if (elem->point(1) > elem->point(5))
        zeta = -zeta_saved;
    }
  // Edge 6
  else if (i < 8 + 7*e)
    {
      i0 = 1;
      i1 = 1;
      i2 = i - 6*e - 6;
      if (elem->point(2) > elem->point(6))
        zeta = -zeta_saved;
    }
  // Edge 7
  else if (i < 8 + 8*e)
    {
      i0 = 0;
      i1 = 1;
      i2 = i - 7*e - 6;
      if (elem->point(3) > elem->point(7))
        zeta = -zeta_saved;
    }
  // Edge 8
  else if (i < 8 + 9*e)
    {
      i0 = i - 8*e - 6;
      i1 = 0;
      i2 = 1;
      if (elem->point(4) > elem->point(5))
        xi = -xi_saved;
    }
  // Edge 9
  else if (i < 8 + 10*e)
    {
      i0 = 1;
      i1 = i - 9*e - 6;
      i2 = 1;
      if (elem->point(5) > elem->point(6))
        eta = -eta_saved;
    }
  // Edge 10
  else if (i < 8 + 11*e)
    {
      i0 = i - 10*e - 6;
      i1 = 1;
      i2 = 1;
      if (elem->point(7) > elem->point(6))
        xi = -xi_saved;
    }
  // Edge 11
  else if (i < 8 + 12*e)
    {
      i0 = 0;
      i1 = i - 11*e - 6;
      i2 = 1;
      if (elem->point(4) > elem->point(7))
        eta = -eta_saved;
    }
  // Face 0
  else if (i < 8 + 12*e + e*e)
    {
      unsigned int basisnum = i - 8 - 12*e;
      i0 = square_number_row[basisnum] + 2;
      i1 = square_number_column[basisnum] + 2;
      i2 = 0;
      const Point min_point = get_min_point(elem, 1, 2, 0, 3);

      if (elem->point(0) == min_point)
        if (elem->point(1) == std::min(elem->point(1), elem->point(3)))
          {
            // Case 1
            xi  = xi_saved;
            eta = eta_saved;
          }
        else
          {
            // Case 2
            xi  = eta_saved;
            eta = xi_saved;
          }

      else if (elem->point(3) == min_point)
        if (elem->point(0) == std::min(elem->point(0), elem->point(2)))
          {
            // Case 3
            xi  = -eta_saved;
            eta = xi_saved;
          }
        else
          {
            // Case 4
            xi  = xi_saved;
            eta = -eta_saved;
          }

      else if (elem->point(2) == min_point)
        if (elem->point(3) == std::min(elem->point(3), elem->point(1)))
          {
            // Case 5
            xi  = -xi_saved;
            eta = -eta_saved;
          }
        else
          {
            // Case 6
            xi  = -eta_saved;
            eta = -xi_saved;
          }

      else if (elem->point(1) == min_point)
        {
          if (elem->point(2) == std::min(elem->point(2), elem->point(0)))
            {
              // Case 7
              xi  = eta_saved;
              eta = -xi_saved;
            }
          else
            {
              // Case 8
              xi  = -xi_saved;
              eta = eta_saved;
            }
        }
    }
  // Face 1
  else if (i < 8 + 12*e + 2*e*e)
    {
      unsigned int basisnum = i - 8 - 12*e - e*e;
      i0 = square_number_row[basisnum] + 2;
      i1 = 0;
      i2 = square_number_column[basisnum] + 2;
      const Point min_point = get_min_point(elem, 0, 1, 5, 4);

      if (elem->point(0) == min_point)
        if (elem->point(1) == std::min(elem->point(1), elem->point(4)))
          {
            // Case 1
            xi   = xi_saved;
            zeta = zeta_saved;
          }
        else
          {
            // Case 2
            xi   = zeta_saved;
            zeta = xi_saved;
          }

      else if (elem->point(1) == min_point)
        if (elem->point(5) == std::min(elem->point(5), elem->point(0)))
          {
            // Case 3
            xi   = zeta_saved;
            zeta = -xi_saved;
          }
        else
          {
            // Case 4
            xi   = -xi_saved;
            zeta = zeta_saved;
          }

      else if (elem->point(5) == min_point)
        if (elem->point(4) == std::min(elem->point(4), elem->point(1)))
          {
            // Case 5
            xi   = -xi_saved;
            zeta = -zeta_saved;
          }
        else
          {
            // Case 6
            xi   = -zeta_saved;
            zeta = -xi_saved;
          }

      else if (elem->point(4) == min_point)
        {
          if (elem->point(0) == std::min(elem->point(0), elem->point(5)))
            {
              // Case 7
              xi   = -xi_saved;
              zeta = zeta_saved;
            }
          else
            {
              // Case 8
              xi   = xi_saved;
              zeta = -zeta_saved;
            }
        }
    }
  // Face 2
  else if (i < 8 + 12*e + 3*e*e)
    {
      unsigned int basisnum = i - 8 - 12*e - 2*e*e;
      i0 = 1;
      i1 = square_number_row[basisnum] + 2;
      i2 = square_number_column[basisnum] + 2;
      const Point min_point = get_min_point(elem, 1, 2, 6, 5);

      if (elem->point(1) == min_point)
        if (elem->point(2) == std::min(elem->point(2), elem->point(5)))
          {
            // Case 1
            eta  = eta_saved;
            zeta = zeta_saved;
          }
        else
          {
            // Case 2
            eta  = zeta_saved;
            zeta = eta_saved;
          }

      else if (elem->point(2) == min_point)
        if (elem->point(6) == std::min(elem->point(6), elem->point(1)))
          {
            // Case 3
            eta  = zeta_saved;
            zeta = -eta_saved;
          }
        else
          {
            // Case 4
            eta  = -eta_saved;
            zeta = zeta_saved;
          }

      else if (elem->point(6) == min_point)
        if (elem->point(5) == std::min(elem->point(5), elem->point(2)))
          {
            // Case 5
            eta  = -eta_saved;
            zeta = -zeta_saved;
          }
        else
          {
            // Case 6
            eta  = -zeta_saved;
            zeta = -eta_saved;
          }

      else if (elem->point(5) == min_point)
        {
          if (elem->point(1) == std::min(elem->point(1), elem->point(6)))
            {
              // Case 7
              eta  = -zeta_saved;
              zeta = eta_saved;
            }
          else
            {
              // Case 8
              eta   = eta_saved;
              zeta = -zeta_saved;
            }
        }
    }
  // Face 3
  else if (i < 8 + 12*e + 4*e*e)
    {
      unsigned int basisnum = i - 8 - 12*e - 3*e*e;
      i0 = square_number_row[basisnum] + 2;
      i1 = 1;
      i2 = square_number_column[basisnum] + 2;
      const Point min_point = get_min_point(elem, 2, 3, 7, 6);

      if (elem->point(3) == min_point)
        if (elem->point(2) == std::min(elem->point(2), elem->point(7)))
          {
            // Case 1
            xi   = xi_saved;
            zeta = zeta_saved;
          }
        else
          {
            // Case 2
            xi   = zeta_saved;
            zeta = xi_saved;
          }

      else if (elem->point(7) == min_point)
        if (elem->point(3) == std::min(elem->point(3), elem->point(6)))
          {
            // Case 3
            xi   = -zeta_saved;
            zeta = xi_saved;
          }
        else
          {
            // Case 4
            xi   = xi_saved;
            zeta = -zeta_saved;
          }

      else if (elem->point(6) == min_point)
        if (elem->point(7) == std::min(elem->point(7), elem->point(2)))
          {
            // Case 5
            xi   = -xi_saved;
            zeta = -zeta_saved;
          }
        else
          {
            // Case 6
            xi   = -zeta_saved;
            zeta = -xi_saved;
          }

      else if (elem->point(2) == min_point)
        {
          if (elem->point(6) == std::min(elem->point(3), elem->point(6)))
            {
              // Case 7
              xi   = zeta_saved;
              zeta = -xi_saved;
            }
          else
            {
              // Case 8
              xi   = -xi_saved;
              zeta = zeta_saved;
            }
        }
    }
  // Face 4
  else if (i < 8 + 12*e + 5*e*e)
    {
      unsigned int basisnum = i - 8 - 12*e - 4*e*e;
      i0 = 0;
      i1 = square_number_row[basisnum] + 2;
      i2 = square_number_column[basisnum] + 2;
      const Point min_point = get_min_point(elem, 3, 0, 4, 7);

      if (elem->point(0) == min_point)
        if (elem->point(3) == std::min(elem->point(3), elem->point(4)))
          {
            // Case 1
            eta  = eta_saved;
            zeta = zeta_saved;
          }
        else
          {
            // Case 2
            eta  = zeta_saved;
            zeta = eta_saved;
          }

      else if (elem->point(4) == min_point)
        if (elem->point(0) == std::min(elem->point(0), elem->point(7)))
          {
            // Case 3
            eta  = -zeta_saved;
            zeta = eta_saved;
          }
        else
          {
            // Case 4
            eta  = eta_saved;
            zeta = -zeta_saved;
          }

      else if (elem->point(7) == min_point)
        if (elem->point(4) == std::min(elem->point(4), elem->point(3)))
          {
            // Case 5
            eta  = -eta_saved;
            zeta = -zeta_saved;
          }
        else
          {
            // Case 6
            eta  = -zeta_saved;
            zeta = -eta_saved;
          }

      else if (elem->point(3) == min_point)
        {
          if (elem->point(7) == std::min(elem->point(7), elem->point(0)))
            {
              // Case 7
              eta   = zeta_saved;
              zeta = -eta_saved;
            }
          else
            {
              // Case 8
              eta  = -eta_saved;
              zeta = zeta_saved;
            }
        }
    }
  // Face 5
  else if (i < 8 + 12*e + 6*e*e)
    {
      unsigned int basisnum = i - 8 - 12*e - 5*e*e;
      i0 = square_number_row[basisnum] + 2;
      i1 = square_number_column[basisnum] + 2;
      i2 = 1;
      const Point min_point = get_min_point(elem, 4, 5, 6, 7);

      if (elem->point(4) == min_point)
        if (elem->point(5) == std::min(elem->point(5), elem->point(7)))
          {
            // Case 1
            xi  = xi_saved;
            eta = eta_saved;
          }
        else
          {
            // Case 2
            xi  = eta_saved;
            eta = xi_saved;
          }

      else if (elem->point(5) == min_point)
        if (elem->point(6) == std::min(elem->point(6), elem->point(4)))
          {
            // Case 3
            xi  = eta_saved;
            eta = -xi_saved;
          }
        else
          {
            // Case 4
            xi  = -xi_saved;
            eta = eta_saved;
          }

      else if (elem->point(6) == min_point)
        if (elem->point(7) == std::min(elem->point(7), elem->point(5)))
          {
            // Case 5
            xi  = -xi_saved;
            eta = -eta_saved;
          }
        else
          {
            // Case 6
            xi  = -eta_saved;
            eta = -xi_saved;
          }

      else if (elem->point(7) == min_point)
        {
          if (elem->point(4) == std::min(elem->point(4), elem->point(6)))
            {
              // Case 7
              xi  = -eta_saved;
              eta = xi_saved;
            }
          else
            {
              // Case 8
              xi  = xi_saved;
              eta = eta_saved;
            }
        }
    }

  // Internal DoFs
  else
    {
      unsigned int basisnum = i - 8 - 12*e - 6*e*e;
      i0 = cube_number_column[basisnum] + 2;
      i1 = cube_number_row[basisnum] + 2;
      i2 = cube_number_page[basisnum] + 2;
    }
}
#endif // LIBMESH_DIM > 2

} // end anonymous namespace



namespace libMesh
{

template <typename RealType>
typename FEShim<3,HIERARCHIC,RealType>::OutputShape FEShim<3,HIERARCHIC,RealType>::shape(const ElemType,
                             const Order,
                             const unsigned int,
                             const Point &)
{
  libmesh_error_msg("Hierarchic shape functions require an Elem for edge/face orientation.");
  return 0.;
}



template <typename RealType>
typename FEShim<3,L2_HIERARCHIC,RealType>::OutputShape FEShim<3,L2_HIERARCHIC,RealType>::shape(const ElemType,
                                const Order,
                                const unsigned int,
                                const Point &)
{
  libmesh_error_msg("Hierarchic shape functions require an Elem for edge/face orientation.");
  return 0.;
}



template <typename RealType>
typename FEShim<3,HIERARCHIC,RealType>::OutputShape FEShim<3,HIERARCHIC,RealType>::shape(const Elem * elem,
                             const Order order,
                             const unsigned int i,
                             const Point & p,
                             const bool add_p_level)
{
  return fe_hierarchic_3D_shape(elem, order, i, p, add_p_level);
}



template <typename RealType>
typename FEShim<3,L2_HIERARCHIC,RealType>::OutputShape FEShim<3,L2_HIERARCHIC,RealType>::shape(const Elem * elem,
                                const Order order,
                                const unsigned int i,
                                const Point & p,
                                const bool add_p_level)
{
  return fe_hierarchic_3D_shape(elem, order, i, p, add_p_level);
}



template <typename RealType>
typename FEShim<3,HIERARCHIC,RealType>::OutputShape FEShim<3,HIERARCHIC,RealType>::shape_deriv(const ElemType,
                                   const Order,
                                   const unsigned int,
                                   const unsigned int,
                                   const Point & )
{
  libmesh_error_msg("Hierarchic shape functions require an Elem for edge/face orientation.");
  return 0.;
}



template <typename RealType>
typename FEShim<3,L2_HIERARCHIC,RealType>::OutputShape FEShim<3,L2_HIERARCHIC,RealType>::shape_deriv(const ElemType,
                                      const Order,
                                      const unsigned int,
                                      const unsigned int,
                                      const Point & )
{
  libmesh_error_msg("Hierarchic shape functions require an Elem for edge/face orientation.");
  return 0.;
}



template <typename RealType>
typename FEShim<3,HIERARCHIC,RealType>::OutputShape FEShim<3,HIERARCHIC,RealType>::shape_deriv(const Elem * elem,
                                   const Order order,
                                   const unsigned int i,
                                   const unsigned int j,
                                   const Point & p,
                                   const bool add_p_level)
{
  return fe_hierarchic_3D_shape_deriv(elem, order, i, j, p, add_p_level);
}


template <typename RealType>
typename FEShim<3,L2_HIERARCHIC,RealType>::OutputShape FEShim<3,L2_HIERARCHIC,RealType>::shape_deriv(const Elem * elem,
                                      const Order order,
                                      const unsigned int i,
                                      const unsigned int j,
                                      const Point & p,
                                      const bool add_p_level)
{
  return fe_hierarchic_3D_shape_deriv(elem, order, i, j, p, add_p_level);
}



#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

template <typename RealType>
typename FEShim<3,HIERARCHIC,RealType>::OutputShape FEShim<3,HIERARCHIC,RealType>::shape_second_deriv(const ElemType,
                                          const Order,
                                          const unsigned int,
                                          const unsigned int,
                                          const Point & )
{
  libmesh_error_msg("Hierarchic shape functions require an Elem for edge/face orientation.");
  return 0.;
}



template <typename RealType>
typename FEShim<3,L2_HIERARCHIC,RealType>::OutputShape FEShim<3,L2_HIERARCHIC,RealType>::shape_second_deriv(const ElemType,
                                             const Order,
                                             const unsigned int,
                                             const unsigned int,
                                             const Point & )
{
  libmesh_error_msg("Hierarchic shape functions require an Elem for edge/face orientation.");
  return 0.;
}



template <typename RealType>
typename FEShim<3,HIERARCHIC,RealType>::OutputShape FEShim<3,HIERARCHIC,RealType>::shape_second_deriv(const Elem * elem,
                                          const Order order,
                                          const unsigned int i,
                                          const unsigned int j,
                                          const Point & p,
                                          const bool add_p_level)
{
  return fe_hierarchic_3D_shape_second_deriv(elem, order, i, j, p, add_p_level);
}



template <typename RealType>
typename FEShim<3,L2_HIERARCHIC,RealType>::OutputShape FEShim<3,L2_HIERARCHIC,RealType>::shape_second_deriv(const Elem * elem,
                                             const Order order,
                                             const unsigned int i,
                                             const unsigned int j,
                                             const Point & p,
                                             const bool add_p_level)
{
  return fe_hierarchic_3D_shape_second_deriv(elem, order, i, j, p, add_p_level);
}

#endif // LIBMESH_ENABLE_SECOND_DERIVATIVES

} // namespace libMesh



namespace
{
using namespace libMesh;

template <typename RealType>
RealType fe_hierarchic_3D_shape(const ElemTempl<RealType> * elem,
                                const Order order,
                                const unsigned int i,
                                const PointTempl<RealType> & p,
                                const bool add_p_level)
{
#if LIBMESH_DIM == 3

  libmesh_assert(elem);
  const ElemType type = elem->type();

  const Order totalorder =
    static_cast<Order>(order+add_p_level*elem->p_level());

  switch (type)
    {
    case HEX8:
    case HEX20:
      libmesh_assert_less (totalorder, 2);
      libmesh_fallthrough();
    case HEX27:
      {
        libmesh_assert_less (i, (totalorder+1u)*(totalorder+1u)*(totalorder+1u));

        // Compute hex shape functions as a tensor-product
        auto xi   = p(0);
        auto eta  = p(1);
        auto zeta = p(2);

        unsigned int i0, i1, i2;

        cube_indices(elem, totalorder, i, xi, eta, zeta, i0, i1, i2);

        return (FEShim<1,HIERARCHIC,RealType>::shape(EDGE3, totalorder, i0, xi)*
                FEShim<1,HIERARCHIC,RealType>::shape(EDGE3, totalorder, i1, eta)*
                FEShim<1,HIERARCHIC,RealType>::shape(EDGE3, totalorder, i2, zeta));
      }

    default:
      libmesh_error_msg("Invalid element type = " << type);
    }

#else // LIBMESH_DIM != 3
  libmesh_ignore(elem, order, i, p, add_p_level);
  libmesh_not_implemented();
#endif
}



template <typename RealType>
RealType fe_hierarchic_3D_shape_deriv(const ElemTempl<RealType> * elem,
                                      const Order order,
                                      const unsigned int i,
                                      const unsigned int j,
                                      const PointTempl<RealType> & p,
                                      const bool add_p_level)
{
  typedef PointTempl<RealType> Point;

#if LIBMESH_DIM == 3
  libmesh_assert(elem);

  libmesh_assert_less (j, 3);

  // cheat by using finite difference approximations:
  const Real eps = 1.e-6;
  Point pp, pm;

  switch (j)
    {
      // d()/dxi
    case 0:
      {
        pp = Point(p(0)+eps, p(1), p(2));
        pm = Point(p(0)-eps, p(1), p(2));
        break;
      }

      // d()/deta
    case 1:
      {
        pp = Point(p(0), p(1)+eps, p(2));
        pm = Point(p(0), p(1)-eps, p(2));
        break;
      }

      // d()/dzeta
    case 2:
      {
        pp = Point(p(0), p(1), p(2)+eps);
        pm = Point(p(0), p(1), p(2)-eps);
        break;
      }

    default:
      libmesh_error_msg("Invalid derivative index j = " << j);
    }

  return (FEShim<3,HIERARCHIC,RealType>::shape(elem, order, i, pp, add_p_level) -
          FEShim<3,HIERARCHIC,RealType>::shape(elem, order, i, pm, add_p_level))/2./eps;
#else // LIBMESH_DIM != 3
  libmesh_ignore(elem, order, i, j, p, add_p_level);
  libmesh_not_implemented();
#endif
}



#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

template <typename RealType>
RealType fe_hierarchic_3D_shape_second_deriv(const ElemTempl<RealType> * elem,
                                             const Order order,
                                             const unsigned int i,
                                             const unsigned int j,
                                             const PointTempl<RealType> & p,
                                             const bool add_p_level)
{
  typedef PointTempl<RealType> Point;

  libmesh_assert(elem);

  const Real eps = 1.e-6;
  Point pp, pm;
  unsigned int prevj = libMesh::invalid_uint;

  switch (j)
    {
      //  d^2()/dxi^2
    case 0:
      {
        pp = Point(p(0)+eps, p(1), p(2));
        pm = Point(p(0)-eps, p(1), p(2));
        prevj = 0;
        break;
      }

      //  d^2()/dxideta
    case 1:
      {
        pp = Point(p(0), p(1)+eps, p(2));
        pm = Point(p(0), p(1)-eps, p(2));
        prevj = 0;
        break;
      }

      //  d^2()/deta^2
    case 2:
      {
        pp = Point(p(0), p(1)+eps, p(2));
        pm = Point(p(0), p(1)-eps, p(2));
        prevj = 1;
        break;
      }

      //  d^2()/dxidzeta
    case 3:
      {
        pp = Point(p(0), p(1), p(2)+eps);
        pm = Point(p(0), p(1), p(2)-eps);
        prevj = 0;
        break;
      }

      //  d^2()/detadzeta
    case 4:
      {
        pp = Point(p(0), p(1), p(2)+eps);
        pm = Point(p(0), p(1), p(2)-eps);
        prevj = 1;
        break;
      }

      //  d^2()/dzeta^2
    case 5:
      {
        pp = Point(p(0), p(1), p(2)+eps);
        pm = Point(p(0), p(1), p(2)-eps);
        prevj = 2;
        break;
      }
    default:
      libmesh_error_msg("Invalid derivative index j = " << j);
    }

  return (FEShim<3,HIERARCHIC,RealType>::shape_deriv(elem, order, i, prevj, pp, add_p_level) -
          FEShim<3,HIERARCHIC,RealType>::shape_deriv(elem, order, i, prevj, pm, add_p_level))
    / 2. / eps;
}

#endif // LIBMESH_ENABLE_SECOND_DERIVATIVES


} // anonymous namespace

#endif // LIBMESH_FE_HIERARCHIC_SHAPE_3D_IMPL_H
