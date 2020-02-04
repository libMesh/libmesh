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

#ifndef LIBMESH_FE_HIERARCHIC_SHAPE_2D_IMPL_H
#define LIBMESH_FE_HIERARCHIC_SHAPE_2D_IMPL_H

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
RealType fe_triangle_helper (const ElemTempl<RealType> & elem,
                         const RealType & edgenumerator,
                         const RealType & crossval,
                         const unsigned int basisorder,
                         const Order totalorder,
                         const unsigned int noden);

template <typename RealType>
RealType fe_hierarchic_2D_shape(const ElemTempl<RealType> * elem,
                            const Order order,
                            const unsigned int i,
                            const PointTempl<RealType> & p,
                            const bool add_p_level);

template <typename RealType>
RealType fe_hierarchic_2D_shape_deriv(const ElemTempl<RealType> * elem,
                                  const Order order,
                                  const unsigned int i,
                                  const unsigned int j,
                                  const PointTempl<RealType> & p,
                                  const bool add_p_level);

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

template <typename RealType>
RealType fe_hierarchic_2D_shape_second_deriv(const ElemTempl<RealType> * elem,
                                         const Order order,
                                         const unsigned int i,
                                         const unsigned int j,
                                         const PointTempl<RealType> & p,
                                         const bool add_p_level);

#endif // LIBMESH_ENABLE_SECOND_DERIVATIVES

} // anonymous namespace



namespace libMesh
{

template <typename RealType>
typename FEShim<2,HIERARCHIC,RealType>::OutputShape FEShim<2,HIERARCHIC,RealType>::shape(const ElemType,
                             const Order,
                             const unsigned int,
                             const Point &)
{
  libmesh_error_msg("Hierarchic shape functions require an Elem for edge orientation.");
  return 0.;
}



template <typename RealType>
typename FEShim<2,L2_HIERARCHIC,RealType>::OutputShape FEShim<2,L2_HIERARCHIC,RealType>::shape(const ElemType,
                                const Order,
                                const unsigned int,
                                const Point &)
{
  libmesh_error_msg("Hierarchic shape functions require an Elem for edge orientation.");
  return 0.;
}



template <typename RealType>
typename FEShim<2,HIERARCHIC,RealType>::OutputShape FEShim<2,HIERARCHIC,RealType>::shape(const ElemTempl<RealType> * elem,
                             const Order order,
                             const unsigned int i,
                             const Point & p,
                             const bool add_p_level)
{
  return fe_hierarchic_2D_shape(elem, order, i, p, add_p_level);
}



template <typename RealType>
typename FEShim<2,L2_HIERARCHIC,RealType>::OutputShape FEShim<2,L2_HIERARCHIC,RealType>::shape(const ElemTempl<RealType> * elem,
                                const Order order,
                                const unsigned int i,
                                const Point & p,
                                const bool add_p_level)
{
  return fe_hierarchic_2D_shape(elem, order, i, p, add_p_level);
}



template <typename RealType>
typename FEShim<2,HIERARCHIC,RealType>::OutputShape FEShim<2,HIERARCHIC,RealType>::shape_deriv(const ElemType,
                                   const Order,
                                   const unsigned int,
                                   const unsigned int,
                                   const Point &)
{
  libmesh_error_msg("Hierarchic shape functions require an Elem for edge orientation.");
  return 0.;
}



template <typename RealType>
typename FEShim<2,L2_HIERARCHIC,RealType>::OutputShape FEShim<2,L2_HIERARCHIC,RealType>::shape_deriv(const ElemType,
                                      const Order,
                                      const unsigned int,
                                      const unsigned int,
                                      const Point &)
{
  libmesh_error_msg("Hierarchic shape functions require an Elem for edge orientation.");
  return 0.;
}



template <typename RealType>
typename FEShim<2,HIERARCHIC,RealType>::OutputShape FEShim<2,HIERARCHIC,RealType>::shape_deriv(const ElemTempl<RealType> * elem,
                                   const Order order,
                                   const unsigned int i,
                                   const unsigned int j,
                                   const Point & p,
                                   const bool add_p_level)
{
  return fe_hierarchic_2D_shape_deriv(elem, order, i, j, p, add_p_level);
}



template <typename RealType>
typename FEShim<2,L2_HIERARCHIC,RealType>::OutputShape FEShim<2,L2_HIERARCHIC,RealType>::shape_deriv(const ElemTempl<RealType> * elem,
                                      const Order order,
                                      const unsigned int i,
                                      const unsigned int j,
                                      const Point & p,
                                      const bool add_p_level)
{
  return fe_hierarchic_2D_shape_deriv(elem, order, i, j, p, add_p_level);
}



#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

template <typename RealType>
typename FEShim<2,HIERARCHIC,RealType>::OutputShape FEShim<2,HIERARCHIC,RealType>::shape_second_deriv(const ElemType,
                                          const Order,
                                          const unsigned int,
                                          const unsigned int,
                                          const Point &)
{
  libmesh_error_msg("Hierarchic shape functions require an Elem for edge orientation.");
  return 0.;
}



template <typename RealType>
typename FEShim<2,L2_HIERARCHIC,RealType>::OutputShape FEShim<2,L2_HIERARCHIC,RealType>::shape_second_deriv(const ElemType,
                                             const Order,
                                             const unsigned int,
                                             const unsigned int,
                                             const Point &)
{
  libmesh_error_msg("Hierarchic shape functions require an Elem for edge orientation.");
  return 0.;
}



template <typename RealType>
typename FEShim<2,HIERARCHIC,RealType>::OutputShape FEShim<2,HIERARCHIC,RealType>::shape_second_deriv(const ElemTempl<RealType> * elem,
                                          const Order order,
                                          const unsigned int i,
                                          const unsigned int j,
                                          const Point & p,
                                          const bool add_p_level)
{
  return fe_hierarchic_2D_shape_second_deriv(elem, order, i, j, p, add_p_level);
}



template <typename RealType>
typename FEShim<2,L2_HIERARCHIC,RealType>::OutputShape FEShim<2,L2_HIERARCHIC,RealType>::shape_second_deriv(const ElemTempl<RealType> * elem,
                                             const Order order,
                                             const unsigned int i,
                                             const unsigned int j,
                                             const Point & p,
                                             const bool add_p_level)
{
  return fe_hierarchic_2D_shape_second_deriv(elem, order, i, j, p, add_p_level);
}

#endif //  LIBMESH_ENABLE_SECOND_DERIVATIVES

} // namespace libMesh



namespace
{
using namespace libMesh;

template <typename RealType>
RealType fe_triangle_helper (const ElemTempl<RealType> & elem,
                         const RealType & edgenumerator,
                         const RealType & crossval,
                         const unsigned int basisorder,
                         const Order totalorder,
                         const unsigned int noden)
{
  // Get factors to account for edge-flipping
  Real flip = 1;
  if (basisorder%2 && (elem.point(noden) > elem.point((noden+1)%3)))
    flip = -1.;

  // Avoid NaN around vertices!
  if (crossval == 0.)
    {
      unsigned int basisfactorial = 1.;
      for (unsigned int n=2; n <= basisorder; ++n)
        basisfactorial *= n;

      return std::pow(edgenumerator, basisorder) / basisfactorial;
    }
  // FIXME - what happens with roundoff when 0 < crossval < O(epsilon)?

  const auto edgeval = edgenumerator / crossval;
  const auto crossfunc = std::pow(crossval, basisorder);

  return flip * crossfunc *
    FEShim<1,HIERARCHIC,RealType>::shape(EDGE3, totalorder,
                            basisorder, edgeval);
}

template <typename RealType>
RealType fe_hierarchic_2D_shape(const ElemTempl<RealType> * elem,
                            const Order order,
                            const unsigned int i,
                            const PointTempl<RealType> & p,
                            const bool add_p_level)
{
  libmesh_assert(elem);

  const Order totalorder =
    static_cast<Order>(order+add_p_level*elem->p_level());
  libmesh_assert_greater (totalorder, 0);

  switch (elem->type())
    {
    case TRI3:
    case TRISHELL3:
    case TRI6:
      {
        const auto zeta1 = p(0);
        const auto zeta2 = p(1);
        const auto zeta0 = 1. - zeta1 - zeta2;

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
            const unsigned int basisorder = i - 1;

            const auto crossval = zeta0 + zeta1;
            const auto edgenumerator = zeta1 - zeta0;

            return fe_triangle_helper(*elem, edgenumerator, crossval,
                                      basisorder, totalorder, 0);
          }
        else if (i < 2u*totalorder + 1)
          {
            const unsigned int basisorder = i - totalorder;

            const auto crossval = zeta2 + zeta1;
            const auto edgenumerator = zeta2 - zeta1;

            return fe_triangle_helper(*elem, edgenumerator, crossval,
                                      basisorder, totalorder, 1);
          }
        else if (i < 3u*totalorder)
          {
            const unsigned int basisorder = i - (2u*totalorder) + 1;

            const auto crossval = zeta0 + zeta2;
            const auto edgenumerator = zeta0 - zeta2;

            return fe_triangle_helper(*elem, edgenumerator, crossval,
                                      basisorder, totalorder, 2);
          }
        // Interior DoFs
        else
          {
            const unsigned int basisnum = i - (3u*totalorder);
            unsigned int exp0 = triangular_number_column[basisnum] + 1;
            unsigned int exp1 = triangular_number_row[basisnum] + 1 -
              triangular_number_column[basisnum];

            RealType returnval = 1;
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

        return f*(FEShim<1,HIERARCHIC,RealType>::shape(EDGE3, totalorder, i0, xi)*
                  FEShim<1,HIERARCHIC,RealType>::shape(EDGE3, totalorder, i1, eta));
      }

    default:
      libmesh_error_msg("ERROR: Unsupported element type = " << elem->type());
    }

  return 0.;
}



template <typename RealType>
RealType fe_hierarchic_2D_shape_deriv(const ElemTempl<RealType> * elem,
                                  const Order order,
                                  const unsigned int i,
                                  const unsigned int j,
                                  const PointTempl<RealType> & p,
                                  const bool add_p_level)
{
  libmesh_assert(elem);

  const ElemType type = elem->type();

  const Order totalorder =
    static_cast<Order>(order+add_p_level*elem->p_level());

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
              const PointTempl<RealType> pp(p(0)+eps, p(1));
              const PointTempl<RealType> pm(p(0)-eps, p(1));

              return (FEShim<2,HIERARCHIC,RealType>::shape(elem, order, i, pp) -
                      FEShim<2,HIERARCHIC,RealType>::shape(elem, order, i, pm))/2./eps;
            }

            // d()/deta
          case 1:
            {
              const PointTempl<RealType> pp(p(0), p(1)+eps);
              const PointTempl<RealType> pm(p(0), p(1)-eps);

              return (FEShim<2,HIERARCHIC,RealType>::shape(elem, order, i, pp) -
                      FEShim<2,HIERARCHIC,RealType>::shape(elem, order, i, pm))/2./eps;
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
        const auto xi  = p(0);
        const auto eta = p(1);

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
            return f*(FEShim<1,HIERARCHIC,RealType>::shape_deriv(EDGE3, totalorder, i0, 0, xi)*
                      FEShim<1,HIERARCHIC,RealType>::shape      (EDGE3, totalorder, i1,    eta));

            // d()/deta
          case 1:
            return f*(FEShim<1,HIERARCHIC,RealType>::shape      (EDGE3, totalorder, i0,    xi)*
                      FEShim<1,HIERARCHIC,RealType>::shape_deriv(EDGE3, totalorder, i1, 0, eta));

          default:
            libmesh_error_msg("Invalid derivative index j = " << j);
          }
      }

    default:
      libmesh_error_msg("ERROR: Unsupported element type = " << type);
    }

  return 0.;
}



#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

template <typename RealType>
RealType fe_hierarchic_2D_shape_second_deriv(const ElemTempl<RealType> * elem,
                                         const Order order,
                                         const unsigned int i,
                                         const unsigned int j,
                                         const PointTempl<RealType> & p,
                                         const bool add_p_level)
{
  typedef PointTempl<RealType> Point;

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

  return (FEShim<2,HIERARCHIC,RealType>::shape_deriv(elem, order, i, prevj, pp, add_p_level) -
          FEShim<2,HIERARCHIC,RealType>::shape_deriv(elem, order, i, prevj, pm, add_p_level)
          )/2./eps;
}

#endif // LIBMESH_ENABLE_SECOND_DERIVATIVES

} // anonymous namespace

#endif // LIBMESH_FE_HIERARCHIC_SHAPE_2D_IMPL_H
