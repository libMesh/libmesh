// The libMesh Finite Element Library.
// Copyright (C) 2002-2021 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


// Local includes
#include "libmesh/fe.h"
#include "libmesh/elem.h"
#include "libmesh/number_lookups.h"
#include "libmesh/enum_to_string.h"

// Anonymous namespace for functions shared by HIERARCHIC and
// L2_HIERARCHIC implementations. Implementations appear at the bottom
// of this file.
namespace
{
using namespace libMesh;

Real fe_triangle_helper (const Elem & elem,
                         const Real edgenumerator,
                         const Real crossval,
                         const unsigned int basisorder,
                         const Order totalorder,
                         const unsigned int noden);

template <FEFamily T>
Real fe_hierarchic_2D_shape(const Elem * elem,
                            const Order order,
                            const unsigned int i,
                            const Point & p,
                            const bool add_p_level);

template <FEFamily T>
Real fe_hierarchic_2D_shape_deriv(const Elem * elem,
                                  const Order order,
                                  const unsigned int i,
                                  const unsigned int j,
                                  const Point & p,
                                  const bool add_p_level);

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

template <FEFamily T>
Real fe_hierarchic_2D_shape_second_deriv(const Elem * elem,
                                         const Order order,
                                         const unsigned int i,
                                         const unsigned int j,
                                         const Point & p,
                                         const bool add_p_level);

#endif // LIBMESH_ENABLE_SECOND_DERIVATIVES

} // anonymous namespace



namespace libMesh
{


LIBMESH_DEFAULT_VECTORIZED_FE(2,HIERARCHIC)
LIBMESH_DEFAULT_VECTORIZED_FE(2,L2_HIERARCHIC)
LIBMESH_DEFAULT_VECTORIZED_FE(2,SIDE_HIERARCHIC)


template <>
Real FE<2,HIERARCHIC>::shape(const ElemType,
                             const Order,
                             const unsigned int,
                             const Point &)
{
  libmesh_error_msg("Hierarchic shape functions require an Elem for edge orientation.");
  return 0.;
}



template <>
Real FE<2,L2_HIERARCHIC>::shape(const ElemType,
                                const Order,
                                const unsigned int,
                                const Point &)
{
  libmesh_error_msg("Hierarchic shape functions require an Elem for edge orientation.");
  return 0.;
}



template <>
Real FE<2,SIDE_HIERARCHIC>::shape(const ElemType,
                                  const Order,
                                  const unsigned int,
                                  const Point &)
{
  libmesh_error_msg("Hierarchic shape functions require an Elem for edge orientation.");
  return 0.;
}



template <>
Real FE<2,HIERARCHIC>::shape(const Elem * elem,
                             const Order order,
                             const unsigned int i,
                             const Point & p,
                             const bool add_p_level)
{
  return fe_hierarchic_2D_shape<HIERARCHIC>(elem, order, i, p, add_p_level);
}



template <>
Real FE<2,HIERARCHIC>::shape(const FEType fet,
                             const Elem * elem,
                             const unsigned int i,
                             const Point & p,
                             const bool add_p_level)
{
  return fe_hierarchic_2D_shape<HIERARCHIC>(elem, fet.order, i, p, add_p_level);
}


template <>
Real FE<2,L2_HIERARCHIC>::shape(const Elem * elem,
                                const Order order,
                                const unsigned int i,
                                const Point & p,
                                const bool add_p_level)
{
  return fe_hierarchic_2D_shape<L2_HIERARCHIC>(elem, order, i, p, add_p_level);
}


template <>
Real FE<2,L2_HIERARCHIC>::shape(const FEType fet,
                                const Elem * elem,
                                const unsigned int i,
                                const Point & p,
                                const bool add_p_level)
{
  return fe_hierarchic_2D_shape<L2_HIERARCHIC>(elem, fet.order, i, p, add_p_level);
}


template <>
Real FE<2,SIDE_HIERARCHIC>::shape(const Elem * elem,
                                  const Order order,
                                  const unsigned int i,
                                  const Point & p,
                                  const bool add_p_level)
{
  libmesh_assert(elem);
  const ElemType type = elem->type();

  const Order totalorder =
    static_cast<Order>(order+add_p_level*elem->p_level());

  const unsigned int dofs_per_side = totalorder+1u;

  switch (type)
    {
    case TRI6:
      {
        libmesh_not_implemented();
        return 0;
      }
    case QUAD8:
    case QUADSHELL8:
    case QUAD9:
      {
        libmesh_assert_less(i, 4*dofs_per_side);

        // Flip odd degree of freedom values if necessary
        // to keep continuity on sides.  We'll flip xi/eta rather than
        // flipping phi, so that we can use this to handle the "nodal"
        // degrees of freedom too.
        Real f = 1.;

        const Real xi = p(0), eta = p(1);
        if (eta < xi)
          {
            if (eta < -xi) // side 0
              {
                if (i >= dofs_per_side)
                  return 0;

                if (totalorder == 0) // special case since raw HIERARCHIC lacks CONSTANTs
                  return 1;

                if ((i < 2 || i % 2) &&
                    elem->point(0) > elem->point(1))
                  f = -1;

                return FE<1,HIERARCHIC>::shape(EDGE3, totalorder, i, f*xi);
              }
            else           // side 1
              {
                if (i < dofs_per_side ||
                    i >= 2*dofs_per_side)
                  return 0;

                if (totalorder == 0) // special case since raw HIERARCHIC lacks CONSTANTs
                  return 1;

                const unsigned int side_i = i - dofs_per_side;

                if ((side_i < 2 || side_i % 2) &&
                    elem->point(1) > elem->point(2))
                  f = -1;

                return FE<1,HIERARCHIC>::shape(EDGE3, totalorder, side_i, f*eta);
              }
          }
        else // xi < eta
          {
            if (eta > -xi)    // side 2
              {
                if (i < 2*dofs_per_side ||
                    i >= 3*dofs_per_side)
                  return 0;

                if (totalorder == 0) // special case since raw HIERARCHIC lacks CONSTANTs
                  return 1;

                const unsigned int side_i = i - 2*dofs_per_side;

                if ((side_i < 2 || side_i % 2) &&
                    elem->point(3) > elem->point(2))
                  f = -1;

                return FE<1,HIERARCHIC>::shape(EDGE3, totalorder, side_i, f*xi);
              }
            else           // side 3
              {
                if (i < 3*dofs_per_side)
                  return 0;

                if (totalorder == 0) // special case since raw HIERARCHIC lacks CONSTANTs
                  return 1;

                const unsigned int side_i = i - 3*dofs_per_side;

                if ((side_i < 2 || side_i % 2) &&
                    elem->point(0) > elem->point(3))
                  f = -1;

                return FE<1,HIERARCHIC>::shape(EDGE3, totalorder, side_i, f*eta);
              }
          }
      }
    default:
      libmesh_error_msg("ERROR: Unsupported element type = " << Utility::enum_to_string(elem->type()));
    }
  return 0;
}


template <>
Real FE<2,SIDE_HIERARCHIC>::shape(const FEType fet,
                                  const Elem * elem,
                                  const unsigned int i,
                                  const Point & p,
                                  const bool add_p_level)
{
  return FE<2,SIDE_HIERARCHIC>::shape(elem, fet.order, i, p, add_p_level);
}


template <>
Real FE<2,HIERARCHIC>::shape_deriv(const ElemType,
                                   const Order,
                                   const unsigned int,
                                   const unsigned int,
                                   const Point &)
{
  libmesh_error_msg("Hierarchic shape functions require an Elem for edge orientation.");
  return 0.;
}



template <>
Real FE<2,L2_HIERARCHIC>::shape_deriv(const ElemType,
                                      const Order,
                                      const unsigned int,
                                      const unsigned int,
                                      const Point &)
{
  libmesh_error_msg("Hierarchic shape functions require an Elem for edge orientation.");
  return 0.;
}



template <>
Real FE<2,SIDE_HIERARCHIC>::shape_deriv(const ElemType,
                                        const Order,
                                        const unsigned int,
                                        const unsigned int,
                                        const Point &)
{
  libmesh_error_msg("Hierarchic shape functions require an Elem for edge orientation.");
  return 0.;
}



template <>
Real FE<2,HIERARCHIC>::shape_deriv(const Elem * elem,
                                   const Order order,
                                   const unsigned int i,
                                   const unsigned int j,
                                   const Point & p,
                                   const bool add_p_level)
{
  return fe_hierarchic_2D_shape_deriv<HIERARCHIC>(elem, order, i, j, p, add_p_level);
}


template <>
Real FE<2,HIERARCHIC>::shape_deriv(const FEType fet,
                                   const Elem * elem,
                                   const unsigned int i,
                                   const unsigned int j,
                                   const Point & p,
                                   const bool add_p_level)
{
  return fe_hierarchic_2D_shape_deriv<HIERARCHIC>(elem, fet.order, i, j, p, add_p_level);
}




template <>
Real FE<2,L2_HIERARCHIC>::shape_deriv(const Elem * elem,
                                      const Order order,
                                      const unsigned int i,
                                      const unsigned int j,
                                      const Point & p,
                                      const bool add_p_level)
{
  return fe_hierarchic_2D_shape_deriv<L2_HIERARCHIC>(elem, order, i, j, p, add_p_level);
}


template <>
Real FE<2,L2_HIERARCHIC>::shape_deriv(const FEType fet,
                                      const Elem * elem,
                                      const unsigned int i,
                                      const unsigned int j,
                                      const Point & p,
                                      const bool add_p_level)
{
  return fe_hierarchic_2D_shape_deriv<L2_HIERARCHIC>(elem, fet.order, i, j, p, add_p_level);
}



template <>
Real FE<2,SIDE_HIERARCHIC>::shape_deriv(const Elem * elem,
                                        const Order order,
                                        const unsigned int i,
                                        const unsigned int j,
                                        const Point & p,
                                        const bool add_p_level)
{
  libmesh_assert(elem);

  const ElemType type = elem->type();

  const Order totalorder =
    static_cast<Order>(order+add_p_level*elem->p_level());

  if (totalorder == 0) // special case since raw HIERARCHIC lacks CONSTANTs
    return 0;

  const unsigned int dofs_per_side = totalorder+1u;

  switch (type)
    {
    case TRI6:
      {
        libmesh_not_implemented();
        return 0;
      }
    case QUAD8:
    case QUADSHELL8:
    case QUAD9:
      {
        libmesh_assert_less(i, 4*dofs_per_side);

        // Flip odd degree of freedom values if necessary
        // to keep continuity on sides.  We'll flip xi/eta rather than
        // flipping phi, so that we can use this to handle the "nodal"
        // degrees of freedom too.
        Real f = 1.;

        const Real xi = p(0), eta = p(1);
        if (eta < xi)
          {
            if (eta < -xi) // side 0
              {
                if (i >= dofs_per_side)
                  return 0;
                if (j != 0)
                  return 0;
                if ((i < 2 || i % 2) &&
                    elem->point(0) > elem->point(1))
                  f = -1;

                return f*FE<1,HIERARCHIC>::shape_deriv(EDGE3, totalorder, i, 0, f*xi);
              }
            else           // side 1
              {
                if (i < dofs_per_side ||
                    i >= 2*dofs_per_side)
                  return 0;
                if (j != 1)
                  return 0;

                const unsigned int side_i = i - dofs_per_side;

                if ((side_i < 2 || side_i % 2) &&
                    elem->point(1) > elem->point(2))
                  f = -1;

                return f*FE<1,HIERARCHIC>::shape_deriv(EDGE3, totalorder, side_i, 0, f*eta);
              }
          }
        else // xi < eta
          {
            if (eta > -xi)    // side 2
              {
                if (i < 2*dofs_per_side ||
                    i >= 3*dofs_per_side)
                  return 0;
                if (j != 0)
                  return 0;

                const unsigned int side_i = i - 2*dofs_per_side;

                if ((side_i < 2 || side_i % 2) &&
                    elem->point(3) > elem->point(2))
                  f = -1;

                return f*FE<1,HIERARCHIC>::shape_deriv(EDGE3, totalorder, side_i, 0, f*xi);
              }
            else           // side 3
              {
                if (i < 3*dofs_per_side)
                  return 0;
                if (j != 1)
                  return 0;

                const unsigned int side_i = i - 3*dofs_per_side;

                if ((side_i < 2 || side_i % 2) &&
                    elem->point(0) > elem->point(3))
                  f = -1;

                return f*FE<1,HIERARCHIC>::shape_deriv(EDGE3, totalorder, side_i, 0, f*eta);
              }
          }
      }
    default:
      libmesh_error_msg("ERROR: Unsupported element type = " << Utility::enum_to_string(elem->type()));
    }
  return 0;
}


template <>
Real FE<2,SIDE_HIERARCHIC>::shape_deriv(const FEType fet,
                                        const Elem * elem,
                                        const unsigned int i,
                                        const unsigned int j,
                                        const Point & p,
                                        const bool add_p_level)
{
  return FE<2,SIDE_HIERARCHIC>::shape_deriv(elem, fet.order, i, j, p, add_p_level);
}


#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

template <>
Real FE<2,HIERARCHIC>::shape_second_deriv(const ElemType,
                                          const Order,
                                          const unsigned int,
                                          const unsigned int,
                                          const Point &)
{
  libmesh_error_msg("Hierarchic shape functions require an Elem for edge orientation.");
  return 0.;
}



template <>
Real FE<2,L2_HIERARCHIC>::shape_second_deriv(const ElemType,
                                             const Order,
                                             const unsigned int,
                                             const unsigned int,
                                             const Point &)
{
  libmesh_error_msg("Hierarchic shape functions require an Elem for edge orientation.");
  return 0.;
}



template <>
Real FE<2,SIDE_HIERARCHIC>::shape_second_deriv(const ElemType,
                                               const Order,
                                               const unsigned int,
                                               const unsigned int,
                                               const Point &)
{
  libmesh_error_msg("Hierarchic shape functions require an Elem for edge orientation.");
  return 0.;
}



template <>
Real FE<2,HIERARCHIC>::shape_second_deriv(const Elem * elem,
                                          const Order order,
                                          const unsigned int i,
                                          const unsigned int j,
                                          const Point & p,
                                          const bool add_p_level)
{
  return fe_hierarchic_2D_shape_second_deriv<HIERARCHIC>(elem, order, i, j, p, add_p_level);
}


template <>
Real FE<2,HIERARCHIC>::shape_second_deriv(const FEType fet,
                                          const Elem * elem,
                                          const unsigned int i,
                                          const unsigned int j,
                                          const Point & p,
                                          const bool add_p_level)
{
  return fe_hierarchic_2D_shape_second_deriv<HIERARCHIC>(elem, fet.order, i, j, p, add_p_level);
}


template <>
Real FE<2,L2_HIERARCHIC>::shape_second_deriv(const Elem * elem,
                                             const Order order,
                                             const unsigned int i,
                                             const unsigned int j,
                                             const Point & p,
                                             const bool add_p_level)
{
  return fe_hierarchic_2D_shape_second_deriv<L2_HIERARCHIC>(elem, order, i, j, p, add_p_level);
}


template <>
Real FE<2,L2_HIERARCHIC>::shape_second_deriv(const FEType fet,
                                             const Elem * elem,
                                             const unsigned int i,
                                             const unsigned int j,
                                             const Point & p,
                                             const bool add_p_level)
{
  return fe_hierarchic_2D_shape_second_deriv<L2_HIERARCHIC>(elem, fet.order, i, j, p, add_p_level);
}


template <>
Real FE<2,SIDE_HIERARCHIC>::shape_second_deriv(const Elem * elem,
                                               const Order order,
                                               const unsigned int i,
                                               const unsigned int j,
                                               const Point & p,
                                               const bool add_p_level)
{
  libmesh_assert(elem);
  const ElemType type = elem->type();

  const Order totalorder =
    static_cast<Order>(order+add_p_level*elem->p_level());

  if (totalorder == 0) // special case since raw HIERARCHIC lacks CONSTANTs
    return 0;

  const unsigned int dofs_per_side = totalorder+1u;

  switch (type)
    {
    case TRI6:
      {
        libmesh_not_implemented();
        return 0;
      }
    case QUAD8:
    case QUADSHELL8:
    case QUAD9:
      {
        libmesh_assert_less(i, 4*dofs_per_side);

        // Flip odd degree of freedom values if necessary
        // to keep continuity on sides.  We'll flip xi/eta rather than
        // flipping phi, so that we can use this to handle the "nodal"
        // degrees of freedom too.
        Real f = 1.;

        const Real xi = p(0), eta = p(1);
        if (eta < xi)
          {
            if (eta < -xi) // side 0
              {
                if (i >= dofs_per_side)
                  return 0;
                if (j != 0)
                  return 0;
                if ((i < 2 || i % 2) &&
                    elem->point(0) > elem->point(1))
                  f = -1;

                return FE<1,HIERARCHIC>::shape_second_deriv(EDGE3, totalorder, i, 0, f*xi);
              }
            else           // side 1
              {
                if (i < dofs_per_side ||
                    i >= 2*dofs_per_side)
                  return 0;
                if (j != 2)
                  return 0;

                const unsigned int side_i = i - dofs_per_side;

                if ((side_i < 2 || side_i % 2) &&
                    elem->point(1) > elem->point(2))
                  f = -1;

                return FE<1,HIERARCHIC>::shape_second_deriv(EDGE3, totalorder, side_i, 0, f*eta);
              }
          }
        else // xi < eta
          {
            if (eta > -xi)    // side 2
              {
                if (i < 2*dofs_per_side ||
                    i >= 3*dofs_per_side)
                  return 0;
                if (j != 0)
                  return 0;

                const unsigned int side_i = i - 2*dofs_per_side;

                if ((side_i < 2 || side_i % 2) &&
                    elem->point(3) > elem->point(2))
                  f = -1;

                return FE<1,HIERARCHIC>::shape_second_deriv(EDGE3, totalorder, side_i, 0, f*xi);
              }
            else           // side 3
              {
                if (i < 3*dofs_per_side)
                  return 0;
                if (j != 2)
                  return 0;

                const unsigned int side_i = i - 3*dofs_per_side;

                if ((side_i < 2 || side_i % 2) &&
                    elem->point(0) > elem->point(3))
                  f = -1;

                return FE<1,HIERARCHIC>::shape_second_deriv(EDGE3, totalorder, side_i, 0, f*eta);
              }
          }
      }
    default:
      libmesh_error_msg("ERROR: Unsupported element type = " << Utility::enum_to_string(elem->type()));
    }
  return 0;
}


template <>
Real FE<2,SIDE_HIERARCHIC>::shape_second_deriv(const FEType fet,
                                               const Elem * elem,
                                               const unsigned int i,
                                               const unsigned int j,
                                               const Point & p,
                                               const bool add_p_level)
{
  return FE<2,SIDE_HIERARCHIC>::shape_second_deriv(elem, fet.order, i, j, p, add_p_level);
}

#endif //  LIBMESH_ENABLE_SECOND_DERIVATIVES

} // namespace libMesh



namespace
{
using namespace libMesh;

Real fe_triangle_helper (const Elem & elem,
                         const Real edgenumerator,
                         const Real crossval,
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

  const Real edgeval = edgenumerator / crossval;
  const Real crossfunc = std::pow(crossval, basisorder);

  return flip * crossfunc *
    FE<1,HIERARCHIC>::shape(EDGE3, totalorder,
                            basisorder, edgeval);
}

template <FEFamily T>
Real fe_hierarchic_2D_shape(const Elem * elem,
                            const Order order,
                            const unsigned int i,
                            const Point & p,
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
        const Real zeta1 = p(0);
        const Real zeta2 = p(1);
        const Real zeta0 = 1. - zeta1 - zeta2;

        libmesh_assert_less (i, (totalorder+1u)*(totalorder+2u)/2);
        libmesh_assert (T == L2_HIERARCHIC || elem->type() == TRI6 || totalorder < 2);

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

            const Real crossval = zeta0 + zeta1;
            const Real edgenumerator = zeta1 - zeta0;

            return fe_triangle_helper(*elem, edgenumerator, crossval,
                                      basisorder, totalorder, 0);
          }
        else if (i < 2u*totalorder + 1)
          {
            const unsigned int basisorder = i - totalorder;

            const Real crossval = zeta2 + zeta1;
            const Real edgenumerator = zeta2 - zeta1;

            return fe_triangle_helper(*elem, edgenumerator, crossval,
                                      basisorder, totalorder, 1);
          }
        else if (i < 3u*totalorder)
          {
            const unsigned int basisorder = i - (2u*totalorder) + 1;

            const Real crossval = zeta0 + zeta2;
            const Real edgenumerator = zeta0 - zeta2;

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
      libmesh_assert (T == L2_HIERARCHIC || totalorder < 2);
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

        return f*(FE<1,T>::shape(EDGE3, totalorder, i0, xi)*
                  FE<1,T>::shape(EDGE3, totalorder, i1, eta));
      }

    default:
      libmesh_error_msg("ERROR: Unsupported element type = " << Utility::enum_to_string(elem->type()));
    }

  return 0.;
}


template <FEFamily T>
Real fe_hierarchic_2D_shape_deriv(const Elem * elem,
                                  const Order order,
                                  const unsigned int i,
                                  const unsigned int j,
                                  const Point & p,
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
              const Point pp(p(0)+eps, p(1));
              const Point pm(p(0)-eps, p(1));

              return (FE<2,T>::shape(elem, order, i, pp) -
                      FE<2,T>::shape(elem, order, i, pm))/2./eps;
            }

            // d()/deta
          case 1:
            {
              const Point pp(p(0), p(1)+eps);
              const Point pm(p(0), p(1)-eps);

              return (FE<2,T>::shape(elem, order, i, pp) -
                      FE<2,T>::shape(elem, order, i, pm))/2./eps;
            }

          default:
            libmesh_error_msg("Invalid derivative index j = " << j);
          }
      }

    case QUAD4:
    case QUADSHELL4:
      libmesh_assert (T == L2_HIERARCHIC || totalorder < 2);
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
            return f*(FE<1,T>::shape_deriv(EDGE3, totalorder, i0, 0, xi)*
                      FE<1,T>::shape      (EDGE3, totalorder, i1,    eta));

            // d()/deta
          case 1:
            return f*(FE<1,T>::shape      (EDGE3, totalorder, i0,    xi)*
                      FE<1,T>::shape_deriv(EDGE3, totalorder, i1, 0, eta));

          default:
            libmesh_error_msg("Invalid derivative index j = " << j);
          }
      }

    default:
      libmesh_error_msg("ERROR: Unsupported element type = " << Utility::enum_to_string(type));
    }

  return 0.;
}



#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

template <FEFamily T>
Real fe_hierarchic_2D_shape_second_deriv(const Elem * elem,
                                         const Order order,
                                         const unsigned int i,
                                         const unsigned int j,
                                         const Point & p,
                                         const bool add_p_level)
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

  return (FE<2,T>::shape_deriv(elem, order, i, prevj, pp, add_p_level) -
          FE<2,T>::shape_deriv(elem, order, i, prevj, pm, add_p_level)
          )/2./eps;
}

#endif // LIBMESH_ENABLE_SECOND_DERIVATIVES

} // anonymous namespace
