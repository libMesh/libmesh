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
#include "libmesh/cell_tet4.h" // We need edge_nodes_map + side_nodes_map

// Anonymous namespace for functions shared by HIERARCHIC and
// L2_HIERARCHIC implementations. Implementations appear at the bottom
// of this file.
namespace
{
using namespace libMesh;

unsigned int cube_side(const Point & p);

Point cube_side_point(unsigned int sidenum, const Point & interior_point);

template <FEFamily T>
Real fe_hierarchic_3D_shape(const Elem * elem,
                            const Order order,
                            const unsigned int i,
                            const Point & p,
                            const bool add_p_level);

template <FEFamily T>
Real fe_hierarchic_3D_shape_deriv(const Elem * elem,
                                  const Order order,
                                  const unsigned int i,
                                  const unsigned int j,
                                  const Point & p,
                                  const bool add_p_level);

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

template <FEFamily T>
Real fe_hierarchic_3D_shape_second_deriv(const Elem * elem,
                                         const Order order,
                                         const unsigned int i,
                                         const unsigned int j,
                                         const Point & p,
                                         const bool add_p_level);

#endif // LIBMESH_ENABLE_SECOND_DERIVATIVES

#if LIBMESH_DIM > 2
Point get_min_point(const Elem * elem,
                    unsigned int a,
                    unsigned int b,
                    unsigned int c,
                    unsigned int d)
{
  return std::min(std::min(elem->point(a),elem->point(b)),
                  std::min(elem->point(c),elem->point(d)));
}

// Remap non-face-nodes based on point ordering
template <unsigned int N_nodes>
unsigned int remap_node(unsigned int n,
                        const Elem & elem,
                        unsigned int nodebegin)
{
  std::array<const Point *, N_nodes> points;

  for (auto i : IntRange<unsigned int>(0, N_nodes))
    points[i] = &elem.point(nodebegin+i);

  std::sort(points.begin(), points.end(),
            [](const Point * a, const Point * b)
            { return *a < *b; });

  const Point * pn = points[n-nodebegin];

  for (auto i : IntRange<unsigned int>(nodebegin, nodebegin+N_nodes))
    if (pn == &elem.point(i))
      return i;

  libmesh_assert(false);
  return libMesh::invalid_uint;
}


void cube_remap(unsigned int & side_i,
                const Elem & side,
                unsigned int totalorder,
                Point & sidep)
{
  // "vertex" nodes are now decoupled from vertices, so we have
  // to order them consistently otherwise
  if (side_i < 4)
    side_i = remap_node<4>(side_i, side, 0);

  // And "edge" nodes are decoupled from edges, so we have to
  // reorder them too!
  else if (side_i < 4u*totalorder)
    {
      unsigned int side_node = (side_i - 4)/(totalorder-1)+4;
      side_node = remap_node<4>(side_node, side, 4);
      side_i = ((side_i - 4) % (totalorder - 1)) // old local edge_i
        + 4 + (side_node-4)*(totalorder-1);
    }

  // Interior dofs in 2D don't care about where xi/eta point in
  // physical space, but here we need them to match from both
  // sides of a face!
  else
    {
      unsigned int min_side_node = remap_node<4>(0, side, 0);
      const bool flip = (side.point(min_side_node) < side.point((min_side_node+1)%4));

      switch (min_side_node) {
      case 0:
        if (flip)
          std::swap(sidep(0), sidep(1));
        break;
      case 1:
        sidep(0) = -sidep(0);
        if (!flip)
          std::swap(sidep(0), sidep(1));
        break;
      case 2:
        sidep(0) = -sidep(0);
        sidep(1) = -sidep(1);
        if (flip)
          std::swap(sidep(0), sidep(1));
        break;
      case 3:
        sidep(1) = -sidep(1);
        if (!flip)
          std::swap(sidep(0), sidep(1));
        break;
      default:
        libmesh_error();
      }
    }
}


void cube_indices(const Elem * elem,
                  const unsigned int totalorder,
                  const unsigned int i,
                  Real & xi, Real & eta, Real & zeta,
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

  Real xi_saved = xi, eta_saved = eta, zeta_saved = zeta;

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


LIBMESH_DEFAULT_VECTORIZED_FE(3,HIERARCHIC)
LIBMESH_DEFAULT_VECTORIZED_FE(3,L2_HIERARCHIC)
LIBMESH_DEFAULT_VECTORIZED_FE(3,SIDE_HIERARCHIC)


template <>
Real FE<3,HIERARCHIC>::shape(const ElemType,
                             const Order,
                             const unsigned int,
                             const Point &)
{
  libmesh_error_msg("Hierarchic shape functions require an Elem for edge/face orientation.");
  return 0.;
}



template <>
Real FE<3,L2_HIERARCHIC>::shape(const ElemType,
                                const Order,
                                const unsigned int,
                                const Point &)
{
  libmesh_error_msg("Hierarchic shape functions require an Elem for edge/face orientation.");
  return 0.;
}



template <>
Real FE<3,SIDE_HIERARCHIC>::shape(const ElemType,
                                  const Order,
                                  const unsigned int,
                                  const Point &)
{
  libmesh_error_msg("Hierarchic shape functions require an Elem for edge/face orientation.");
  return 0.;
}



template <>
Real FE<3,HIERARCHIC>::shape(const Elem * elem,
                             const Order order,
                             const unsigned int i,
                             const Point & p,
                             const bool add_p_level)
{
  return fe_hierarchic_3D_shape<HIERARCHIC>(elem, order, i, p, add_p_level);
}


template <>
Real FE<3,HIERARCHIC>::shape(const FEType fet,
                             const Elem * elem,
                             const unsigned int i,
                             const Point & p,
                             const bool add_p_level)
{
  return fe_hierarchic_3D_shape<HIERARCHIC>(elem, fet.order, i, p, add_p_level);
}




template <>
Real FE<3,L2_HIERARCHIC>::shape(const Elem * elem,
                                const Order order,
                                const unsigned int i,
                                const Point & p,
                                const bool add_p_level)
{
  return fe_hierarchic_3D_shape<L2_HIERARCHIC>(elem, order, i, p, add_p_level);
}


template <>
Real FE<3,L2_HIERARCHIC>::shape(const FEType fet,
                                const Elem * elem,
                                const unsigned int i,
                                const Point & p,
                                const bool add_p_level)
{
  return fe_hierarchic_3D_shape<L2_HIERARCHIC>(elem, fet.order, i, p, add_p_level);
}



template <>
Real FE<3,SIDE_HIERARCHIC>::shape(const Elem * elem,
                                  const Order order,
                                  const unsigned int i,
                                  const Point & p,
                                  const bool add_p_level)
{
#if LIBMESH_DIM == 3
  libmesh_assert(elem);
  const ElemType type = elem->type();

  const Order totalorder =
    static_cast<Order>(order+add_p_level*elem->p_level());

  switch (type)
    {
    case HEX27:
      {
        const unsigned int dofs_per_side = (totalorder+1u)*(totalorder+1u);
        libmesh_assert_less(i, 6*dofs_per_side);

        const unsigned int sidenum = cube_side(p);
        const unsigned int dof_offset = sidenum * dofs_per_side;

        if (i < dof_offset) // i is on a previous side
          return 0;

        if (i >= dof_offset + dofs_per_side) // i is on a later side
          return 0;

        if (totalorder == 0) // special case since raw HIERARCHIC lacks CONSTANTs
          return 1;

        unsigned int side_i = i - dof_offset;

        std::unique_ptr<const Elem> side = elem->build_side_ptr(sidenum);

        Point sidep = cube_side_point(sidenum, p);

        cube_remap(side_i, *side, totalorder, sidep);

        return FE<2,HIERARCHIC>::shape(side.get(), order, side_i, sidep, add_p_level);
      }

    default:
      libmesh_error_msg("Invalid element type = " << Utility::enum_to_string(type));
    }

#else // LIBMESH_DIM != 3
  libmesh_ignore(elem, order, i, p, add_p_level);
  libmesh_not_implemented();
#endif
}


template <>
Real FE<3,SIDE_HIERARCHIC>::shape(const FEType fet,
                                  const Elem * elem,
                                  const unsigned int i,
                                  const Point & p,
                                  const bool add_p_level)
{
  return FE<3,SIDE_HIERARCHIC>::shape(elem,fet.order, i, p, add_p_level);
}


template <>
Real FE<3,HIERARCHIC>::shape_deriv(const ElemType,
                                   const Order,
                                   const unsigned int,
                                   const unsigned int,
                                   const Point & )
{
  libmesh_error_msg("Hierarchic shape functions require an Elem for edge/face orientation.");
  return 0.;
}



template <>
Real FE<3,L2_HIERARCHIC>::shape_deriv(const ElemType,
                                      const Order,
                                      const unsigned int,
                                      const unsigned int,
                                      const Point & )
{
  libmesh_error_msg("Hierarchic shape functions require an Elem for edge/face orientation.");
  return 0.;
}



template <>
Real FE<3,SIDE_HIERARCHIC>::shape_deriv(const ElemType,
                                        const Order,
                                        const unsigned int,
                                        const unsigned int,
                                        const Point & )
{
  libmesh_error_msg("Hierarchic shape functions require an Elem for edge/face orientation.");
  return 0.;
}



template <>
Real FE<3,HIERARCHIC>::shape_deriv(const Elem * elem,
                                   const Order order,
                                   const unsigned int i,
                                   const unsigned int j,
                                   const Point & p,
                                   const bool add_p_level)
{
  return fe_hierarchic_3D_shape_deriv<HIERARCHIC>(elem, order, i, j, p, add_p_level);
}


template <>
Real FE<3,HIERARCHIC>::shape_deriv(const FEType fet,
                                   const Elem * elem,
                                   const unsigned int i,
                                   const unsigned int j,
                                   const Point & p,
                                   const bool add_p_level)
{
  return fe_hierarchic_3D_shape_deriv<HIERARCHIC>(elem, fet.order, i, j, p, add_p_level);
}



template <>
Real FE<3,L2_HIERARCHIC>::shape_deriv(const Elem * elem,
                                      const Order order,
                                      const unsigned int i,
                                      const unsigned int j,
                                      const Point & p,
                                      const bool add_p_level)
{
  return fe_hierarchic_3D_shape_deriv<L2_HIERARCHIC>(elem, order, i, j, p, add_p_level);
}


template <>
Real FE<3,L2_HIERARCHIC>::shape_deriv(const FEType fet,
                                      const Elem * elem,
                                      const unsigned int i,
                                      const unsigned int j,
                                      const Point & p,
                                      const bool add_p_level)
{
  return fe_hierarchic_3D_shape_deriv<L2_HIERARCHIC>(elem, fet.order, i, j, p, add_p_level);
}



template <>
Real FE<3,SIDE_HIERARCHIC>::shape_deriv(const Elem * elem,
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
    static_cast<Order>(order+add_p_level*elem->p_level());

  if (totalorder == 0) // special case since raw HIERARCHIC lacks CONSTANTs
    return 0; // constants have zero derivative

  switch (type)
    {
    case HEX27:
      {
        const unsigned int dofs_per_side = (totalorder+1u)*(totalorder+1u);
        libmesh_assert_less(i, 6*dofs_per_side);

        const unsigned int sidenum = cube_side(p);
        const unsigned int dof_offset = sidenum * dofs_per_side;

        if (i < dof_offset) // i is on a previous side
          return 0;

        if (i >= dof_offset + dofs_per_side) // i is on a later side
          return 0;

        unsigned int side_i = i - dof_offset;

        std::unique_ptr<const Elem> side = elem->build_side_ptr(sidenum);

        Point sidep = cube_side_point(sidenum, p);

        cube_remap(side_i, *side, totalorder, sidep);

        // What direction on the side corresponds to the derivative
        // direction we want?
        unsigned int sidej = 100;

        // Do we need a -1 here to flip that direction?
        Real f = 1.;

        switch (j)
          {
          case 0: // d()/dxi
            {
              switch (sidenum)
                {
                case 0:
                  sidej = 1;
                  break;
                case 1:
                  sidej = 0;
                  break;
                case 2:
                  return 0;
                case 3:
                  sidej = 0;
                  f = -1;
                  break;
                case 4:
                  return 0;
                case 5:
                  sidej = 0;
                  break;
                default:
                  libmesh_error();
                }
              break;
            }
          case 1: // d()/deta
            {
              switch (sidenum)
                {
                case 0:
                  sidej = 0;
                  break;
                case 1:
                  return 0;
                case 2:
                  sidej = 0;
                  break;
                case 3:
                  return 0;
                case 4:
                  sidej = 0;
                  f = -1;
                  break;
                case 5:
                  sidej = 1;
                  break;
                default:
                  libmesh_error();
                }
              break;
            }
          case 2: // d()/dzeta
            {
              switch (sidenum)
                {
                case 0:
                  return 0;
                case 1:
                case 2:
                case 3:
                case 4:
                  sidej = 1;
                  break;
                case 5:
                  return 0;
                default:
                  libmesh_error();
                }
              break;
            }

          default:
            libmesh_error_msg("Invalid derivative index j = " << j);
          }

        return f * FE<2,HIERARCHIC>::shape_deriv(side.get(), order,
                                                 side_i, sidej, sidep,
                                                 add_p_level);
      }

    default:
      libmesh_error_msg("Invalid element type = " << Utility::enum_to_string(type));
    }

#else // LIBMESH_DIM != 3
  libmesh_ignore(elem, order, i, j, p, add_p_level);
  libmesh_not_implemented();
#endif
}


template <>
Real FE<3,SIDE_HIERARCHIC>::shape_deriv(const FEType fet,
                                        const Elem * elem,
                                        const unsigned int i,
                                        const unsigned int j,
                                        const Point & p,
                                        const bool add_p_level)
{
  return FE<3,SIDE_HIERARCHIC>::shape_deriv(elem, fet.order, i, j, p, add_p_level);
}


#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

template <>
Real FE<3,HIERARCHIC>::shape_second_deriv(const ElemType,
                                          const Order,
                                          const unsigned int,
                                          const unsigned int,
                                          const Point & )
{
  libmesh_error_msg("Hierarchic shape functions require an Elem for edge/face orientation.");
  return 0.;
}



template <>
Real FE<3,L2_HIERARCHIC>::shape_second_deriv(const ElemType,
                                             const Order,
                                             const unsigned int,
                                             const unsigned int,
                                             const Point & )
{
  libmesh_error_msg("Hierarchic shape functions require an Elem for edge/face orientation.");
  return 0.;
}



template <>
Real FE<3,SIDE_HIERARCHIC>::shape_second_deriv(const ElemType,
                                               const Order,
                                               const unsigned int,
                                               const unsigned int,
                                               const Point & )
{
  libmesh_error_msg("Hierarchic shape functions require an Elem for edge/face orientation.");
  return 0.;
}



template <>
Real FE<3,HIERARCHIC>::shape_second_deriv(const Elem * elem,
                                          const Order order,
                                          const unsigned int i,
                                          const unsigned int j,
                                          const Point & p,
                                          const bool add_p_level)
{
  return fe_hierarchic_3D_shape_second_deriv<HIERARCHIC>(elem, order, i, j, p, add_p_level);
}



template <>
Real FE<3,HIERARCHIC>::shape_second_deriv(const FEType fet,
                                          const Elem * elem,
                                          const unsigned int i,
                                          const unsigned int j,
                                          const Point & p,
                                          const bool add_p_level)
{
  return fe_hierarchic_3D_shape_second_deriv<HIERARCHIC>(elem, fet.order, i, j, p, add_p_level);
}



template <>
Real FE<3,L2_HIERARCHIC>::shape_second_deriv(const Elem * elem,
                                             const Order order,
                                             const unsigned int i,
                                             const unsigned int j,
                                             const Point & p,
                                             const bool add_p_level)
{
  return fe_hierarchic_3D_shape_second_deriv<L2_HIERARCHIC>(elem, order, i, j, p, add_p_level);
}


template <>
Real FE<3,L2_HIERARCHIC>::shape_second_deriv(const FEType fet,
                                             const Elem * elem,
                                             const unsigned int i,
                                             const unsigned int j,
                                             const Point & p,
                                             const bool add_p_level)
{
  return fe_hierarchic_3D_shape_second_deriv<L2_HIERARCHIC>(elem, fet.order, i, j, p, add_p_level);
}


template <>
Real FE<3,SIDE_HIERARCHIC>::shape_second_deriv(const Elem * elem,
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
    static_cast<Order>(order+add_p_level*elem->p_level());

  if (totalorder == 0) // special case since raw HIERARCHIC lacks CONSTANTs
    return 0; // constants have zero derivative

  switch (type)
    {
    case HEX27:
      {
        const unsigned int dofs_per_side = (totalorder+1u)*(totalorder+1u);
        libmesh_assert_less(i, 6*dofs_per_side);

        const unsigned int sidenum = cube_side(p);
        const unsigned int dof_offset = sidenum * dofs_per_side;

        if (i < dof_offset) // i is on a previous side
          return 0;

        if (i >= dof_offset + dofs_per_side) // i is on a later side
          return 0;

        unsigned int side_i = i - dof_offset;

        std::unique_ptr<const Elem> side = elem->build_side_ptr(sidenum);

        Point sidep = cube_side_point(sidenum, p);

        cube_remap(side_i, *side, totalorder, sidep);

        // What second derivative or mixed derivative on the side
        // corresponds to the xi/eta/zeta mix we were asked for?
        unsigned int sidej = 100;

        // Do we need a -1 here to flip the final derivative value?
        Real f = 1.;

        switch (j)
          {
          case 0: // d^2()/dxi^2
            {
              switch (sidenum)
                {
                case 0:
                  sidej = 2;
                  break;
                case 1:
                  sidej = 0;
                  break;
                case 2:
                  return 0;
                case 3:
                  sidej = 0;
                  break;
                case 4:
                  return 0;
                case 5:
                  sidej = 0;
                  break;
                default:
                  libmesh_error();
                }
              break;
            }
          case 1: // d^2()/dxideta
            {
              switch (sidenum)
                {
                case 0:
                  sidej = 1;
                  break;
                case 1:
                case 2:
                case 3:
                case 4:
                  return 0;
                case 5:
                  sidej = 1;
                  break;
                default:
                  libmesh_error();
                }
              break;
            }
          case 2: // d^2()/deta^2
            {
              switch (sidenum)
                {
                case 0:
                  sidej = 0;
                  break;
                case 1:
                  return 0;
                case 2:
                  sidej = 0;
                  break;
                case 3:
                  return 0;
                case 4:
                  sidej = 0;
                  break;
                case 5:
                  sidej = 2;
                  break;
                default:
                  libmesh_error();
                }
              break;
            }
          case 3: // d^2()/dxidzeta
            {
              switch (sidenum)
                {
                case 0:
                  return 0;
                case 1:
                  sidej = 1;
                  break;
                case 2:
                  return 0;
                case 3:
                  sidej = 1;
                  f = -1;
                  break;
                case 4:
                case 5:
                  return 0;
                default:
                  libmesh_error();
                }
              break;
            }
          case 4: // d^2()/detadzeta
            {
              switch (sidenum)
                {
                case 0:
                case 1:
                  return 0;
                case 2:
                  sidej = 1;
                  break;
                case 3:
                  return 0;
                case 4:
                  sidej = 1;
                  f = -1;
                  break;
                case 5:
                  return 0;
                default:
                  libmesh_error();
                }
              break;
            }
          case 5: // d^2()/dzeta^2
            {
              switch (sidenum)
                {
                case 0:
                  return 0;
                case 1:
                case 2:
                case 3:
                case 4:
                  sidej = 2;
                  break;
                case 5:
                  return 0;
                default:
                  libmesh_error();
                }
              break;
            }

          default:
            libmesh_error_msg("Invalid derivative index j = " << j);
          }

        return f * FE<2,HIERARCHIC>::shape_second_deriv(side.get(),
                                                        order, side_i,
                                                        sidej, sidep,
                                                        add_p_level);
      }

    default:
      libmesh_error_msg("Invalid element type = " << Utility::enum_to_string(type));
    }

#else // LIBMESH_DIM != 3
  libmesh_ignore(elem, order, i, j, p, add_p_level);
  libmesh_not_implemented();
#endif
}


template <>
Real FE<3,SIDE_HIERARCHIC>::shape_second_deriv(const FEType fet,
                                               const Elem * elem,
                                               const unsigned int i,
                                               const unsigned int j,
                                               const Point & p,
                                               const bool add_p_level)
{
  return FE<3,SIDE_HIERARCHIC>::shape_second_deriv(elem, fet.order, i, j, p, add_p_level);
}

#endif // LIBMESH_ENABLE_SECOND_DERIVATIVES

} // namespace libMesh



namespace
{
using namespace libMesh;


unsigned int cube_side (const Point & p)
{
  const Real xi = p(0), eta = p(1), zeta = p(2);
  const Real absxi   = std::abs(xi),
             abseta  = std::abs(eta),
             abszeta = std::abs(zeta);
  const Real maxabs_xi_eta   = std::max(absxi, abseta),
             maxabs_xi_zeta  = std::max(absxi, abszeta),
             maxabs_eta_zeta = std::max(abseta, abszeta);

  if (zeta < -maxabs_xi_eta)
    return 0;
  else if (eta < -maxabs_xi_zeta)
    return 1;
  else if (xi > maxabs_eta_zeta)
    return 2;
  else if (eta > maxabs_xi_zeta)
    return 3;
  else if (xi < -maxabs_eta_zeta)
    return 4;
  else if (zeta > maxabs_xi_eta)
    return 5;
  else
    libmesh_error_msg("Cannot determine side to evaluate");
}



Point cube_side_point(unsigned int sidenum, const Point & p)
{
  Point sidep;

  switch (sidenum)
    {
    case 0:
      sidep(0) = p(1);
      sidep(1) = p(0);
      break;
    case 1:
      sidep(0) = p(0);
      sidep(1) = p(2);
      break;
    case 2:
      sidep(0) = p(1);
      sidep(1) = p(2);
      break;
    case 3:
      sidep(0) = -p(0);
      sidep(1) = p(2);
      break;
    case 4:
      sidep(0) = -p(1);
      sidep(1) = p(2);
      break;
    case 5:
      sidep(0) = p(0);
      sidep(1) = p(1);
      break;
    default:
      libmesh_error();
    }

  return sidep;
}

template <FEFamily T>
Real fe_hierarchic_3D_shape(const Elem * elem,
                            const Order order,
                            const unsigned int i,
                            const Point & p,
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
      libmesh_assert (T == L2_HIERARCHIC || totalorder < 2);
      libmesh_fallthrough();
    case HEX27:
      {
        libmesh_assert_less (i, (totalorder+1u)*(totalorder+1u)*(totalorder+1u));

        // Compute hex shape functions as a tensor-product
        Real xi   = p(0);
        Real eta  = p(1);
        Real zeta = p(2);

        unsigned int i0, i1, i2;

        cube_indices(elem, totalorder, i, xi, eta, zeta, i0, i1, i2);

        return (FE<1,T>::shape(EDGE3, totalorder, i0, xi)*
                FE<1,T>::shape(EDGE3, totalorder, i1, eta)*
                FE<1,T>::shape(EDGE3, totalorder, i2, zeta));
      }

    case TET4:
    case TET10:
    case TET14:
      {
        const Real zeta[4] = { 1. - p(0) - p(1) - p(2), p(0), p(1), p(2) };

        // Nodal DoFs
        if (i < 4)
          return zeta[i];

        // Edge DoFs
        else if (i < 6u*totalorder - 2u)
          {
            const unsigned int edge_num = (i - 4) / (totalorder - 1u);
            // const int edge_node = edge_num + 4;
            const unsigned int basisorder = i - 2 - ((totalorder - 1u) * edge_num);

            const unsigned int edgevertex0 = Tet4::edge_nodes_map[edge_num][0],
                               edgevertex1 = Tet4::edge_nodes_map[edge_num][1];

            // Get factors to account for edge-flipping
            Real flip = 1;
            if (basisorder%2 &&
                (elem->point(edgevertex0) >
                 elem->point(edgevertex1)))
              flip = -1;

            const Real crossval = zeta[edgevertex0] + zeta[edgevertex1];
            const Real edgenumerator = zeta[edgevertex1] - zeta[edgevertex0];

            if (crossval == 0.) // Yes, exact comparison; we seem numerically stable otherwise
              {
                unsigned int basisfactorial = 1.;
                for (unsigned int n=2; n <= basisorder; ++n)
                  basisfactorial *= n;

                return std::pow(edgenumerator, basisorder) / basisfactorial;
              }

            const Real edgeval = edgenumerator / crossval;
            const Real crossfunc = std::pow(crossval, basisorder);

            return flip * crossfunc *
              FE<1,HIERARCHIC>::shape(EDGE3, totalorder,
                                      basisorder, edgeval);
          }

        // Face DoFs
        else if (i < 2u*totalorder*totalorder + 2u)
          {
            const int dofs_per_face = (totalorder - 1u) * (totalorder - 2u) / 2;
            const int face_num = (i - (6u*totalorder - 2u)) / dofs_per_face;

            // Reorient nodes to account for flipping and rotation.
            // We could try to identify indices with symmetric shape
            // functions, to skip this in those cases, if we really
            // need to optimize later.
            unsigned int facevertex0 = Tet4::side_nodes_map[face_num][0],
                         facevertex1 = Tet4::side_nodes_map[face_num][1],
                         facevertex2 = Tet4::side_nodes_map[face_num][2];

            // With only 3 items, we should bubble sort!
            // Programming-for-MechE's class pays off!
            bool lastcheck = true;
            if (elem->point(facevertex0) > elem->point(facevertex1))
              {
                std::swap(facevertex0, facevertex1);
                lastcheck = true;
              }
            if (elem->point(facevertex1) > elem->point(facevertex2))
              std::swap(facevertex1, facevertex2);
            if (lastcheck && elem->point(facevertex0) > elem->point(facevertex1))
              std::swap(facevertex0, facevertex1);

            const Real zeta0 = zeta[facevertex0],
                       zeta1 = zeta[facevertex1],
                       zeta2 = zeta[facevertex2];

            const unsigned int basisnum =
              i - 4 -
              (totalorder - 1u) * /*n_edges*/6 -
              (dofs_per_face * face_num);

            const unsigned int exp0 = triangular_number_column[basisnum] + 1;
            const unsigned int exp1 = triangular_number_row[basisnum] + 1 -
              triangular_number_column[basisnum];

            Real returnval = 1;
            for (unsigned int n = 0; n != exp0; ++n)
              returnval *= zeta0;
            for (unsigned int n = 0; n != exp1; ++n)
              returnval *= zeta1;
            returnval *= zeta2;
            return returnval;
          }

        // Interior DoFs
        else
          {
            const unsigned int basisnum = i - 2u*totalorder*totalorder - 2u;
            const unsigned int exp0 = tetrahedral_number_column[basisnum] + 1;
            const unsigned int exp1 = tetrahedral_number_row[basisnum] + 1 -
                                      tetrahedral_number_column[basisnum] -
                                      tetrahedral_number_page[basisnum];
            const unsigned int exp2 = tetrahedral_number_page[basisnum] + 1;

            Real returnval = 1;
            for (unsigned int n = 0; n != exp0; ++n)
              returnval *= zeta[0];
            for (unsigned int n = 0; n != exp1; ++n)
              returnval *= zeta[1];
            for (unsigned int n = 0; n != exp2; ++n)
              returnval *= zeta[2];
            returnval *= zeta[3];
            return returnval;
          }

        libmesh_error();
      }

    default:
      libmesh_error_msg("Invalid element type = " << Utility::enum_to_string(type));
    }

#else // LIBMESH_DIM != 3
  libmesh_ignore(elem, order, i, p, add_p_level);
  libmesh_not_implemented();
#endif
}



template <FEFamily T>
Real fe_hierarchic_3D_shape_deriv(const Elem * elem,
                                  const Order order,
                                  const unsigned int i,
                                  const unsigned int j,
                                  const Point & p,
                                  const bool add_p_level)
{
  return fe_fdm_deriv(elem, order, i, j, p, add_p_level, FE<3,T>::shape);
}



#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

template <FEFamily T>
Real fe_hierarchic_3D_shape_second_deriv(const Elem * elem,
                                         const Order order,
                                         const unsigned int i,
                                         const unsigned int j,
                                         const Point & p,
                                         const bool add_p_level)
{
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

  return (FE<3,T>::shape_deriv(elem, order, i, prevj, pp, add_p_level) -
          FE<3,T>::shape_deriv(elem, order, i, prevj, pm, add_p_level))
    / 2. / eps;
}

#endif // LIBMESH_ENABLE_SECOND_DERIVATIVES


} // anonymous namespace
