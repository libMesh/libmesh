// The libMesh Finite Element Library.
// Copyright (C) 2002-2026 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/quadrature_grid.h"
#include "libmesh/enum_to_string.h"

namespace libMesh
{


void QGrid::init_3D()
{
#if LIBMESH_DIM == 3

  switch (_type)
    {
      //---------------------------------------------
      // Hex quadrature rules
    case HEX8:
    case HEX20:
    case HEX27:
      {
        // We compute the 3D quadrature rule as a tensor
        // product of the 1D quadrature rule.
        //
        // We ignore p_level in QGrid, since the "order" here is a
        // user-requested number of grid points, nothing to do with
        // the order of exactly-integrated polynomial spaces
        QGrid q1D(1,_order);

        tensor_product_hex( q1D );

        return;
      }



      //---------------------------------------------
      // Tetrahedral quadrature rules
    case TET4:
    case TET10:
    case TET14:
      {
        // This is the composite centroid rule on a uniform refinement
        // of the reference tet with order+1 intervals per edge: one
        // point at the centroid of each of the (order+1)^3 subtets.
        // The refinement consists of upright tets, inverted tets,
        // and octahedra, and we split each octahedron into four tets
        // about the same diagonal.  Every subtet has the same volume
        // regardless of which diagonal is chosen, so the weights are
        // equal.  The rule is exact for linears and has O(dx^2) error
        // on smooth functions, like the midpoint grid on hexes.
        const unsigned int np = (_order+1)*(_order+1)*(_order+1);
        // Master tet has 1x1 triangle base, height 1, so volume = 1/6
        const Real weight = Real(1)/Real(6)/np;
        const Real dx = Real(1)/static_cast<int>(_order+1);
        _points.resize(np);
        _weights.resize(np);

        unsigned int pt = 0;
        auto add_point = [this, &pt, weight, dx](Real x, Real y, Real z)
        {
          _points[pt] = Point(x*dx, y*dx, z*dx);
          _weights[pt] = weight;
          pt++;
        };

        for (int i = 0; i != _order + 1; ++i)
          {
            for (int j = 0; j != _order + 1 - i; ++j)
              {
                for (int k = 0; k != _order + 1 - i - j; ++k)
                  {
                    // Centroid of the upright tet with vertices
                    // (i,j,k), (i+1,j,k), (i,j+1,k), (i,j,k+1)
                    add_point(i+Real(1)/4, j+Real(1)/4, k+Real(1)/4);

                    // The octahedron with vertices (i+1,j,k),
                    // (i,j+1,k), (i,j,k+1), (i+1,j+1,k),
                    // (i+1,j,k+1), (i,j+1,k+1) lies inside the
                    // reference tet only if i+j+k < order.  We split
                    // it into four tets about the diagonal from
                    // (i+1,j,k) to (i,j+1,k+1); their centroids are
                    // offset from the octahedron center by
                    // (-1/4,0,0), (1/4,0,0), (0,1/4,-1/4), and
                    // (0,-1/4,1/4).
                    if (i + j + k < _order)
                      {
                        const Real ci = i+Real(1)/2;
                        const Real cj = j+Real(1)/2;
                        const Real ck = k+Real(1)/2;
                        add_point(ci-Real(1)/4, cj, ck);
                        add_point(ci+Real(1)/4, cj, ck);
                        add_point(ci, cj+Real(1)/4, ck-Real(1)/4);
                        add_point(ci, cj-Real(1)/4, ck+Real(1)/4);
                      }

                    // Centroid of the inverted tet with vertices
                    // (i+1,j+1,k), (i+1,j,k+1), (i,j+1,k+1),
                    // (i+1,j+1,k+1), which lies inside the reference
                    // tet only if i+j+k < order-1
                    if (i + j + k + 1 < _order)
                      add_point(i+Real(3)/4, j+Real(3)/4, k+Real(3)/4);
                  }
              }
          }
        libmesh_assert_equal_to(pt, np);
        return;
      }


      // Prism quadrature rules
    case PRISM6:
    case PRISM15:
    case PRISM18:
    case PRISM20:
    case PRISM21:
      {
        // We compute the 3D quadrature rule as a tensor
        // product of the 1D quadrature rule and a 2D
        // triangle quadrature rule
        //
        // We ignore p_level in QGrid, since the "order" here is a
        // user-requested number of grid points, nothing to do with
        // the order of exactly-integrated polynomial spaces
        QGrid q1D(1,_order);
        QGrid q2D(2,_order);

        // Initialize the 2D rule (1D is pre-initialized)
        q2D.init(TRI3, 0, /*simple_type_only=*/true);

        tensor_product_prism(q1D, q2D);

        return;
      }



      //---------------------------------------------
      // Pyramid
    case PYRAMID5:
    case PYRAMID13:
    case PYRAMID14:
    case PYRAMID18:
      {
        // This is the composite centroid rule on a uniform refinement
        // of the reference pyramid with order+1 intervals per edge.
        // Layer k of the refinement, between z = k*dz and z =
        // (k+1)*dz, has an m x m grid of cells on its bottom face and
        // an (m-1) x (m-1) grid on its top face, where m = order+1-k.
        // It consists of
        // * m^2 upright pyramids, one on each bottom cell,
        // * (m-1)^2 inverted pyramids, one under each top cell, and
        // * 2m(m-1) tets, one per interior edge of the bottom grid,
        //   between the two upright pyramids sharing that edge.
        // Each pyramid has twice the volume of each tet, so we split
        // every pyramid into two tets along a base diagonal, giving
        // 2(order+1)^3 tets of equal volume.  We place one point at
        // the centroid of each tet.  The rule is exact for linears
        // and has O(dz^2) error on smooth functions.
        const unsigned int np = 2*(_order+1)*(_order+1)*(_order+1);
        _points.resize(np);
        _weights.resize(np);
        // Master pyramid has 2x2 base, height 1, so volume = 4/3
        const Real weight = Real(4)/Real(3)/np;
        const Real dz = Real(1)/static_cast<int>(_order+1);

        // Lattice point (i,j) on the bottom of layer k is at
        // ((2i+k)*dz - 1, (2j+k)*dz - 1, k*dz).  This map is affine,
        // so each centroid below is computed as the average of its
        // tet's vertices in (i,j,k) lattice coordinates.
        unsigned int pt = 0;
        auto add_point = [this, &pt, weight, dz](Real i, Real j, Real k)
        {
          _points[pt] = Point((2*i+k)*dz - 1, (2*j+k)*dz - 1, k*dz);
          _weights[pt] = weight;
          pt++;
        };

        for (int k = 0; k != _order + 1; ++k)
          {
            const int m = _order + 1 - k;
            for (int i = 0; i != m; ++i)
              {
                for (int j = 0; j != m; ++j)
                  {
                    // The upright pyramid with base vertices
                    // (i,j,k), (i+1,j,k), (i+1,j+1,k), (i,j+1,k) and
                    // apex (i,j,k+1), split along the diagonal from
                    // (i,j,k) to (i+1,j+1,k)
                    add_point(i+Real(1)/2, j+Real(1)/4, k+Real(1)/4);
                    add_point(i+Real(1)/4, j+Real(1)/2, k+Real(1)/4);

                    // The inverted pyramid with base vertices
                    // (i,j,k+1), (i+1,j,k+1), (i+1,j+1,k+1),
                    // (i,j+1,k+1) and apex (i+1,j+1,k), split along
                    // the diagonal from (i,j,k+1) to (i+1,j+1,k+1)
                    if (i + 1 < m && j + 1 < m)
                      {
                        add_point(i+Real(3)/4, j+Real(1)/2, k+Real(3)/4);
                        add_point(i+Real(1)/2, j+Real(3)/4, k+Real(3)/4);
                      }

                    // The tet with vertices (i+1,j,k), (i+1,j+1,k),
                    // (i,j,k+1), (i+1,j,k+1), between the upright
                    // pyramids on cells (i,j) and (i+1,j)
                    if (i + 1 < m)
                      add_point(i+Real(3)/4, j+Real(1)/4, k+Real(1)/2);

                    // The tet with vertices (i,j+1,k), (i+1,j+1,k),
                    // (i,j,k+1), (i,j+1,k+1), between the upright
                    // pyramids on cells (i,j) and (i,j+1)
                    if (j + 1 < m)
                      add_point(i+Real(1)/4, j+Real(3)/4, k+Real(1)/2);
                  }
              }
          }
        libmesh_assert_equal_to(pt, np);
        return;
      }



      //---------------------------------------------
      // Unsupported type
    default:
      libmesh_error_msg("ERROR: Unsupported type: " << Utility::enum_to_string(_type));
    }
#endif
}

} // namespace libMesh
