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

#ifndef LIBMESH_CELL_HEX_IMPL_H
#define LIBMESH_CELL_HEX_IMPL_H

// C++ includes
#include <algorithm> // for std::min, std::max

// Local includes
#include "libmesh/cell_hex.h"
#include "libmesh/cell_hex8.h"
#include "libmesh/face_quad4.h"
#include "libmesh/enum_elem_quality.h"
#include "libmesh/tensor_value.h"

#include <array>

namespace libMesh
{

// ------------------------------------------------------------
// Hex class member functions
template <typename RealType>
dof_id_type HexTempl<RealType>::key (const unsigned int s) const
{
  libmesh_assert_less (s, this->n_sides());

  return this->compute_key(this->node_id(Hex8::side_nodes_map[s][0]),
                           this->node_id(Hex8::side_nodes_map[s][1]),
                           this->node_id(Hex8::side_nodes_map[s][2]),
                           this->node_id(Hex8::side_nodes_map[s][3]));
}



template <typename RealType>
unsigned int HexTempl<RealType>::which_node_am_i(unsigned int side,
                                  unsigned int side_node) const
{
  libmesh_assert_less (side, this->n_sides());
  libmesh_assert_less (side_node, 4);

  return Hex8::side_nodes_map[side][side_node];
}



template <typename RealType>
std::unique_ptr<ElemTempl<RealType>> HexTempl<RealType>::side_ptr (const unsigned int i)
{
  libmesh_assert_less (i, this->n_sides());

  std::unique_ptr<Elem> face = libmesh_make_unique<Quad4>();

  for (auto n : face->node_index_range())
    face->set_node(n) = this->node_ptr(Hex8::side_nodes_map[i][n]);

  return face;
}



template <typename RealType>
void HexTempl<RealType>::side_ptr (std::unique_ptr<Elem> & side,
                    const unsigned int i)
{
  this->template simple_side_ptr<Hex,Hex8>(side, i, QUAD4);
}



template <typename RealType>
bool HexTempl<RealType>::is_child_on_side(const unsigned int c,
                           const unsigned int s) const
{
  libmesh_assert_less (c, this->n_children());
  libmesh_assert_less (s, this->n_sides());

  // This array maps the Hex8 node numbering to the Hex8 child
  // numbering.  I.e.
  //   node 6 touches child 7, and
  //   node 7 touches child 6, etc.
  const unsigned int node_child_map[8] = { 0, 1, 3, 2, 4, 5, 7, 6 };

  for (unsigned int i = 0; i != 4; ++i)
    if (node_child_map[Hex8::side_nodes_map[s][i]] == c)
      return true;

  return false;
}



template <typename RealType>
bool HexTempl<RealType>::is_edge_on_side(const unsigned int e,
                          const unsigned int s) const
{
  libmesh_assert_less (e, this->n_edges());
  libmesh_assert_less (s, this->n_sides());

  return (this->is_node_on_side(Hex8::edge_nodes_map[e][0],s) &&
          this->is_node_on_side(Hex8::edge_nodes_map[e][1],s));
}



template <typename RealType>
unsigned int HexTempl<RealType>::opposite_side(const unsigned int side_in) const
{
  libmesh_assert_less (side_in, 6);
  static const unsigned char hex_opposites[6] = {5, 3, 4, 1, 2, 0};
  return hex_opposites[side_in];
}



template <typename RealType>
unsigned int HexTempl<RealType>::opposite_node(const unsigned int node_in,
                                const unsigned int side_in) const
{
  libmesh_assert_less (node_in, 26);
  libmesh_assert_less (node_in, this->n_nodes());
  libmesh_assert_less (side_in, this->n_sides());
  libmesh_assert(this->is_node_on_side(node_in, side_in));

  static const unsigned char side05_nodes_map[] =
    {4, 5, 6, 7, 0, 1, 2, 3, 16, 17, 18, 19, 255, 255, 255, 255, 8, 9, 10, 11, 25, 255, 255, 255, 255, 20};
  static const unsigned char side13_nodes_map[] =
    {3, 2, 1, 0, 7, 6, 5, 4, 10, 255, 8, 255, 15, 14, 13, 12, 18, 255, 16, 255, 255, 23, 255, 21, 255, 255};
  static const unsigned char side24_nodes_map[] =
    {1, 0, 3, 2, 5, 4, 7, 6, 255, 11, 255, 9, 13, 12, 15, 14, 255, 19, 255, 17, 255, 255, 24, 255, 22, 255};

  switch (side_in)
    {
    case 0:
    case 5:
      return side05_nodes_map[node_in];
    case 1:
    case 3:
      return side13_nodes_map[node_in];
    case 2:
    case 4:
      return side24_nodes_map[node_in];
    default:
      libmesh_error_msg("Unsupported side_in = " << side_in);
    }
}



template <typename RealType>
Real HexTempl<RealType>::quality (const ElemQuality q) const
{
  switch (q)
    {

#if LIBMESH_DIM >= 3
      /**
       * Compute the min/max diagonal ratio.
       * Source: CUBIT User's Manual.
       */
    case DIAGONAL:
      {
        // Diagonal between node 0 and node 6
        const Real d06 = this->length(0,6);

        // Diagonal between node 3 and node 5
        const Real d35 = this->length(3,5);

        // Diagonal between node 1 and node 7
        const Real d17 = this->length(1,7);

        // Diagonal between node 2 and node 4
        const Real d24 = this->length(2,4);

        // Find the biggest and smallest diagonals
        const Real min = std::min(d06, std::min(d35, std::min(d17, d24)));
        const Real max = std::max(d06, std::max(d35, std::max(d17, d24)));

        libmesh_assert_not_equal_to (max, 0.0);

        return min / max;

        break;
      }

      /**
       * Minimum ratio of lengths derived from opposite edges.
       * Source: CUBIT User's Manual.
       */
    case TAPER:
      {

        /**
         * Compute the side lengths.
         */
        const Real d01 = this->length(0,1);
        const Real d12 = this->length(1,2);
        const Real d23 = this->length(2,3);
        const Real d03 = this->length(0,3);
        const Real d45 = this->length(4,5);
        const Real d56 = this->length(5,6);
        const Real d67 = this->length(6,7);
        const Real d47 = this->length(4,7);
        const Real d04 = this->length(0,4);
        const Real d15 = this->length(1,5);
        const Real d37 = this->length(3,7);
        const Real d26 = this->length(2,6);

        std::vector<Real> edge_ratios(12);
        // Front
        edge_ratios[0] = std::min(d01, d45) / std::max(d01, d45);
        edge_ratios[1] = std::min(d04, d15) / std::max(d04, d15);

        // Right
        edge_ratios[2] = std::min(d15, d26) / std::max(d15, d26);
        edge_ratios[3] = std::min(d12, d56) / std::max(d12, d56);

        // Back
        edge_ratios[4] = std::min(d67, d23) / std::max(d67, d23);
        edge_ratios[5] = std::min(d26, d37) / std::max(d26, d37);

        // Left
        edge_ratios[6] = std::min(d04, d37) / std::max(d04, d37);
        edge_ratios[7] = std::min(d03, d47) / std::max(d03, d47);

        // Bottom
        edge_ratios[8] = std::min(d01, d23) / std::max(d01, d23);
        edge_ratios[9] = std::min(d03, d12) / std::max(d03, d12);

        // Top
        edge_ratios[10] = std::min(d45, d67) / std::max(d45, d67);
        edge_ratios[11] = std::min(d56, d47) / std::max(d56, d47);

        return *(std::min_element(edge_ratios.begin(), edge_ratios.end())) ;

        break;
      }


      /**
       * Minimum edge length divided by max diagonal length.
       * Source: CUBIT User's Manual.
       */
    case STRETCH:
      {
        const Real sqrt3 = 1.73205080756888;

        /**
         * Compute the maximum diagonal.
         */
        const Real d06 = this->length(0,6);
        const Real d17 = this->length(1,7);
        const Real d35 = this->length(3,5);
        const Real d24 = this->length(2,4);
        const Real max_diag = std::max(d06, std::max(d17, std::max(d35, d24)));

        libmesh_assert_not_equal_to ( max_diag, 0.0 );

        /**
         * Compute the minimum edge length.
         */
        std::vector<Real> edges(12);
        edges[0]  = this->length(0,1);
        edges[1]  = this->length(1,2);
        edges[2]  = this->length(2,3);
        edges[3]  = this->length(0,3);
        edges[4]  = this->length(4,5);
        edges[5]  = this->length(5,6);
        edges[6]  = this->length(6,7);
        edges[7]  = this->length(4,7);
        edges[8]  = this->length(0,4);
        edges[9]  = this->length(1,5);
        edges[10] = this->length(2,6);
        edges[11] = this->length(3,7);

        const Real min_edge = *(std::min_element(edges.begin(), edges.end()));
        return sqrt3 * min_edge / max_diag ;
      }


    case SHAPE:
    case SKEW:
      {
        // From: P. Knupp, "Algebraic mesh quality metrics for
        // unstructured initial meshes," Finite Elements in Analysis
        // and Design 39, 2003, p. 217-241, Sections 6.2 and 6.3.

        // Make local copies of points, we will access these several
        // times below.
        const Point
          x0 = this->point(0), x1 = this->point(1), x2 = this->point(2), x3 = this->point(3),
          x4 = this->point(4), x5 = this->point(5), x6 = this->point(6), x7 = this->point(7);

        // The columns of the Jacobian matrices are:
        // \vec{x}_{\xi}   = \vec{a1}*eta*zeta + \vec{b1}*eta + \vec{c1}*zeta + \vec{d1}
        // \vec{x}_{\eta}  = \vec{a2}*xi*zeta  + \vec{b2}*xi  + \vec{c2}*zeta + \vec{d2}
        // \vec{x}_{\zeta} = \vec{a3}*xi*eta   + \vec{b3}*xi  + \vec{c3}*eta  + \vec{d3}
        // where the ai, bi, ci, and di are constants defined below.
        const Point a1 = -x0 + x1 - x2 + x3 + x4 - x5 + x6 - x7;
        const Point b1 =  x0 - x1 + x2 - x3 + x4 - x5 + x6 - x7;
        const Point c1 =  x0 - x1 - x2 + x3 - x4 + x5 + x6 - x7;
        const Point d1 = -x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7;

        const Point a2 =  a1;
        const Point b2 =  b1;
        const Point c2 =  x0 + x1 - x2 - x3 - x4 - x5 + x6 + x7;
        const Point d2 = -x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7;

        const Point a3 =  a1;
        const Point b3 =  c1;
        const Point c3 =  c2;
        const Point d3 = -x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7;

        // Form the nodal Jacobians. These were computed using a
        // Python script and the formulas above. Note that we are
        // actually computing the Jacobian _columns_ and passing them
        // to the RealTensor constructor which expects _rows_, but
        // it's OK because we are only interested in determinants and
        // products which are not affected by taking the transpose.
        std::array<RealTensor, 8> A =
          {
            RealTensor(d1, d2, d3),
            RealTensor(d1, b2 + d2, b3 + d3),
            RealTensor(b1 + d1, b2 + d2, a3 + b3 + c3 + d3),
            RealTensor(b1 + d1, d2, c3 + d3),
            RealTensor(c1 + d1, c2 + d2, d3),
            RealTensor(c1 + d1, a2 + b2 + c2 + d2, b3 + d3),
            RealTensor(a1 + b1 + c1 + d1, a2 + b2 + c2 + d2, a3 + b3 + c3 + d3),
            RealTensor(a1 + b1 + c1 + d1, c2 + d2, c3 + d3)
          };

        // Compute Nodal areas, alpha_k = det(A_k).
        // If any of these are zero or negative, we return zero
        // (lowest possible value) for the quality, since the formulas
        // below require positive nodal areas.
        std::array<Real, 8> alpha;
        for (unsigned int k=0; k<alpha.size(); ++k)
          {
            alpha[k] = A[k].det();
            if (alpha[k] <= 0.)
              return 0.;
          }

        // Compute metric tensors, T_k = A_k^T * A_k.
        std::array<RealTensor, 8> T;
        for (unsigned int k=0; k<T.size(); ++k)
          T[k] = A[k] * A[k].transpose();

        // Compute and return the shape metric. These only use the
        // diagonal entries of the T_k.
        Real den = 0.;
        if (q == SHAPE)
          {
            for (unsigned int k=0; k<T.size(); ++k)
              den += T[k].tr() / std::pow(alpha[k], 2./3.);
            return (den == 0.) ? 0 : (24. / den);
          }
        else
          {
            for (unsigned int k=0; k<T.size(); ++k)
              den += std::pow(std::sqrt(T[k](0,0) * T[k](1,1) * T[k](2,2)) / alpha[k], 2./3.);
            return (den == 0.) ? 0 : (8. / den);
          }
      }
#endif // LIBMESH_DIM >= 3

      /**
       * I don't know what to do for this metric.
       * Maybe the base class knows...
       */
    default:
      return ElemTempl<RealType>::quality(q);
    }
}



template <typename RealType>
std::pair<Real, Real> HexTempl<RealType>::qual_bounds (const ElemQuality q) const
{
  std::pair<Real, Real> bounds;

  switch (q)
    {

    case ASPECT_RATIO:
      bounds.first  = 1.;
      bounds.second = 4.;
      break;

    case SKEW:
      bounds.first  = 0.;
      bounds.second = 0.5;
      break;

    case SHEAR:
    case SHAPE:
      bounds.first  = 0.3;
      bounds.second = 1.;
      break;

    case CONDITION:
      bounds.first  = 1.;
      bounds.second = 8.;
      break;

    case JACOBIAN:
      bounds.first  = 0.5;
      bounds.second = 1.;
      break;

    case DISTORTION:
      bounds.first  = 0.6;
      bounds.second = 1.;
      break;

    case TAPER:
      bounds.first  = 0.;
      bounds.second = 0.4;
      break;

    case STRETCH:
      bounds.first  = 0.25;
      bounds.second = 1.;
      break;

    case DIAGONAL:
      bounds.first  = 0.65;
      bounds.second = 1.;
      break;

    case SIZE:
      bounds.first  = 0.5;
      bounds.second = 1.;
      break;

    default:
      libMesh::out << "Warning: Invalid quality measure chosen." << std::endl;
      bounds.first  = -1;
      bounds.second = -1;
    }

  return bounds;
}

} // namespace libMesh

#endif // LIBMESH_CELL_HEX_IMPL_H
