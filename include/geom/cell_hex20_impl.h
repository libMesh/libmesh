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

#ifndef LIBMESH_CELL_HEX20_IMPL_H
#define LIBMESH_CELL_HEX20_IMPL_H

// Local includes
#include "libmesh/side.h"
#include "libmesh/cell_hex20.h"
#include "libmesh/edge_edge3.h"
#include "libmesh/face_quad8.h"
#include "libmesh/enum_io_package.h"
#include "libmesh/enum_order.h"

namespace libMesh
{
// ------------------------------------------------------------
// Hex20 class member functions

template <typename RealType>
bool Hex20Templ<RealType>::is_vertex(const unsigned int i) const
{
  if (i < 8)
    return true;
  return false;
}

template <typename RealType>
bool Hex20Templ<RealType>::is_edge(const unsigned int i) const
{
  if (i > 7)
    return true;
  return false;
}

template <typename RealType>
bool Hex20Templ<RealType>::is_face(const unsigned int) const
{
  return false;
}

template <typename RealType>
bool Hex20Templ<RealType>::is_node_on_side(const unsigned int n,
                            const unsigned int s) const
{
  libmesh_assert_less (s, this->n_sides());
  return std::find(std::begin(side_nodes_map[s]),
                   std::end(side_nodes_map[s]),
                   n) != std::end(side_nodes_map[s]);
}

template <typename RealType>
std::vector<unsigned>
Hex20Templ<RealType>::nodes_on_side(const unsigned int s) const
{
  libmesh_assert_less(s, this->n_sides());
  return {std::begin(side_nodes_map[s]), std::end(side_nodes_map[s])};
}

template <typename RealType>
bool Hex20Templ<RealType>::is_node_on_edge(const unsigned int n,
                            const unsigned int e) const
{
  libmesh_assert_less (e, this->n_edges());
  return std::find(std::begin(edge_nodes_map[e]),
                   std::end(edge_nodes_map[e]),
                   n) != std::end(edge_nodes_map[e]);
}



template <typename RealType>
bool Hex20Templ<RealType>::has_affine_map() const
{
  // Make sure x-edge endpoints are affine
  Point v = this->point(1) - this->point(0);
  if (!v.relative_fuzzy_equals(this->point(2) - this->point(3)) ||
      !v.relative_fuzzy_equals(this->point(5) - this->point(4)) ||
      !v.relative_fuzzy_equals(this->point(6) - this->point(7)))
    return false;
  // Make sure x-edges are straight
  v /= 2;
  if (!v.relative_fuzzy_equals(this->point(8) - this->point(0)) ||
      !v.relative_fuzzy_equals(this->point(10) - this->point(3)) ||
      !v.relative_fuzzy_equals(this->point(16) - this->point(4)) ||
      !v.relative_fuzzy_equals(this->point(18) - this->point(7)))
    return false;
  // Make sure xz-faces are identical parallelograms
  v = this->point(4) - this->point(0);
  if (!v.relative_fuzzy_equals(this->point(7) - this->point(3)))
    return false;
  v /= 2;
  if (!v.relative_fuzzy_equals(this->point(12) - this->point(0)) ||
      !v.relative_fuzzy_equals(this->point(13) - this->point(1)) ||
      !v.relative_fuzzy_equals(this->point(14) - this->point(2)) ||
      !v.relative_fuzzy_equals(this->point(15) - this->point(3)))
    return false;
  // Make sure y-edges are straight
  v = (this->point(3) - this->point(0))/2;
  if (!v.relative_fuzzy_equals(this->point(11) - this->point(0)) ||
      !v.relative_fuzzy_equals(this->point(9) - this->point(1)) ||
      !v.relative_fuzzy_equals(this->point(17) - this->point(5)) ||
      !v.relative_fuzzy_equals(this->point(19) - this->point(4)))
    return false;
  // If all the above checks out, the map is affine
  return true;
}



template <typename RealType>
Order Hex20Templ<RealType>::default_order() const
{
  return SECOND;
}



template <typename RealType>
std::unique_ptr<ElemTempl<RealType>> Hex20Templ<RealType>::build_side_ptr (const unsigned int i,
                                             bool proxy )
{
  libmesh_assert_less (i, this->n_sides());

  if (proxy)
    return libmesh_make_unique<Side<Quad8,Hex20>>(this,i);

  else
    {
      std::unique_ptr<Elem> face = libmesh_make_unique<Quad8>();
      face->subdomain_id() = this->subdomain_id();

      for (auto n : face->node_index_range())
        face->set_node(n) = this->node_ptr(Hex20::side_nodes_map[i][n]);

      return face;
    }
}



template <typename RealType>
void Hex20Templ<RealType>::build_side_ptr (std::unique_ptr<Elem> & side,
                            const unsigned int i)
{
  this->template simple_build_side_ptr<Hex20>(side, i, QUAD8);
}



template <typename RealType>
unsigned int Hex20Templ<RealType>::which_node_am_i(unsigned int side,
                                    unsigned int side_node) const
{
  libmesh_assert_less (side, this->n_sides());
  libmesh_assert_less (side_node, Hex20::nodes_per_side);

  return Hex20::side_nodes_map[side][side_node];
}



template <typename RealType>
std::unique_ptr<ElemTempl<RealType>> Hex20Templ<RealType>::build_edge_ptr (const unsigned int i)
{
  libmesh_assert_less (i, this->n_edges());

  return libmesh_make_unique<SideEdge<Edge3,Hex20>>(this,i);
}



template <typename RealType>
void Hex20Templ<RealType>::connectivity(const unsigned int sc,
                         const IOPackage iop,
                         std::vector<dof_id_type> & conn) const
{
  libmesh_assert(this->_nodes);
  libmesh_assert_less (sc, this->n_sub_elem());
  libmesh_assert_not_equal_to (iop, INVALID_IO_PACKAGE);


  switch (iop)
    {
    case TECPLOT:
      {
        switch (sc)
          {
          case 0:
            conn.resize(8);
            conn[0] = this->node_id(0)+1;
            conn[1] = this->node_id(1)+1;
            conn[2] = this->node_id(2)+1;
            conn[3] = this->node_id(3)+1;
            conn[4] = this->node_id(4)+1;
            conn[5] = this->node_id(5)+1;
            conn[6] = this->node_id(6)+1;
            conn[7] = this->node_id(7)+1;

            return;

          default:
            libmesh_error_msg("Unknown sc = " << sc);
          }
      }

    case VTK:
      {
        switch (sc)
          {
          case 0:
            conn.resize(20);
            conn[0] = this->node_id(0);
            conn[1] = this->node_id(1);
            conn[2] = this->node_id(2);
            conn[3] = this->node_id(3);
            conn[4] = this->node_id(4);
            conn[5] = this->node_id(5);
            conn[6] = this->node_id(6);
            conn[7] = this->node_id(7);
            conn[8] = this->node_id(8);
            conn[9] = this->node_id(9);
            conn[10] = this->node_id(10);
            conn[11] = this->node_id(11);
            conn[12] = this->node_id(16);
            conn[13] = this->node_id(17);
            conn[14] = this->node_id(18);
            conn[15] = this->node_id(19);
            conn[16] = this->node_id(12);
            conn[17] = this->node_id(13);
            conn[18] = this->node_id(14);
            conn[19] = this->node_id(15);
            return;

          default:
            libmesh_error_msg("Unknown sc = " << sc);
          }
      }

    default:
      libmesh_error_msg("Unsupported IO package " << iop);
    }
}




template <typename RealType>
unsigned short int Hex20Templ<RealType>::second_order_adjacent_vertex (const unsigned int n,
                                                        const unsigned int v) const
{
  libmesh_assert_greater_equal (n, this->n_vertices());
  libmesh_assert_less (n, this->n_nodes());
  libmesh_assert_less (v, 2);
  /*
   * the _second_order_adjacent_vertices matrix is
   * stored in cell_hex.C, since this matrix is identical
   * for Hex20 and Hex27 (for the first 12 higher-order nodes)
   */
  return this->_second_order_adjacent_vertices[n-this->n_vertices()][v];
}



template <typename RealType>
std::pair<unsigned short int, unsigned short int>
Hex20Templ<RealType>::second_order_child_vertex (const unsigned int n) const
{
  libmesh_assert_greater_equal (n, this->n_vertices());
  libmesh_assert_less (n, this->n_nodes());
  /*
   * the _second_order_vertex_child_* vectors are
   * stored in cell_hex.C, since they are identical
   * for Hex20 and Hex27 (for the first 12 higher-order nodes)
   */
  return std::pair<unsigned short int, unsigned short int>
    (this->_second_order_vertex_child_number[n],
     this->_second_order_vertex_child_index[n]);
}



template <typename RealType>
RealType Hex20Templ<RealType>::volume () const
{
  // This specialization is good for Lagrange mappings only
  if (this->mapping_type() != LAGRANGE_MAP)
    return this->Elem::volume();

  // Make copies of our points.  It makes the subsequent calculations a bit
  // shorter and avoids dereferencing the same pointer multiple times.
  Point
    x0 = this->point(0),   x1 = this->point(1),   x2 = this->point(2),   x3 = this->point(3),   x4 = this->point(4),
    x5 = this->point(5),   x6 = this->point(6),   x7 = this->point(7),   x8 = this->point(8),   x9 = this->point(9),
    x10 = this->point(10), x11 = this->point(11), x12 = this->point(12), x13 = this->point(13), x14 = this->point(14),
    x15 = this->point(15), x16 = this->point(16), x17 = this->point(17), x18 = this->point(18), x19 = this->point(19);

  // The constant components of the dx/dxi vector,
  // dx/dxi = \vec{a000} + \vec{a001}*zeta + \vec{a002}*zeta^2 + ...
  // These were copied directly from the output of a Python script.
  // There are at most 17 terms with total degree <=3, but only 12
  // of them are non-zero for each direction.
  Point dx_dxi[17] =
    {
      x0/8 - x1/8 - x11/4 - x12/4 + x13/4 + x14/4 - x15/4 + x17/4 - x19/4 - x2/8 + x3/8 + x4/8 - x5/8 - x6/8 + x7/8 + x9/4,
      x11/4 + x17/4 - x19/4 - x9/4,
      -x0/8 + x1/8 + x12/4 - x13/4 - x14/4 + x15/4 + x2/8 - x3/8 - x4/8 + x5/8 + x6/8 - x7/8,
      x12/4 - x13/4 + x14/4 - x15/4,
      -x0/8 + x1/8 - x2/8 + x3/8 + x4/8 - x5/8 + x6/8 - x7/8,
      x0/8 - x1/8 - x12/4 + x13/4 - x14/4 + x15/4 + x2/8 - x3/8 + x4/8 - x5/8 + x6/8 - x7/8,
      -x0/8 + x1/8 + x11/4 - x17/4 + x19/4 + x2/8 - x3/8 - x4/8 + x5/8 + x6/8 - x7/8 - x9/4,
      x0/8 - x1/8 - x11/4 - x17/4 + x19/4 - x2/8 + x3/8 - x4/8 + x5/8 + x6/8 - x7/8 + x9/4,
      x0/4 + x1/4 - x10/2 - x16/2 - x18/2 + x2/4 + x3/4 + x4/4 + x5/4 + x6/4 + x7/4 - x8/2,
      -x0/4 - x1/4 + x10/2 - x16/2 - x18/2 - x2/4 - x3/4 + x4/4 + x5/4 + x6/4 + x7/4 + x8/2,
      Point(0,0,0),
      -x0/4 - x1/4 - x10/2 + x16/2 - x18/2 + x2/4 + x3/4 - x4/4 - x5/4 + x6/4 + x7/4 + x8/2,
      x0/4 + x1/4 + x10/2 + x16/2 - x18/2 - x2/4 - x3/4 - x4/4 - x5/4 + x6/4 + x7/4 - x8/2,
      Point(0,0,0),
      Point(0,0,0),
      Point(0,0,0),
      Point(0,0,0)
    };

  // The constant components of the dx/deta vector. These were copied
  // directly from the output of a Python script.  There are at most
  // 17 terms with total degree <=3, but only 12 of them are non-zero
  // for each direction.
  Point dx_deta[17] =
    {
      x0/8 + x1/8 + x10/4 - x12/4 - x13/4 + x14/4 + x15/4 - x16/4 + x18/4 - x2/8 - x3/8 + x4/8 + x5/8 - x6/8 - x7/8 - x8/4,
      -x10/4 - x16/4 + x18/4 + x8/4,
      -x0/8 - x1/8 + x12/4 + x13/4 - x14/4 - x15/4 + x2/8 + x3/8 - x4/8 - x5/8 + x6/8 + x7/8,
      x0/4 + x1/4 - x11/2 - x17/2 - x19/2 + x2/4 + x3/4 + x4/4 + x5/4 + x6/4 + x7/4 - x9/2,
      -x0/4 - x1/4 + x11/2 - x17/2 - x19/2 - x2/4 - x3/4 + x4/4 + x5/4 + x6/4 + x7/4 + x9/2,
      Point(0,0,0),
      Point(0,0,0),
      Point(0,0,0),
      x12/4 - x13/4 + x14/4 - x15/4,
      -x0/8 + x1/8 - x2/8 + x3/8 + x4/8 - x5/8 + x6/8 - x7/8,
      x0/8 - x1/8 - x12/4 + x13/4 - x14/4 + x15/4 + x2/8 - x3/8 + x4/8 - x5/8 + x6/8 - x7/8,
      -x0/4 + x1/4 + x11/2 - x17/2 + x19/2 + x2/4 - x3/4 - x4/4 + x5/4 + x6/4 - x7/4 - x9/2,
      x0/4 - x1/4 - x11/2 - x17/2 + x19/2 - x2/4 + x3/4 - x4/4 + x5/4 + x6/4 - x7/4 + x9/2,
      Point(0,0,0),
      -x0/8 - x1/8 - x10/4 + x16/4 - x18/4 + x2/8 + x3/8 - x4/8 - x5/8 + x6/8 + x7/8 + x8/4,
      x0/8 + x1/8 + x10/4 + x16/4 - x18/4 - x2/8 - x3/8 - x4/8 - x5/8 + x6/8 + x7/8 - x8/4,
      Point(0,0,0)
    };

  // The constant components of the dx/dzeta vector. These were copied
  // directly from the output of a Python script.  There are at most
  // 17 terms with total degree <=3, but only 12 of them are non-zero
  // for each direction.
  Point dx_dzeta[17] =
    {
      x0/8 + x1/8 - x10/4 - x11/4 + x16/4 + x17/4 + x18/4 + x19/4 + x2/8 + x3/8 - x4/8 - x5/8 - x6/8 - x7/8 - x8/4 - x9/4,
      x0/4 + x1/4 - x12/2 - x13/2 - x14/2 - x15/2 + x2/4 + x3/4 + x4/4 + x5/4 + x6/4 + x7/4,
      Point(0,0,0),
      -x10/4 - x16/4 + x18/4 + x8/4,
      -x0/4 - x1/4 + x12/2 + x13/2 - x14/2 - x15/2 + x2/4 + x3/4 - x4/4 - x5/4 + x6/4 + x7/4,
      Point(0,0,0),
      -x0/8 - x1/8 + x11/4 - x17/4 - x19/4 - x2/8 - x3/8 + x4/8 + x5/8 + x6/8 + x7/8 + x9/4,
      Point(0,0,0),
      x11/4 + x17/4 - x19/4 - x9/4,
      -x0/4 + x1/4 + x12/2 - x13/2 - x14/2 + x15/2 + x2/4 - x3/4 - x4/4 + x5/4 + x6/4 - x7/4,
      Point(0,0,0),
      -x0/8 + x1/8 - x2/8 + x3/8 + x4/8 - x5/8 + x6/8 - x7/8,
      x0/4 - x1/4 - x12/2 + x13/2 - x14/2 + x15/2 + x2/4 - x3/4 + x4/4 - x5/4 + x6/4 - x7/4,
      x0/8 - x1/8 - x11/4 - x17/4 + x19/4 - x2/8 + x3/8 - x4/8 + x5/8 + x6/8 - x7/8 + x9/4,
      -x0/8 - x1/8 + x10/4 - x16/4 - x18/4 - x2/8 - x3/8 + x4/8 + x5/8 + x6/8 + x7/8 + x8/4,
      Point(0,0,0),
      x0/8 + x1/8 + x10/4 + x16/4 - x18/4 - x2/8 - x3/8 - x4/8 - x5/8 + x6/8 + x7/8 - x8/4,
    };

  // The integer exponents for each term.
  static const int exponents[17][3] =
    {
      {0, 0, 0},
      {0, 0, 1},
      {0, 0, 2},
      {0, 1, 0},
      {0, 1, 1},
      {0, 1, 2},
      {0, 2, 0},
      {0, 2, 1},
      {1, 0, 0},
      {1, 0, 1},
      {1, 0, 2},
      {1, 1, 0},
      {1, 1, 1},
      {1, 2, 0},
      {2, 0, 0},
      {2, 0, 1},
      {2, 1, 0}
    };


  // 3x3 quadrature, exact for bi-quintics
  const int N = 3;
  const Real w[N] = {5./9, 8./9, 5./9};

  // Quadrature point locations raised to powers.  q[0][2] is
  // quadrature point 0, squared, q[1][1] is quadrature point 1 to the
  // first power, etc.
  const Real q[N][N] =
    {
      //^0   ^1                 ^2
      {  1., -std::sqrt(15)/5., 15./25},
      {  1., 0.,                0.},
      {  1., std::sqrt(15)/5.,  15./25}
    };


  RealType vol = 0.;
  for (int i=0; i<N; ++i)
    for (int j=0; j<N; ++j)
      for (int k=0; k<N; ++k)
        {
          // Compute dx_dxi, dx_deta, dx_dzeta at the current quadrature point.
          Point dx_dxi_q, dx_deta_q, dx_dzeta_q;
          for (int c=0; c<17; ++c)
            {
              Real coeff =
                q[i][exponents[c][0]] *
                q[j][exponents[c][1]] *
                q[k][exponents[c][2]];

              dx_dxi_q   += coeff * dx_dxi[c];
              dx_deta_q  += coeff * dx_deta[c];
              dx_dzeta_q += coeff * dx_dzeta[c];
            }

          // Compute scalar triple product, multiply by weight, and accumulate volume.
          vol += w[i] * w[j] * w[k] * triple_product(dx_dxi_q, dx_deta_q, dx_dzeta_q);
        }

  return vol;
}

} // namespace libMesh

#endif // LIBMESH_CELL_HEX20_IMPL_H
