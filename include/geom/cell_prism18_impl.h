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

#ifndef LIBMESH_CELL_PRISM18_IMPL_H
#define LIBMESH_CELL_PRISM18_IMPL_H

// Local includes
#include "libmesh/side.h"
#include "libmesh/cell_prism18.h"
#include "libmesh/edge_edge3.h"
#include "libmesh/face_quad9.h"
#include "libmesh/face_tri6.h"
#include "libmesh/enum_io_package.h"
#include "libmesh/enum_order.h"
#include "libmesh/int_range.h"

namespace libMesh
{

// ------------------------------------------------------------
// Prism18 class member functions

template <typename RealType>
bool Prism18Templ<RealType>::is_vertex(const unsigned int i) const
{
  if (i < 6)
    return true;
  return false;
}

template <typename RealType>
bool Prism18Templ<RealType>::is_edge(const unsigned int i) const
{
  if (i < 6)
    return false;
  if (i > 14)
    return false;
  return true;
}

template <typename RealType>
bool Prism18Templ<RealType>::is_face(const unsigned int i) const
{
  if (i > 14)
    return true;
  return false;
}

template <typename RealType>
bool Prism18Templ<RealType>::is_node_on_side(const unsigned int n,
                              const unsigned int s) const
{
  libmesh_assert_less (s, this->n_sides());
  return std::find(std::begin(side_nodes_map[s]),
                   std::end(side_nodes_map[s]),
                   n) != std::end(side_nodes_map[s]);
}

template <typename RealType>
std::vector<unsigned>
Prism18Templ<RealType>::nodes_on_side(const unsigned int s) const
{
  libmesh_assert_less(s, this->n_sides());
  auto trim = (s > 0 && s < 4) ? 0 : 3;
  return {std::begin(side_nodes_map[s]), std::end(side_nodes_map[s]) - trim};
}

template <typename RealType>
bool Prism18Templ<RealType>::is_node_on_edge(const unsigned int n,
                              const unsigned int e) const
{
  libmesh_assert_less (e, this->n_edges());
  return std::find(std::begin(edge_nodes_map[e]),
                   std::end(edge_nodes_map[e]),
                   n) != std::end(edge_nodes_map[e]);
}



template <typename RealType>
bool Prism18Templ<RealType>::has_affine_map() const
{
  // Make sure z edges are affine
  Point v = this->point(3) - this->point(0);
  if (!v.relative_fuzzy_equals(this->point(4) - this->point(1)) ||
      !v.relative_fuzzy_equals(this->point(5) - this->point(2)))
    return false;
  // Make sure edges are straight
  v /= 2;
  if (!v.relative_fuzzy_equals(this->point(9) - this->point(0)) ||
      !v.relative_fuzzy_equals(this->point(10) - this->point(1)) ||
      !v.relative_fuzzy_equals(this->point(11) - this->point(2)) ||
      !v.relative_fuzzy_equals(this->point(15) - this->point(6)) ||
      !v.relative_fuzzy_equals(this->point(16) - this->point(7)) ||
      !v.relative_fuzzy_equals(this->point(17) - this->point(8)))
    return false;
  v = (this->point(1) - this->point(0))/2;
  if (!v.relative_fuzzy_equals(this->point(6) - this->point(0)) ||
      !v.relative_fuzzy_equals(this->point(12) - this->point(3)))
    return false;
  v = (this->point(2) - this->point(0))/2;
  if (!v.relative_fuzzy_equals(this->point(8) - this->point(0)) ||
      !v.relative_fuzzy_equals(this->point(14) - this->point(3)))
    return false;
  v = (this->point(2) - this->point(1))/2;
  if (!v.relative_fuzzy_equals(this->point(7) - this->point(1)) ||
      !v.relative_fuzzy_equals(this->point(13) - this->point(4)))
    return false;
  return true;
}



template <typename RealType>
Order Prism18Templ<RealType>::default_order() const
{
  return SECOND;
}

template <typename RealType>
dof_id_type Prism18Templ<RealType>::key (const unsigned int s) const
{
  libmesh_assert_less (s, this->n_sides());

  switch (s)
    {
    case 0:  // the triangular face at z=0
      {
        return Prism::key(0);
      }
    case 1:  // the quad face at y=0
      {
        return Elem::compute_key (this->node_id(15));
      }
    case 2:  // the other quad face
      {
        return Elem::compute_key (this->node_id(16));
      }
    case 3: // the quad face at x=0
      {
        return Elem::compute_key (this->node_id(17));
      }
    case 4: // the triangular face at z=1
      {
        return Prism::key(4);
      }
    default:
      libmesh_error_msg("Invalid side " << s);
    }
}



template <typename RealType>
unsigned int Prism18Templ<RealType>::which_node_am_i(unsigned int side,
                                      unsigned int side_node) const
{
  libmesh_assert_less (side, this->n_sides());

  // Never more than 9 nodes per side.
  libmesh_assert_less(side_node, 9);

  // Some sides have 6 nodes.
  libmesh_assert(!(side==0 || side==4) || side_node < 6);

  return Prism18::side_nodes_map[side][side_node];
}



template <typename RealType>
std::unique_ptr<ElemTempl<RealType>> Prism18Templ<RealType>::build_side_ptr (const unsigned int i,
                                               bool proxy)
{
  libmesh_assert_less (i, this->n_sides());

  if (proxy)
    {
      switch(i)
        {
        case 0:
        case 4:
          return libmesh_make_unique<Side<Tri6Templ<Real>,Prism18>>(this,i);

        case 1:
        case 2:
        case 3:
          return libmesh_make_unique<Side<Quad9Templ<Real>,Prism18>>(this,i);

        default:
          libmesh_error_msg("Invalid side i = " << i);
        }
    }

  else
    {
      // Return value
      std::unique_ptr<Elem> face;

      switch (i)
        {
        case 0: // the triangular face at z=-1
        case 4: // the triangular face at z=1
          {
            face = libmesh_make_unique<Tri6Templ<Real>>();
            break;
          }
        case 1: // the quad face at y=0
        case 2: // the other quad face
        case 3: // the quad face at x=0
          {
            face = libmesh_make_unique<Quad9Templ<Real>>();
            break;
          }
        default:
          libmesh_error_msg("Invalid side i = " << i);
        }

      face->subdomain_id() = this->subdomain_id();

      // Set the nodes
      for (auto n : face->node_index_range())
        face->set_node(n) = this->node_ptr(Prism18::side_nodes_map[i][n]);

      return face;
    }
}



template <typename RealType>
void Prism18Templ<RealType>::build_side_ptr (std::unique_ptr<Elem> & side,
                              const unsigned int i)
{
  libmesh_assert_less (i, this->n_sides());

  switch (i)
    {
    case 0: // the triangular face at z=-1
    case 4: // the triangular face at z=1
      {
        if (!side.get() || side->type() != TRI6)
          {
            side = this->build_side_ptr(i, false);
            return;
          }
        break;
      }

    case 1: // the quad face at y=0
    case 2: // the other quad face
    case 3: // the quad face at x=0
      {
        if (!side.get() || side->type() != QUAD9)
          {
            side = this->build_side_ptr(i, false);
            return;
          }
        break;
      }

    default:
      libmesh_error_msg("Invalid side i = " << i);
    }

  side->subdomain_id() = this->subdomain_id();

  // Set the nodes
  for (auto n : side->node_index_range())
    side->set_node(n) = this->node_ptr(Prism18::side_nodes_map[i][n]);
}



template <typename RealType>
std::unique_ptr<ElemTempl<RealType>> Prism18Templ<RealType>::build_edge_ptr (const unsigned int i)
{
  libmesh_assert_less (i, this->n_edges());

  return libmesh_make_unique<SideEdge<Edge3Templ<Real>,Prism18>>(this,i);
}



template <typename RealType>
void Prism18Templ<RealType>::connectivity(const unsigned int sc,
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
        conn.resize(8);
        switch (sc)
          {

          case 0:
            {
              conn[0] = this->node_id(0)+1;
              conn[1] = this->node_id(6)+1;
              conn[2] = this->node_id(8)+1;
              conn[3] = this->node_id(8)+1;
              conn[4] = this->node_id(9)+1;
              conn[5] = this->node_id(15)+1;
              conn[6] = this->node_id(17)+1;
              conn[7] = this->node_id(17)+1;

              return;
            }

          case 1:
            {
              conn[0] = this->node_id(6)+1;
              conn[1] = this->node_id(1)+1;
              conn[2] = this->node_id(7)+1;
              conn[3] = this->node_id(7)+1;
              conn[4] = this->node_id(15)+1;
              conn[5] = this->node_id(10)+1;
              conn[6] = this->node_id(16)+1;
              conn[7] = this->node_id(16)+1;

              return;
            }

          case 2:
            {
              conn[0] = this->node_id(8)+1;
              conn[1] = this->node_id(7)+1;
              conn[2] = this->node_id(2)+1;
              conn[3] = this->node_id(2)+1;
              conn[4] = this->node_id(17)+1;
              conn[5] = this->node_id(16)+1;
              conn[6] = this->node_id(11)+1;
              conn[7] = this->node_id(11)+1;

              return;
            }

          case 3:
            {
              conn[0] = this->node_id(6)+1;
              conn[1] = this->node_id(7)+1;
              conn[2] = this->node_id(8)+1;
              conn[3] = this->node_id(8)+1;
              conn[4] = this->node_id(15)+1;
              conn[5] = this->node_id(16)+1;
              conn[6] = this->node_id(17)+1;
              conn[7] = this->node_id(17)+1;

              return;
            }

          case 4:
            {
              conn[0] = this->node_id(9)+1;
              conn[1] = this->node_id(15)+1;
              conn[2] = this->node_id(17)+1;
              conn[3] = this->node_id(17)+1;
              conn[4] = this->node_id(3)+1;
              conn[5] = this->node_id(12)+1;
              conn[6] = this->node_id(14)+1;
              conn[7] = this->node_id(14)+1;

              return;
            }

          case 5:
            {
              conn[0] = this->node_id(15)+1;
              conn[1] = this->node_id(10)+1;
              conn[2] = this->node_id(16)+1;
              conn[3] = this->node_id(16)+1;
              conn[4] = this->node_id(12)+1;
              conn[5] = this->node_id(4)+1;
              conn[6] = this->node_id(13)+1;
              conn[7] = this->node_id(13)+1;

              return;
            }

          case 6:
            {
              conn[0] = this->node_id(17)+1;
              conn[1] = this->node_id(16)+1;
              conn[2] = this->node_id(11)+1;
              conn[3] = this->node_id(11)+1;
              conn[4] = this->node_id(14)+1;
              conn[5] = this->node_id(13)+1;
              conn[6] = this->node_id(5)+1;
              conn[7] = this->node_id(5)+1;

              return;
            }

          case 7:
            {
              conn[0] = this->node_id(15)+1;
              conn[1] = this->node_id(16)+1;
              conn[2] = this->node_id(17)+1;
              conn[3] = this->node_id(17)+1;
              conn[4] = this->node_id(12)+1;
              conn[5] = this->node_id(13)+1;
              conn[6] = this->node_id(14)+1;
              conn[7] = this->node_id(14)+1;

              return;
            }

          default:
            libmesh_error_msg("Invalid sc = " << sc);
          }

      }

    case VTK:
      {
        // VTK now supports VTK_BIQUADRATIC_QUADRATIC_WEDGE directly
        const unsigned int conn_size = 18;
        conn.resize(conn_size);

        // VTK's VTK_BIQUADRATIC_QUADRATIC_WEDGE first 9 (vertex) and
        // last 3 (mid-face) nodes match.  The middle and top layers
        // of mid-edge nodes are reversed from LibMesh's.
        for (auto i : index_range(conn))
          conn[i] = this->node_id(i);

        // top "ring" of mid-edge nodes
        conn[9]  = this->node_id(12);
        conn[10] = this->node_id(13);
        conn[11] = this->node_id(14);

        // middle "ring" of mid-edge nodes
        conn[12] = this->node_id(9);
        conn[13] = this->node_id(10);
        conn[14] = this->node_id(11);

        return;
      }

    default:
      libmesh_error_msg("Unsupported IO package " << iop);
    }
}




template <typename RealType>
unsigned int Prism18Templ<RealType>::n_second_order_adjacent_vertices (const unsigned int n) const
{
  switch (n)
    {
    case 6:
    case 7:
    case 8:
    case 9:
    case 10:
    case 11:
    case 12:
    case 13:
    case 14:
      return 2;

    case 15:
    case 16:
    case 17:
      return 4;

    default:
      libmesh_error_msg("Invalid node n = " << n);
    }
}





template <typename RealType>
unsigned short int Prism18Templ<RealType>::second_order_adjacent_vertex (const unsigned int n,
                                                          const unsigned int v) const
{
  libmesh_assert_greater_equal (n, this->n_vertices());
  libmesh_assert_less (n, this->n_nodes());

  switch (n)
    {
      /*
       * These nodes are unique to \p Prism18,
       * let our _remaining_... matrix handle
       * this.
       */
    case 15:
    case 16:
    case 17:
      {
        libmesh_assert_less (v, 4);
        return _remaining_second_order_adjacent_vertices[n-15][v];
      }

      /*
       * All other second-order nodes (6,...,14) are
       * identical with Prism15 and are therefore
       * delegated to the _second_order matrix of
       * \p Prism
       */
    default:
      {
        libmesh_assert_less (v, 2);
        return this->_second_order_adjacent_vertices[n-this->n_vertices()][v];
      }

    }

  libmesh_error_msg("We'll never ge here!");
  return static_cast<unsigned short int>(-1);
}



template <typename RealType>
std::pair<unsigned short int, unsigned short int>
Prism18Templ<RealType>::second_order_child_vertex (const unsigned int n) const
{
  libmesh_assert_greater_equal (n, this->n_vertices());
  libmesh_assert_less (n, this->n_nodes());

  return std::pair<unsigned short int, unsigned short int>
    (this->_second_order_vertex_child_number[n],
     this->_second_order_vertex_child_index[n]);
}



template <typename RealType>
RealType Prism18Templ<RealType>::volume () const
{
  // This specialization is good for Lagrange mappings only
  if (this->mapping_type() != LAGRANGE_MAP)
    return this->Elem::volume();

  // Make copies of our points.  It makes the subsequent calculations a bit
  // shorter and avoids dereferencing the same pointer multiple times.
  Point
    x0 = this->point(0),   x1 = this->point(1),   x2 = this->point(2),   x3 = this->point(3),   x4 = this->point(4), x5 = this->point(5),
    x6 = this->point(6),   x7 = this->point(7),   x8 = this->point(8),   x9 = this->point(9),  x10 = this->point(10), x11 = this->point(11),
    x12 = this->point(12), x13 = this->point(13), x14 = this->point(14), x15 = this->point(15), x16 = this->point(16), x17 = this->point(17);

  // The number of components in the dx_dxi, dx_deta, and dx_dzeta arrays.
  const int n_components = 16;

  // Terms are copied directly from a Python script.
  Point dx_dxi[n_components] =
    {
      -x10 + 4*x15 - 3*x9,
      3*x0/2 + x1/2 + 2*x12 - 3*x3/2 - x4/2 - 2*x6,
      -3*x0/2 - x1/2 + x10 + 2*x12 - 4*x15 - 3*x3/2 - x4/2 + 2*x6 + 3*x9,
      -4*x15 + 4*x16 - 4*x17 + 4*x9,
      -2*x0 - 2*x12 + 2*x13 - 2*x14 + 2*x3 + 2*x6 - 2*x7 + 2*x8,
      2*x0 - 2*x12 + 2*x13 - 2*x14 + 4*x15 - 4*x16 + 4*x17 + 2*x3 - 2*x6 + 2*x7 - 2*x8 - 4*x9,
      Point(0,0,0),
      Point(0,0,0),
      4*x10 - 8*x15 + 4*x9,
      -2*x0 - 2*x1 - 4*x12 + 2*x3 + 2*x4 + 4*x6,
      2*x0 + 2*x1 - 4*x10 - 4*x12 + 8*x15 + 2*x3 + 2*x4 - 4*x6 - 4*x9,
      Point(0,0,0),
      Point(0,0,0),
      Point(0,0,0),
      Point(0,0,0),
      Point(0,0,0)
    };

  Point dx_deta[n_components] =
    {
      -x11 + 4*x17 - 3*x9,
      3*x0/2 + 2*x14 + x2/2 - 3*x3/2 - x5/2 - 2*x8,
      -3*x0/2 + x11 + 2*x14 - 4*x17 - x2/2 - 3*x3/2 - x5/2 + 2*x8 + 3*x9,
      4*x11 - 8*x17 + 4*x9,
      -2*x0 - 4*x14 - 2*x2 + 2*x3 + 2*x5 + 4*x8,
      2*x0 - 4*x11 - 4*x14 + 8*x17 + 2*x2 + 2*x3 + 2*x5 - 4*x8 - 4*x9,
      Point(0,0,0),
      Point(0,0,0),
      -4*x15 + 4*x16 - 4*x17 + 4*x9,
      -2*x0 - 2*x12 + 2*x13 - 2*x14 + 2*x3 + 2*x6 - 2*x7 + 2*x8,
      2*x0 - 2*x12 + 2*x13 - 2*x14 + 4*x15 - 4*x16 + 4*x17 + 2*x3 - 2*x6 + 2*x7 - 2*x8 - 4*x9,
      Point(0,0,0),
      Point(0,0,0),
      Point(0,0,0),
      Point(0,0,0),
      Point(0,0,0)
    };

  Point dx_dzeta[n_components] =
    {
      -x0/2 + x3/2,
      x0 + x3 - 2*x9,
      Point(0,0,0),
      3*x0/2 + 2*x14 + x2/2 - 3*x3/2 - x5/2 - 2*x8,
      -3*x0 + 2*x11 + 4*x14 - 8*x17 - x2 - 3*x3 - x5 + 4*x8 + 6*x9,
      Point(0,0,0),
      -x0 - 2*x14 - x2 + x3 + x5 + 2*x8,
      2*x0 - 4*x11 - 4*x14 + 8*x17 + 2*x2 + 2*x3 + 2*x5 - 4*x8 - 4*x9,
      3*x0/2 + x1/2 + 2*x12 - 3*x3/2 - x4/2 - 2*x6,
      -3*x0 - x1 + 2*x10 + 4*x12 - 8*x15 - 3*x3 - x4 + 4*x6 + 6*x9,
      Point(0,0,0),
      -2*x0 - 2*x12 + 2*x13 - 2*x14 + 2*x3 + 2*x6 - 2*x7 + 2*x8,
      4*x0 - 4*x12 + 4*x13 - 4*x14 + 8*x15 - 8*x16 + 8*x17 + 4*x3 - 4*x6 + 4*x7 - 4*x8 - 8*x9,
      Point(0,0,0),
      -x0 - x1 - 2*x12 + x3 + x4 + 2*x6,
      2*x0 + 2*x1 - 4*x10 - 4*x12 + 8*x15 + 2*x3 + 2*x4 - 4*x6 - 4*x9
    };

  // The quadrature rule for the Prism18 is a tensor product between a
  // FOURTH-order TRI rule (in xi, eta) and a FIFTH-order EDGE rule
  // in zeta.

  // Number of points in the 2D quadrature rule.
  const int N2D = 6;

  // Parameters of the 2D rule
  static const Real
    w1 = Real(1.1169079483900573284750350421656140e-01L),
    w2 = Real(5.4975871827660933819163162450105264e-02L),
    a1 = Real(4.4594849091596488631832925388305199e-01L),
    a2 = Real(9.1576213509770743459571463402201508e-02L);

  // Points and weights of the 2D rule
  static const Real w2D[N2D] = {w1, w1, w1, w2, w2, w2};

  // Quadrature point locations raised to powers.  xi[0][2] is
  // quadrature point 0, squared, xi[1][1] is quadrature point 1 to the
  // first power, etc.  This lets us avoid calling std::pow inside the
  // loops below.
  static const Real xi[N2D][3] =
    {
      // ^0   ^1      ^2
      {   1., a1,     a1*a1},
      {   1., 1-2*a1, (1-2*a1)*(1-2*a1)},
      {   1., a1,     a1*a1},
      {   1., a2,     a2*a2},
      {   1., 1-2*a2, (1-2*a2)*(1-2*a2)},
      {   1., a2,     a2*a2}
    };

  static const Real eta[N2D][3] =
    {
      // ^0   ^1      ^2
      {   1., a1,     a1*a1},
      {   1., a1,     a1*a1},
      {   1., 1-2*a1, (1-2*a1)*(1-2*a1)},
      {   1., a2,     a2*a2},
      {   1., a2,     a2*a2},
      {   1., 1-2*a2, (1-2*a2)*(1-2*a2)}
    };

  // Number of points in the 1D quadrature rule.
  const int N1D = 3;

  // Points and weights of the 1D quadrature rule.
  static const Real w1D[N1D] = {5./9, 8./9, 5./9};

  const Real zeta[N1D][3] =
    {
      //^0   ^1                 ^2
      {  1., -std::sqrt(15)/5., 15./25},
      {  1., 0.,                0.},
      {  1., std::sqrt(15)/5.,  15./25}
    };

  // The integer exponents for each term.
  static const int exponents[n_components][3] =
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
      {2, 0, 1}
    };

  RealType vol = 0.;
  for (int i=0; i<N2D; ++i)
    for (int j=0; j<N1D; ++j)
      {
        // Compute dx_dxi, dx_deta, dx_dzeta at the current quadrature point.
        Point dx_dxi_q, dx_deta_q, dx_dzeta_q;
        for (int c=0; c<n_components; ++c)
          {
            Real coeff =
              xi[i][exponents[c][0]]*
              eta[i][exponents[c][1]]*
              zeta[j][exponents[c][2]];

            dx_dxi_q   += coeff * dx_dxi[c];
            dx_deta_q  += coeff * dx_deta[c];
            dx_dzeta_q += coeff * dx_dzeta[c];
          }

        // Compute scalar triple product, multiply by weight, and accumulate volume.
        vol += w2D[i] * w1D[j] * triple_product(dx_dxi_q, dx_deta_q, dx_dzeta_q);
      }

  return vol;
}


} // namespace libMesh

#endif // LIBMESH_CELL_PRISM18_IMPL_H
