// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

// Local includes
#include "libmesh/side.h"
#include "libmesh/cell_prism18.h"
#include "libmesh/edge_edge3.h"
#include "libmesh/face_quad9.h"
#include "libmesh/face_tri6.h"

namespace libMesh
{



// ------------------------------------------------------------
// Prism18 class static member initializations
const unsigned int Prism18::side_nodes_map[5][9] =
  {
    {0, 2, 1,  8,  7,  6, 99, 99, 99}, // Side 0
    {0, 1, 4,  3,  6, 10, 12,  9, 15}, // Side 1
    {1, 2, 5,  4,  7, 11, 13, 10, 16}, // Side 2
    {2, 0, 3,  5,  8,  9, 14, 11, 17}, // Side 3
    {3, 4, 5, 12, 13, 14, 99, 99, 99}  // Side 4
  };

const unsigned int Prism18::edge_nodes_map[9][3] =
  {
    {0, 1, 6},  // Side 0
    {1, 2, 7},  // Side 1
    {0, 2, 8},  // Side 2
    {0, 3, 9},  // Side 3
    {1, 4, 10}, // Side 4
    {2, 5, 11}, // Side 5
    {3, 4, 12}, // Side 6
    {4, 5, 13}, // Side 7
    {3, 5, 14}  // Side 8
  };


// ------------------------------------------------------------
// Prism18 class member functions

bool Prism18::is_vertex(const unsigned int i) const
{
  if (i < 6)
    return true;
  return false;
}

bool Prism18::is_edge(const unsigned int i) const
{
  if (i < 6)
    return false;
  if (i > 14)
    return false;
  return true;
}

bool Prism18::is_face(const unsigned int i) const
{
  if (i > 14)
    return true;
  return false;
}

bool Prism18::is_node_on_side(const unsigned int n,
                              const unsigned int s) const
{
  libmesh_assert_less (s, n_sides());
  for (unsigned int i = 0; i != 9; ++i)
    if (side_nodes_map[s][i] == n)
      return true;
  return false;
}

bool Prism18::is_node_on_edge(const unsigned int n,
                              const unsigned int e) const
{
  libmesh_assert_less (e, n_edges());
  for (unsigned int i = 0; i != 3; ++i)
    if (edge_nodes_map[e][i] == n)
      return true;
  return false;
}



bool Prism18::has_affine_map() const
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



dof_id_type Prism18::key (const unsigned int s) const
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

  libmesh_error_msg("We'll never get here!");
  return 0;
}





UniquePtr<Elem> Prism18::build_side_ptr (const unsigned int i,
                                         bool proxy)
{
  libmesh_assert_less (i, this->n_sides());

  if (proxy)
    {
      switch(i)
        {
        case 0:
        case 4:
          return UniquePtr<Elem>(new Side<Tri6,Prism18>(this,i));

        case 1:
        case 2:
        case 3:
          return UniquePtr<Elem>(new Side<Quad9,Prism18>(this,i));

        default:
          libmesh_error_msg("Invalid side i = " << i);
        }
    }

  else
    {
      // Create NULL pointer to be initialized, returned later.
      Elem * face = libmesh_nullptr;

      switch (i)
        {
        case 0: // the triangular face at z=-1
        case 4: // the triangular face at z=1
          {
            face = new Tri6;
            break;
          }
        case 1: // the quad face at y=0
        case 2: // the other quad face
        case 3: // the quad face at x=0
          {
            face = new Quad9;
            break;
          }
        default:
          libmesh_error_msg("Invalid side i = " << i);
        }

      face->subdomain_id() = this->subdomain_id();

      // Set the nodes
      for (unsigned n=0; n<face->n_nodes(); ++n)
        face->set_node(n) = this->node_ptr(Prism18::side_nodes_map[i][n]);

      return UniquePtr<Elem>(face);
    }

  libmesh_error_msg("We'll never get here!");
  return UniquePtr<Elem>();
}



UniquePtr<Elem> Prism18::build_edge_ptr (const unsigned int i)
{
  libmesh_assert_less (i, this->n_edges());

  return UniquePtr<Elem>(new SideEdge<Edge3,Prism18>(this,i));
}



void Prism18::connectivity(const unsigned int sc,
                           const IOPackage iop,
                           std::vector<dof_id_type> & conn) const
{
  libmesh_assert(_nodes);
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
        conn.resize(18);

        // VTK's VTK_BIQUADRATIC_QUADRATIC_WEDGE first 9 (vertex) and
        // last 3 (mid-face) nodes match.  The middle and top layers
        // of mid-edge nodes are reversed from LibMesh's.
        for (std::size_t i=0; i<conn.size(); ++i)
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

        /*
          conn.resize(6);
          switch (sc)
          {

          case 0:
          {
          conn[0] = this->node_id(0);
          conn[1] = this->node_id(6);
          conn[2] = this->node_id(8);
          conn[3] = this->node_id(9);
          conn[4] = this->node_id(15);
          conn[5] = this->node_id(17);

          return;
          }

          case 1:
          {
          conn[0] = this->node_id(6);
          conn[1] = this->node_id(1);
          conn[2] = this->node_id(7);
          conn[3] = this->node_id(15);
          conn[4] = this->node_id(10);
          conn[5] = this->node_id(16);

          return;
          }

          case 2:
          {
          conn[0] = this->node_id(8);
          conn[1] = this->node_id(7);
          conn[2] = this->node_id(2);
          conn[3] = this->node_id(17);
          conn[4] = this->node_id(16);
          conn[5] = this->node_id(11);

          return;
          }

          case 3:
          {
          conn[0] = this->node_id(6);
          conn[1] = this->node_id(7);
          conn[2] = this->node_id(8);
          conn[3] = this->node_id(15);
          conn[4] = this->node_id(16);
          conn[5] = this->node_id(17);

          return;
          }

          case 4:
          {
          conn[0] = this->node_id(9);
          conn[1] = this->node_id(15);
          conn[2] = this->node_id(17);
          conn[3] = this->node_id(3);
          conn[4] = this->node_id(12);
          conn[5] = this->node_id(14);

          return;
          }

          case 5:
          {
          conn[0] = this->node_id(15);
          conn[1] = this->node_id(10);
          conn[2] = this->node_id(16);
          conn[3] = this->node_id(12);
          conn[4] = this->node_id(4);
          conn[5] = this->node_id(13);

          return;
          }

          case 6:
          {
          conn[0] = this->node_id(17);
          conn[1] = this->node_id(16);
          conn[2] = this->node_id(11);
          conn[3] = this->node_id(14);
          conn[4] = this->node_id(13);
          conn[5] = this->node_id(5);

          return;
          }

          case 7:
          {
          conn[0] = this->node_id(15);
          conn[1] = this->node_id(16);
          conn[2] = this->node_id(17);
          conn[3] = this->node_id(12);
          conn[4] = this->node_id(13);
          conn[5] = this->node_id(14);

          return;
          }

          default:
          libmesh_error_msg("Invalid sc = " << sc);
          }
        */
      }

    default:
      libmesh_error_msg("Unsupported IO package " << iop);
    }
}




unsigned int Prism18::n_second_order_adjacent_vertices (const unsigned int n) const
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

  libmesh_error_msg("We'll never get here!");
  return libMesh::invalid_uint;
}





unsigned short int Prism18::second_order_adjacent_vertex (const unsigned int n,
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
        return _second_order_adjacent_vertices[n-this->n_vertices()][v];
      }

    }

  libmesh_error_msg("We'll never ge here!");
  return static_cast<unsigned short int>(-1);
}



const unsigned short int Prism18::_remaining_second_order_adjacent_vertices[3][4] =
  {
    { 0,  1,  3,  4}, // vertices adjacent to node 15
    { 1,  2,  4,  5}, // vertices adjacent to node 16
    { 0,  2,  3,  5}  // vertices adjacent to node 17
  };



std::pair<unsigned short int, unsigned short int>
Prism18::second_order_child_vertex (const unsigned int n) const
{
  libmesh_assert_greater_equal (n, this->n_vertices());
  libmesh_assert_less (n, this->n_nodes());

  return std::pair<unsigned short int, unsigned short int>
    (_second_order_vertex_child_number[n],
     _second_order_vertex_child_index[n]);
}



Real Prism18::volume () const
{
  // Make copies of our points.  It makes the subsequent calculations a bit
  // shorter and avoids dereferencing the same pointer multiple times.
  Point
    x0 = point(0),   x1 = point(1),   x2 = point(2),   x3 = point(3),   x4 = point(4), x5 = point(5),
    x6 = point(6),   x7 = point(7),   x8 = point(8),   x9 = point(9),  x10 = point(10), x11 = point(11),
    x12 = point(12), x13 = point(13), x14 = point(14), x15 = point(15), x16 = point(16), x17 = point(17);

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
    w1 = 1.1169079483900573284750350421656140e-01L,
    w2 = 5.4975871827660933819163162450105264e-02L,
    a1 = 4.4594849091596488631832925388305199e-01L,
    a2 = 9.1576213509770743459571463402201508e-02L;

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

  Real vol = 0.;
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




#ifdef LIBMESH_ENABLE_AMR

const float Prism18::_embedding_matrix[8][18][18] =
  {
    // embedding matrix for child 0
    {
      //       0         1         2         3         4         5         6         7         8         9        10        11        12        13        14        15        16        17
      {       1.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 0
      {       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 1
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 2
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 3
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.}, // 4
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.}, // 5
      {    0.375,   -0.125,       0.,       0.,       0.,       0.,     0.75,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 6
      {       0.,   -0.125,   -0.125,       0.,       0.,       0.,      0.5,     0.25,      0.5,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 7
      {    0.375,       0.,   -0.125,       0.,       0.,       0.,       0.,       0.,     0.75,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 8
      {    0.375,       0.,       0.,   -0.125,       0.,       0.,       0.,       0.,       0.,     0.75,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 9
      {       0.,       0.,       0.,       0.,       0.,       0.,    0.375,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,       0.,     0.75,       0.,       0.}, // 10
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,    0.375,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,       0.,     0.75}, // 11
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,    0.375,   -0.125,       0.,       0.,       0.,       0.,     0.75,       0.,       0.}, // 12
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,   -0.125,       0.,       0.,       0.,      0.5,     0.25,      0.5}, // 13
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,    0.375,       0.,   -0.125,       0.,       0.,       0.,       0.,       0.,     0.75}, // 14
      { 0.140625,-0.046875,       0.,-0.046875, 0.015625,       0.,  0.28125,       0.,       0.,  0.28125, -0.09375,       0., -0.09375,       0.,       0.,   0.5625,       0.,       0.}, // 15
      {       0.,-0.046875,-0.046875,       0., 0.015625, 0.015625,   0.1875,  0.09375,   0.1875,       0., -0.09375, -0.09375,  -0.0625, -0.03125,  -0.0625,    0.375,   0.1875,    0.375}, // 16
      { 0.140625,       0.,-0.046875,-0.046875,       0., 0.015625,       0.,       0.,  0.28125,  0.28125,       0., -0.09375,       0.,       0., -0.09375,       0.,       0.,   0.5625}  // 17
    },

    // embedding matrix for child 1
    {
      //       0         1         2         3         4         5         6         7         8         9        10        11        12        13        14        15        16        17
      {       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 0
      {       0.,       1.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 1
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 2
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.}, // 3
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 4
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.}, // 5
      {   -0.125,    0.375,       0.,       0.,       0.,       0.,     0.75,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 6
      {       0.,    0.375,   -0.125,       0.,       0.,       0.,       0.,     0.75,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 7
      {   -0.125,       0.,   -0.125,       0.,       0.,       0.,      0.5,      0.5,     0.25,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 8
      {       0.,       0.,       0.,       0.,       0.,       0.,    0.375,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,       0.,     0.75,       0.,       0.}, // 9
      {       0.,    0.375,       0.,       0.,   -0.125,       0.,       0.,       0.,       0.,       0.,     0.75,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 10
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,    0.375,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,       0.,     0.75,       0.}, // 11
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,    0.375,       0.,       0.,       0.,       0.,     0.75,       0.,       0.}, // 12
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,    0.375,   -0.125,       0.,       0.,       0.,       0.,     0.75,       0.}, // 13
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,   -0.125,       0.,       0.,       0.,      0.5,      0.5,     0.25}, // 14
      {-0.046875, 0.140625,       0., 0.015625,-0.046875,       0.,  0.28125,       0.,       0., -0.09375,  0.28125,       0., -0.09375,       0.,       0.,   0.5625,       0.,       0.}, // 15
      {       0., 0.140625,-0.046875,       0.,-0.046875, 0.015625,       0.,  0.28125,       0.,       0.,  0.28125, -0.09375,       0., -0.09375,       0.,       0.,   0.5625,       0.}, // 16
      {-0.046875,       0.,-0.046875, 0.015625,       0., 0.015625,   0.1875,   0.1875,  0.09375, -0.09375,       0., -0.09375,  -0.0625,  -0.0625, -0.03125,    0.375,    0.375,   0.1875}  // 17
    },

    // embedding matrix for child 2
    {
      //       0         1         2         3         4         5         6         7         8         9        10        11        12        13        14        15        16        17
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 0
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 1
      {       0.,       0.,       1.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 2
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.}, // 3
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.}, // 4
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.,       0.,       0.}, // 5
      {   -0.125,   -0.125,       0.,       0.,       0.,       0.,     0.25,      0.5,      0.5,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 6
      {       0.,   -0.125,    0.375,       0.,       0.,       0.,       0.,     0.75,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 7
      {   -0.125,       0.,    0.375,       0.,       0.,       0.,       0.,       0.,     0.75,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 8
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,    0.375,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,       0.,     0.75}, // 9
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,    0.375,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,       0.,     0.75,       0.}, // 10
      {       0.,       0.,    0.375,       0.,       0.,   -0.125,       0.,       0.,       0.,       0.,       0.,     0.75,       0.,       0.,       0.,       0.,       0.,       0.}, // 11
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,   -0.125,       0.,       0.,       0.,       0.,     0.25,      0.5,      0.5}, // 12
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,    0.375,       0.,       0.,       0.,       0.,     0.75,       0.}, // 13
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,    0.375,       0.,       0.,       0.,       0.,       0.,     0.75}, // 14
      {-0.046875,-0.046875,       0., 0.015625, 0.015625,       0.,  0.09375,   0.1875,   0.1875, -0.09375, -0.09375,       0., -0.03125,  -0.0625,  -0.0625,   0.1875,    0.375,    0.375}, // 15
      {       0.,-0.046875, 0.140625,       0., 0.015625,-0.046875,       0.,  0.28125,       0.,       0., -0.09375,  0.28125,       0., -0.09375,       0.,       0.,   0.5625,       0.}, // 16
      {-0.046875,       0., 0.140625, 0.015625,       0.,-0.046875,       0.,       0.,  0.28125, -0.09375,       0.,  0.28125,       0.,       0., -0.09375,       0.,       0.,   0.5625}  // 17
    },

    // embedding matrix for child 3
    {
      //       0         1         2         3         4         5         6         7         8         9        10        11        12        13        14        15        16        17
      {       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 0
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 1
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 2
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.}, // 3
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.}, // 4
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.}, // 5
      {   -0.125,       0.,   -0.125,       0.,       0.,       0.,      0.5,      0.5,     0.25,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 6
      {   -0.125,   -0.125,       0.,       0.,       0.,       0.,     0.25,      0.5,      0.5,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 7
      {       0.,   -0.125,   -0.125,       0.,       0.,       0.,      0.5,     0.25,      0.5,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 8
      {       0.,       0.,       0.,       0.,       0.,       0.,    0.375,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,       0.,     0.75,       0.,       0.}, // 9
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,    0.375,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,       0.,     0.75,       0.}, // 10
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,    0.375,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,       0.,     0.75}, // 11
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,   -0.125,       0.,       0.,       0.,      0.5,      0.5,     0.25}, // 12
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,   -0.125,       0.,       0.,       0.,       0.,     0.25,      0.5,      0.5}, // 13
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,   -0.125,       0.,       0.,       0.,      0.5,     0.25,      0.5}, // 14
      {-0.046875,       0.,-0.046875, 0.015625,       0., 0.015625,   0.1875,   0.1875,  0.09375, -0.09375,       0., -0.09375,  -0.0625,  -0.0625, -0.03125,    0.375,    0.375,   0.1875}, // 15
      {-0.046875,-0.046875,       0., 0.015625, 0.015625,       0.,  0.09375,   0.1875,   0.1875, -0.09375, -0.09375,       0., -0.03125,  -0.0625,  -0.0625,   0.1875,    0.375,    0.375}, // 16
      {       0.,-0.046875,-0.046875,       0., 0.015625, 0.015625,   0.1875,  0.09375,   0.1875,       0., -0.09375, -0.09375,  -0.0625, -0.03125,  -0.0625,    0.375,   0.1875,    0.375}  // 17
    },

    // embedding matrix for child 4
    {
      //       0         1         2         3         4         5         6         7         8         9        10        11        12        13        14        15        16        17
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 0
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.}, // 1
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.}, // 2
      {       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 3
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.,       0.}, // 4
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.}, // 5
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,    0.375,   -0.125,       0.,       0.,       0.,       0.,     0.75,       0.,       0.}, // 6
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,   -0.125,       0.,       0.,       0.,      0.5,     0.25,      0.5}, // 7
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,    0.375,       0.,   -0.125,       0.,       0.,       0.,       0.,       0.,     0.75}, // 8
      {   -0.125,       0.,       0.,    0.375,       0.,       0.,       0.,       0.,       0.,     0.75,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 9
      {       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,       0.,       0.,       0.,       0.,    0.375,       0.,       0.,     0.75,       0.,       0.}, // 10
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,       0.,       0.,       0.,       0.,    0.375,       0.,       0.,     0.75}, // 11
      {       0.,       0.,       0.,    0.375,   -0.125,       0.,       0.,       0.,       0.,       0.,       0.,       0.,     0.75,       0.,       0.,       0.,       0.,       0.}, // 12
      {       0.,       0.,       0.,       0.,   -0.125,   -0.125,       0.,       0.,       0.,       0.,       0.,       0.,      0.5,     0.25,      0.5,       0.,       0.,       0.}, // 13
      {       0.,       0.,       0.,    0.375,       0.,   -0.125,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,     0.75,       0.,       0.,       0.}, // 14
      {-0.046875, 0.015625,       0., 0.140625,-0.046875,       0., -0.09375,       0.,       0.,  0.28125, -0.09375,       0.,  0.28125,       0.,       0.,   0.5625,       0.,       0.}, // 15
      {       0., 0.015625, 0.015625,       0.,-0.046875,-0.046875,  -0.0625, -0.03125,  -0.0625,       0., -0.09375, -0.09375,   0.1875,  0.09375,   0.1875,    0.375,   0.1875,    0.375}, // 16
      {-0.046875,       0., 0.015625, 0.140625,       0.,-0.046875,       0.,       0., -0.09375,  0.28125,       0., -0.09375,       0.,       0.,  0.28125,       0.,       0.,   0.5625}  // 17
    },

    // embedding matrix for child 5
    {
      //       0         1         2         3         4         5         6         7         8         9        10        11        12        13        14        15        16        17
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.}, // 0
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 1
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.}, // 2
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.,       0.}, // 3
      {       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 4
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.}, // 5
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,    0.375,       0.,       0.,       0.,       0.,     0.75,       0.,       0.}, // 6
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,    0.375,   -0.125,       0.,       0.,       0.,       0.,     0.75,       0.}, // 7
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,   -0.125,       0.,       0.,       0.,      0.5,      0.5,     0.25}, // 8
      {       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,       0.,       0.,       0.,       0.,    0.375,       0.,       0.,     0.75,       0.,       0.}, // 9
      {       0.,   -0.125,       0.,       0.,    0.375,       0.,       0.,       0.,       0.,       0.,     0.75,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 10
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,       0.,       0.,       0.,       0.,    0.375,       0.,       0.,     0.75,       0.}, // 11
      {       0.,       0.,       0.,   -0.125,    0.375,       0.,       0.,       0.,       0.,       0.,       0.,       0.,     0.75,       0.,       0.,       0.,       0.,       0.}, // 12
      {       0.,       0.,       0.,       0.,    0.375,   -0.125,       0.,       0.,       0.,       0.,       0.,       0.,       0.,     0.75,       0.,       0.,       0.,       0.}, // 13
      {       0.,       0.,       0.,   -0.125,       0.,   -0.125,       0.,       0.,       0.,       0.,       0.,       0.,      0.5,      0.5,     0.25,       0.,       0.,       0.}, // 14
      { 0.015625,-0.046875,       0.,-0.046875, 0.140625,       0., -0.09375,       0.,       0., -0.09375,  0.28125,       0.,  0.28125,       0.,       0.,   0.5625,       0.,       0.}, // 15
      {       0.,-0.046875, 0.015625,       0., 0.140625,-0.046875,       0., -0.09375,       0.,       0.,  0.28125, -0.09375,       0.,  0.28125,       0.,       0.,   0.5625,       0.}, // 16
      { 0.015625,       0., 0.015625,-0.046875,       0.,-0.046875,  -0.0625,  -0.0625, -0.03125, -0.09375,       0., -0.09375,   0.1875,   0.1875,  0.09375,    0.375,    0.375,   0.1875}  // 17
    },

    // embedding matrix for child 6
    {
      //       0         1         2         3         4         5         6         7         8         9        10        11        12        13        14        15        16        17
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.}, // 0
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.}, // 1
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.,       0.,       0.}, // 2
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.}, // 3
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.}, // 4
      {       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 5
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,   -0.125,       0.,       0.,       0.,       0.,     0.25,      0.5,      0.5}, // 6
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,    0.375,       0.,       0.,       0.,       0.,     0.75,       0.}, // 7
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,    0.375,       0.,       0.,       0.,       0.,       0.,     0.75}, // 8
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,       0.,       0.,       0.,       0.,    0.375,       0.,       0.,     0.75}, // 9
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,       0.,       0.,       0.,       0.,    0.375,       0.,       0.,     0.75,       0.}, // 10
      {       0.,       0.,   -0.125,       0.,       0.,    0.375,       0.,       0.,       0.,       0.,       0.,     0.75,       0.,       0.,       0.,       0.,       0.,       0.}, // 11
      {       0.,       0.,       0.,   -0.125,   -0.125,       0.,       0.,       0.,       0.,       0.,       0.,       0.,     0.25,      0.5,      0.5,       0.,       0.,       0.}, // 12
      {       0.,       0.,       0.,       0.,   -0.125,    0.375,       0.,       0.,       0.,       0.,       0.,       0.,       0.,     0.75,       0.,       0.,       0.,       0.}, // 13
      {       0.,       0.,       0.,   -0.125,       0.,    0.375,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,     0.75,       0.,       0.,       0.}, // 14
      { 0.015625, 0.015625,       0.,-0.046875,-0.046875,       0., -0.03125,  -0.0625,  -0.0625, -0.09375, -0.09375,       0.,  0.09375,   0.1875,   0.1875,   0.1875,    0.375,    0.375}, // 15
      {       0., 0.015625,-0.046875,       0.,-0.046875, 0.140625,       0., -0.09375,       0.,       0., -0.09375,  0.28125,       0.,  0.28125,       0.,       0.,   0.5625,       0.}, // 16
      { 0.015625,       0.,-0.046875,-0.046875,       0., 0.140625,       0.,       0., -0.09375, -0.09375,       0.,  0.28125,       0.,       0.,  0.28125,       0.,       0.,   0.5625}  // 17
    },

    // embedding matrix for child 7
    {
      //       0         1         2         3         4         5         6         7         8         9        10        11        12        13        14        15        16        17
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.}, // 0
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.}, // 1
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.}, // 2
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.,       0.}, // 3
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.}, // 4
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.}, // 5
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,   -0.125,       0.,       0.,       0.,      0.5,      0.5,     0.25}, // 6
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,   -0.125,       0.,       0.,       0.,       0.,     0.25,      0.5,      0.5}, // 7
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,   -0.125,       0.,       0.,       0.,      0.5,     0.25,      0.5}, // 8
      {       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,       0.,       0.,       0.,       0.,    0.375,       0.,       0.,     0.75,       0.,       0.}, // 9
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,       0.,       0.,       0.,       0.,    0.375,       0.,       0.,     0.75,       0.}, // 10
      {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,       0.,       0.,       0.,       0.,    0.375,       0.,       0.,     0.75}, // 11
      {       0.,       0.,       0.,   -0.125,       0.,   -0.125,       0.,       0.,       0.,       0.,       0.,       0.,      0.5,      0.5,     0.25,       0.,       0.,       0.}, // 12
      {       0.,       0.,       0.,   -0.125,   -0.125,       0.,       0.,       0.,       0.,       0.,       0.,       0.,     0.25,      0.5,      0.5,       0.,       0.,       0.}, // 13
      {       0.,       0.,       0.,       0.,   -0.125,   -0.125,       0.,       0.,       0.,       0.,       0.,       0.,      0.5,     0.25,      0.5,       0.,       0.,       0.}, // 14
      { 0.015625,       0., 0.015625,-0.046875,       0.,-0.046875,  -0.0625,  -0.0625, -0.03125, -0.09375,       0., -0.09375,   0.1875,   0.1875,  0.09375,    0.375,    0.375,   0.1875}, // 15
      { 0.015625, 0.015625,       0.,-0.046875,-0.046875,       0., -0.03125,  -0.0625,  -0.0625, -0.09375, -0.09375,       0.,  0.09375,   0.1875,   0.1875,   0.1875,    0.375,    0.375}, // 16
      {       0., 0.015625, 0.015625,       0.,-0.046875,-0.046875,  -0.0625, -0.03125,  -0.0625,       0., -0.09375, -0.09375,   0.1875,  0.09375,   0.1875,    0.375,   0.1875,    0.375}  // 17
    }
  };

#endif

} // namespace libMesh
