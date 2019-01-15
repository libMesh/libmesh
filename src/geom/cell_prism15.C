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


// Local includes
#include "libmesh/side.h"
#include "libmesh/cell_prism15.h"
#include "libmesh/edge_edge3.h"
#include "libmesh/face_quad8.h"
#include "libmesh/face_tri6.h"
#include "libmesh/enum_io_package.h"
#include "libmesh/enum_order.h"

namespace libMesh
{



// ------------------------------------------------------------
// Prism15 class static member initializations
const int Prism15::num_nodes;
const int Prism15::num_sides;
const int Prism15::num_edges;
const int Prism15::num_children;
const int Prism15::nodes_per_side;
const int Prism15::nodes_per_edge;

const unsigned int Prism15::side_nodes_map[Prism15::num_sides][Prism15::nodes_per_side] =
  {
    {0, 2, 1,  8,  7,  6, 99, 99}, // Side 0
    {0, 1, 4,  3,  6, 10, 12,  9}, // Side 1
    {1, 2, 5,  4,  7, 11, 13, 10}, // Side 2
    {2, 0, 3,  5,  8,  9, 14, 11}, // Side 3
    {3, 4, 5, 12, 13, 14, 99, 99}  // Side 4
  };

const unsigned int Prism15::edge_nodes_map[Prism15::num_edges][Prism15::nodes_per_edge] =
  {
    {0, 1,  6}, // Edge 0
    {1, 2,  7}, // Edge 1
    {0, 2,  8}, // Edge 2
    {0, 3,  9}, // Edge 3
    {1, 4, 10}, // Edge 4
    {2, 5, 11}, // Edge 5
    {3, 4, 12}, // Edge 6
    {4, 5, 13}, // Edge 7
    {3, 5, 14}  // Edge 8
  };


// ------------------------------------------------------------
// Prism15 class member functions

bool Prism15::is_vertex(const unsigned int i) const
{
  if (i < 6)
    return true;
  return false;
}

bool Prism15::is_edge(const unsigned int i) const
{
  if (i < 6)
    return false;
  return true;
}

bool Prism15::is_face(const unsigned int) const
{
  return false;
}

bool Prism15::is_node_on_side(const unsigned int n,
                              const unsigned int s) const
{
  libmesh_assert_less (s, n_sides());
  return std::find(std::begin(side_nodes_map[s]),
                   std::end(side_nodes_map[s]),
                   n) != std::end(side_nodes_map[s]);
}

std::vector<unsigned int>
Prism15::nodes_on_side(const unsigned int s) const
{
  libmesh_assert_less(s, n_sides());
  auto trim = (s > 0 && s < 4) ? 0 : 2;
  return {std::begin(side_nodes_map[s]), std::end(side_nodes_map[s]) - trim};
}

bool Prism15::is_node_on_edge(const unsigned int n,
                              const unsigned int e) const
{
  libmesh_assert_less (e, n_edges());
  return std::find(std::begin(edge_nodes_map[e]),
                   std::end(edge_nodes_map[e]),
                   n) != std::end(edge_nodes_map[e]);
}



bool Prism15::has_affine_map() const
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
      !v.relative_fuzzy_equals(this->point(11) - this->point(2)))
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



Order Prism15::default_order() const
{
  return SECOND;
}



unsigned int Prism15::which_node_am_i(unsigned int side,
                                      unsigned int side_node) const
{
  libmesh_assert_less (side, this->n_sides());

  // Never more than 8 nodes per side.
  libmesh_assert_less(side_node, Prism15::nodes_per_side);

  // Some sides have 6 nodes.
  libmesh_assert(!(side==0 || side==4) || side_node < 6);

  return Prism15::side_nodes_map[side][side_node];
}



std::unique_ptr<Elem> Prism15::build_side_ptr (const unsigned int i,
                                               bool proxy)
{
  libmesh_assert_less (i, this->n_sides());

  if (proxy)
    {
      switch (i)
        {
        case 0:  // the triangular face at z=-1
        case 4:
          return libmesh_make_unique<Side<Tri6,Prism15>>(this,i);

        case 1:
        case 2:
        case 3:
          return libmesh_make_unique<Side<Quad8,Prism15>>(this,i);

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
            face = libmesh_make_unique<Tri6>();
            break;
          }
        case 1: // the quad face at y=0
        case 2: // the other quad face
        case 3: // the quad face at x=0
          {
            face = libmesh_make_unique<Quad8>();
            break;
          }
        default:
          libmesh_error_msg("Invalid side i = " << i);
        }

      face->subdomain_id() = this->subdomain_id();

      // Set the nodes
      for (unsigned n=0; n<face->n_nodes(); ++n)
        face->set_node(n) = this->node_ptr(Prism15::side_nodes_map[i][n]);

      return face;
    }
}


void Prism15::build_side_ptr (std::unique_ptr<Elem> & side,
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
        if (!side.get() || side->type() != QUAD8)
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
    side->set_node(n) = this->node_ptr(Prism15::side_nodes_map[i][n]);
}



std::unique_ptr<Elem> Prism15::build_edge_ptr (const unsigned int i)
{
  libmesh_assert_less (i, this->n_edges());

  return libmesh_make_unique<SideEdge<Edge3,Prism15>>(this,i);
}


void Prism15::connectivity(const unsigned int libmesh_dbg_var(sc),
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
        conn[0] = this->node_id(0)+1;
        conn[1] = this->node_id(1)+1;
        conn[2] = this->node_id(2)+1;
        conn[3] = this->node_id(2)+1;
        conn[4] = this->node_id(3)+1;
        conn[5] = this->node_id(4)+1;
        conn[6] = this->node_id(5)+1;
        conn[7] = this->node_id(5)+1;
        return;
      }

    case VTK:
      {
        /*
          conn.resize(6);
          conn[0] = this->node_id(0);
          conn[1] = this->node_id(2);
          conn[2] = this->node_id(1);
          conn[3] = this->node_id(3);
          conn[4] = this->node_id(5);
          conn[5] = this->node_id(4);
        */

        // VTK's VTK_QUADRATIC_WEDGE first 9 nodes match, then their
        // middle and top layers of mid-edge nodes are reversed from
        // LibMesh's.
        conn.resize(15);
        for (unsigned i=0; i<9; ++i)
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




unsigned short int Prism15::second_order_adjacent_vertex (const unsigned int n,
                                                          const unsigned int v) const
{
  libmesh_assert_greater_equal (n, this->n_vertices());
  libmesh_assert_less (n, this->n_nodes());
  libmesh_assert_less (v, 2);
  return _second_order_adjacent_vertices[n-this->n_vertices()][v];
}



std::pair<unsigned short int, unsigned short int>
Prism15::second_order_child_vertex (const unsigned int n) const
{
  libmesh_assert_greater_equal (n, this->n_vertices());
  libmesh_assert_less (n, this->n_nodes());

  return std::pair<unsigned short int, unsigned short int>
    (_second_order_vertex_child_number[n],
     _second_order_vertex_child_index[n]);
}



Real Prism15::volume () const
{
  // Make copies of our points.  It makes the subsequent calculations a bit
  // shorter and avoids dereferencing the same pointer multiple times.
  Point
    x0 = point(0),   x1 = point(1),   x2 = point(2),   x3 = point(3),   x4 = point(4),
    x5 = point(5),   x6 = point(6),   x7 = point(7),   x8 = point(8),   x9 = point(9),
    x10 = point(10), x11 = point(11), x12 = point(12), x13 = point(13), x14 = point(14);

  // Terms are copied directly from a Python script.
  Point dx_dxi[10] =
    {
      -x0 - x1 + x10 + 2*x12 - x3 - x4 + 2*x6 - x9,
      3*x0/2 + x1/2 + 2*x12 - 3*x3/2 - x4/2 - 2*x6,
      -x0/2 + x1/2 - x10 - x3/2 + x4/2 + x9,
      2*x0 - 2*x12 + 2*x13 - 2*x14 + 2*x3 - 2*x6 + 2*x7 - 2*x8,
      -2*x0 - 2*x12 + 2*x13 - 2*x14 + 2*x3 + 2*x6 - 2*x7 + 2*x8,
      Point(0,0,0),
      2*x0 + 2*x1 - 4*x12 + 2*x3 + 2*x4 - 4*x6,
      -2*x0 - 2*x1 - 4*x12 + 2*x3 + 2*x4 + 4*x6,
      Point(0,0,0),
      Point(0,0,0)
    };

  Point dx_deta[10] =
    {
      -x0 + x11 + 2*x14 - x2 - x3 - x5 + 2*x8 - x9,
      3*x0/2 + 2*x14 + x2/2 - 3*x3/2 - x5/2 - 2*x8,
      -x0/2 - x11 + x2/2 - x3/2 + x5/2 + x9,
      2*x0 - 4*x14 + 2*x2 + 2*x3 + 2*x5 - 4*x8,
      -2*x0 - 4*x14 - 2*x2 + 2*x3 + 2*x5 + 4*x8,
      Point(0,0,0),
      2*x0 - 2*x12 + 2*x13 - 2*x14 + 2*x3 - 2*x6 + 2*x7 - 2*x8,
      -2*x0 - 2*x12 + 2*x13 - 2*x14 + 2*x3 + 2*x6 - 2*x7 + 2*x8,
      Point(0,0,0),
      Point(0,0,0)
    };

  Point dx_dzeta[10] =
    {
      -x0/2 + x3/2,
      x0 + x3 - 2*x9,
      Point(0,0,0),
      3*x0/2 + 2*x14 + x2/2 - 3*x3/2 - x5/2 - 2*x8,
      -x0 - 2*x11 + x2 - x3 + x5 + 2*x9,
      -x0 - 2*x14 - x2 + x3 + x5 + 2*x8,
      3*x0/2 + x1/2 + 2*x12 - 3*x3/2 - x4/2 - 2*x6,
      -x0 + x1 - 2*x10 - x3 + x4 + 2*x9,
      -2*x0 - 2*x12 + 2*x13 - 2*x14 + 2*x3 + 2*x6 - 2*x7 + 2*x8,
      -x0 - x1 - 2*x12 + x3 + x4 + 2*x6
    };

  // The quadrature rule for the Prism15 is a tensor product between a
  // FOURTH-order TRI3 rule (in xi, eta) and a FIFTH-order EDGE2 rule
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
  static const int exponents[10][3] =
    {
      {0, 0, 0},
      {0, 0, 1},
      {0, 0, 2},
      {0, 1, 0},
      {0, 1, 1},
      {0, 2, 0},
      {1, 0, 0},
      {1, 0, 1},
      {1, 1, 0},
      {2, 0, 0}
    };

  Real vol = 0.;
  for (int i=0; i<N2D; ++i)
    for (int j=0; j<N1D; ++j)
      {
        // Compute dx_dxi, dx_deta, dx_dzeta at the current quadrature point.
        Point dx_dxi_q, dx_deta_q, dx_dzeta_q;
        for (int c=0; c<10; ++c)
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

const float Prism15::_embedding_matrix[Prism15::num_children][Prism15::num_nodes][Prism15::num_nodes] =
  {
    // Embedding matrix for child 0
    {
      //       0        1        2        3        4        5        6        7        8        9       10       11       12       13       14
      {       1,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0 }, //  0
      {       0,       0,       0,       0,       0,       0,       1,       0,       0,       0,       0,       0,       0,       0,       0 }, //  1
      {       0,       0,       0,       0,       0,       0,       0,       0,       1,       0,       0,       0,       0,       0,       0 }, //  2
      {       0,       0,       0,       0,       0,       0,       0,       0,       0,       1,       0,       0,       0,       0,       0 }, //  3
      {   -0.25,   -0.25,       0,   -0.25,   -0.25,       0,     0.5,       0,       0,     0.5,     0.5,       0,     0.5,       0,       0 }, //  4
      {   -0.25,       0,   -0.25,   -0.25,       0,   -0.25,       0,       0,     0.5,     0.5,       0,     0.5,       0,       0,     0.5 }, //  5
      {   0.375,  -0.125,       0,       0,       0,       0,    0.75,       0,       0,       0,       0,       0,       0,       0,       0 }, //  6
      {       0,  -0.125,  -0.125,       0,       0,       0,     0.5,    0.25,     0.5,       0,       0,       0,       0,       0,       0 }, //  7
      {   0.375,       0,  -0.125,       0,       0,       0,       0,       0,    0.75,       0,       0,       0,       0,       0,       0 }, //  8
      {   0.375,       0,       0,  -0.125,       0,       0,       0,       0,       0,    0.75,       0,       0,       0,       0,       0 }, //  9
      { -0.1875, -0.1875,       0, -0.1875, -0.1875,       0,    0.75,       0,       0,   0.375,   0.375,       0,    0.25,       0,       0 }, // 10
      { -0.1875,       0, -0.1875, -0.1875,       0, -0.1875,       0,       0,    0.75,   0.375,       0,   0.375,       0,       0,    0.25 }, // 11
      { -0.1875, -0.1875,       0, -0.1875, -0.1875,       0,   0.375,       0,       0,    0.75,    0.25,       0,   0.375,       0,       0 }, // 12
      {   -0.25, -0.1875, -0.1875,   -0.25, -0.1875, -0.1875,    0.25,   0.125,    0.25,     0.5,    0.25,    0.25,    0.25,   0.125,    0.25 }, // 13
      { -0.1875,       0, -0.1875, -0.1875,       0, -0.1875,       0,       0,   0.375,    0.75,       0,    0.25,       0,       0,   0.375 }  // 14
    },

    // Embedding matrix for child 1
    {
      //       0        1        2        3        4        5        6        7        8        9       10       11       12       13       14
      {       0,       0,       0,       0,       0,       0,       1,       0,       0,       0,       0,       0,       0,       0,       0 }, //  0
      {       0,       1,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0 }, //  1
      {       0,       0,       0,       0,       0,       0,       0,       1,       0,       0,       0,       0,       0,       0,       0 }, //  2
      {   -0.25,   -0.25,       0,   -0.25,   -0.25,       0,     0.5,       0,       0,     0.5,     0.5,       0,     0.5,       0,       0 }, //  3
      {       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       1,       0,       0,       0,       0 }, //  4
      {       0,   -0.25,   -0.25,       0,   -0.25,   -0.25,       0,     0.5,       0,       0,     0.5,     0.5,       0,     0.5,       0 }, //  5
      {  -0.125,   0.375,       0,       0,       0,       0,    0.75,       0,       0,       0,       0,       0,       0,       0,       0 }, //  6
      {       0,   0.375,  -0.125,       0,       0,       0,       0,    0.75,       0,       0,       0,       0,       0,       0,       0 }, //  7
      {  -0.125,       0,  -0.125,       0,       0,       0,     0.5,     0.5,    0.25,       0,       0,       0,       0,       0,       0 }, //  8
      { -0.1875, -0.1875,       0, -0.1875, -0.1875,       0,    0.75,       0,       0,   0.375,   0.375,       0,    0.25,       0,       0 }, //  9
      {       0,   0.375,       0,       0,  -0.125,       0,       0,       0,       0,       0,    0.75,       0,       0,       0,       0 }, // 10
      {       0, -0.1875, -0.1875,       0, -0.1875, -0.1875,       0,    0.75,       0,       0,   0.375,   0.375,       0,    0.25,       0 }, // 11
      { -0.1875, -0.1875,       0, -0.1875, -0.1875,       0,   0.375,       0,       0,    0.25,    0.75,       0,   0.375,       0,       0 }, // 12
      {       0, -0.1875, -0.1875,       0, -0.1875, -0.1875,       0,   0.375,       0,       0,    0.75,    0.25,       0,   0.375,       0 }, // 13
      { -0.1875,   -0.25, -0.1875, -0.1875,   -0.25, -0.1875,    0.25,    0.25,   0.125,    0.25,     0.5,    0.25,    0.25,    0.25,   0.125 }  // 14
    },

    // Embedding matrix for child 2
    {
      //       0        1        2        3        4        5        6        7        8        9       10       11       12       13       14
      {       0,       0,       0,       0,       0,       0,       0,       0,       1,       0,       0,       0,       0,       0,       0 }, //  0
      {       0,       0,       0,       0,       0,       0,       0,       1,       0,       0,       0,       0,       0,       0,       0 }, //  1
      {       0,       0,       1,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0 }, //  2
      {   -0.25,       0,   -0.25,   -0.25,       0,   -0.25,       0,       0,     0.5,     0.5,       0,     0.5,       0,       0,     0.5 }, //  3
      {       0,   -0.25,   -0.25,       0,   -0.25,   -0.25,       0,     0.5,       0,       0,     0.5,     0.5,       0,     0.5,       0 }, //  4
      {       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       1,       0,       0,       0 }, //  5
      {  -0.125,  -0.125,       0,       0,       0,       0,    0.25,     0.5,     0.5,       0,       0,       0,       0,       0,       0 }, //  6
      {       0,  -0.125,   0.375,       0,       0,       0,       0,    0.75,       0,       0,       0,       0,       0,       0,       0 }, //  7
      {  -0.125,       0,   0.375,       0,       0,       0,       0,       0,    0.75,       0,       0,       0,       0,       0,       0 }, //  8
      { -0.1875,       0, -0.1875, -0.1875,       0, -0.1875,       0,       0,    0.75,   0.375,       0,   0.375,       0,       0,    0.25 }, //  9
      {       0, -0.1875, -0.1875,       0, -0.1875, -0.1875,       0,    0.75,       0,       0,   0.375,   0.375,       0,    0.25,       0 }, // 10
      {       0,       0,   0.375,       0,       0,  -0.125,       0,       0,       0,       0,       0,    0.75,       0,       0,       0 }, // 11
      { -0.1875, -0.1875,   -0.25, -0.1875, -0.1875,   -0.25,   0.125,    0.25,    0.25,    0.25,    0.25,     0.5,   0.125,    0.25,    0.25 }, // 12
      {       0, -0.1875, -0.1875,       0, -0.1875, -0.1875,       0,   0.375,       0,       0,    0.25,    0.75,       0,   0.375,       0 }, // 13
      { -0.1875,       0, -0.1875, -0.1875,       0, -0.1875,       0,       0,   0.375,    0.25,       0,    0.75,       0,       0,   0.375 }  // 14
    },

    // Embedding matrix for child 3
    {
      //       0        1        2        3        4        5        6        7        8        9       10       11       12       13       14
      {       0,       0,       0,       0,       0,       0,       1,       0,       0,       0,       0,       0,       0,       0,       0 }, //  0
      {       0,       0,       0,       0,       0,       0,       0,       1,       0,       0,       0,       0,       0,       0,       0 }, //  1
      {       0,       0,       0,       0,       0,       0,       0,       0,       1,       0,       0,       0,       0,       0,       0 }, //  2
      {   -0.25,   -0.25,       0,   -0.25,   -0.25,       0,     0.5,       0,       0,     0.5,     0.5,       0,     0.5,       0,       0 }, //  3
      {       0,   -0.25,   -0.25,       0,   -0.25,   -0.25,       0,     0.5,       0,       0,     0.5,     0.5,       0,     0.5,       0 }, //  4
      {   -0.25,       0,   -0.25,   -0.25,       0,   -0.25,       0,       0,     0.5,     0.5,       0,     0.5,       0,       0,     0.5 }, //  5
      {  -0.125,       0,  -0.125,       0,       0,       0,     0.5,     0.5,    0.25,       0,       0,       0,       0,       0,       0 }, //  6
      {  -0.125,  -0.125,       0,       0,       0,       0,    0.25,     0.5,     0.5,       0,       0,       0,       0,       0,       0 }, //  7
      {       0,  -0.125,  -0.125,       0,       0,       0,     0.5,    0.25,     0.5,       0,       0,       0,       0,       0,       0 }, //  8
      { -0.1875, -0.1875,       0, -0.1875, -0.1875,       0,    0.75,       0,       0,   0.375,   0.375,       0,    0.25,       0,       0 }, //  9
      {       0, -0.1875, -0.1875,       0, -0.1875, -0.1875,       0,    0.75,       0,       0,   0.375,   0.375,       0,    0.25,       0 }, // 10
      { -0.1875,       0, -0.1875, -0.1875,       0, -0.1875,       0,       0,    0.75,   0.375,       0,   0.375,       0,       0,    0.25 }, // 11
      { -0.1875,   -0.25, -0.1875, -0.1875,   -0.25, -0.1875,    0.25,    0.25,   0.125,    0.25,     0.5,    0.25,    0.25,    0.25,   0.125 }, // 12
      { -0.1875, -0.1875,   -0.25, -0.1875, -0.1875,   -0.25,   0.125,    0.25,    0.25,    0.25,    0.25,     0.5,   0.125,    0.25,    0.25 }, // 13
      {   -0.25, -0.1875, -0.1875,   -0.25, -0.1875, -0.1875,    0.25,   0.125,    0.25,     0.5,    0.25,    0.25,    0.25,   0.125,    0.25 }  // 14
    },

    // Embedding matrix for child 4
    {
      //       0        1        2        3        4        5        6        7        8        9       10       11       12       13       14
      {       0,       0,       0,       0,       0,       0,       0,       0,       0,       1,       0,       0,       0,       0,       0 }, //  0
      {   -0.25,   -0.25,       0,   -0.25,   -0.25,       0,     0.5,       0,       0,     0.5,     0.5,       0,     0.5,       0,       0 }, //  1
      {   -0.25,       0,   -0.25,   -0.25,       0,   -0.25,       0,       0,     0.5,     0.5,       0,     0.5,       0,       0,     0.5 }, //  2
      {       0,       0,       0,       1,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0 }, //  3
      {       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       1,       0,       0 }, //  4
      {       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       1 }, //  5
      { -0.1875, -0.1875,       0, -0.1875, -0.1875,       0,   0.375,       0,       0,    0.75,    0.25,       0,   0.375,       0,       0 }, //  6
      {   -0.25, -0.1875, -0.1875,   -0.25, -0.1875, -0.1875,    0.25,   0.125,    0.25,     0.5,    0.25,    0.25,    0.25,   0.125,    0.25 }, //  7
      { -0.1875,       0, -0.1875, -0.1875,       0, -0.1875,       0,       0,   0.375,    0.75,       0,    0.25,       0,       0,   0.375 }, //  8
      {  -0.125,       0,       0,   0.375,       0,       0,       0,       0,       0,    0.75,       0,       0,       0,       0,       0 }, //  9
      { -0.1875, -0.1875,       0, -0.1875, -0.1875,       0,    0.25,       0,       0,   0.375,   0.375,       0,    0.75,       0,       0 }, // 10
      { -0.1875,       0, -0.1875, -0.1875,       0, -0.1875,       0,       0,    0.25,   0.375,       0,   0.375,       0,       0,    0.75 }, // 11
      {       0,       0,       0,   0.375,  -0.125,       0,       0,       0,       0,       0,       0,       0,    0.75,       0,       0 }, // 12
      {       0,       0,       0,       0,  -0.125,  -0.125,       0,       0,       0,       0,       0,       0,     0.5,    0.25,     0.5 }, // 13
      {       0,       0,       0,   0.375,       0,  -0.125,       0,       0,       0,       0,       0,       0,       0,       0,    0.75 }  // 14
    },

    // Embedding matrix for child 5
    {
      //       0        1        2        3        4        5        6        7        8        9       10       11       12       13       14
      {   -0.25,   -0.25,       0,   -0.25,   -0.25,       0,     0.5,       0,       0,     0.5,     0.5,       0,     0.5,       0,       0 }, //  0
      {       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       1,       0,       0,       0,       0 }, //  1
      {       0,   -0.25,   -0.25,       0,   -0.25,   -0.25,       0,     0.5,       0,       0,     0.5,     0.5,       0,     0.5,       0 }, //  2
      {       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       1,       0,       0 }, //  3
      {       0,       0,       0,       0,       1,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0 }, //  4
      {       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       1,       0 }, //  5
      { -0.1875, -0.1875,       0, -0.1875, -0.1875,       0,   0.375,       0,       0,    0.25,    0.75,       0,   0.375,       0,       0 }, //  6
      {       0, -0.1875, -0.1875,       0, -0.1875, -0.1875,       0,   0.375,       0,       0,    0.75,    0.25,       0,   0.375,       0 }, //  7
      { -0.1875,   -0.25, -0.1875, -0.1875,   -0.25, -0.1875,    0.25,    0.25,   0.125,    0.25,     0.5,    0.25,    0.25,    0.25,   0.125 }, //  8
      { -0.1875, -0.1875,       0, -0.1875, -0.1875,       0,    0.25,       0,       0,   0.375,   0.375,       0,    0.75,       0,       0 }, //  9
      {       0,  -0.125,       0,       0,   0.375,       0,       0,       0,       0,       0,    0.75,       0,       0,       0,       0 }, // 10
      {       0, -0.1875, -0.1875,       0, -0.1875, -0.1875,       0,    0.25,       0,       0,   0.375,   0.375,       0,    0.75,       0 }, // 11
      {       0,       0,       0,  -0.125,   0.375,       0,       0,       0,       0,       0,       0,       0,    0.75,       0,       0 }, // 12
      {       0,       0,       0,       0,   0.375,  -0.125,       0,       0,       0,       0,       0,       0,       0,    0.75,       0 }, // 13
      {       0,       0,       0,  -0.125,       0,  -0.125,       0,       0,       0,       0,       0,       0,     0.5,     0.5,    0.25 }  // 14
    },

    // Embedding matrix for child 6
    {
      //       0        1        2        3        4        5        6        7        8        9       10       11       12       13       14
      {   -0.25,       0,   -0.25,   -0.25,       0,   -0.25,       0,       0,     0.5,     0.5,       0,     0.5,       0,       0,     0.5 }, //  0
      {       0,   -0.25,   -0.25,       0,   -0.25,   -0.25,       0,     0.5,       0,       0,     0.5,     0.5,       0,     0.5,       0 }, //  1
      {       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       1,       0,       0,       0 }, //  2
      {       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       1 }, //  3
      {       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       1,       0 }, //  4
      {       0,       0,       0,       0,       0,       1,       0,       0,       0,       0,       0,       0,       0,       0,       0 }, //  5
      { -0.1875, -0.1875,   -0.25, -0.1875, -0.1875,   -0.25,   0.125,    0.25,    0.25,    0.25,    0.25,     0.5,   0.125,    0.25,    0.25 }, //  6
      {       0, -0.1875, -0.1875,       0, -0.1875, -0.1875,       0,   0.375,       0,       0,    0.25,    0.75,       0,   0.375,       0 }, //  7
      { -0.1875,       0, -0.1875, -0.1875,       0, -0.1875,       0,       0,   0.375,    0.25,       0,    0.75,       0,       0,   0.375 }, //  8
      { -0.1875,       0, -0.1875, -0.1875,       0, -0.1875,       0,       0,    0.25,   0.375,       0,   0.375,       0,       0,    0.75 }, //  9
      {       0, -0.1875, -0.1875,       0, -0.1875, -0.1875,       0,    0.25,       0,       0,   0.375,   0.375,       0,    0.75,       0 }, // 10
      {       0,       0,  -0.125,       0,       0,   0.375,       0,       0,       0,       0,       0,    0.75,       0,       0,       0 }, // 11
      {       0,       0,       0,  -0.125,  -0.125,       0,       0,       0,       0,       0,       0,       0,    0.25,     0.5,     0.5 }, // 12
      {       0,       0,       0,       0,  -0.125,   0.375,       0,       0,       0,       0,       0,       0,       0,    0.75,       0 }, // 13
      {       0,       0,       0,  -0.125,       0,   0.375,       0,       0,       0,       0,       0,       0,       0,       0,    0.75 }  // 14
    },

    // Embedding matrix for child 7
    {
      //       0        1        2        3        4        5        6        7        8        9       10       11       12       13       14
      {   -0.25,   -0.25,       0,   -0.25,   -0.25,       0,     0.5,       0,       0,     0.5,     0.5,       0,     0.5,       0,       0 }, //  0
      {       0,   -0.25,   -0.25,       0,   -0.25,   -0.25,       0,     0.5,       0,       0,     0.5,     0.5,       0,     0.5,       0 }, //  1
      {   -0.25,       0,   -0.25,   -0.25,       0,   -0.25,       0,       0,     0.5,     0.5,       0,     0.5,       0,       0,     0.5 }, //  2
      {       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       1,       0,       0 }, //  3
      {       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       1,       0 }, //  4
      {       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       1 }, //  5
      { -0.1875,   -0.25, -0.1875, -0.1875,   -0.25, -0.1875,    0.25,    0.25,   0.125,    0.25,     0.5,    0.25,    0.25,    0.25,   0.125 }, //  6
      { -0.1875, -0.1875,   -0.25, -0.1875, -0.1875,   -0.25,   0.125,    0.25,    0.25,    0.25,    0.25,     0.5,   0.125,    0.25,    0.25 }, //  7
      {   -0.25, -0.1875, -0.1875,   -0.25, -0.1875, -0.1875,    0.25,   0.125,    0.25,     0.5,    0.25,    0.25,    0.25,   0.125,    0.25 }, //  8
      { -0.1875, -0.1875,       0, -0.1875, -0.1875,       0,    0.25,       0,       0,   0.375,   0.375,       0,    0.75,       0,       0 }, //  9
      {       0, -0.1875, -0.1875,       0, -0.1875, -0.1875,       0,    0.25,       0,       0,   0.375,   0.375,       0,    0.75,       0 }, // 10
      { -0.1875,       0, -0.1875, -0.1875,       0, -0.1875,       0,       0,    0.25,   0.375,       0,   0.375,       0,       0,    0.75 }, // 11
      {       0,       0,       0,  -0.125,       0,  -0.125,       0,       0,       0,       0,       0,       0,     0.5,     0.5,    0.25 }, // 12
      {       0,       0,       0,  -0.125,  -0.125,       0,       0,       0,       0,       0,       0,       0,    0.25,     0.5,     0.5 }, // 13
      {       0,       0,       0,       0,  -0.125,  -0.125,       0,       0,       0,       0,       0,       0,     0.5,    0.25,     0.5 }  // 14
    }
  };

#endif

} // namespace libMesh
