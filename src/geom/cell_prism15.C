// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

std::vector<unsigned>
Prism15::nodes_on_edge(const unsigned int e) const
{
  libmesh_assert_less(e, n_edges());
  return {std::begin(edge_nodes_map[e]), std::end(edge_nodes_map[e])};
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
  if (!v.relative_fuzzy_equals(this->point(4) - this->point(1), affine_tol) ||
      !v.relative_fuzzy_equals(this->point(5) - this->point(2), affine_tol))
    return false;
  // Make sure edges are straight
  v /= 2;
  if (!v.relative_fuzzy_equals(this->point(9) - this->point(0), affine_tol) ||
      !v.relative_fuzzy_equals(this->point(10) - this->point(1), affine_tol) ||
      !v.relative_fuzzy_equals(this->point(11) - this->point(2), affine_tol))
    return false;
  v = (this->point(1) - this->point(0))/2;
  if (!v.relative_fuzzy_equals(this->point(6) - this->point(0), affine_tol) ||
      !v.relative_fuzzy_equals(this->point(12) - this->point(3), affine_tol))
    return false;
  v = (this->point(2) - this->point(0))/2;
  if (!v.relative_fuzzy_equals(this->point(8) - this->point(0), affine_tol) ||
      !v.relative_fuzzy_equals(this->point(14) - this->point(3), affine_tol))
    return false;
  v = (this->point(2) - this->point(1))/2;
  if (!v.relative_fuzzy_equals(this->point(7) - this->point(1), affine_tol) ||
      !v.relative_fuzzy_equals(this->point(13) - this->point(4), affine_tol))
    return false;
  return true;
}



Order Prism15::default_order() const
{
  return SECOND;
}



unsigned int Prism15::local_side_node(unsigned int side,
                                      unsigned int side_node) const
{
  libmesh_assert_less (side, this->n_sides());

  // Never more than 8 nodes per side.
  libmesh_assert_less(side_node, Prism15::nodes_per_side);

  // Some sides have 6 nodes.
  libmesh_assert(!(side==0 || side==4) || side_node < 6);

  return Prism15::side_nodes_map[side][side_node];
}



unsigned int Prism15::local_edge_node(unsigned int edge,
                                      unsigned int edge_node) const
{
  libmesh_assert_less(edge, this->n_edges());
  libmesh_assert_less(edge_node, Prism15::nodes_per_edge);

  return Prism15::edge_nodes_map[edge][edge_node];
}



std::unique_ptr<Elem> Prism15::build_side_ptr (const unsigned int i)
{
  libmesh_assert_less (i, this->n_sides());

  std::unique_ptr<Elem> face;

  switch (i)
    {
    case 0: // the triangular face at z=-1
    case 4: // the triangular face at z=1
      {
        face = std::make_unique<Tri6>();
        break;
      }
    case 1: // the quad face at y=0
    case 2: // the other quad face
    case 3: // the quad face at x=0
      {
        face = std::make_unique<Quad8>();
        break;
      }
    default:
      libmesh_error_msg("Invalid side i = " << i);
    }

  // Set the nodes
  for (auto n : face->node_index_range())
    face->set_node(n, this->node_ptr(Prism15::side_nodes_map[i][n]));

  face->set_interior_parent(this);
  face->inherit_data_from(*this);

  return face;
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
            side = this->build_side_ptr(i);
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
            side = this->build_side_ptr(i);
            return;
          }
        break;
      }

    default:
      libmesh_error_msg("Invalid side i = " << i);
    }

  side->inherit_data_from(*this);

  // Set the nodes
  for (auto n : side->node_index_range())
    side->set_node(n, this->node_ptr(Prism15::side_nodes_map[i][n]));
}



std::unique_ptr<Elem> Prism15::build_edge_ptr (const unsigned int i)
{
  return this->simple_build_edge_ptr<Edge3,Prism15>(i);
}



void Prism15::build_edge_ptr (std::unique_ptr<Elem> & edge, const unsigned int i)
{
  this->simple_build_edge_ptr<Prism15>(edge, i, EDGE3);
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
  // This specialization is good for Lagrange mappings only
  if (this->mapping_type() != LAGRANGE_MAP)
    return this->Elem::volume();

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
    w1 = 1.1169079483900573284750350421656140e-01_R,
    w2 = 5.4975871827660933819163162450105264e-02_R,
    a1 = 4.4594849091596488631832925388305199e-01_R,
    a2 = 9.1576213509770743459571463402201508e-02_R;

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

const Real Prism15::_embedding_matrix[Prism15::num_children][Prism15::num_nodes][Prism15::num_nodes] =
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


void
Prism15::permute(unsigned int perm_num)
{
  libmesh_assert_less (perm_num, 6);
  const unsigned int side = perm_num % 2;
  const unsigned int rotate = perm_num / 2;

  for (unsigned int i = 0; i != rotate; ++i)
    {
      swap3nodes(0,1,2);
      swap3nodes(3,4,5);
      swap3nodes(6,7,8);
      swap3nodes(9,10,11);
      swap3nodes(12,13,14);
      swap3neighbors(1,2,3);
    }

  switch (side) {
  case 0:
    break;
  case 1:
    swap2nodes(1,3);
    swap2nodes(0,4);
    swap2nodes(2,5);
    swap2nodes(6,12);
    swap2nodes(9,10);
    swap2nodes(7,14);
    swap2nodes(8,13);
    swap2neighbors(0,4);
    swap2neighbors(2,3);
    break;
  default:
    libmesh_error();
  }

}


void
Prism15::flip(BoundaryInfo * boundary_info)
{
  libmesh_assert(boundary_info);

  swap2nodes(0,1);
  swap2nodes(3,4);
  swap2nodes(7,8);
  swap2nodes(9,10);
  swap2nodes(13,14);
  swap2neighbors(2,3);
  swap2boundarysides(2,3,boundary_info);
  swap2boundaryedges(0,1,boundary_info);
  swap2boundaryedges(3,4,boundary_info);
  swap2boundaryedges(7,8,boundary_info);
}


ElemType
Prism15::side_type (const unsigned int s) const
{
  libmesh_assert_less (s, 5);
  if (s == 0 || s == 4)
    return TRI6;
  return QUAD8;
}

} // namespace libMesh
