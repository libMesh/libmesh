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
#include "libmesh/cell_pyramid14.h"
#include "libmesh/edge_edge3.h"
#include "libmesh/face_tri6.h"
#include "libmesh/face_quad9.h"
#include "libmesh/enum_io_package.h"
#include "libmesh/enum_order.h"

namespace libMesh
{




// ------------------------------------------------------------
// Pyramid14 class static member initializations
const int Pyramid14::num_nodes;
const int Pyramid14::nodes_per_side;
const int Pyramid14::nodes_per_edge;

const unsigned int Pyramid14::side_nodes_map[Pyramid14::num_sides][Pyramid14::nodes_per_side] =
  {
    {0, 1, 4, 5, 10,  9, 99, 99, 99}, // Side 0 (front)
    {1, 2, 4, 6, 11, 10, 99, 99, 99}, // Side 1 (right)
    {2, 3, 4, 7, 12, 11, 99, 99, 99}, // Side 2 (back)
    {3, 0, 4, 8,  9, 12, 99, 99, 99}, // Side 3 (left)
    {0, 3, 2, 1,  8,  7,  6,  5, 13}  // Side 4 (base)
  };

const unsigned int Pyramid14::edge_nodes_map[Pyramid14::num_edges][Pyramid14::nodes_per_edge] =
  {
    {0, 1,  5}, // Edge 0
    {1, 2,  6}, // Edge 1
    {2, 3,  7}, // Edge 2
    {0, 3,  8}, // Edge 3
    {0, 4,  9}, // Edge 4
    {1, 4, 10}, // Edge 5
    {2, 4, 11}, // Edge 6
    {3, 4, 12}  // Edge 7
  };

// ------------------------------------------------------------
// Pyramid14 class member functions

bool Pyramid14::is_vertex(const unsigned int i) const
{
  if (i < 5)
    return true;
  return false;
}



bool Pyramid14::is_edge(const unsigned int i) const
{
  if (i < 5)
    return false;
  if (i == 13)
    return false;
  return true;
}



bool Pyramid14::is_face(const unsigned int i) const
{
  if (i == 13)
    return true;
  return false;
}



bool Pyramid14::is_node_on_side(const unsigned int n,
                                const unsigned int s) const
{
  libmesh_assert_less (s, n_sides());
  return std::find(std::begin(side_nodes_map[s]),
                   std::end(side_nodes_map[s]),
                   n) != std::end(side_nodes_map[s]);
}

std::vector<unsigned>
Pyramid14::nodes_on_side(const unsigned int s) const
{
  libmesh_assert_less(s, n_sides());
  auto trim = (s == 4) ? 0 : 3;
  return {std::begin(side_nodes_map[s]), std::end(side_nodes_map[s]) - trim};
}

std::vector<unsigned>
Pyramid14::nodes_on_edge(const unsigned int e) const
{
  libmesh_assert_less(e, n_edges());
  return {std::begin(edge_nodes_map[e]), std::end(edge_nodes_map[e])};
}

bool Pyramid14::is_node_on_edge(const unsigned int n,
                                const unsigned int e) const
{
  libmesh_assert_less (e, n_edges());
  return std::find(std::begin(edge_nodes_map[e]),
                   std::end(edge_nodes_map[e]),
                   n) != std::end(edge_nodes_map[e]);
}

bool Pyramid14::has_affine_map() const
{
  // TODO: If the base is a parallelogram and all the triangular faces are planar,
  // the map should be linear, but I need to test this theory...
  return false;
}



Order Pyramid14::default_order() const
{
  return SECOND;
}



dof_id_type Pyramid14::key (const unsigned int s) const
{
  libmesh_assert_less (s, this->n_sides());

  switch (s)
    {
    case 0: // triangular face 1
    case 1: // triangular face 2
    case 2: // triangular face 3
    case 3: // triangular face 4
      return Pyramid::key(s);

    case 4:  // the quad face at z=0
      return this->compute_key (this->node_id(13));

    default:
      libmesh_error_msg("Invalid side s = " << s);
    }
}



unsigned int Pyramid14::local_side_node(unsigned int side,
                                        unsigned int side_node) const
{
  libmesh_assert_less (side, this->n_sides());

  // Never more than 9 nodes per side.
  libmesh_assert_less(side_node, Pyramid14::nodes_per_side);

  // Some sides have 6 nodes.
  libmesh_assert(side == 4 || side_node < 6);

  return Pyramid14::side_nodes_map[side][side_node];
}



unsigned int Pyramid14::local_edge_node(unsigned int edge,
                                        unsigned int edge_node) const
{
  libmesh_assert_less(edge, this->n_edges());
  libmesh_assert_less(edge_node, Pyramid14::nodes_per_edge);

  return Pyramid14::edge_nodes_map[edge][edge_node];
}



std::unique_ptr<Elem> Pyramid14::build_side_ptr (const unsigned int i)
{
  libmesh_assert_less (i, this->n_sides());

  std::unique_ptr<Elem> face;

  switch (i)
    {
    case 0: // triangular face 1
    case 1: // triangular face 2
    case 2: // triangular face 3
    case 3: // triangular face 4
      {
        face = std::make_unique<Tri6>();
        break;
      }
    case 4: // the quad face at z=0
      {
        face = std::make_unique<Quad9>();
        break;
      }
    default:
      libmesh_error_msg("Invalid side i = " << i);
    }

  // Set the nodes
  for (auto n : face->node_index_range())
    face->set_node(n, this->node_ptr(Pyramid14::side_nodes_map[i][n]));

  face->set_interior_parent(this);
  face->inherit_data_from(*this);

  return face;
}



void Pyramid14::build_side_ptr (std::unique_ptr<Elem> & side,
                                const unsigned int i)
{
  libmesh_assert_less (i, this->n_sides());

  switch (i)
    {
    case 0: // triangular face 1
    case 1: // triangular face 2
    case 2: // triangular face 3
    case 3: // triangular face 4
      {
        if (!side.get() || side->type() != TRI6)
          {
            side = this->build_side_ptr(i);
            return;
          }
        break;
      }
    case 4:  // the quad face at z=0
      {
        if (!side.get() || side->type() != QUAD9)
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
    side->set_node(n, this->node_ptr(Pyramid14::side_nodes_map[i][n]));
}



std::unique_ptr<Elem> Pyramid14::build_edge_ptr (const unsigned int i)
{
  return this->simple_build_edge_ptr<Edge3,Pyramid14>(i);
}



void Pyramid14::build_edge_ptr (std::unique_ptr<Elem> & edge, const unsigned int i)
{
  this->simple_build_edge_ptr<Pyramid14>(edge, i, EDGE3);
}



void Pyramid14::connectivity(const unsigned int libmesh_dbg_var(sc),
                             const IOPackage iop,
                             std::vector<dof_id_type> & /*conn*/) const
{
  libmesh_assert(_nodes);
  libmesh_assert_less (sc, this->n_sub_elem());
  libmesh_assert_not_equal_to (iop, INVALID_IO_PACKAGE);

  switch (iop)
    {
    case TECPLOT:
      {
        // TODO
        libmesh_not_implemented();
      }

    case VTK:
      {
        // TODO
        libmesh_not_implemented();
      }

    default:
      libmesh_error_msg("Unsupported IO package " << iop);
    }
}



unsigned int Pyramid14::n_second_order_adjacent_vertices (const unsigned int n) const
{
  switch (n)
    {
    case 5:
    case 6:
    case 7:
    case 8:
    case 9:
    case 10:
    case 11:
    case 12:
      return 2;

    case 13:
      return 4;

    default:
      libmesh_error_msg("Invalid node n = " << n);
    }
}


unsigned short int Pyramid14::second_order_adjacent_vertex (const unsigned int n,
                                                            const unsigned int v) const
{
  libmesh_assert_greater_equal (n, this->n_vertices());
  libmesh_assert_less (n, this->n_nodes());

  switch (n)
    {
    case 5:
    case 6:
    case 7:
    case 8:
    case 9:
    case 10:
    case 11:
    case 12:
      {
        libmesh_assert_less (v, 2);

        // This is the analog of the static, const arrays
        // {Hex,Prism,Tet10}::_second_order_adjacent_vertices
        // defined in the respective source files... possibly treat
        // this similarly once the Pyramid13 has been added?
        unsigned short node_list[8][2] =
          {
            {0,1},
            {1,2},
            {2,3},
            {0,3},
            {0,4},
            {1,4},
            {2,4},
            {3,4}
          };

        return node_list[n-5][v];
      }

      // mid-face node on bottom
    case 13:
      {
        libmesh_assert_less (v, 4);

        // The vertex nodes surrounding node 13 are 0, 1, 2, and 3.
        // Thus, the v'th node is simply = v.
        return cast_int<unsigned short>(v);
      }

    default:
      libmesh_error_msg("Invalid n = " << n);

    }
}



Real Pyramid14::volume () const
{
  // This specialization is good for Lagrange mappings only
  if (this->mapping_type() != LAGRANGE_MAP)
    return this->Elem::volume();

  // Make copies of our points.  It makes the subsequent calculations a bit
  // shorter and avoids dereferencing the same pointer multiple times.
  Point
    x0 = point(0), x1 = point(1), x2 = point(2), x3 = point(3),   x4 = point(4),   x5 = point(5),   x6 = point(6),
    x7 = point(7), x8 = point(8), x9 = point(9), x10 = point(10), x11 = point(11), x12 = point(12), x13 = point(13);

  // dx/dxi and dx/deta have 15 components while dx/dzeta has 20.
  // These are copied directly from the output of a Python script.
  Point dx_dxi[15] =
    {
      x6/2 - x8/2,
      x0/4 - x1/4 + x10 + x11 - x12 - x2/4 + x3/4 - 3*x6/2 + 3*x8/2 - x9,
      -x0/2 + x1/2 - 2*x10 - 2*x11 + 2*x12 + x2/2 - x3/2 + 3*x6/2 - 3*x8/2 + 2*x9,
      x0/4 - x1/4 + x10 + x11 - x12 - x2/4 + x3/4 - x6/2 + x8/2 - x9,
      x0/4 - x1/4 + x2/4 - x3/4,
      -3*x0/4 + 3*x1/4 - x10 + x11 - x12 - 3*x2/4 + 3*x3/4 + x9,
      x0/2 - x1/2 + x10 - x11 + x12 + x2/2 - x3/2 - x9,
      -x0/4 + x1/4 + x2/4 - x3/4 - x6/2 + x8/2,
      x0/4 - x1/4 - x2/4 + x3/4 + x6/2 - x8/2,
      -2*x13 + x6 + x8,
      4*x13 - 2*x6 - 2*x8,
      -2*x13 + x6 + x8,
      -x0/2 - x1/2 + x2/2 + x3/2 + x5 - x7,
      x0/2 + x1/2 - x2/2 - x3/2 - x5 + x7,
      x0/2 + x1/2 + 2*x13 + x2/2 + x3/2 - x5 - x6 - x7 - x8
    };

  // dx/dxi and dx/deta have 15 components while dx/dzeta has 20.
  // These are copied directly from the output of a Python script.
  Point dx_deta[15] =
    {
      -x5/2 + x7/2,
      x0/4 + x1/4 - x10 + x11 + x12 - x2/4 - x3/4 + 3*x5/2 - 3*x7/2 - x9,
      -x0/2 - x1/2 + 2*x10 - 2*x11 - 2*x12 + x2/2 + x3/2 - 3*x5/2 + 3*x7/2 + 2*x9,
      x0/4 + x1/4 - x10 + x11 + x12 - x2/4 - x3/4 + x5/2 - x7/2 - x9,
      -2*x13 + x5 + x7,
      4*x13 - 2*x5 - 2*x7,
      -2*x13 + x5 + x7,
      x0/4 - x1/4 + x2/4 - x3/4,
      -3*x0/4 + 3*x1/4 - x10 + x11 - x12 - 3*x2/4 + 3*x3/4 + x9,
      x0/2 - x1/2 + x10 - x11 + x12 + x2/2 - x3/2 - x9,
      -x0/2 + x1/2 + x2/2 - x3/2 - x6 + x8,
      x0/2 - x1/2 - x2/2 + x3/2 + x6 - x8,
      -x0/4 - x1/4 + x2/4 + x3/4 + x5/2 - x7/2,
      x0/4 + x1/4 - x2/4 - x3/4 - x5/2 + x7/2,
      x0/2 + x1/2 + 2*x13 + x2/2 + x3/2 - x5 - x6 - x7 - x8
    };

  // dx/dxi and dx/deta have 15 components while dx/dzeta has 20.
  // These are copied directly from the output of a Python script.
  Point dx_dzeta[20] =
    {
      -x0/4 - x1/4 + x10 + x11 + x12 - 2*x13 - x2/4 - x3/4 - x4 + x9,
      5*x0/4 + 5*x1/4 - 5*x10 - 5*x11 - 5*x12 + 8*x13 + 5*x2/4 + 5*x3/4 + 7*x4 - 5*x9,
      -9*x0/4 - 9*x1/4 + 9*x10 + 9*x11 + 9*x12 - 12*x13 - 9*x2/4 - 9*x3/4 - 15*x4 + 9*x9,
      7*x0/4 + 7*x1/4 - 7*x10 - 7*x11 - 7*x12 + 8*x13 + 7*x2/4 + 7*x3/4 + 13*x4 - 7*x9,
      -x0/2 - x1/2 + 2*x10 + 2*x11 + 2*x12 - 2*x13 - x2/2 - x3/2 - 4*x4 + 2*x9,
      x0/4 + x1/4 - x10 + x11 + x12 - x2/4 - x3/4 + x5/2 - x7/2 - x9,
      -3*x0/4 - 3*x1/4 + 3*x10 - 3*x11 - 3*x12 + 3*x2/4 + 3*x3/4 - 3*x5/2 + 3*x7/2 + 3*x9,
      3*x0/4 + 3*x1/4 - 3*x10 + 3*x11 + 3*x12 - 3*x2/4 - 3*x3/4 + 3*x5/2 - 3*x7/2 - 3*x9,
      -x0/4 - x1/4 + x10 - x11 - x12 + x2/4 + x3/4 - x5/2 + x7/2 + x9,
      x0/4 - x1/4 + x10 + x11 - x12 - x2/4 + x3/4 - x6/2 + x8/2 - x9,
      -3*x0/4 + 3*x1/4 - 3*x10 - 3*x11 + 3*x12 + 3*x2/4 - 3*x3/4 + 3*x6/2 - 3*x8/2 + 3*x9,
      3*x0/4 - 3*x1/4 + 3*x10 + 3*x11 - 3*x12 - 3*x2/4 + 3*x3/4 - 3*x6/2 + 3*x8/2 - 3*x9,
      -x0/4 + x1/4 - x10 - x11 + x12 + x2/4 - x3/4 + x6/2 - x8/2 + x9,
      -x0/4 + x1/4 - x10 + x11 - x12 - x2/4 + x3/4 + x9,
      x0/4 - x1/4 + x10 - x11 + x12 + x2/4 - x3/4 - x9,
      -x0/4 + x1/4 + x2/4 - x3/4 - x6/2 + x8/2,
      x0/4 - x1/4 - x2/4 + x3/4 + x6/2 - x8/2,
      -x0/4 - x1/4 + x2/4 + x3/4 + x5/2 - x7/2,
      x0/4 + x1/4 - x2/4 - x3/4 - x5/2 + x7/2,
      x0/2 + x1/2 + 2*x13 + x2/2 + x3/2 - x5 - x6 - x7 - x8
    };

  // The (xi, eta, zeta) exponents for each of the dx_dxi terms
  static const int dx_dxi_exponents[15][3] =
    {
      {0, 0, 0},
      {0, 0, 1},
      {0, 0, 2},
      {0, 0, 3},
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
      {1, 2, 0}
    };

  // The (xi, eta, zeta) exponents for each of the dx_deta terms
  static const int dx_deta_exponents[15][3] =
    {
      {0, 0, 0},
      {0, 0, 1},
      {0, 0, 2},
      {0, 0, 3},
      {0, 1, 0},
      {0, 1, 1},
      {0, 1, 2},
      {1, 0, 0},
      {1, 0, 1},
      {1, 0, 2},
      {1, 1, 0},
      {1, 1, 1},
      {2, 0, 0},
      {2, 0, 1},
      {2, 1, 0}
    };

  // The (xi, eta, zeta) exponents for each of the dx_dzeta terms
  static const int dx_dzeta_exponents[20][3] =
    {
      {0, 0, 0},
      {0, 0, 1},
      {0, 0, 2},
      {0, 0, 3},
      {0, 0, 4},
      {0, 1, 0},
      {0, 1, 1},
      {0, 1, 2},
      {0, 1, 3},
      {1, 0, 0},
      {1, 0, 1},
      {1, 0, 2},
      {1, 0, 3},
      {1, 1, 0},
      {1, 1, 1},
      {1, 2, 0},
      {1, 2, 1},
      {2, 1, 0},
      {2, 1, 1},
      {2, 2, 0},
    };

  // Number of points in the quadrature rule
  const int N = 27;

  // Parameters of the quadrature rule
  static const Real
    // Parameters used for (xi, eta) quadrature points.
    a1 = -7.1805574131988893873307823958101e-01_R,
    a2 = -5.0580870785392503961340276902425e-01_R,
    a3 = -2.2850430565396735359574351631314e-01_R,
    // Parameters used for zeta quadrature points.
    b1 = 7.2994024073149732155837979012003e-02_R,
    b2 = 3.4700376603835188472176354340395e-01_R,
    b3 = 7.0500220988849838312239847758405e-01_R,
    // There are 9 unique weight values since there are three
    // for each of the three unique zeta values.
    w1 = 4.8498876871878584357513834016440e-02_R,
    w2 = 4.5137737425884574692441981593901e-02_R,
    w3 = 9.2440441384508327195915094925393e-03_R,
    w4 = 7.7598202995005734972022134426305e-02_R,
    w5 = 7.2220379881415319507907170550242e-02_R,
    w6 = 1.4790470621521332351346415188063e-02_R,
    w7 = 1.2415712479200917595523541508209e-01_R,
    w8 = 1.1555260781026451121265147288039e-01_R,
    w9 = 2.3664752994434131762154264300901e-02_R;

  // The points and weights of the 3x3x3 quadrature rule
  static const Real xi[N][3] =
    {// ^0   ^1  ^2
      { 1.,  a1, a1*a1},
      { 1.,  a2, a2*a2},
      { 1.,  a3, a3*a3},
      { 1.,  a1, a1*a1},
      { 1.,  a2, a2*a2},
      { 1.,  a3, a3*a3},
      { 1.,  a1, a1*a1},
      { 1.,  a2, a2*a2},
      { 1.,  a3, a3*a3},
      { 1.,  0., 0.   },
      { 1.,  0., 0.   },
      { 1.,  0., 0.   },
      { 1.,  0., 0.   },
      { 1.,  0., 0.   },
      { 1.,  0., 0.   },
      { 1.,  0., 0.   },
      { 1.,  0., 0.   },
      { 1.,  0., 0.   },
      { 1., -a1, a1*a1},
      { 1., -a2, a2*a2},
      { 1., -a3, a3*a3},
      { 1., -a1, a1*a1},
      { 1., -a2, a2*a2},
      { 1., -a3, a3*a3},
      { 1., -a1, a1*a1},
      { 1., -a2, a2*a2},
      { 1., -a3, a3*a3}
    };

  static const Real eta[N][3] =
    {// ^0   ^1  ^2
      { 1.,  a1, a1*a1},
      { 1.,  a2, a2*a2},
      { 1.,  a3, a3*a3},
      { 1.,  0., 0.   },
      { 1.,  0., 0.   },
      { 1.,  0., 0.   },
      { 1., -a1, a1*a1},
      { 1., -a2, a2*a2},
      { 1., -a3, a3*a3},
      { 1.,  a1, a1*a1},
      { 1.,  a2, a2*a2},
      { 1.,  a3, a3*a3},
      { 1.,  0., 0.   },
      { 1.,  0., 0.   },
      { 1.,  0., 0.   },
      { 1., -a1, a1*a1},
      { 1., -a2, a2*a2},
      { 1., -a3, a3*a3},
      { 1.,  a1, a1*a1},
      { 1.,  a2, a2*a2},
      { 1.,  a3, a3*a3},
      { 1.,  0., 0.   },
      { 1.,  0., 0.   },
      { 1.,  0., 0.   },
      { 1., -a1, a1*a1},
      { 1., -a2, a2*a2},
      { 1., -a3, a3*a3}
    };

  static const Real zeta[N][5] =
    {// ^0  ^1  ^2     ^3        ^4
      { 1., b1, b1*b1, b1*b1*b1, b1*b1*b1*b1},
      { 1., b2, b2*b2, b2*b2*b2, b2*b2*b2*b2},
      { 1., b3, b3*b3, b3*b3*b3, b3*b3*b3*b3},
      { 1., b1, b1*b1, b1*b1*b1, b1*b1*b1*b1},
      { 1., b2, b2*b2, b2*b2*b2, b2*b2*b2*b2},
      { 1., b3, b3*b3, b3*b3*b3, b3*b3*b3*b3},
      { 1., b1, b1*b1, b1*b1*b1, b1*b1*b1*b1},
      { 1., b2, b2*b2, b2*b2*b2, b2*b2*b2*b2},
      { 1., b3, b3*b3, b3*b3*b3, b3*b3*b3*b3},
      { 1., b1, b1*b1, b1*b1*b1, b1*b1*b1*b1},
      { 1., b2, b2*b2, b2*b2*b2, b2*b2*b2*b2},
      { 1., b3, b3*b3, b3*b3*b3, b3*b3*b3*b3},
      { 1., b1, b1*b1, b1*b1*b1, b1*b1*b1*b1},
      { 1., b2, b2*b2, b2*b2*b2, b2*b2*b2*b2},
      { 1., b3, b3*b3, b3*b3*b3, b3*b3*b3*b3},
      { 1., b1, b1*b1, b1*b1*b1, b1*b1*b1*b1},
      { 1., b2, b2*b2, b2*b2*b2, b2*b2*b2*b2},
      { 1., b3, b3*b3, b3*b3*b3, b3*b3*b3*b3},
      { 1., b1, b1*b1, b1*b1*b1, b1*b1*b1*b1},
      { 1., b2, b2*b2, b2*b2*b2, b2*b2*b2*b2},
      { 1., b3, b3*b3, b3*b3*b3, b3*b3*b3*b3},
      { 1., b1, b1*b1, b1*b1*b1, b1*b1*b1*b1},
      { 1., b2, b2*b2, b2*b2*b2, b2*b2*b2*b2},
      { 1., b3, b3*b3, b3*b3*b3, b3*b3*b3*b3},
      { 1., b1, b1*b1, b1*b1*b1, b1*b1*b1*b1},
      { 1., b2, b2*b2, b2*b2*b2, b2*b2*b2*b2},
      { 1., b3, b3*b3, b3*b3*b3, b3*b3*b3*b3}
    };

  static const Real w[N] = {w1, w2, w3, w4, w5, w6, // 0-5
                            w1, w2, w3, w4, w5, w6, // 6-11
                            w7, w8, w9, w4, w5, w6, // 12-17
                            w1, w2, w3, w4, w5, w6, // 18-23
                            w1, w2, w3};            // 24-26

  Real vol = 0.;
  for (int q=0; q<N; ++q)
    {
      // Compute denominators for the current q.
      Real
        den2 = (1. - zeta[q][1])*(1. - zeta[q][1]),
        den3 = den2*(1. - zeta[q][1]);

      // Compute dx/dxi and dx/deta at the current q.
      Point dx_dxi_q, dx_deta_q;
      for (int c=0; c<15; ++c)
        {
          dx_dxi_q +=
            xi[q][dx_dxi_exponents[c][0]]*
            eta[q][dx_dxi_exponents[c][1]]*
            zeta[q][dx_dxi_exponents[c][2]]*dx_dxi[c];

          dx_deta_q +=
            xi[q][dx_deta_exponents[c][0]]*
            eta[q][dx_deta_exponents[c][1]]*
            zeta[q][dx_deta_exponents[c][2]]*dx_deta[c];
        }

      // Compute dx/dzeta at the current q.
      Point dx_dzeta_q;
      for (int c=0; c<20; ++c)
        {
          dx_dzeta_q +=
            xi[q][dx_dzeta_exponents[c][0]]*
            eta[q][dx_dzeta_exponents[c][1]]*
            zeta[q][dx_dzeta_exponents[c][2]]*dx_dzeta[c];
        }

      // Scale everything appropriately
      dx_dxi_q /= den2;
      dx_deta_q /= den2;
      dx_dzeta_q /= den3;

      // Compute scalar triple product, multiply by weight, and accumulate volume.
      vol += w[q] * triple_product(dx_dxi_q, dx_deta_q, dx_dzeta_q);
    }

  return vol;
}


void Pyramid14::permute(unsigned int perm_num)
{
  libmesh_assert_less (perm_num, 4);

  for (unsigned int i = 0; i != perm_num; ++i)
    {
      swap4nodes(0,1,2,3);
      swap4nodes(5,6,7,8);
      swap4nodes(9,10,11,12);
      swap4neighbors(0,1,2,3);
    }
}


void Pyramid14::flip(BoundaryInfo * boundary_info)
{
  libmesh_assert(boundary_info);

  swap2nodes(0,1);
  swap2nodes(2,3);
  swap2nodes(6,8);
  swap2nodes(9,10);
  swap2nodes(11,12);
  swap2neighbors(1,3);
  swap2boundarysides(1,3,boundary_info);
  swap2boundaryedges(1,3,boundary_info);
  swap2boundaryedges(4,5,boundary_info);
  swap2boundaryedges(6,7,boundary_info);
}


unsigned int Pyramid14::center_node_on_side(const unsigned short side) const
{
  libmesh_assert_less (side, Pyramid14::num_sides);
  return side == 4 ? 13 : invalid_uint;
}


ElemType Pyramid14::side_type (const unsigned int s) const
{
  libmesh_assert_less (s, 5);
  if (s < 4)
    return TRI6;
  return QUAD9;
}


} // namespace libMesh
