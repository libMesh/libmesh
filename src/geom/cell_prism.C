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


// C++ includes

// Local includes
#include "libmesh/cell_prism.h"
#include "libmesh/cell_prism6.h"
#include "libmesh/face_quad4.h"
#include "libmesh/face_tri3.h"

namespace libMesh
{


// ------------------------------------------------------------
// Prism class static member initializations
const int Prism::num_sides;
const int Prism::num_edges;
const int Prism::num_children;

const Real Prism::_master_points[18][3] =
  {
    {0, 0, -1},
    {1, 0, -1},
    {0, 1, -1},
    {0, 0, 1},
    {1, 0, 1},
    {0, 1, 1},
    {0.5, 0, -1},
    {0.5, 0.5, -1},
    {0, 0.5, -1},
    {0, 0, 0},
    {1, 0, 0},
    {0, 1, 0},
    {0.5, 0, 1},
    {0.5, 0.5, 1},
    {0, 0.5, 1},
    {0.5, 0, 0},
    {0.5, 0.5, 0},
    {0, 0.5, 0}
  };

const unsigned int Prism::edge_sides_map[9][2] =
  {
    {0, 1}, // Edge 0
    {0, 2}, // Edge 1
    {0, 3}, // Edge 2
    {1, 3}, // Edge 3
    {1, 2}, // Edge 4
    {2, 3}, // Edge 5
    {1, 4}, // Edge 6
    {2, 4}, // Edge 7
    {3, 4}  // Edge 8
  };

const unsigned int Prism::adjacent_edges_map[/*num_vertices*/6][/*n_adjacent_edges*/3] =
  {
    {0, 2, 3},  // Edges adjacent to node 0
    {0, 1, 4},  // Edges adjacent to node 1
    {1, 2, 5},  // Edges adjacent to node 2
    {3, 6, 8},  // Edges adjacent to node 3
    {4, 6, 7},  // Edges adjacent to node 4
    {5, 7, 8},  // Edges adjacent to node 5
  };


// ------------------------------------------------------------
// Prism class member functions
dof_id_type Prism::key (const unsigned int s) const
{
  libmesh_assert_less (s, this->n_sides());

  switch (s)
    {
    case 0: // the triangular face at z=0
    case 4: // the triangular face at z=1
      return this->compute_key (this->node_id(Prism6::side_nodes_map[s][0]),
                                this->node_id(Prism6::side_nodes_map[s][1]),
                                this->node_id(Prism6::side_nodes_map[s][2]));

    case 1: // the quad face at y=0
    case 2: // the other quad face
    case 3: // the quad face at x=0
      return this->compute_key (this->node_id(Prism6::side_nodes_map[s][0]),
                                this->node_id(Prism6::side_nodes_map[s][1]),
                                this->node_id(Prism6::side_nodes_map[s][2]),
                                this->node_id(Prism6::side_nodes_map[s][3]));

    default:
      libmesh_error_msg("Invalid side " << s);
    }
}



dof_id_type Prism::low_order_key (const unsigned int s) const
{
  libmesh_assert_less (s, this->n_sides());

  switch (s)
    {
    case 0: // the triangular face at z=0
    case 4: // the triangular face at z=1
      return this->compute_key (this->node_id(Prism6::side_nodes_map[s][0]),
                                this->node_id(Prism6::side_nodes_map[s][1]),
                                this->node_id(Prism6::side_nodes_map[s][2]));

    case 1: // the quad face at y=0
    case 2: // the other quad face
    case 3: // the quad face at x=0
      return this->compute_key (this->node_id(Prism6::side_nodes_map[s][0]),
                                this->node_id(Prism6::side_nodes_map[s][1]),
                                this->node_id(Prism6::side_nodes_map[s][2]),
                                this->node_id(Prism6::side_nodes_map[s][3]));

    default:
      libmesh_error_msg("Invalid side " << s);
    }
}



unsigned int Prism::local_side_node(unsigned int side,
                                    unsigned int side_node) const
{
  libmesh_assert_less (side, this->n_sides());

  // Never more than 4 nodes per side.
  libmesh_assert_less(side_node, Prism6::nodes_per_side);

  // Some sides have 3 nodes.
  libmesh_assert(!(side==0 || side==4) || side_node < 3);

  return Prism6::side_nodes_map[side][side_node];
}



unsigned int Prism::local_edge_node(unsigned int edge,
                                    unsigned int edge_node) const
{
  libmesh_assert_less(edge, this->n_edges());
  libmesh_assert_less(edge_node, Prism6::nodes_per_edge);

  return Prism6::edge_nodes_map[edge][edge_node];
}



std::unique_ptr<Elem> Prism::side_ptr (const unsigned int i)
{
  libmesh_assert_less (i, this->n_sides());

  std::unique_ptr<Elem> face;

  // Set up the type of element
  switch (i)
    {
    case 0: // the triangular face at z=0
    case 4: // the triangular face at z=1
      {
        face = std::make_unique<Tri3>();
        break;
      }
    case 1: // the quad face at y=0
    case 2: // the other quad face
    case 3: // the quad face at x=0
      {
        face = std::make_unique<Quad4>();
        break;
      }
    default:
      libmesh_error_msg("Invalid side i = " << i);
    }

  // Set the nodes
  for (auto n : face->node_index_range())
    face->set_node(n, this->node_ptr(Prism6::side_nodes_map[i][n]));

  return face;
}



void Prism::side_ptr (std::unique_ptr<Elem> & side,
                      const unsigned int i)
{
  libmesh_assert_less (i, this->n_sides());

  switch (i)
    {
      // the base face
    case 0: // the triangular face at z=0
    case 4: // the triangular face at z=1
      {
        if (!side.get() || side->type() != TRI3)
          {
            side = this->side_ptr(i);
            return;
          }
        break;
      }

    case 1: // the quad face at y=0
    case 2: // the other quad face
    case 3: // the quad face at x=0
      {
        if (!side.get() || side->type() != QUAD4)
          {
            side = this->side_ptr(i);
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
    side->set_node(n, this->node_ptr(Prism6::side_nodes_map[i][n]));
}



bool Prism::is_child_on_side(const unsigned int c,
                             const unsigned int s) const
{
  libmesh_assert_less (c, this->n_children());
  libmesh_assert_less (s, this->n_sides());

  for (unsigned int i = 0; i != 4; ++i)
    if (Prism6::side_elems_map[s][i] == c)
      return true;
  return false;
}



bool Prism::is_edge_on_side(const unsigned int e,
                            const unsigned int s) const
{
  libmesh_assert_less (e, this->n_edges());
  libmesh_assert_less (s, this->n_sides());

  return (edge_sides_map[e][0] == s || edge_sides_map[e][1] == s);
}



std::vector<unsigned int> Prism::sides_on_edge(const unsigned int e) const
{
  libmesh_assert_less(e, this->n_edges());

  return {edge_sides_map[e][0], edge_sides_map[e][1]};
}



unsigned int Prism::opposite_side(const unsigned int side_in) const
{
  libmesh_assert_less (side_in, 5);
  static const unsigned int prism_opposites[5] =
    {4, invalid_uint, invalid_uint, invalid_uint, 0};
  return prism_opposites[side_in];
}



bool
Prism::is_flipped() const
{
  return (triple_product(this->point(1)-this->point(0),
                        this->point(2)-this->point(0),
                        this->point(3)-this->point(0)) < 0);
}


std::vector<unsigned int>
Prism::edges_adjacent_to_node(const unsigned int n) const
{
  libmesh_assert_less(n, this->n_nodes());

  // For vertices, we use the Prism::adjacent_sides_map, otherwise each
  // of the mid-edge nodes is adjacent only to the edge it is on, and
  // face/internal nodes are not adjacent to any edge.
  if (this->is_vertex(n))
    return {std::begin(adjacent_edges_map[n]), std::end(adjacent_edges_map[n])};
  else if (this->is_edge(n))
    return {n - this->n_vertices()};

  libmesh_assert(this->is_face(n) || this->is_internal(n));
  return {};
}



bool Prism::on_reference_element(const Point & p,
                                 const Real eps) const
{
  const Real & xi = p(0);
  const Real & eta = p(1);
  const Real & zeta = p(2);

  // inside the reference triangle with zeta in [-1,1]
  return ((xi   >=  0.-eps) &&
          (eta  >=  0.-eps) &&
          (zeta >= -1.-eps) &&
          (zeta <=  1.+eps) &&
          ((xi + eta) <= 1.+eps));
}



const unsigned short int Prism::_second_order_vertex_child_number[18] =
  {
    99,99,99,99,99,99, // Vertices
    0,1,0,0,1,2,3,4,3, // Edges
    0,1,0              // Faces
  };



const unsigned short int Prism::_second_order_vertex_child_index[18] =
  {
    99,99,99,99,99,99, // Vertices
    1,2,2,3,4,5,4,5,5, // Edges
    4,5,5              // Faces
  };


const unsigned short int Prism::_second_order_adjacent_vertices[9][2] =
  {
    { 0,  1}, // vertices adjacent to node 6
    { 1,  2}, // vertices adjacent to node 7
    { 0,  2}, // vertices adjacent to node 8

    { 0,  3}, // vertices adjacent to node 9
    { 1,  4}, // vertices adjacent to node 10
    { 2,  5}, // vertices adjacent to node 11

    { 3,  4}, // vertices adjacent to node 12
    { 4,  5}, // vertices adjacent to node 13
    { 3,  5}  // vertices adjacent to node 14
  };

} // namespace libMesh
