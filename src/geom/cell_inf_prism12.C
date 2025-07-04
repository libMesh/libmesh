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

#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

// Local includes
#include "libmesh/cell_inf_prism12.h"
#include "libmesh/edge_edge3.h"
#include "libmesh/edge_inf_edge2.h"
#include "libmesh/face_tri6.h"
#include "libmesh/face_inf_quad6.h"
#include "libmesh/enum_io_package.h"
#include "libmesh/enum_order.h"

namespace libMesh
{


// ------------------------------------------------------------
// InfPrism12 class static member initializations
const int InfPrism12::num_nodes;
const int InfPrism12::nodes_per_side;
const int InfPrism12::nodes_per_edge;

const unsigned int InfPrism12::side_nodes_map[InfPrism12::num_sides][InfPrism12::nodes_per_side] =
  {
    { 0, 1, 2, 6, 7, 8},  // Side 0
    { 0, 1, 3, 4, 6, 9},  // Side 1
    { 1, 2, 4, 5, 7, 10}, // Side 2
    { 2, 0, 5, 3, 8, 11}  // Side 3
  };

const unsigned int InfPrism12::edge_nodes_map[InfPrism12::num_edges][InfPrism12::nodes_per_edge] =
  {
    {0, 1,  6}, // Side 0
    {1, 2,  7}, // Side 1
    {0, 2,  8}, // Side 2
    {0, 3, 99}, // Side 3
    {1, 4, 99}, // Side 4
    {2, 5, 99}  // Side 5
  };

// ------------------------------------------------------------
// InfPrism12 class member functions

bool InfPrism12::is_node_on_side(const unsigned int n,
                                 const unsigned int s) const
{
  libmesh_assert_less (s, n_sides());
  return std::find(std::begin(side_nodes_map[s]),
                   std::end(side_nodes_map[s]),
                   n) != std::end(side_nodes_map[s]);
}

std::vector<unsigned>
InfPrism12::nodes_on_side(const unsigned int s) const
{
  libmesh_assert_less(s, n_sides());
  return {std::begin(side_nodes_map[s]), std::end(side_nodes_map[s])};
}

std::vector<unsigned>
InfPrism12::nodes_on_edge(const unsigned int e) const
{
  libmesh_assert_less(e, n_edges());
  auto trim = (e < 3) ? 0 : 1;
  return {std::begin(edge_nodes_map[e]), std::end(edge_nodes_map[e]) - trim};
}

bool InfPrism12::is_node_on_edge(const unsigned int n,
                                 const unsigned int e) const
{
  libmesh_assert_less (e, n_edges());
  return std::find(std::begin(edge_nodes_map[e]),
                   std::end(edge_nodes_map[e]),
                   n) != std::end(edge_nodes_map[e]);
}

Order InfPrism12::default_order() const
{
  return SECOND;
}

unsigned int InfPrism12::local_side_node(unsigned int side,
                                         unsigned int side_node) const
{
  libmesh_assert_less (side, this->n_sides());

  // Never more than 6 nodes per side.
  libmesh_assert_less(side_node, InfPrism12::nodes_per_side);

  return InfPrism12::side_nodes_map[side][side_node];
}



unsigned int InfPrism12::local_edge_node(unsigned int edge,
                                         unsigned int edge_node) const
{
  libmesh_assert_less (edge, this->n_edges());

  // Never more than 3 nodes per edge.
  libmesh_assert_less(edge_node, InfPrism12::nodes_per_edge);

  // Some edges only have 2 nodes.
  libmesh_assert(edge < 3 || edge_node < 2);

  return InfPrism12::edge_nodes_map[edge][edge_node];
}



std::unique_ptr<Elem> InfPrism12::build_side_ptr (const unsigned int i)
{
  libmesh_assert_less (i, this->n_sides());

  std::unique_ptr<Elem> face;

  switch (i)
    {
    case 0: // the triangular face at z=-1, base face
      {
        face = std::make_unique<Tri6>();
        break;
      }

    case 1: // the quad face at y=0
    case 2: // the other quad face
    case 3: // the quad face at x=0
      {
        face = std::make_unique<InfQuad6>();
        break;
      }

    default:
      libmesh_error_msg("Invalid side i = " << i);
    }

  // Set the nodes
  for (auto n : face->node_index_range())
    face->set_node(n, this->node_ptr(InfPrism12::side_nodes_map[i][n]));

  face->set_interior_parent(this);
  face->inherit_data_from(*this);

  return face;
}


void InfPrism12::build_side_ptr (std::unique_ptr<Elem> & side,
                                 const unsigned int i)
{
  libmesh_assert_less (i, this->n_sides());

  switch (i)
    {
    case 0: // the triangular face at z=-1, base face
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
        if (!side.get() || side->type() != INFQUAD6)
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
    side->set_node(n, this->node_ptr(InfPrism12::side_nodes_map[i][n]));
}



std::unique_ptr<Elem> InfPrism12::build_edge_ptr (const unsigned int i)
{
  if (i < 3) // base edges
    return this->simple_build_edge_ptr<Edge3,InfPrism12>(i);

  // infinite edges
  return this->simple_build_edge_ptr<InfEdge2,InfPrism12>(i);
}



void InfPrism12::build_edge_ptr (std::unique_ptr<Elem> & edge,
                                 const unsigned int i)
{
  libmesh_assert_less (i, this->n_edges());

  switch (i)
    {
      // the base edges
    case 0:
    case 1:
    case 2:
      {
        if (!edge.get() || edge->type() != EDGE3)
          {
            edge = this->build_edge_ptr(i);
            return;
          }
        break;
      }

      // the infinite edges
    case 3:
    case 4:
    case 5:
      {
        if (!edge.get() || edge->type() != INFEDGE2)
          {
            edge = this->build_edge_ptr(i);
            return;
          }
        break;
      }

    default:
      libmesh_error_msg("Invalid edge i = " << i);
    }

  edge->inherit_data_from(*this);

  // Set the nodes
  for (auto n : edge->node_index_range())
    edge->set_node(n, this->node_ptr(InfPrism12::edge_nodes_map[i][n]));
}



void InfPrism12::connectivity(const unsigned int sc,
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

            // guess this is a collapsed hex8
            conn[0] = this->node_id(0)+1;
            conn[1] = this->node_id(6)+1;
            conn[2] = this->node_id(8)+1;
            conn[3] = this->node_id(8)+1;
            conn[4] = this->node_id(3)+1;
            conn[5] = this->node_id(9)+1;
            conn[6] = this->node_id(11)+1;
            conn[7] = this->node_id(11)+1;

            return;

          case 1:

            conn[0] = this->node_id(6)+1;
            conn[1] = this->node_id(7)+1;
            conn[2] = this->node_id(8)+1;
            conn[3] = this->node_id(8)+1;
            conn[4] = this->node_id(9)+1;
            conn[5] = this->node_id(10)+1;
            conn[6] = this->node_id(11)+1;
            conn[7] = this->node_id(11)+1;

            return;

          case 2:

            conn[0] = this->node_id(6)+1;
            conn[1] = this->node_id(1)+1;
            conn[2] = this->node_id(7)+1;
            conn[3] = this->node_id(7)+1;
            conn[4] = this->node_id(9)+1;
            conn[5] = this->node_id(4)+1;
            conn[6] = this->node_id(10)+1;
            conn[7] = this->node_id(10)+1;

            return;

          case 3:

            conn[0] = this->node_id(8)+1;
            conn[1] = this->node_id(7)+1;
            conn[2] = this->node_id(2)+1;
            conn[3] = this->node_id(2)+1;
            conn[4] = this->node_id(11)+1;
            conn[5] = this->node_id(10)+1;
            conn[6] = this->node_id(5)+1;
            conn[7] = this->node_id(5)+1;

            return;

          default:
            libmesh_error_msg("Invalid sc = " << sc);
          }
      }

    default:
      libmesh_error_msg("Unsupported IO package " << iop);
    }
}





unsigned short int InfPrism12::second_order_adjacent_vertex (const unsigned int n,
                                                             const unsigned int v) const
{
  libmesh_assert_greater_equal (n, this->n_vertices());
  libmesh_assert_less (n, this->n_nodes());
  libmesh_assert_less (v, 2);
  return _second_order_adjacent_vertices[n-this->n_vertices()][v];
}



const unsigned short int InfPrism12::_second_order_adjacent_vertices[InfPrism12::num_edges][2] =
  {
    { 0,  1}, // vertices adjacent to node 6
    { 1,  2}, // vertices adjacent to node 7
    { 0,  2}, // vertices adjacent to node 8

    { 3,  4}, // vertices adjacent to node 9
    { 4,  5}, // vertices adjacent to node 10
    { 3,  5}  // vertices adjacent to node 11
  };



std::pair<unsigned short int, unsigned short int>
InfPrism12::second_order_child_vertex (const unsigned int n) const
{
  libmesh_assert_greater_equal (n, this->n_vertices());
  libmesh_assert_less (n, this->n_nodes());

  return std::pair<unsigned short int, unsigned short int>
    (_second_order_vertex_child_number[n],
     _second_order_vertex_child_index[n]);
}



const unsigned short int InfPrism12::_second_order_vertex_child_number[InfPrism12::num_nodes] =
  {
    99,99,99,99,99,99, // Vertices
    0,1,0,             // Edges
    0,1,0              // Faces
  };



const unsigned short int InfPrism12::_second_order_vertex_child_index[InfPrism12::num_nodes] =
  {
    99,99,99,99,99,99, // Vertices
    1,2,2,             // Edges
    4,5,5              // Faces
  };



#ifdef LIBMESH_ENABLE_AMR

const Real InfPrism12::_embedding_matrix[InfPrism12::num_children][InfPrism12::num_nodes][InfPrism12::num_nodes] =
  {
    // embedding matrix for child 0
    {
      //          0           1           2           3           4           5           6           7           8           9          10          11 th parent Node
      {         1.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 0th child N.
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 1
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0}, // 2
      {         0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 3
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0}, // 4
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0}, // 5
      {       0.375,     -0.125,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0,        0.0,        0.0,        0.0}, // 6
      {         0.0,     -0.125,     -0.125,        0.0,        0.0,        0.0,        0.5,       0.25,        0.5,        0.0,        0.0,        0.0}, // 7
      {       0.375,        0.0,     -0.125,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0,        0.0}, // 8
      {         0.0,        0.0,        0.0,      0.375,     -0.125,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0}, // 9
      {         0.0,        0.0,        0.0,        0.0,     -0.125,     -0.125,        0.0,        0.0,        0.0,        0.5,       0.25,        0.5}, // 10
      {         0.0,        0.0,        0.0,      0.375,        0.0,     -0.125,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75}  // 11
    },

    // embedding matrix for child 1
    {
      //          0           1           2           3           4           5           6           7           8           9          10          11 th parent Node
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 0th child N.
      {         0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 1
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0}, // 2
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0}, // 3
      {         0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 4
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0}, // 5
      {      -0.125,      0.375,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0,        0.0,        0.0,        0.0}, // 6
      {         0.0,      0.375,     -0.125,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0,        0.0,        0.0}, // 7
      {      -0.125,        0.0,     -0.125,        0.0,        0.0,        0.0,        0.5,        0.5,       0.25,        0.0,        0.0,        0.0}, // 8
      {         0.0,        0.0,        0.0,     -0.125,      0.375,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0}, // 9
      {         0.0,        0.0,        0.0,        0.0,      0.375,     -0.125,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0}, // 10
      {         0.0,        0.0,        0.0,     -0.125,        0.0,     -0.125,        0.0,        0.0,        0.0,        0.5,        0.5,       0.25}  // 11
    },

    // embedding matrix for child 2
    {
      //          0           1           2           3           4           5           6           7           8           9          10          11 th parent Node
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0}, // 0th child N.
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0}, // 1
      {         0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 2
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0}, // 3
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0}, // 4
      {         0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 5
      {      -0.125,     -0.125,        0.0,        0.0,        0.0,        0.0,       0.25,        0.5,        0.5,        0.0,        0.0,        0.0}, // 6
      {         0.0,     -0.125,      0.375,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0,        0.0,        0.0}, // 7
      {      -0.125,        0.0,      0.375,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0,        0.0}, // 8
      {         0.0,        0.0,        0.0,     -0.125,     -0.125,        0.0,        0.0,        0.0,        0.0,       0.25,        0.5,        0.5}, // 9
      {         0.0,        0.0,        0.0,        0.0,     -0.125,      0.375,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0}, // 10
      {         0.0,        0.0,        0.0,     -0.125,        0.0,      0.375,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75}  // 11
    },

    // embedding matrix for child 3
    {
      //          0           1           2           3           4           5           6           7           8           9          10          11 th parent Node
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 0th child N.
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0}, // 1
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0}, // 2
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0}, // 3
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0}, // 4
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0}, // 5
      {      -0.125,        0.0,     -0.125,        0.0,        0.0,        0.0,        0.5,        0.5,       0.25,        0.0,        0.0,        0.0}, // 6
      {      -0.125,     -0.125,        0.0,        0.0,        0.0,        0.0,       0.25,        0.5,        0.5,        0.0,        0.0,        0.0}, // 7
      {         0.0,     -0.125,     -0.125,        0.0,        0.0,        0.0,        0.5,       0.25,        0.5,        0.0,        0.0,        0.0}, // 8
      {         0.0,        0.0,        0.0,     -0.125,        0.0,     -0.125,        0.0,        0.0,        0.0,        0.5,        0.5,       0.25}, // 9
      {         0.0,        0.0,        0.0,     -0.125,     -0.125,        0.0,        0.0,        0.0,        0.0,       0.25,        0.5,        0.5}, // 10
      {         0.0,        0.0,        0.0,        0.0,     -0.125,     -0.125,        0.0,        0.0,        0.0,        0.5,       0.25,        0.5}  // 11
    }

  };


#endif


void
InfPrism12::permute(unsigned int perm_num)
{
  libmesh_assert_less (perm_num, 3);

  for (unsigned int i = 0; i != perm_num; ++i)
    {
      swap3nodes(0,1,2);
      swap3nodes(3,4,5);
      swap3nodes(6,7,8);
      swap3nodes(9,10,11);
      swap3neighbors(1,2,3);
    }
}


void
InfPrism12::flip(BoundaryInfo * boundary_info)
{
  libmesh_assert(boundary_info);

  swap2nodes(0,1);
  swap2nodes(3,4);
  swap2nodes(7,8);
  swap2nodes(10,11);
  swap2neighbors(2,3);
  swap2boundarysides(2,3,boundary_info);
  swap2boundaryedges(3,4,boundary_info);
}


ElemType
InfPrism12::side_type (const unsigned int s) const
{
  libmesh_assert_less (s, 4);
  if (s == 0)
    return TRI6;
  return INFQUAD6;
}


} // namespace libMesh

#endif // ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
