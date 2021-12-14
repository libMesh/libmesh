// The libMesh Finite Element Library.
// Copyright (C) 2002-2021 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

// Local includes cont'd
#include "libmesh/cell_inf_hex16.h"
#include "libmesh/edge_edge3.h"
#include "libmesh/edge_inf_edge2.h"
#include "libmesh/face_quad8.h"
#include "libmesh/face_inf_quad6.h"
#include "libmesh/side.h"
#include "libmesh/enum_io_package.h"
#include "libmesh/enum_order.h"

namespace libMesh
{


// ------------------------------------------------------------
// InfHex16 class static member initializations
const int InfHex16::num_nodes;
const int InfHex16::num_sides;
const int InfHex16::num_edges;
const int InfHex16::num_children;
const int InfHex16::nodes_per_side;
const int InfHex16::nodes_per_edge;

const unsigned int InfHex16::side_nodes_map[InfHex16::num_sides][InfHex16::nodes_per_side] =
  {
    { 0, 1, 2, 3, 8, 9, 10, 11},   // Side 0
    { 0, 1, 4, 5, 8, 12, 99, 99},  // Side 1
    { 1, 2, 5, 6, 9, 13, 99, 99},  // Side 2
    { 2, 3, 6, 7, 10, 14, 99, 99}, // Side 3
    { 3, 0, 7, 4, 11, 15, 99, 99}  // Side 4
  };

const unsigned int InfHex16::edge_nodes_map[InfHex16::num_edges][InfHex16::nodes_per_edge] =
  {
    {0, 1,  8}, // Edge 0
    {1, 2,  9}, // Edge 1
    {2, 3, 10}, // Edge 2
    {0, 3, 11}, // Edge 3
    {0, 4, 99}, // Edge 4
    {1, 5, 99}, // Edge 5
    {2, 6, 99}, // Edge 6
    {3, 7, 99}  // Edge 7
  };

// ------------------------------------------------------------
// InfHex16 class member functions

bool InfHex16::is_vertex(const unsigned int i) const
{
  if (i < 4)
    return true;
  return false;
}

bool InfHex16::is_edge(const unsigned int i) const
{
  if (i < 4)
    return false;
  if (i > 11)
    return false;
  return true;
}

bool InfHex16::is_face(const unsigned int i) const
{
  if (i > 11)
    return true;
  return false;
}

bool InfHex16::is_node_on_side(const unsigned int n,
                               const unsigned int s) const
{
  libmesh_assert_less (s, n_sides());
  return std::find(std::begin(side_nodes_map[s]),
                   std::end(side_nodes_map[s]),
                   n) != std::end(side_nodes_map[s]);
}

std::vector<unsigned>
InfHex16::nodes_on_side(const unsigned int s) const
{
  libmesh_assert_less(s, n_sides());
  auto trim = (s == 0) ? 0 : 2;
  return {std::begin(side_nodes_map[s]), std::end(side_nodes_map[s]) - trim};
}

std::vector<unsigned>
InfHex16::nodes_on_edge(const unsigned int e) const
{
  libmesh_assert_less(e, n_edges());
  auto trim = (e < 4) ? 0 : 1;
  return {std::begin(edge_nodes_map[e]), std::end(edge_nodes_map[e]) - trim};
}

bool InfHex16::is_node_on_edge(const unsigned int n,
                               const unsigned int e) const
{
  libmesh_assert_less (e, n_edges());
  return std::find(std::begin(edge_nodes_map[e]),
                   std::end(edge_nodes_map[e]),
                   n) != std::end(edge_nodes_map[e]);
}



Order InfHex16::default_order() const
{
  return SECOND;
}



unsigned int InfHex16::local_side_node(unsigned int side,
                                       unsigned int side_node) const
{
  libmesh_assert_less (side, this->n_sides());

  // Never more than 8 nodes per side.
  libmesh_assert_less (side_node, InfHex16::nodes_per_side);

  // Some sides have 6 nodes.
  libmesh_assert(side == 0 || side_node < 6);

  return InfHex16::side_nodes_map[side][side_node];
}



unsigned int InfHex16::local_edge_node(unsigned int edge,
                                       unsigned int edge_node) const
{
  libmesh_assert_less (edge, this->n_edges());

  // Never more than 3 nodes per edge.
  libmesh_assert_less (edge_node, InfHex16::nodes_per_edge);

  // Some edges only have 2 nodes.
  libmesh_assert(edge < 4 || edge_node < 2);

  return InfHex16::edge_nodes_map[edge][edge_node];
}



std::unique_ptr<Elem> InfHex16::build_side_ptr (const unsigned int i,
                                                bool proxy)
{
  libmesh_assert_less (i, this->n_sides());

  std::unique_ptr<Elem> face;
  if (proxy)
    {
#ifdef LIBMESH_ENABLE_DEPRECATED
      libmesh_deprecated();
      switch (i)
        {
          // base
        case 0:
          {
            face = libmesh_make_unique<Side<Quad8,InfHex16>>(this,i);
            break;
          }

          // ifem sides
        case 1:
        case 2:
        case 3:
        case 4:
          {
            face = libmesh_make_unique<Side<InfQuad6,InfHex16>>(this,i);
            break;
          }

        default:
          libmesh_error_msg("Invalid side i = " << i);
        }
#else
      libmesh_error();
#endif // LIBMESH_ENABLE_DEPRECATED
    }
  else
    {
      // Think of a unit cube: (-1,1) x (-1,1) x (1,1)
      switch (i)
        {
          // the base face
        case 0:
          {
            face = libmesh_make_unique<Quad8>(this);
            break;
          }

          // connecting to another infinite element
        case 1:
        case 2:
        case 3:
        case 4:
          {
            face = libmesh_make_unique<InfQuad6>(this);
            break;
          }

        default:
          libmesh_error_msg("Invalid side i = " << i);
        }

      // Set the nodes
      for (auto n : face->node_index_range())
        face->set_node(n) = this->node_ptr(InfHex16::side_nodes_map[i][n]);
    }

#ifdef LIBMESH_ENABLE_DEPRECATED
  if (!proxy) // proxy sides used to leave parent() set
#endif
    face->set_parent(nullptr);
  face->set_interior_parent(this);

  face->subdomain_id() = this->subdomain_id();
  face->set_mapping_type(this->mapping_type());
#ifdef LIBMESH_ENABLE_AMR
  face->set_p_level(this->p_level());
#endif

  return face;
}



void InfHex16::build_side_ptr (std::unique_ptr<Elem> & side,
                               const unsigned int i)
{
  libmesh_assert_less (i, this->n_sides());

  // Think of a unit cube: (-1,1) x (-1,1) x (1,1)
  switch (i)
    {
      // the base face
    case 0:
      {
        if (!side.get() || side->type() != QUAD8)
          {
            side = this->build_side_ptr(i, false);
            return;
          }
        break;
      }

      // connecting to another infinite element
    case 1:
    case 2:
    case 3:
    case 4:
      {
        if (!side.get() || side->type() != INFQUAD6)
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
  side->set_mapping_type(this->mapping_type());

  // Set the nodes
  for (auto n : side->node_index_range())
    side->set_node(n) = this->node_ptr(InfHex16::side_nodes_map[i][n]);
}



std::unique_ptr<Elem> InfHex16::build_edge_ptr (const unsigned int i)
{
  if (i < 4) // base edges
    return this->simple_build_edge_ptr<Edge3,InfHex16>(i);

  // infinite edges
  return this->simple_build_edge_ptr<InfEdge2,InfHex16>(i);
}



void InfHex16::build_edge_ptr (std::unique_ptr<Elem> & edge,
                               const unsigned int i)
{
  libmesh_assert_less (i, this->n_edges());

  switch (i)
    {
      // the base edges
    case 0:
    case 1:
    case 2:
    case 3:
      {
        if (!edge.get() || edge->type() != EDGE3)
          {
            edge = this->build_edge_ptr(i);
            return;
          }
        break;
      }

      // the infinite edges
    case 4:
    case 5:
    case 6:
    case 7:
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

  edge->subdomain_id() = this->subdomain_id();
  edge->set_mapping_type(this->mapping_type());
#ifdef LIBMESH_ENABLE_AMR
  edge->set_p_level(this->p_level());
#endif

  // Set the nodes
  for (auto n : edge->node_index_range())
    edge->set_node(n) = this->node_ptr(InfHex16::edge_nodes_map[i][n]);
}



void InfHex16::connectivity(const unsigned int sc,
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
        switch (sc)
          {
          case 0:

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
            libmesh_error_msg("Invalid sc = " << sc);
          }
      }

    default:
      libmesh_error_msg("Unsupported IO package " << iop);
    }
}




unsigned short int InfHex16::second_order_adjacent_vertex (const unsigned int n,
                                                           const unsigned int v) const
{
  libmesh_assert_greater_equal (n, this->n_vertices());
  libmesh_assert_less (n, this->n_nodes());
  libmesh_assert_less (v, 2);
  // note that the _second_order_adjacent_vertices matrix is
  // stored in \p InfHex
  return _second_order_adjacent_vertices[n-this->n_vertices()][v];
}



std::pair<unsigned short int, unsigned short int>
InfHex16::second_order_child_vertex (const unsigned int n) const
{
  libmesh_assert_greater_equal (n, this->n_vertices());
  libmesh_assert_less (n, this->n_nodes());
  /*
   * the _second_order_vertex_child_* vectors are
   * stored in cell_inf_hex.C, since they are identical
   * for InfHex16 and InfHex18
   */
  return std::pair<unsigned short int, unsigned short int>
    (_second_order_vertex_child_number[n],
     _second_order_vertex_child_index[n]);
}





#ifdef LIBMESH_ENABLE_AMR

const Real InfHex16::_embedding_matrix[InfHex16::num_children][InfHex16::num_nodes][InfHex16::num_nodes] =
  {
    // embedding matrix for child 0
    {
      //          0           1           2           3           4           5           6           7           8           9          10          11          12          13          14          15 th parent Node
      {         1.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 0th child N.
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 1
      {       -0.25,      -0.25,      -0.25,      -0.25,        0.0,        0.0,        0.0,        0.0,        0.5,        0.5,        0.5,        0.5,        0.0,        0.0,        0.0,        0.0}, // 2
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0}, // 3
      {         0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 4
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0}, // 5
      {         0.0,        0.0,        0.0,        0.0,      -0.25,      -0.25,      -0.25,      -0.25,        0.0,        0.0,        0.0,        0.0,        0.5,        0.5,        0.5,        0.5}, // 6
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0}, // 7
      {       0.375,     -0.125,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 8
      {     -0.1875,    -0.1875,    -0.1875,    -0.1875,        0.0,        0.0,        0.0,        0.0,       0.75,      0.375,       0.25,      0.375,        0.0,        0.0,        0.0,        0.0}, // 9
      {     -0.1875,    -0.1875,    -0.1875,    -0.1875,        0.0,        0.0,        0.0,        0.0,      0.375,       0.25,      0.375,       0.75,        0.0,        0.0,        0.0,        0.0}, // 10
      {       0.375,        0.0,        0.0,     -0.125,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0,        0.0,        0.0}, // 11
      {         0.0,        0.0,        0.0,        0.0,      0.375,     -0.125,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0,        0.0}, // 12
      {         0.0,        0.0,        0.0,        0.0,    -0.1875,    -0.1875,    -0.1875,    -0.1875,        0.0,        0.0,        0.0,        0.0,       0.75,      0.375,       0.25,      0.375}, // 13
      {         0.0,        0.0,        0.0,        0.0,    -0.1875,    -0.1875,    -0.1875,    -0.1875,        0.0,        0.0,        0.0,        0.0,      0.375,       0.25,      0.375,       0.75}, // 14
      {         0.0,        0.0,        0.0,        0.0,      0.375,        0.0,        0.0,     -0.125,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75}  // 15
    },

    // embedding matrix for child 1
    {
      //          0           1           2           3           4           5           6           7           8           9          10          11          12          13          14          15 th parent Node
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 0th child N.
      {         0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 1
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 2
      {       -0.25,      -0.25,      -0.25,      -0.25,        0.0,        0.0,        0.0,        0.0,        0.5,        0.5,        0.5,        0.5,        0.0,        0.0,        0.0,        0.0}, // 3
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0}, // 4
      {         0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 5
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0}, // 6
      {         0.0,        0.0,        0.0,        0.0,      -0.25,      -0.25,      -0.25,      -0.25,        0.0,        0.0,        0.0,        0.0,        0.5,        0.5,        0.5,        0.5}, // 7
      {      -0.125,      0.375,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 8
      {         0.0,      0.375,     -0.125,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 9
      {     -0.1875,    -0.1875,    -0.1875,    -0.1875,        0.0,        0.0,        0.0,        0.0,      0.375,       0.75,      0.375,       0.25,        0.0,        0.0,        0.0,        0.0}, // 10
      {     -0.1875,    -0.1875,    -0.1875,    -0.1875,        0.0,        0.0,        0.0,        0.0,       0.75,      0.375,       0.25,      0.375,        0.0,        0.0,        0.0,        0.0}, // 11
      {         0.0,        0.0,        0.0,        0.0,     -0.125,      0.375,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0,        0.0}, // 12
      {         0.0,        0.0,        0.0,        0.0,        0.0,      0.375,     -0.125,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0}, // 13
      {         0.0,        0.0,        0.0,        0.0,    -0.1875,    -0.1875,    -0.1875,    -0.1875,        0.0,        0.0,        0.0,        0.0,      0.375,       0.75,      0.375,       0.25}, // 14
      {         0.0,        0.0,        0.0,        0.0,    -0.1875,    -0.1875,    -0.1875,    -0.1875,        0.0,        0.0,        0.0,        0.0,       0.75,      0.375,       0.25,      0.375}  // 15
    },

    // embedding matrix for child 2
    {
      //          0           1           2           3           4           5           6           7           8           9          10          11          12          13          14          15 th parent Node
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0}, // 0th child N.
      {       -0.25,      -0.25,      -0.25,      -0.25,        0.0,        0.0,        0.0,        0.0,        0.5,        0.5,        0.5,        0.5,        0.0,        0.0,        0.0,        0.0}, // 1
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 2
      {         0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 3
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0}, // 4
      {         0.0,        0.0,        0.0,        0.0,      -0.25,      -0.25,      -0.25,      -0.25,        0.0,        0.0,        0.0,        0.0,        0.5,        0.5,        0.5,        0.5}, // 5
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0}, // 6
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 7
      {     -0.1875,    -0.1875,    -0.1875,    -0.1875,        0.0,        0.0,        0.0,        0.0,      0.375,       0.25,      0.375,       0.75,        0.0,        0.0,        0.0,        0.0}, // 8
      {     -0.1875,    -0.1875,    -0.1875,    -0.1875,        0.0,        0.0,        0.0,        0.0,       0.25,      0.375,       0.75,      0.375,        0.0,        0.0,        0.0,        0.0}, // 9
      {         0.0,        0.0,     -0.125,      0.375,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0,        0.0,        0.0,        0.0}, // 10
      {      -0.125,        0.0,        0.0,      0.375,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0,        0.0,        0.0}, // 11
      {         0.0,        0.0,        0.0,        0.0,    -0.1875,    -0.1875,    -0.1875,    -0.1875,        0.0,        0.0,        0.0,        0.0,      0.375,       0.25,      0.375,       0.75}, // 12
      {         0.0,        0.0,        0.0,        0.0,    -0.1875,    -0.1875,    -0.1875,    -0.1875,        0.0,        0.0,        0.0,        0.0,       0.25,      0.375,       0.75,      0.375}, // 13
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,     -0.125,      0.375,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0}, // 14
      {         0.0,        0.0,        0.0,        0.0,     -0.125,        0.0,        0.0,      0.375,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75}  // 15
    },

    // embedding matrix for child 3
    {
      //          0           1           2           3           4           5           6           7           8           9          10          11          12          13          14          15 th parent Node
      {       -0.25,      -0.25,      -0.25,      -0.25,        0.0,        0.0,        0.0,        0.0,        0.5,        0.5,        0.5,        0.5,        0.0,        0.0,        0.0,        0.0}, // 0th child N.
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 1
      {         0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 2
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 3
      {         0.0,        0.0,        0.0,        0.0,      -0.25,      -0.25,      -0.25,      -0.25,        0.0,        0.0,        0.0,        0.0,        0.5,        0.5,        0.5,        0.5}, // 4
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0}, // 5
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 6
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0}, // 7
      {     -0.1875,    -0.1875,    -0.1875,    -0.1875,        0.0,        0.0,        0.0,        0.0,      0.375,       0.75,      0.375,       0.25,        0.0,        0.0,        0.0,        0.0}, // 8
      {         0.0,     -0.125,      0.375,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 9
      {         0.0,        0.0,      0.375,     -0.125,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0,        0.0,        0.0,        0.0}, // 10
      {     -0.1875,    -0.1875,    -0.1875,    -0.1875,        0.0,        0.0,        0.0,        0.0,       0.25,      0.375,       0.75,      0.375,        0.0,        0.0,        0.0,        0.0}, // 11
      {         0.0,        0.0,        0.0,        0.0,    -0.1875,    -0.1875,    -0.1875,    -0.1875,        0.0,        0.0,        0.0,        0.0,      0.375,       0.75,      0.375,       0.25}, // 12
      {         0.0,        0.0,        0.0,        0.0,        0.0,     -0.125,      0.375,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0}, // 13
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,      0.375,     -0.125,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0}, // 14
      {         0.0,        0.0,        0.0,        0.0,    -0.1875,    -0.1875,    -0.1875,    -0.1875,        0.0,        0.0,        0.0,        0.0,       0.25,      0.375,       0.75,      0.375}  // 15
    }
  };


#endif


void
InfHex16::permute(unsigned int perm_num)
{
  libmesh_assert_less (perm_num, 4);

  for (unsigned int i = 0; i != perm_num; ++i)
    {
      swap4nodes(0,1,2,3);
      swap4nodes(4,5,6,7);
      swap4nodes(8,9,10,11);
      swap4nodes(12,13,14,15);
      swap4neighbors(1,2,3,4);
    }
}


ElemType
InfHex16::side_type (const unsigned int s) const
{
  libmesh_assert_less (s, 5);
  if (s == 0)
    return QUAD8;
  return INFQUAD6;
}


} // namespace libMesh

#endif  // ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
