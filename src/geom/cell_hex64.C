// The libMesh Finite Element Library.
// Copyright (C) 2002-2026 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/cell_hex64.h"
#include "libmesh/edge_edge4.h"
#include "libmesh/face_quad16.h"
#include "libmesh/enum_io_package.h"
#include "libmesh/enum_order.h"

namespace libMesh
{



// ------------------------------------------------------------
// Hex64 class static member initializations
const int Hex64::num_nodes;
const int Hex64::nodes_per_side;
const int Hex64::nodes_per_edge;

// ------------------------------------------------------------
// Hex64 class member functions

bool Hex64::is_vertex(const unsigned int i) const
{
  if (i < 8)
    return true;
  return false;
}

bool Hex64::is_edge(const unsigned int i) const
{
  if (i < 8)
    return false;
  if (i > 31)
    return false;
  return true;
}

bool Hex64::is_face(const unsigned int i) const
{
  if (i > 55)
    return false;
  if (i > 31)
    return true;
  return false;
}

bool Hex64::is_node_on_side(const unsigned int n,
                            const unsigned int s) const
{
  libmesh_assert_less (s, n_sides());
  return std::find(std::begin(side_nodes_map[s]),
                   std::end(side_nodes_map[s]),
                   n) != std::end(side_nodes_map[s]);
}

std::vector<unsigned>
Hex64::nodes_on_side(const unsigned int s) const
{
  libmesh_assert_less(s, n_sides());
  return {std::begin(side_nodes_map[s]), std::end(side_nodes_map[s])};
}

std::vector<unsigned>
Hex64::nodes_on_edge(const unsigned int e) const
{
  libmesh_assert_less(e, n_edges());
  return {std::begin(edge_nodes_map[e]), std::end(edge_nodes_map[e])};
}

bool Hex64::is_node_on_edge(const unsigned int n,
                            const unsigned int e) const
{
  libmesh_assert_less (e, n_edges());
  return std::find(std::begin(edge_nodes_map[e]),
                   std::end(edge_nodes_map[e]),
                   n) != std::end(edge_nodes_map[e]);
}




Order Hex64::default_order() const
{
  return THIRD;
}


// TODO
dof_id_type Hex64::key (const unsigned int s) const
{
  libmesh_assert_less (s, this->n_sides());

  // Think of a unit cube: (-1,1) x (-1,1) x (1,1)
  switch (s)
    {
    case 0:  // the face at z=0

      return
        this->compute_key (this->node_id(20));

    case 1:  // the face at y = 0

      return
        this->compute_key (this->node_id(21));

    case 2:  // the face at x=1

      return
        this->compute_key (this->node_id(22));

    case 3: // the face at y=1

      return
        this->compute_key (this->node_id(23));

    case 4: // the face at x=0

      return
        this->compute_key (this->node_id(24));

    case 5: // the face at z=1

      return
        this->compute_key (this->node_id(25));

    default:
      libmesh_error_msg("Invalid side " << s);
    }
}



unsigned int Hex64::local_side_node(unsigned int side,
                                    unsigned int side_node) const
{
  libmesh_assert_less (side, this->n_sides());
  libmesh_assert_less (side_node, Hex64::nodes_per_side);

  return Hex64::side_nodes_map[side][side_node];
}



unsigned int Hex64::local_edge_node(unsigned int edge,
                                    unsigned int edge_node) const
{
  libmesh_assert_less (edge, this->n_edges());
  libmesh_assert_less (edge_node, Hex64::nodes_per_edge);

  return Hex64::edge_nodes_map[edge][edge_node];
}



std::unique_ptr<Elem> Hex64::build_side_ptr (const unsigned int i)
{
  return this->simple_build_side_ptr<Quad16, Hex64>(i);
}



void Hex64::build_side_ptr (std::unique_ptr<Elem> & side,
                            const unsigned int i)
{
  this->simple_build_side_ptr<Hex64>(side, i, QUAD16);
}



std::unique_ptr<Elem> Hex64::build_edge_ptr (const unsigned int i)
{
  return this->simple_build_edge_ptr<Edge4,Hex64>(i);
}



void Hex64::build_edge_ptr (std::unique_ptr<Elem> & edge, const unsigned int i)
{
  this->simple_build_edge_ptr<Hex64>(edge, i, EDGE4);
}



void Hex64::connectivity(const unsigned int sc,
                         const IOPackage iop,
                         std::vector<dof_id_type> &) const
{
  libmesh_assert(_nodes);
  libmesh_assert_less (sc, this->n_sub_elem());
  libmesh_assert_not_equal_to (iop, INVALID_IO_PACKAGE);
  libmesh_error_msg("Unsupported IO package for Hex64: " << iop);
}




//TODO
unsigned int Hex64::n_second_order_adjacent_vertices (const unsigned int n) const
{
  switch (n)
    {
    case 8:
    case 9:
    case 10:
    case 11:
    case 12:
    case 13:
    case 14:
    case 15:
    case 16:
    case 17:
    case 18:
    case 19:
      return 2;

    case 20:
    case 21:
    case 22:
    case 23:
    case 24:
    case 25:
      return 4;

    case 26:
      return 8;

    default:
      libmesh_error_msg("Invalid node number n = " << n);
    }
}


//TODO
unsigned short int Hex64::second_order_adjacent_vertex (const unsigned int n,
                                                        const unsigned int v) const
{
  libmesh_assert_greater_equal (n, this->n_vertices());
  libmesh_assert_less (n, this->n_nodes());

  switch (n)
    {
      /*
       * these are all nodes that are unique to Hex64,
       * use our _remaining.... matrix
       */
    case 20:
    case 21:
    case 22:
    case 23:
    case 24:
    case 25:
      {
        libmesh_assert_less (v, 4);
        return _remaining_second_order_adjacent_vertices[n-20][v];
      }

      /*
       * for the bubble node the return value is simply v.
       * Why? -- the user asks for the v-th adjacent vertex,
       * from \p n_second_order_adjacent_vertices() there
       * are 8 adjacent vertices, and these happen to be
       * 0..7
       */
    case 26:
      {
        libmesh_assert_less (v, 8);
        return static_cast<unsigned short int>(v);
      }

      /*
       * nodes 8..19:
       * these are all nodes that are identical for
       * Hex20 and Hex64.  Therefore use the
       * matrix stored in cell_hex.C
       */
    default:
      {
        libmesh_assert_less (v, 2);
        return _second_order_adjacent_vertices[n-this->n_vertices()][v];
      }
    }
}


//TODO
const unsigned short int Hex64::_remaining_second_order_adjacent_vertices[6][4] =
  {
    { 0,  1,  2,  3}, // vertices adjacent to node 20   face nodes
    { 0,  1,  4,  5}, // vertices adjacent to node 21
    { 1,  2,  5,  6}, // vertices adjacent to node 22
    { 2,  3,  6,  7}, // vertices adjacent to node 23
    { 0,  3,  4,  7}, // vertices adjacent to node 24
    { 4,  5,  6,  7}, // vertices adjacent to node 25
  };



std::pair<unsigned short int, unsigned short int>
Hex64::second_order_child_vertex (const unsigned int n) const
{
  libmesh_assert_greater_equal (n, this->n_vertices());
  libmesh_assert_less (n, this->n_nodes());
  /*
   * the _second_order_vertex_child_* vectors are
   * stored in cell_hex.C, since they are identical
   * for Hex20 and Hex64 (for the first 12 higher-order nodes)
   */
  return std::pair<unsigned short int, unsigned short int>
    (_second_order_vertex_child_number[n],
     _second_order_vertex_child_index[n]);
}

//TODO
void
Hex64::permute(unsigned int perm_num)
{
  libmesh_assert_less (perm_num, 24);
  const unsigned int side = perm_num % 6;
  const unsigned int rotate = perm_num / 6;

  for (unsigned int i = 0; i != rotate; ++i)
    {
      swap4nodes(0,1,2,3);
      swap4nodes(4,5,6,7);
      swap4nodes(8,9,10,11);
      swap4nodes(12,13,14,15);
      swap4nodes(16,17,18,19);
      swap4nodes(21,22,23,24);
      swap4neighbors(1,2,3,4);
    }

  switch (side) {
  case 0:
    break;
  case 1:
    swap4nodes(3,7,4,0);
    swap4nodes(11,15,19,12);
    swap4nodes(10,18,16,8);
    swap4nodes(2,6,5,1);
    swap4nodes(9,14,17,13);
    swap4nodes(20,23,25,21);
    swap4neighbors(0,3,5,1);
    break;
  case 2:
    swap4nodes(0,4,5,1);
    swap4nodes(8,12,16,13);
    swap4nodes(3,7,6,2);
    swap4nodes(10,15,18,14);
    swap4nodes(11,19,17,9);
    swap4nodes(20,24,25,22);
    swap4neighbors(0,4,5,2);
    break;
  case 3:
    swap4nodes(0,4,7,3);
    swap4nodes(12,19,15,11);
    swap4nodes(8,16,18,10);
    swap4nodes(1,5,6,2);
    swap4nodes(13,17,14,9);
    swap4nodes(20,21,25,23);
    swap4neighbors(0,1,5,3);
    break;
  case 4:
    swap4nodes(1,5,4,0);
    swap4nodes(8,13,16,12);
    swap4nodes(9,17,19,11);
    swap4nodes(2,6,7,3);
    swap4nodes(10,14,18,15);
    swap4nodes(20,22,25,24);
    swap4neighbors(0,2,5,4);
    break;
  case 5:
    swap2nodes(0,7);
    swap2nodes(8,18);
    swap2nodes(1,6);
    swap2nodes(2,5);
    swap2nodes(10,16);
    swap2nodes(3,4);
    swap2nodes(11,19);
    swap2nodes(12,15);
    swap2nodes(9,17);
    swap2nodes(13,14);
    swap2nodes(20,25);
    swap2nodes(21,23);
    swap2neighbors(0,5);
    swap2neighbors(1,3);
    break;
  default:
    libmesh_error();
  }
}

//TODO
void
Hex64::flip(BoundaryInfo * boundary_info)
{
  libmesh_assert(boundary_info);

  swap2nodes(0,1);
  swap2nodes(2,3);
  swap2nodes(4,5);
  swap2nodes(6,7);
  swap2nodes(9,11);
  swap2nodes(12,13);
  swap2nodes(14,15);
  swap2nodes(17,19);
  swap2nodes(22,24);
  swap2neighbors(2,4);
  swap2boundarysides(2,4,boundary_info);
  swap2boundaryedges(1,3,boundary_info);
  swap2boundaryedges(4,5,boundary_info);
  swap2boundaryedges(6,7,boundary_info);
  swap2boundaryedges(9,11,boundary_info);
}


ElemType
Hex64::side_type (const unsigned int libmesh_dbg_var(s)) const
{
  libmesh_assert_less (s, 6);
  return QUAD16;
}

const Real Hex64::_master_points[64][3] =
  {
    // vertices
    {   -1.,    -1.,    -1.}, // 0
    {    1.,    -1.,    -1.}, // 1
    {    1.,     1.,    -1.}, // 2
    {   -1.,     1.,    -1.}, // 3
    {   -1.,    -1.,     1.}, // 4
    {    1.,    -1.,     1.}, // 5
    {    1.,     1.,     1.}, // 6
    {   -1.,     1.,     1.}, // 7
    // edge nodes (2 per edge, from the edge's first vertex to its second)
    {-1./3.,    -1.,    -1.}, // 8
    { 1./3.,    -1.,    -1.}, // 9
    {    1., -1./3.,    -1.}, // 10
    {    1.,  1./3.,    -1.}, // 11
    { 1./3.,     1.,    -1.}, // 12
    {-1./3.,     1.,    -1.}, // 13
    {   -1.,  1./3.,    -1.}, // 14
    {   -1., -1./3.,    -1.}, // 15
    {   -1.,    -1., -1./3.}, // 16
    {   -1.,    -1.,  1./3.}, // 17
    {    1.,    -1., -1./3.}, // 18
    {    1.,    -1.,  1./3.}, // 19
    {    1.,     1., -1./3.}, // 20
    {    1.,     1.,  1./3.}, // 21
    {   -1.,     1., -1./3.}, // 22
    {   -1.,     1.,  1./3.}, // 23
    {-1./3.,    -1.,     1.}, // 24
    { 1./3.,    -1.,     1.}, // 25
    {    1., -1./3.,     1.}, // 26
    {    1.,  1./3.,     1.}, // 27
    { 1./3.,     1.,     1.}, // 28
    {-1./3.,     1.,     1.}, // 29
    {   -1.,  1./3.,     1.}, // 30
    {   -1., -1./3.,     1.}, // 31
    // face nodes (4 per side)
    {-1./3., -1./3.,    -1.}, // 32
    {-1./3.,  1./3.,    -1.}, // 33
    { 1./3.,  1./3.,    -1.}, // 34
    { 1./3., -1./3.,    -1.}, // 35
    {-1./3.,    -1., -1./3.}, // 36
    { 1./3.,    -1., -1./3.}, // 37
    { 1./3.,    -1.,  1./3.}, // 38
    {-1./3.,    -1.,  1./3.}, // 39
    {    1., -1./3., -1./3.}, // 40
    {    1.,  1./3., -1./3.}, // 41
    {    1.,  1./3.,  1./3.}, // 42
    {    1., -1./3.,  1./3.}, // 43
    { 1./3.,     1., -1./3.}, // 44
    {-1./3.,     1., -1./3.}, // 45
    {-1./3.,     1.,  1./3.}, // 46
    { 1./3.,     1.,  1./3.}, // 47
    {   -1.,  1./3., -1./3.}, // 48
    {   -1., -1./3., -1./3.}, // 49
    {   -1., -1./3.,  1./3.}, // 50
    {   -1.,  1./3.,  1./3.}, // 51
    {-1./3., -1./3.,     1.}, // 52
    { 1./3., -1./3.,     1.}, // 53
    { 1./3.,  1./3.,     1.}, // 54
    {-1./3.,  1./3.,     1.}, // 55
    // interior nodes (56+v is nearest vertex v)
    {-1./3., -1./3., -1./3.}, // 56
    { 1./3., -1./3., -1./3.}, // 57
    { 1./3.,  1./3., -1./3.}, // 58
    {-1./3.,  1./3., -1./3.}, // 59
    {-1./3., -1./3.,  1./3.}, // 60
    { 1./3., -1./3.,  1./3.}, // 61
    { 1./3.,  1./3.,  1./3.}, // 62
    {-1./3.,  1./3.,  1./3.}  // 63
};


} // namespace libMesh
