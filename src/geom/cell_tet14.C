// The libMesh Finite Element Library.
// Copyright (C) 2002-2022 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/cell_tet14.h"
#include "libmesh/edge_edge3.h"
#include "libmesh/face_tri7.h"
#include "libmesh/enum_io_package.h"
#include "libmesh/enum_order.h"

namespace {
  constexpr libMesh::Real r18 = 18;
  constexpr libMesh::Real r64 = 64;
  constexpr libMesh::Real r72 = 72;
}

namespace libMesh
{



// ------------------------------------------------------------
// Tet14 class static member initializations
const int Tet14::num_nodes;
const int Tet14::num_sides;
const int Tet14::num_edges;
const int Tet14::num_children;
const int Tet14::nodes_per_side;
const int Tet14::nodes_per_edge;

const unsigned int Tet14::side_nodes_map[Tet14::num_sides][Tet14::nodes_per_side] =
  {
    {0, 2, 1, 6, 5, 4, 10}, // Side 0
    {0, 1, 3, 4, 8, 7, 11}, // Side 1
    {1, 2, 3, 5, 9, 8, 12}, // Side 2
    {2, 0, 3, 6, 7, 9, 13}  // Side 3
  };

const unsigned int Tet14::edge_nodes_map[Tet14::num_edges][Tet14::nodes_per_edge] =
  {
    {0, 1, 4}, // Edge 0
    {1, 2, 5}, // Edge 1
    {0, 2, 6}, // Edge 2
    {0, 3, 7}, // Edge 3
    {1, 3, 8}, // Edge 4
    {2, 3, 9}  // Edge 5
  };

// ------------------------------------------------------------
// Tet14 class member functions

bool Tet14::is_vertex(const unsigned int i) const
{
  if (i < 4)
    return true;
  return false;
}

bool Tet14::is_edge(const unsigned int i) const
{
  if (i < 4 || i > 9)
    return false;
  return true;
}

bool Tet14::is_face(const unsigned int i) const
{
  if (i > 9)
    return true;
  return false;
}

bool Tet14::is_node_on_side(const unsigned int n,
                            const unsigned int s) const
{
  libmesh_assert_less (s, n_sides());
  return std::find(std::begin(side_nodes_map[s]),
                   std::end(side_nodes_map[s]),
                   n) != std::end(side_nodes_map[s]);
}

std::vector<unsigned>
Tet14::nodes_on_side(const unsigned int s) const
{
  libmesh_assert_less(s, n_sides());
  return {std::begin(side_nodes_map[s]), std::end(side_nodes_map[s])};
}

std::vector<unsigned>
Tet14::nodes_on_edge(const unsigned int e) const
{
  libmesh_assert_less(e, n_edges());
  return {std::begin(edge_nodes_map[e]), std::end(edge_nodes_map[e])};
}

bool Tet14::is_node_on_edge(const unsigned int n,
                            const unsigned int e) const
{
  libmesh_assert_less (e, n_edges());
  return std::find(std::begin(edge_nodes_map[e]),
                   std::end(edge_nodes_map[e]),
                   n) != std::end(edge_nodes_map[e]);
}


#ifdef LIBMESH_ENABLE_AMR

// This function only works if LIBMESH_ENABLE_AMR...
bool Tet14::is_child_on_side(const unsigned int c,
                             const unsigned int s) const
{
  // Table of local IDs for the midege nodes on the side opposite a given node.
  // See the ASCII art in the header file for this class to confirm this.
  const unsigned int midedge_nodes_opposite[4][3] =
    {
      {5,8,9}, // midedge nodes opposite node 0
      {6,7,9}, // midedge nodes opposite node 1
      {4,7,8}, // midedge nodes opposite node 2
      {4,5,6}  // midedge nodes opposite node 3
    };

  // Call the base class helper function
  return Tet::is_child_on_side_helper(c, s, midedge_nodes_opposite);
}

#else

bool Tet14::is_child_on_side(const unsigned int /*c*/,
                             const unsigned int /*s*/) const
{
  libmesh_not_implemented();
  return false;
}

#endif //LIBMESH_ENABLE_AMR



bool Tet14::has_affine_map() const
{
  // Make sure edges are straight
  Point v10 = this->point(1) - this->point(0);
  if (!v10.relative_fuzzy_equals
      ((this->point(4) - this->point(0))*2))
    return false;
  Point v21 = this->point(2) - this->point(1);
  if (!v21.relative_fuzzy_equals
      ((this->point(5) - this->point(1))*2))
    return false;
  Point v20 = this->point(2) - this->point(0);
  if (!v20.relative_fuzzy_equals
      ((this->point(6) - this->point(0))*2))
    return false;
  Point v30 = this->point(3) - this->point(0);
  if (!v30.relative_fuzzy_equals
      ((this->point(7) - this->point(0))*2))
    return false;
  Point v31 = this->point(3) - this->point(1);
  if (!v31.relative_fuzzy_equals
      ((this->point(8) - this->point(1))*2))
    return false;
  Point v32 = this->point(3) - this->point(2);
  if (!v32.relative_fuzzy_equals
      ((this->point(9) - this->point(2))*2))
    return false;

  // Make sure midface nodes are midface
  if (!(v20+v10).relative_fuzzy_equals
      ((this->point(10) - this->point(0))*3))
    return false;
  if (!(v30+v10).relative_fuzzy_equals
      ((this->point(11) - this->point(0))*3))
    return false;
  if (!(v31+v21).relative_fuzzy_equals
      ((this->point(12) - this->point(1))*3))
    return false;
  if (!(v30+v20).relative_fuzzy_equals
      ((this->point(13) - this->point(0))*3))
    return false;

  return true;
}



Order Tet14::default_order() const
{
  return THIRD;
}



unsigned int Tet14::local_side_node(unsigned int side,
                                    unsigned int side_node) const
{
  libmesh_assert_less (side, this->n_sides());
  libmesh_assert_less (side_node, Tet14::nodes_per_side);

  return Tet14::side_nodes_map[side][side_node];
}



unsigned int Tet14::local_edge_node(unsigned int edge,
                                    unsigned int edge_node) const
{
  libmesh_assert_less (edge, this->n_edges());
  libmesh_assert_less (edge_node, Tet14::nodes_per_edge);

  return Tet14::edge_nodes_map[edge][edge_node];
}



std::unique_ptr<Elem> Tet14::build_side_ptr (const unsigned int i,
                                             bool proxy)
{
  return this->simple_build_side_ptr<Tri7, Tet14>(i, proxy);
}



void Tet14::build_side_ptr (std::unique_ptr<Elem> & side,
                            const unsigned int i)
{
  this->simple_build_side_ptr<Tet14>(side, i, TRI7);
}



std::unique_ptr<Elem> Tet14::build_edge_ptr (const unsigned int i)
{
  return this->simple_build_edge_ptr<Edge3,Tet14>(i);
}



void Tet14::build_edge_ptr (std::unique_ptr<Elem> & edge, const unsigned int i)
{
  this->simple_build_edge_ptr<Tet14>(edge, i, EDGE3);
}



void Tet14::connectivity(const unsigned int sc,
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


            // Linear sub-tet 0
          case 0:

            conn[0] = this->node_id(0)+1;
            conn[1] = this->node_id(4)+1;
            conn[2] = this->node_id(6)+1;
            conn[3] = this->node_id(6)+1;
            conn[4] = this->node_id(7)+1;
            conn[5] = this->node_id(7)+1;
            conn[6] = this->node_id(7)+1;
            conn[7] = this->node_id(7)+1;

            return;

            // Linear sub-tet 1
          case 1:

            conn[0] = this->node_id(4)+1;
            conn[1] = this->node_id(1)+1;
            conn[2] = this->node_id(5)+1;
            conn[3] = this->node_id(5)+1;
            conn[4] = this->node_id(8)+1;
            conn[5] = this->node_id(8)+1;
            conn[6] = this->node_id(8)+1;
            conn[7] = this->node_id(8)+1;

            return;

            // Linear sub-tet 2
          case 2:

            conn[0] = this->node_id(5)+1;
            conn[1] = this->node_id(2)+1;
            conn[2] = this->node_id(6)+1;
            conn[3] = this->node_id(6)+1;
            conn[4] = this->node_id(9)+1;
            conn[5] = this->node_id(9)+1;
            conn[6] = this->node_id(9)+1;
            conn[7] = this->node_id(9)+1;

            return;

            // Linear sub-tet 3
          case 3:

            conn[0] = this->node_id(7)+1;
            conn[1] = this->node_id(8)+1;
            conn[2] = this->node_id(9)+1;
            conn[3] = this->node_id(9)+1;
            conn[4] = this->node_id(3)+1;
            conn[5] = this->node_id(3)+1;
            conn[6] = this->node_id(3)+1;
            conn[7] = this->node_id(3)+1;

            return;

            // Linear sub-tet 4
          case 4:

            conn[0] = this->node_id(4)+1;
            conn[1] = this->node_id(8)+1;
            conn[2] = this->node_id(6)+1;
            conn[3] = this->node_id(6)+1;
            conn[4] = this->node_id(7)+1;
            conn[5] = this->node_id(7)+1;
            conn[6] = this->node_id(7)+1;
            conn[7] = this->node_id(7)+1;

            return;

            // Linear sub-tet 5
          case 5:

            conn[0] = this->node_id(4)+1;
            conn[1] = this->node_id(5)+1;
            conn[2] = this->node_id(6)+1;
            conn[3] = this->node_id(6)+1;
            conn[4] = this->node_id(8)+1;
            conn[5] = this->node_id(8)+1;
            conn[6] = this->node_id(8)+1;
            conn[7] = this->node_id(8)+1;

            return;

            // Linear sub-tet 6
          case 6:

            conn[0] = this->node_id(5)+1;
            conn[1] = this->node_id(9)+1;
            conn[2] = this->node_id(6)+1;
            conn[3] = this->node_id(6)+1;
            conn[4] = this->node_id(8)+1;
            conn[5] = this->node_id(8)+1;
            conn[6] = this->node_id(8)+1;
            conn[7] = this->node_id(8)+1;

            return;

            // Linear sub-tet 7
          case 7:

            conn[0] = this->node_id(7)+1;
            conn[1] = this->node_id(6)+1;
            conn[2] = this->node_id(9)+1;
            conn[3] = this->node_id(9)+1;
            conn[4] = this->node_id(8)+1;
            conn[5] = this->node_id(8)+1;
            conn[6] = this->node_id(8)+1;
            conn[7] = this->node_id(8)+1;

            return;


          default:
            libmesh_error_msg("Invalid sc = " << sc);
          }
      }

    case VTK:
      {
        conn.resize(10);
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
        return;
        /*
          conn.resize(4);
          switch (sc)
          {
          // Linear sub-tet 0
          case 0:

          conn[0] = this->node_id(0);
          conn[1] = this->node_id(4);
          conn[2] = this->node_id(6);
          conn[3] = this->node_id(7);

          return;

          // Linear sub-tet 1
          case 1:

          conn[0] = this->node_id(4);
          conn[1] = this->node_id(1);
          conn[2] = this->node_id(5);
          conn[3] = this->node_id(8);

          return;

          // Linear sub-tet 2
          case 2:

          conn[0] = this->node_id(5);
          conn[1] = this->node_id(2);
          conn[2] = this->node_id(6);
          conn[3] = this->node_id(9);

          return;

          // Linear sub-tet 3
          case 3:

          conn[0] = this->node_id(7);
          conn[1] = this->node_id(8);
          conn[2] = this->node_id(9);
          conn[3] = this->node_id(3);

          return;

          // Linear sub-tet 4
          case 4:

          conn[0] = this->node_id(4);
          conn[1] = this->node_id(8);
          conn[2] = this->node_id(6);
          conn[3] = this->node_id(7);

          return;

          // Linear sub-tet 5
          case 5:

          conn[0] = this->node_id(4);
          conn[1] = this->node_id(5);
          conn[2] = this->node_id(6);
          conn[3] = this->node_id(8);

          return;

          // Linear sub-tet 6
          case 6:

          conn[0] = this->node_id(5);
          conn[1] = this->node_id(9);
          conn[2] = this->node_id(6);
          conn[3] = this->node_id(8);

          return;

          // Linear sub-tet 7
          case 7:

          conn[0] = this->node_id(7);
          conn[1] = this->node_id(6);
          conn[2] = this->node_id(9);
          conn[3] = this->node_id(8);

          return;


          default:

          libmesh_error_msg("Invalid sc = " << sc);
          }
        */
      }

    default:
      libmesh_error_msg("Unsupported IO package " << iop);
    }
}



unsigned int Tet14::n_second_order_adjacent_vertices (const unsigned int n) const
{
  switch (n)
    {
    case 4:
    case 5:
    case 6:
    case 7:
    case 8:
    case 9:
      return 2;

    case 10:
    case 11:
    case 12:
    case 13:
      return 3;

    default:
      libmesh_error_msg("Invalid n = " << n);
    }
}




const unsigned short int Tet14::_second_order_vertex_child_number[14] =
  {
    99,99,99,99, // Vertices
    0,1,0,0,1,2, // Edges
    5,4,6,7      // Faces
  };



const unsigned short int Tet14::_second_order_vertex_child_index[14] =
  {
    99,99,99,99, // Vertices
    1,2,2,3,3,3, // Edges
    10,13,12,13  // Faces
  };



std::pair<unsigned short int, unsigned short int>
Tet14::second_order_child_vertex (const unsigned int n) const
{
  libmesh_assert_greater_equal (n, this->n_vertices());
  libmesh_assert_less (n, this->n_nodes());
  return std::pair<unsigned short int, unsigned short int>
    (_second_order_vertex_child_number[n],
     _second_order_vertex_child_index[n]);
}



unsigned short int Tet14::second_order_adjacent_vertex (const unsigned int n,
                                                        const unsigned int v) const
{
  libmesh_assert_greater_equal (n, this->n_vertices());
  libmesh_assert_less (n, this->n_nodes());
  libmesh_assert_less (v, 3);
  libmesh_assert (n > 9 || v < 2); // Only face nodes have multiple adjacencies
  return _second_order_adjacent_vertices[n-this->n_vertices()][v];
}



const unsigned short int Tet14::_second_order_adjacent_vertices[10][3] =
  {
    {0, 1, 99}, // vertices adjacent to node 4
    {1, 2, 99}, // vertices adjacent to node 5
    {0, 2, 99}, // vertices adjacent to node 6
    {0, 3, 99}, // vertices adjacent to node 7
    {1, 3, 99}, // vertices adjacent to node 8
    {2, 3, 99}, // vertices adjacent to node 9
    {0, 1,  2}, // vertices adjacent to node 10
    {0, 1,  3}, // vertices adjacent to node 11
    {1, 2,  3}, // vertices adjacent to node 12
    {0, 2,  3}, // vertices adjacent to node 13
  };





#ifdef LIBMESH_ENABLE_AMR

const Real Tet14::_embedding_matrix[Tet14::num_children][Tet14::num_nodes][Tet14::num_nodes] =
  {
    // embedding matrix for child 0
    {
      //    0      1      2      3      4      5      6      7      8      9     10     11     12     13
      {    1.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.}, // 0
      {    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.}, // 1
      {    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.,    0.,    0.,    0.,    0.}, // 2
      {    0.,    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.,    0.,    0.,    0.}, // 3
      { 0.375,-0.125,    0.,    0.,  0.75,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.}, // 4
      {.09375,-.03125,-.03125,  0., 0.125,-0.125, 0.125,    0.,    0.,    0.,.84375,    0.,    0.,    0.}, // 5
      { 0.375,    0.,-0.125,    0.,    0.,    0.,  0.75,    0.,    0.,    0.,    0.,    0.,    0.,    0.}, // 6
      { 0.375,    0.,    0.,-0.125,    0.,    0.,    0.,  0.75,    0.,    0.,    0.,    0.,    0.,    0.}, // 7
      {.09375,-.03125,   0.,-.03125,0.125,    0.,    0., 0.125,-0.125,    0.,    0.,.84375,    0.,    0.}, // 8
      {.09375,    0.,-.03125,-.03125,  0.,    0., 0.125, 0.125,    0.,-0.125,    0.,    0.,    0.,.84375}, // 9
      { 5/r18,-1/r18,-1/r18,    0., 4/r18,-2/r18, 4/r18,    0.,    0.,    0.,   0.5,    0.,    0.,    0.}, // 10
      { 5/r18,-1/r18,    0.,-1/r18, 4/r18,    0.,    0., 4/r18,-2/r18,    0.,    0.,   0.5,    0.,    0.}, // 11
      { 0.125,-1/r72,-1/r72,-1/r72,    0.,-2/r18,    0.,    0.,-2/r18,-2/r18, 0.375, 0.375, 0.125, 0.375}, // 12
      { 5/r18,    0.,-1/r18,-1/r18,    0.,    0., 4/r18, 4/r18,    0.,-2/r18,    0.,    0.,    0.,   0.5}  // 13
    },

    // embedding matrix for child 1?
    {
      //    0      1      2      3      4      5      6      7      8      9     10     11     12     13
      {    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.}, // 0
      {    0.,    1.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.}, // 1
      {    0.,    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.}, // 2
      {    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.,    0.,    0.}, // 3
      {-0.125, 0.375,    0.,    0.,  0.75,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.}, // 4
      {    0., 0.375,-0.125,    0.,    0.,  0.75,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.}, // 5
      {-.03125,.09375,-.03125,  0., 0.125, 0.125,-0.125,    0.,    0.,    0.,.84375,    0.,    0.,    0.}, // 6
      {-.03125,.09375,   0.,-.03125,0.125,    0.,    0.,-0.125, 0.125,    0.,    0.,.84375,    0.,    0.}, // 7
      {    0., 0.375,    0.,-0.125,    0.,    0.,    0.,    0.,  0.75,    0.,    0.,    0.,    0.,    0.}, // 8
      {    0.,.09375,-0.03125,-0.03125,0., 0.125,    0.,    0., 0.125,-0.125,    0.,    0.,.84375,    0.}, // 9
      {-1/r18, 5/r18,-1/r18,    0., 4/r18, 4/r18,-2/r18,    0.,    0.,    0.,   0.5,    0.,    0.,    0.}, // 10
      {-1/r18, 5/r18,    0.,-1/r18, 4/r18,    0.,    0.,-2/r18, 4/r18,    0.,    0.,   0.5,    0.,    0.}, // 11
      {    0., 5/r18,-1/r18,-1/r18,    0., 4/r18,    0.,    0., 4/r18,-2/r18,    0.,    0.,   0.5,    0.}, // 12
      {-1/r72, 0.125,-1/r72,-1/r72,    0.,    0.,-2/r18,-2/r18,    0.,-2/r18, 0.375, 0.375, 0.375, 0.125}  // 13
    },

    // embedding matrix for child 2
    {
      //    0      1      2      3      4      5      6      7      8      9     10     11     12     13
      {    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.,    0.,    0.,    0.,    0.}, // 0
      {    0.,    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.}, // 1
      {    0.,    0.,    1.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.}, // 2
      {    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.,    0.}, // 3
      {-.03125,-.03125,.09375,  0.,-0.125, 0.125, 0.125,    0.,    0.,    0.,.84375,    0.,    0.,    0.}, // 4
      {    0.,-0.125, 0.375,    0.,    0.,  0.75,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.}, // 5
      {-0.125,    0., 0.375,    0.,    0.,    0.,  0.75,    0.,    0.,    0.,    0.,    0.,    0.,    0.}, // 6
      {-.03125,   0.,.09375,-.03125,   0.,    0., 0.125,-0.125,    0., 0.125,    0.,    0.,    0.,.84375}, // 7
      {    0.,-.03125,.09375,-.03125,  0., 0.125,    0.,    0.,-0.125, 0.125,    0.,    0.,.84375,    0.}, // 8
      {    0.,    0., 0.375,-0.125,    0.,    0.,    0.,    0.,    0.,  0.75,    0.,    0.,    0.,    0.}, // 9
      {-1/r18,-1/r18, 5/r18,    0.,-2/r18, 4/r18, 4/r18,    0.,    0.,    0.,   0.5,    0.,    0.,    0.}, // 10
      {-1/r72,-1/r72, 0.125,-1/r72,-2/r18,    0.,    0.,-2/r18,-2/r18,    0., 0.375, 0.125, 0.375, 0.375}, // 11
      {    0.,-1/r18, 5/r18,-1/r18,    0., 4/r18,    0.,    0.,-2/r18, 4/r18,    0.,    0.,   0.5,    0.}, // 12
      {-1/r18,    0., 5/r18,-1/r18,    0.,    0., 4/r18,-2/r18,    0., 4/r18,    0.,    0.,    0.,   0.5}  // 13
    },

    // embedding matrix for child 3
    {
      //    0      1      2      3      4      5      6      7      8      9     10     11     12     13
      {    0.,    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.,    0.,    0.,    0.}, // 0
      {    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.,    0.,    0.}, // 1
      {    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.,    0.}, // 2
      {    0.,    0.,    0.,    1.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.}, // 3
      {-.03125,-.03125,  0.,.09375,-0.125,    0.,    0., 0.125, 0.125,    0.,    0.,.84375,    0.,    0.}, // 4
      {    0.,-.03125,-.03125,.09375,  0.,-0.125,    0.,    0., 0.125, 0.125,    0.,    0.,.84375,    0.}, // 5
      {-.03125,    0.,-.03125,.09375,  0.,    0.,-0.125, 0.125,    0., 0.125,    0.,    0.,    0.,.84375}, // 6
      {-0.125,    0.,    0., 0.375,    0.,    0.,    0.,  0.75,    0.,    0.,    0.,    0.,    0.,    0.}, // 7
      {    0.,-0.125,    0., 0.375,    0.,    0.,    0.,    0.,  0.75,    0.,    0.,    0.,    0.,    0.}, // 8
      {    0.,    0.,-0.125, 0.375,    0.,    0.,    0.,    0.,    0.,  0.75,    0.,    0.,    0.,    0.}, // 9
      {-1/r72,-1/r72,-1/r72, 0.125,-2/r18,-2/r18,-2/r18,    0.,    0.,    0., 0.125, 0.375, 0.375, 0.375}, // 10
      {-1/r18,-1/r18,    0., 5/r18,-2/r18,    0.,    0., 4/r18, 4/r18,    0.,    0.,   0.5,    0.,    0.}, // 11
      {    0.,-1/r18,-1/r18, 5/r18,    0.,-2/r18,    0.,    0., 4/r18, 4/r18,    0.,    0.,   0.5,    0.}, // 12
      {-1/r18,    0.,-1/r18, 5/r18,    0.,    0.,-2/r18, 4/r18,    0., 4/r18,    0.,    0.,    0.,   0.5}  // 13
    },

    // embedding matrix for child 4
    {
      //    0      1      2      3      4      5      6      7      8      9     10     11     12     13
      {    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.}, // 0
      {    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.,    0.,    0.}, // 1
      {    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.,    0.,    0.,    0.,    0.}, // 2
      {    0.,    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.,    0.,    0.,    0.}, // 3
      {-.03125,.09375,   0.,-.03125,0.125,    0.,    0.,-0.125, 0.125,    0.,    0.,.84375,    0.,    0.}, // 4
      {-.03125,.09375,   0.,-.03125,0.125,    0.,    0.,-0.125, 0.125,    0.,    0.,.84375,    0.,    0.}, // 5
      {.09375,-.03125,-.03125,  0., 0.125,-0.125, 0.125,    0.,    0.,    0.,.84375,    0.,    0.,    0.}, // 6
      {.09375,-.03125,   0.,-.03125,0.125,    0.,    0., 0.125,-0.125,    0.,    0.,.84375,    0.,    0.}, // 7
      {-.03125,-.03125,  0.,.09375,-0.125,    0.,    0., 0.125, 0.125,    0.,    0.,.84375,    0.,    0.}, // 8
      {.09375,    0.,-.03125,-.03125,  0.,    0., 0.125, 0.125,    0.,-0.125,    0.,    0.,    0.,.84375}, // 9
      { 2/r72, 2/r72,    0.,    0.,    0.,-2/r18,-2/r18,-2/r18,-2/r18,-2/r18,   0.5,   0.5,  0.25,  0.25}, // 10
      {    0.,    0,     0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.,    0.}, // 11
      { 2/r72,    0.,    0., 2/r72,-2/r18,-2/r18,-2/r18,    0.,-2/r18,-2/r18,  0.25,   0.5,  0.25,   0.5}, // 12
      { 5/r18,    0.,-1/r18,-1/r18,    0.,    0., 4/r18, 4/r18,    0.,-2/r18,    0.,    0.,    0.,   0.5}  // 13
    },

    // embedding matrix for child 5?
    {
      //    0      1      2      3      4      5      6      7      8      9     10     11     12     13
      {    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.}, // 0
      {    0.,    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.}, // 1
      {    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.,    0.,    0.,    0.,    0.}, // 2
      {    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.,    0.,    0.}, // 3
      {-.03125,.09375,-.03125,  0., 0.125, 0.125,-0.125,    0.,    0.,    0.,.84375,    0.,    0.,    0.}, // 4
      {-.03125,-.03125,.09375,  0.,-0.125, 0.125, 0.125,    0.,    0.,    0.,.84375,    0.,    0.,    0.}, // 5
      {.09375,-.03125,-.03125,  0., 0.125,-0.125, 0.125,    0.,    0.,    0.,.84375,    0.,    0.,    0.}, // 6
      {-.03125,.09375,   0.,-.03125,0.125,    0.,    0.,-0.125, 0.125,    0.,    0.,.84375,    0.,    0.}, // 7
      {    0.,.09375,-.03125,-.03125,  0., 0.125,    0.,    0., 0.125,-0.125,    0.,    0.,.84375,    0.}, // 8
      { 1/r64, 1/r64, 1/r64, 1/r64,-0.125,-0.125,-0.125,-0.125,-0.125,-0.125,27/r64,27/r64,27/r64,27/r64}, // 9
      {    0.,    0,     0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.}, // 10
      {-1/r72, 0.125,-1/r72,-1/r72,    0.,    0.,-2/r18,-2/r18,    0.,-2/r18, 0.375, 0.375, 0.375, 0.125}, // 11
      {    0., 2/r72, 2/r72,    0.,-2/r18,    0.,-2/r18,-2/r18,-2/r18,-2/r18,   0.5,  0.25,   0.5,  0.25}, // 12
      { 2/r72, 2/r72,    0.,    0.,    0.,-2/r18,-2/r18,-2/r18,-2/r18,-2/r18,   0.5,   0.5,  0.25,  0.25}  // 13
    },

    // embedding matrix for child 6?
    {
      //    0      1      2      3      4      5      6      7      8      9     10     11     12     13
      {    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.,    0.,    0.,    0.,    0.}, // 0
      {    0.,    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.}, // 1
      {    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.,    0.}, // 2
      {    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.,    0.,    0.}, // 3
      {-.03125,-.03125,.09375,  0.,-0.125, 0.125, 0.125,    0.,    0.,    0.,.84375,    0.,    0.,    0.}, // 4
      {    0.,-.03125,.09375,-.03125,  0., 0.125,    0.,    0.,-0.125, 0.125,    0.,    0.,.84375,    0.}, // 5
      {-.03125,    0.,.09375,-.03125,  0.,    0., 0.125,-0.125,    0., 0.125,    0.,    0.,    0.,.84375}, // 6
      { 1/r64, 1/r64, 1/r64, 1/r64,-0.125,-0.125,-0.125,-0.125,-0.125,-0.125,27/r64,27/r64,27/r64,27/r64}, // 7
      {    0.,.09375,-.03125,-.03125,  0., 0.125,    0.,    0., 0.125,-0.125,    0.,    0.,.84375,    0.}, // 8
      {    0.,-.03125,-.03125,.09375,  0.,-0.125,    0.,    0., 0.125, 0.125,    0.,    0.,.84375,    0.}, // 9
      {-1/r72,-1/r72, 0.125,-1/r72,-2/r18,    0.,    0.,-2/r18,-2/r18,    0., 0.375, 0.125, 0.375, 0.375}, // 10
      {    0., 2/r72, 2/r72,    0.,-2/r18,    0.,-2/r18,-2/r18,-2/r18,-2/r18,   0.5,  0.25,   0.5,  0.25}, // 11
      {    0.,    0,     0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.}, // 12
      {    0.,    0., 2/r72, 2/r72,-2/r18,-2/r18,-2/r18,-2/r18,-2/r18,    0.,  0.25,  0.25,   0.5,   0.5}  // 13
    },

    // embedding matrix for child 7?
    {
      //    0      1      2      3      4      5      6      7      8      9     10     11     12     13
      {    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.,    0.,    0.,    0.,    0.}, // 0
      {    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.,    0.,    0.}, // 1
      {    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.,    0.}, // 2
      {    0.,    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.,    0.,    0.,    0.}, // 3
      { 1/r64, 1/r64, 1/r64, 1/r64,-0.125,-0.125,-0.125,-0.125,-0.125,-0.125,27/r64,27/r64,27/r64,27/r64}, // 4
      {-.03125,-.03125,  0.,.09375,-0.125,    0.,    0., 0.125, 0.125,    0.,    0.,.84375,    0.,    0.}, // 5
      {.09375,    0.,-.03125,-.03125,  0.,    0., 0.125, 0.125,    0.,-0.125,    0.,    0.,    0.,.84375}, // 6
      {-.03125,    0.,.09375,-.03125,  0.,    0., 0.125,-0.125,    0., 0.125,    0.,    0.,    0.,.84375}, // 7
      {    0.,-.03125,-.03125,.09375,  0.,-0.125,    0.,    0., 0.125, 0.125,    0.,    0.,.84375,    0.}, // 8
      {-.03125,    0.,-.03125,.09375,  0.,    0.,-0.125, 0.125,    0., 0.125,    0.,    0.,    0.,.84375}, // 9
      {    0.,    0., 2/r72, 2/r72,-2/r18,-2/r18,-2/r18,-2/r18,-2/r18,    0.,  0.25,  0.25,   0.5,   0.5}, // 10
      { 2/r72,    0.,    0., 2/r72,-2/r18,-2/r18,-2/r18,    0.,-2/r18,-2/r18,  0.25,   0.5,  0.25,   0.5}, // 11
      {-1/r72,-1/r72,-1/r72, 0.125,-2/r18,-2/r18,-2/r18,    0.,    0.,    0., 0.125, 0.375, 0.375, 0.375}, // 12
      {    0.,    0,     0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    1.}, // 13
    }
  };



Real Tet14::embedding_matrix (const unsigned int i,
                              const unsigned int j,
                              const unsigned int k) const
{
  // Choose an optimal diagonal, if one has not already been selected
  this->choose_diagonal();

  // Permuted j and k indices
  unsigned int
    jp=j,
    kp=k;

  if ((i>3) && (this->_diagonal_selection!=DIAG_02_13))
    {
      // Just the enum value cast to an unsigned int...
      const unsigned ds = static_cast<unsigned int>(this->_diagonal_selection); // == 1 or 2

      // Instead of doing a lot of arithmetic, use these
      // straightforward arrays for the permutations.  Note that 3 ->
      // 3 and 10 -> 10, and the first array consists of "forward"
      // permutations of the sets {0,1,2}, {4,5,6}, {7,8,9}, and
      // {11,12,13} while the second array consists of "reverse"
      // permutations of the same sets.

      const unsigned int perms[2][14] =
        {
          {1, 2, 0, 3, 5, 6, 4, 8, 9, 7, 10, 12, 13, 11},
          {2, 0, 1, 3, 6, 4, 5, 9, 7, 8, 10, 13, 11, 12}
        };

      // Permute j
      jp = perms[ds-1][j];
      //      if (jp<3)
      //        jp = (jp+ds)%3;
      //      else if (jp>3)
      //        jp = (jp-1+ds)%3 + 1 + 3*((jp-1)/3);

      // Permute k
      kp = perms[ds-1][k];
      //      if (kp<3)
      //        kp = (kp+ds)%3;
      //      else if (kp>3)
      //        kp = (kp-1+ds)%3 + 1 + 3*((kp-1)/3);
    }

  // Debugging:
  // libMesh::err << "Selected diagonal " << _diagonal_selection << std::endl;
  // libMesh::err << "j=" << j << std::endl;
  // libMesh::err << "k=" << k << std::endl;
  // libMesh::err << "jp=" << jp << std::endl;
  // libMesh::err << "kp=" << kp << std::endl;

  // Call embedding matrix with permuted indices
  return this->_embedding_matrix[i][jp][kp];
}

const std::vector<std::pair<unsigned char, unsigned char>>
  Tet14::_parent_bracketing_nodes[Tet14::num_children][Tet14::num_nodes] =
  {
    // Node 0,      1,      2,      3,      4,      5,      6,      7,      8,      9,
    //             10,             11,             12,             13
    // Child 0
    {    {{}},{{0,1}},{{0,2}},{{0,3}},{{0,4}},{{4,6}},{{0,6}},{{0,7}},{{4,7}},{{6,7}},{{0,10}},{{0,11}},{{0,12}},{{0,13}} },
    // Child 1
    { {{0,1}},   {{}},{{1,2}},{{1,3}},{{1,4}},{{1,5}},{{4,5}},{{4,8}},{{1,8}},{{5,8}},{{1,10}},{{1,11}},{{1,12}},{{1,13}} },
    // Child 2
    { {{0,2}},{{1,2}},   {{}},{{2,3}},{{5,6}},{{2,5}},{{2,6}},{{6,9}},{{5,9}},{{2,9}},{{2,10}},{{2,11}},{{2,12}},{{2,13}} },
    // Child 3
    { {{0,3}},{{1,3}},{{2,3}},   {{}},{{7,8}},{{8,9}},{{7,9}},{{3,7}},{{3,8}},{{3,9}},{{3,10}},{{3,11}},{{3,12}},{{3,13}} },
    // Child 4
    { {{0,1}},{{1,3}},{{0,2}},{{0,3}},{{4,8}},{{4,8}},{{4,6}},{{4,7}},{{7,8}},{{6,7}},
            {{10,11}},{{0,8},{1,7},{3,4}},  {{11,13}},       {{0,13}} },
    // Child 5
    { {{0,1}},{{1,2}},{{0,2}},{{1,3}},{{4,5}},{{5,6}},{{4,6}},{{4,8}},{{5,8}},{{6,8}},
  {{0,5},{1,6},{2,4}},       {{1,13}},      {{10,12}},      {{10,11}} },
    // Child 6
    { {{0,2}},{{1,2}},{{2,3}},{{1,3}},{{5,6}},{{5,9}},{{6,9}},{{6,8}},{{5,8}},{{8,9}},
             {{2,11}},      {{10,12}},{{1,9},{2,8},{3,5}},  {{12,13}} },
    // Child 7
    { {{0,2}},{{1,3}},{{2,3}},{{0,3}},{{6,8}},{{7,8}},{{6,7}},{{6,9}},{{8,9}},{{7,9}},
            {{12,13}},      {{11,13}},       {{3,10}},{{0,9},{2,7},{3,6}} }
  };

#endif // #ifdef LIBMESH_ENABLE_AMR



void Tet14::permute(unsigned int perm_num)
{
  libmesh_assert_less (perm_num, 12);

  const unsigned int side = perm_num % 4;
  const unsigned int rotate = perm_num / 4;

  for (unsigned int i = 0; i != rotate; ++i)
    {
      swap3nodes(0,1,2);
      swap3nodes(4,5,6);
      swap3nodes(7,8,9);
      swap3nodes(11,12,13);
      swap3neighbors(1,2,3);
    }

  switch (side) {
  case 0:
    break;
  case 1:
    swap3nodes(0,2,3);
    swap3nodes(4,5,8);
    swap3nodes(6,9,7);
    swap3nodes(10,12,11);
    swap3neighbors(0,2,1);
    break;
  case 2:
    swap3nodes(2,0,3);
    swap3nodes(5,4,8);
    swap3nodes(6,7,9);
    swap3nodes(10,11,12);
    swap3neighbors(0,1,2);
    break;
  case 3:
    swap3nodes(2,1,3);
    swap3nodes(5,8,9);
    swap3nodes(6,4,7);
    swap3nodes(10,11,13);
    swap3neighbors(0,1,3);
    break;
  default:
    libmesh_error();
  }
}


ElemType Tet14::side_type (const unsigned int libmesh_dbg_var(s)) const
{
  libmesh_assert_less (s, 4);
  return TRI7;
}


} // namespace libMesh
