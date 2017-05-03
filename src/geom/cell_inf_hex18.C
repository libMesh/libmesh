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

// Local includes
#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

// C++ includes

// Local includes cont'd
#include "libmesh/cell_inf_hex18.h"
#include "libmesh/edge_edge3.h"
#include "libmesh/edge_inf_edge2.h"
#include "libmesh/face_quad9.h"
#include "libmesh/face_inf_quad6.h"
#include "libmesh/side.h"

namespace libMesh
{


// ------------------------------------------------------------
// InfHex18 class static member initializations
const unsigned int InfHex18::side_nodes_map[5][9] =
  {
    { 0, 1, 2, 3, 8, 9, 10, 11, 16},   // Side 0
    { 0, 1, 4, 5, 8, 12, 99, 99, 99},  // Side 1
    { 1, 2, 5, 6, 9, 13, 99, 99, 99},  // Side 2
    { 2, 3, 6, 7, 10, 14, 99, 99, 99}, // Side 3
    { 3, 0, 7, 4, 11, 15, 99, 99, 99}  // Side 4
  };

const unsigned int InfHex18::edge_nodes_map[8][3] =
  {
    { 0, 1, 8},  // Side 0
    { 1, 2, 9},  // Side 1
    { 2, 3, 10}, // Side 2
    { 0, 3, 11}, // Side 3
    { 0, 4, 99}, // Side 4
    { 1, 5, 99}, // Side 5
    { 2, 6, 99}, // Side 6
    { 3, 7, 99}  // Side 7
  };

// ------------------------------------------------------------
// InfHex18 class member functions

bool InfHex18::is_vertex(const unsigned int i) const
{
  if (i < 4)
    return true;
  return false;
}

bool InfHex18::is_edge(const unsigned int i) const
{
  if (i < 4)
    return false;
  if (i > 11)
    return false;
  return true;
}

bool InfHex18::is_face(const unsigned int i) const
{
  if (i > 11)
    return true;
  return false;
}

bool InfHex18::is_node_on_side(const unsigned int n,
                               const unsigned int s) const
{
  libmesh_assert_less (s, n_sides());
  for (unsigned int i = 0; i != 9; ++i)
    if (side_nodes_map[s][i] == n)
      return true;
  return false;
}

bool InfHex18::is_node_on_edge(const unsigned int n,
                               const unsigned int e) const
{
  libmesh_assert_less (e, n_edges());
  for (unsigned int i = 0; i != 3; ++i)
    if (edge_nodes_map[e][i] == n)
      return true;
  return false;
}



dof_id_type InfHex18::key (const unsigned int s) const
{
  libmesh_assert_less (s, this->n_sides());

  // Think of a unit cube: (-1,1) x (-1,1) x (1,1)
  switch (s)
    {
    case 0: // the base face
      return this->compute_key (this->node_id(16));

    case 1: // the face at y = -1
    case 2: // the face at x = 1
    case 3: // the face at y = 1
    case 4: // the face at x = -1
      return InfHex::key(s);

    default:
      libmesh_error_msg("Invalid side s = " << s);
    }

  libmesh_error_msg("We'll never get here!");
  return 0;
}



unsigned int InfHex18::which_node_am_i(unsigned int side,
                                       unsigned int side_node) const
{
  libmesh_assert_less (side, this->n_sides());

  // Never more than 9 nodes per side.
  libmesh_assert_less (side_node, 9);

  // Some sides have 6 nodes.
  libmesh_assert(side == 0 || side_node < 6);

  return InfHex18::side_nodes_map[side][side_node];
}



UniquePtr<Elem> InfHex18::build_side_ptr (const unsigned int i,
                                          bool proxy)
{
  libmesh_assert_less (i, this->n_sides());

  if (proxy)
    {
      switch (i)
        {
          // base
        case 0:
          return UniquePtr<Elem>(new Side<Quad9,InfHex18>(this,i));

          // ifem sides
        case 1:
        case 2:
        case 3:
        case 4:
          return UniquePtr<Elem>(new Side<InfQuad6,InfHex18>(this,i));

        default:
          libmesh_error_msg("Invalid side i = " << i);
        }
    }

  else
    {
      // Create NULL pointer to be initialized, returned later.
      Elem * face = libmesh_nullptr;

      // Think of a unit cube: (-1,1) x (-1,1) x (1,1)
      switch (i)
        {
          // the base face
        case 0:
          {
            face = new Quad9;
            break;
          }

          // connecting to another infinite element
        case 1:
        case 2:
        case 3:
        case 4:
          {
            face = new InfQuad6;
            break;
          }

        default:
          libmesh_error_msg("Invalid side i = " << i);
        }

      face->subdomain_id() = this->subdomain_id();

      // Set the nodes
      for (unsigned n=0; n<face->n_nodes(); ++n)
        face->set_node(n) = this->node_ptr(InfHex18::side_nodes_map[i][n]);

      return UniquePtr<Elem>(face);
    }

  libmesh_error_msg("We'll never get here!");
  return UniquePtr<Elem>();
}



UniquePtr<Elem> InfHex18::build_edge_ptr (const unsigned int i)
{
  libmesh_assert_less (i, this->n_edges());

  if (i < 4) // base edges
    return UniquePtr<Elem>(new SideEdge<Edge3,InfHex18>(this,i));
  // infinite edges
  return UniquePtr<Elem>(new SideEdge<InfEdge2,InfHex18>(this,i));
}



void InfHex18::connectivity(const unsigned int sc,
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
            conn[1] = this->node_id(8)+1;
            conn[2] = this->node_id(16)+1;
            conn[3] = this->node_id(11)+1;
            conn[4] = this->node_id(4)+1;
            conn[5] = this->node_id(12)+1;
            conn[6] = this->node_id(17)+1;
            conn[7] = this->node_id(15)+1;

            return;

          case 1:

            conn[0] = this->node_id(8)+1;
            conn[1] = this->node_id(1)+1;
            conn[2] = this->node_id(9)+1;
            conn[3] = this->node_id(16)+1;
            conn[4] = this->node_id(12)+1;
            conn[5] = this->node_id(5)+1;
            conn[6] = this->node_id(13)+1;
            conn[7] = this->node_id(17)+1;

            return;

          case 2:

            conn[0] = this->node_id(11)+1;
            conn[1] = this->node_id(16)+1;
            conn[2] = this->node_id(10)+1;
            conn[3] = this->node_id(3)+1;
            conn[4] = this->node_id(15)+1;
            conn[5] = this->node_id(17)+1;
            conn[6] = this->node_id(14)+1;
            conn[7] = this->node_id(7)+1;

            return;

          case 3:

            conn[0] = this->node_id(16)+1;
            conn[1] = this->node_id(9)+1;
            conn[2] = this->node_id(2)+1;
            conn[3] = this->node_id(10)+1;
            conn[4] = this->node_id(17)+1;
            conn[5] = this->node_id(13)+1;
            conn[6] = this->node_id(6)+1;
            conn[7] = this->node_id(14)+1;

            return;

          default:
            libmesh_error_msg("Invalid sc = " << sc);
          }
      }

    default:
      libmesh_error_msg("Unsupported IO package " << iop);
    }
}




unsigned int InfHex18::n_second_order_adjacent_vertices (const unsigned int n) const
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
      return 2;

    case 16:
    case 17:
      return 4;

    default:
      libmesh_error_msg("Invalid node n = " << n);
    }

  libmesh_error_msg("We'll never get here!");
  return libMesh::invalid_uint;
}



unsigned short int InfHex18::second_order_adjacent_vertex (const unsigned int n,
                                                           const unsigned int v) const
{
  libmesh_assert_greater_equal (n, this->n_vertices());
  libmesh_assert_less (n, this->n_nodes());
  libmesh_assert_less (v, this->n_second_order_adjacent_vertices(n));

  if (n == 16)
    /*
     * for the bubble node in the base the return value is
     * simply v.  Why? -- the user asks for the v-th
     * adjacent vertex, from \p n_second_order_adjacent_vertices()
     * there are 4 adjacent vertices, and these happen to be
     * 0..3
     */
    return static_cast<unsigned short int>(v);
  else if (n == 17)
    /*
     * for the bubble node further out similar reasoning works,
     * but v must be shifted to the further-out nodes:
     * simply add 4
     */
    return static_cast<unsigned short int>(v+4);

  else
    /*
     * all others are stored in the vertices matrix -- note
     * that this matrix is kept in \p InfHex to foster
     * code-reuse
     */
    return _second_order_adjacent_vertices[n-this->n_vertices()][v];
}



std::pair<unsigned short int, unsigned short int>
InfHex18::second_order_child_vertex (const unsigned int n) const
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

const float InfHex18::_embedding_matrix[4][18][18] =
  {
    // embedding matrix for child 0
    {
      //          0           1           2           3           4           5           6           7           8           9          10          11          12          13          14          15          16          17 th parent Node
      {         1.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 0th child N.
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 1
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0}, // 2
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 3
      {         0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 4
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 5
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0}, // 6
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0}, // 7
      {       0.375,     -0.125,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 8
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,      0.375,        0.0,     -0.125,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0}, // 9
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,     -0.125,        0.0,      0.375,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0}, // 10
      {       0.375,        0.0,        0.0,     -0.125,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 11
      {         0.0,        0.0,        0.0,        0.0,      0.375,     -0.125,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0,        0.0,        0.0,        0.0}, // 12
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,      0.375,        0.0,     -0.125,        0.0,        0.0,       0.75}, // 13
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,     -0.125,        0.0,      0.375,        0.0,       0.75}, // 14
      {         0.0,        0.0,        0.0,        0.0,      0.375,        0.0,        0.0,     -0.125,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0}, // 15
      {    0.140625,  -0.046875,   0.015625,  -0.046875,        0.0,        0.0,        0.0,        0.0,    0.28125,   -0.09375,   -0.09375,    0.28125,        0.0,        0.0,        0.0,        0.0,     0.5625,        0.0}, // 16
      {         0.0,        0.0,        0.0,        0.0,   0.140625,  -0.046875,   0.015625,  -0.046875,        0.0,        0.0,        0.0,        0.0,    0.28125,   -0.09375,   -0.09375,    0.28125,        0.0,     0.5625}  // 17
    },

    // embedding matrix for child 1
    {
      //          0           1           2           3           4           5           6           7           8           9          10          11          12          13          14          15          16          17 th parent Node
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 0th child N.
      {         0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 1
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 2
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0}, // 3
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 4
      {         0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 5
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0}, // 6
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0}, // 7
      {      -0.125,      0.375,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 8
      {         0.0,      0.375,     -0.125,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 9
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,      0.375,        0.0,     -0.125,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0}, // 10
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,      0.375,        0.0,     -0.125,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0}, // 11
      {         0.0,        0.0,        0.0,        0.0,     -0.125,      0.375,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0,        0.0,        0.0,        0.0}, // 12
      {         0.0,        0.0,        0.0,        0.0,        0.0,      0.375,     -0.125,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0,        0.0,        0.0}, // 13
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,      0.375,        0.0,     -0.125,        0.0,       0.75}, // 14
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,      0.375,        0.0,     -0.125,        0.0,        0.0,       0.75}, // 15
      {   -0.046875,   0.140625,  -0.046875,   0.015625,        0.0,        0.0,        0.0,        0.0,    0.28125,    0.28125,   -0.09375,   -0.09375,        0.0,        0.0,        0.0,        0.0,     0.5625,        0.0}, // 16
      {         0.0,        0.0,        0.0,        0.0,  -0.046875,   0.140625,  -0.046875,   0.015625,        0.0,        0.0,        0.0,        0.0,    0.28125,    0.28125,   -0.09375,   -0.09375,        0.0,     0.5625}  // 17
    },

    // embedding matrix for child 2
    {
      //          0           1           2           3           4           5           6           7           8           9          10          11          12          13          14          15          16          17 th parent Node
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 0th child N.
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0}, // 1
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 2
      {         0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 3
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0}, // 4
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0}, // 5
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0}, // 6
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 7
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,     -0.125,        0.0,      0.375,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0}, // 8
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,     -0.125,        0.0,      0.375,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0}, // 9
      {         0.0,        0.0,     -0.125,      0.375,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 10
      {      -0.125,        0.0,        0.0,      0.375,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 11
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,     -0.125,        0.0,      0.375,        0.0,       0.75}, // 12
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,     -0.125,        0.0,      0.375,        0.0,        0.0,       0.75}, // 13
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,     -0.125,      0.375,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0,        0.0}, // 14
      {         0.0,        0.0,        0.0,        0.0,     -0.125,        0.0,        0.0,      0.375,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0}, // 15
      {   -0.046875,   0.015625,  -0.046875,   0.140625,        0.0,        0.0,        0.0,        0.0,   -0.09375,   -0.09375,    0.28125,    0.28125,        0.0,        0.0,        0.0,        0.0,     0.5625,        0.0}, // 16
      {         0.0,        0.0,        0.0,        0.0,  -0.046875,   0.015625,  -0.046875,   0.140625,        0.0,        0.0,        0.0,        0.0,   -0.09375,   -0.09375,    0.28125,    0.28125,        0.0,     0.5625}  // 17
    },

    // embedding matrix for child 3
    {
      //          0           1           2           3           4           5           6           7           8           9          10          11          12          13          14          15          16          17 th parent Node
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0}, // 0th child N.
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 1
      {         0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 2
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 3
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0}, // 4
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0}, // 5
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 6
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0}, // 7
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,      0.375,        0.0,     -0.125,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0}, // 8
      {         0.0,     -0.125,      0.375,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 9
      {         0.0,        0.0,      0.375,     -0.125,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 10
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,     -0.125,        0.0,      0.375,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0}, // 11
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,      0.375,        0.0,     -0.125,        0.0,       0.75}, // 12
      {         0.0,        0.0,        0.0,        0.0,        0.0,     -0.125,      0.375,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0,        0.0,        0.0}, // 13
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,      0.375,     -0.125,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0,        0.0}, // 14
      {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,     -0.125,        0.0,      0.375,        0.0,        0.0,       0.75}, // 15
      {    0.015625,  -0.046875,   0.140625,  -0.046875,        0.0,        0.0,        0.0,        0.0,   -0.09375,    0.28125,    0.28125,   -0.09375,        0.0,        0.0,        0.0,        0.0,     0.5625,        0.0}, // 16
      {         0.0,        0.0,        0.0,        0.0,   0.015625,  -0.046875,   0.140625,  -0.046875,        0.0,        0.0,        0.0,        0.0,   -0.09375,    0.28125,    0.28125,   -0.09375,        0.0,     0.5625}  // 17
    }
  };




#endif

} // namespace libMesh

#endif // ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
