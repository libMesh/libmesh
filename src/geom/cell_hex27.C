// The libMesh Finite Element Library.
// Copyright (C) 2002-2014 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/cell_hex27.h"
#include "libmesh/edge_edge3.h"
#include "libmesh/face_quad9.h"

namespace libMesh
{



// ------------------------------------------------------------
// Hex27 class static member initializations
const unsigned int Hex27::side_nodes_map[6][9] =
  {
    {0, 3, 2, 1, 11, 10,  9,  8, 20}, // Side 0
    {0, 1, 5, 4,  8, 13, 16, 12, 21}, // Side 1
    {1, 2, 6, 5,  9, 14, 17, 13, 22}, // Side 2
    {2, 3, 7, 6, 10, 15, 18, 14, 23}, // Side 3
    {3, 0, 4, 7, 11, 12, 19, 15, 24}, // Side 4
    {4, 5, 6, 7, 16, 17, 18, 19, 25}  // Side 5
  };

const unsigned int Hex27::edge_nodes_map[12][3] =
  {
    {0, 1, 8},  // Side 0
    {1, 2, 9},  // Side 1
    {2, 3, 10}, // Side 2
    {0, 3, 11}, // Side 3
    {0, 4, 12}, // Side 4
    {1, 5, 13}, // Side 5
    {2, 6, 14}, // Side 6
    {3, 7, 15}, // Side 7
    {4, 5, 16}, // Side 8
    {5, 6, 17}, // Side 9
    {6, 7, 18}, // Side 10
    {4, 7, 19}  // Side 11
  };



// ------------------------------------------------------------
// Hex27 class member functions

bool Hex27::is_vertex(const unsigned int i) const
{
  if (i < 8)
    return true;
  return false;
}

bool Hex27::is_edge(const unsigned int i) const
{
  if (i < 8)
    return false;
  if (i > 19)
    return false;
  return true;
}

bool Hex27::is_face(const unsigned int i) const
{
  if (i == 26)
    return false;
  if (i > 19)
    return true;
  return false;
}

bool Hex27::is_node_on_side(const unsigned int n,
                            const unsigned int s) const
{
  libmesh_assert_less (s, n_sides());
  for (unsigned int i = 0; i != 9; ++i)
    if (side_nodes_map[s][i] == n)
      return true;
  return false;
}

bool Hex27::is_node_on_edge(const unsigned int n,
                            const unsigned int e) const
{
  libmesh_assert_less (e, n_edges());
  for (unsigned int i = 0; i != 3; ++i)
    if (edge_nodes_map[e][i] == n)
      return true;
  return false;
}



bool Hex27::has_affine_map() const
{
  // Make sure x-edge endpoints are affine
  Point v = this->point(1) - this->point(0);
  if (!v.relative_fuzzy_equals(this->point(2) - this->point(3)) ||
      !v.relative_fuzzy_equals(this->point(5) - this->point(4)) ||
      !v.relative_fuzzy_equals(this->point(6) - this->point(7)))
    return false;
  // Make sure x-edges are straight
  // and x-face and center points are centered
  v /= 2;
  if (!v.relative_fuzzy_equals(this->point(8) - this->point(0)) ||
      !v.relative_fuzzy_equals(this->point(10) - this->point(3)) ||
      !v.relative_fuzzy_equals(this->point(16) - this->point(4)) ||
      !v.relative_fuzzy_equals(this->point(18) - this->point(7)) ||
      !v.relative_fuzzy_equals(this->point(20) - this->point(11)) ||
      !v.relative_fuzzy_equals(this->point(21) - this->point(12)) ||
      !v.relative_fuzzy_equals(this->point(23) - this->point(15)) ||
      !v.relative_fuzzy_equals(this->point(25) - this->point(19)) ||
      !v.relative_fuzzy_equals(this->point(26) - this->point(24)))
    return false;
  // Make sure xz-faces are identical parallelograms
  v = this->point(4) - this->point(0);
  if (!v.relative_fuzzy_equals(this->point(7) - this->point(3)))
    return false;
  v /= 2;
  if (!v.relative_fuzzy_equals(this->point(12) - this->point(0)) ||
      !v.relative_fuzzy_equals(this->point(13) - this->point(1)) ||
      !v.relative_fuzzy_equals(this->point(14) - this->point(2)) ||
      !v.relative_fuzzy_equals(this->point(15) - this->point(3)) ||
      !v.relative_fuzzy_equals(this->point(22) - this->point(9)) ||
      !v.relative_fuzzy_equals(this->point(24) - this->point(11)))
    return false;
  // Make sure y-edges are straight
  v = (this->point(3) - this->point(0))/2;
  if (!v.relative_fuzzy_equals(this->point(11) - this->point(0)) ||
      !v.relative_fuzzy_equals(this->point(9) - this->point(1)) ||
      !v.relative_fuzzy_equals(this->point(17) - this->point(5)) ||
      !v.relative_fuzzy_equals(this->point(19) - this->point(4)))
    return false;
  // If all the above checks out, the map is affine
  return true;
}



dof_id_type Hex27::key (const unsigned int s) const
{
  libmesh_assert_less (s, this->n_sides());

  // Think of a unit cube: (-1,1) x (-1,1) x (1,1)
  switch (s)
    {
    case 0:  // the face at z=0

      return
        this->compute_key (this->node(20));

    case 1:  // the face at y = 0

      return
        this->compute_key (this->node(21));

    case 2:  // the face at x=1

      return
        this->compute_key (this->node(22));

    case 3: // the face at y=1

      return
        this->compute_key (this->node(23));

    case 4: // the face at x=0

      return
        this->compute_key (this->node(24));

    case 5: // the face at z=1

      return
        this->compute_key (this->node(25));

    default:
      libmesh_error_msg("Invalid side " << s);
    }

  libmesh_error_msg("We'll never get here!");
  return 0;
}



AutoPtr<Elem> Hex27::build_side (const unsigned int i,
                                 bool proxy) const
{
  libmesh_assert_less (i, this->n_sides());

  if (proxy)
    return AutoPtr<Elem>(new Side<Quad9,Hex27>(this,i));

  else
    {
      Elem* face = new Quad9;
      face->subdomain_id() = this->subdomain_id();

      // Think of a unit cube: (-1,1) x (-1,1) x (1,1)
      switch (i)
        {
        case 0:  // the face at z=0
          {
            face->set_node(0) = this->get_node(0);
            face->set_node(1) = this->get_node(3);
            face->set_node(2) = this->get_node(2);
            face->set_node(3) = this->get_node(1);
            face->set_node(4) = this->get_node(11);
            face->set_node(5) = this->get_node(10);
            face->set_node(6) = this->get_node(9);
            face->set_node(7) = this->get_node(8);
            face->set_node(8) = this->get_node(20);
            break;
          }
        case 1:  // the face at y = 0
          {
            face->set_node(0) = this->get_node(0);
            face->set_node(1) = this->get_node(1);
            face->set_node(2) = this->get_node(5);
            face->set_node(3) = this->get_node(4);
            face->set_node(4) = this->get_node(8);
            face->set_node(5) = this->get_node(13);
            face->set_node(6) = this->get_node(16);
            face->set_node(7) = this->get_node(12);
            face->set_node(8) = this->get_node(21);
            break;
          }
        case 2:  // the face at x=1
          {
            face->set_node(0) = this->get_node(1);
            face->set_node(1) = this->get_node(2);
            face->set_node(2) = this->get_node(6);
            face->set_node(3) = this->get_node(5);
            face->set_node(4) = this->get_node(9);
            face->set_node(5) = this->get_node(14);
            face->set_node(6) = this->get_node(17);
            face->set_node(7) = this->get_node(13);
            face->set_node(8) = this->get_node(22);
            break;
          }
        case 3: // the face at y=1
          {
            face->set_node(0) = this->get_node(2);
            face->set_node(1) = this->get_node(3);
            face->set_node(2) = this->get_node(7);
            face->set_node(3) = this->get_node(6);
            face->set_node(4) = this->get_node(10);
            face->set_node(5) = this->get_node(15);
            face->set_node(6) = this->get_node(18);
            face->set_node(7) = this->get_node(14);
            face->set_node(8) = this->get_node(23);
            break;
          }
        case 4: // the face at x=0
          {
            face->set_node(0) = this->get_node(3);
            face->set_node(1) = this->get_node(0);
            face->set_node(2) = this->get_node(4);
            face->set_node(3) = this->get_node(7);
            face->set_node(4) = this->get_node(11);
            face->set_node(5) = this->get_node(12);
            face->set_node(6) = this->get_node(19);
            face->set_node(7) = this->get_node(15);
            face->set_node(8) = this->get_node(24);
            break;
          }
        case 5: // the face at z=1
          {
            face->set_node(0) = this->get_node(4);
            face->set_node(1) = this->get_node(5);
            face->set_node(2) = this->get_node(6);
            face->set_node(3) = this->get_node(7);
            face->set_node(4) = this->get_node(16);
            face->set_node(5) = this->get_node(17);
            face->set_node(6) = this->get_node(18);
            face->set_node(7) = this->get_node(19);
            face->set_node(8) = this->get_node(25);
            break;
          }
        default:
          libmesh_error_msg("Invalid side i = " << i);
        }

      return AutoPtr<Elem>(face);
    }

  libmesh_error_msg("We'll never get here!");
  return AutoPtr<Elem>();
}



AutoPtr<Elem> Hex27::build_edge (const unsigned int i) const
{
  libmesh_assert_less (i, this->n_edges());

  return AutoPtr<Elem>(new SideEdge<Edge3,Hex27>(this,i));
}



void Hex27::connectivity(const unsigned int sc,
                         const IOPackage iop,
                         std::vector<dof_id_type>& conn) const
{
  libmesh_assert(_nodes);
  libmesh_assert_less (sc, this->n_sub_elem());
  libmesh_assert_not_equal_to (iop, INVALID_IO_PACKAGE);

  conn.resize(8);

  switch (iop)
    {
    case TECPLOT:
      {
        switch (sc)
          {
          case 0:

            conn[0] = this->node(0)+1;
            conn[1] = this->node(8)+1;
            conn[2] = this->node(20)+1;
            conn[3] = this->node(11)+1;
            conn[4] = this->node(12)+1;
            conn[5] = this->node(21)+1;
            conn[6] = this->node(26)+1;
            conn[7] = this->node(24)+1;

            return;

          case 1:

            conn[0] = this->node(8)+1;
            conn[1] = this->node(1)+1;
            conn[2] = this->node(9)+1;
            conn[3] = this->node(20)+1;
            conn[4] = this->node(21)+1;
            conn[5] = this->node(13)+1;
            conn[6] = this->node(22)+1;
            conn[7] = this->node(26)+1;

            return;

          case 2:

            conn[0] = this->node(11)+1;
            conn[1] = this->node(20)+1;
            conn[2] = this->node(10)+1;
            conn[3] = this->node(3)+1;
            conn[4] = this->node(24)+1;
            conn[5] = this->node(26)+1;
            conn[6] = this->node(23)+1;
            conn[7] = this->node(15)+1;

            return;

          case 3:

            conn[0] = this->node(20)+1;
            conn[1] = this->node(9)+1;
            conn[2] = this->node(2)+1;
            conn[3] = this->node(10)+1;
            conn[4] = this->node(26)+1;
            conn[5] = this->node(22)+1;
            conn[6] = this->node(14)+1;
            conn[7] = this->node(23)+1;

            return;

          case 4:

            conn[0] = this->node(12)+1;
            conn[1] = this->node(21)+1;
            conn[2] = this->node(26)+1;
            conn[3] = this->node(24)+1;
            conn[4] = this->node(4)+1;
            conn[5] = this->node(16)+1;
            conn[6] = this->node(25)+1;
            conn[7] = this->node(19)+1;

            return;

          case 5:

            conn[0] = this->node(21)+1;
            conn[1] = this->node(13)+1;
            conn[2] = this->node(22)+1;
            conn[3] = this->node(26)+1;
            conn[4] = this->node(16)+1;
            conn[5] = this->node(5)+1;
            conn[6] = this->node(17)+1;
            conn[7] = this->node(25)+1;

            return;

          case 6:

            conn[0] = this->node(24)+1;
            conn[1] = this->node(26)+1;
            conn[2] = this->node(23)+1;
            conn[3] = this->node(15)+1;
            conn[4] = this->node(19)+1;
            conn[5] = this->node(25)+1;
            conn[6] = this->node(18)+1;
            conn[7] = this->node(7)+1;

            return;

          case 7:

            conn[0] = this->node(26)+1;
            conn[1] = this->node(22)+1;
            conn[2] = this->node(14)+1;
            conn[3] = this->node(23)+1;
            conn[4] = this->node(25)+1;
            conn[5] = this->node(17)+1;
            conn[6] = this->node(6)+1;
            conn[7] = this->node(18)+1;

            return;

          default:
            libmesh_error_msg("Invalid sc = " << sc);
          }
      }

    case VTK:
      {
        // VTK now supports VTK_TRIQUADRATIC_HEXAHEDRON directly
        conn.resize(27);

        conn[0] = this->node(0);
        conn[1] = this->node(1);
        conn[2] = this->node(2);
        conn[3] = this->node(3);
        conn[4] = this->node(4);
        conn[5] = this->node(5);
        conn[6] = this->node(6);
        conn[7] = this->node(7);
        conn[8] = this->node(8);
        conn[9] = this->node(9);
        conn[10] = this->node(10);
        conn[11] = this->node(11); //
        conn[12] = this->node(16);
        conn[13] = this->node(17);
        conn[14] = this->node(18);
        conn[15] = this->node(19);
        conn[16] = this->node(12);
        conn[17] = this->node(13); //
        conn[18] = this->node(14);
        conn[19] = this->node(15);
        conn[20] = this->node(24);
        conn[21] = this->node(22);
        conn[22] = this->node(21);
        conn[23] = this->node(23);
        conn[24] = this->node(20);
        conn[25] = this->node(25);
        conn[26] = this->node(26);

        return;

        /*
          switch (sc)
          {
          case 0:

          conn[0] = this->node(0);
          conn[1] = this->node(8);
          conn[2] = this->node(20);
          conn[3] = this->node(11);
          conn[4] = this->node(12);
          conn[5] = this->node(21);
          conn[6] = this->node(26);
          conn[7] = this->node(24);

          return;

          case 1:

          conn[0] = this->node(8);
          conn[1] = this->node(1);
          conn[2] = this->node(9);
          conn[3] = this->node(20);
          conn[4] = this->node(21);
          conn[5] = this->node(13);
          conn[6] = this->node(22);
          conn[7] = this->node(26);

          return;

          case 2:

          conn[0] = this->node(11);
          conn[1] = this->node(20);
          conn[2] = this->node(10);
          conn[3] = this->node(3);
          conn[4] = this->node(24);
          conn[5] = this->node(26);
          conn[6] = this->node(23);
          conn[7] = this->node(15);

          return;

          case 3:

          conn[0] = this->node(20);
          conn[1] = this->node(9);
          conn[2] = this->node(2);
          conn[3] = this->node(10);
          conn[4] = this->node(26);
          conn[5] = this->node(22);
          conn[6] = this->node(14);
          conn[7] = this->node(23);

          return;

          case 4:

          conn[0] = this->node(12);
          conn[1] = this->node(21);
          conn[2] = this->node(26);
          conn[3] = this->node(24);
          conn[4] = this->node(4);
          conn[5] = this->node(16);
          conn[6] = this->node(25);
          conn[7] = this->node(19);

          return;

          case 5:

          conn[0] = this->node(21);
          conn[1] = this->node(13);
          conn[2] = this->node(22);
          conn[3] = this->node(26);
          conn[4] = this->node(16);
          conn[5] = this->node(5);
          conn[6] = this->node(17);
          conn[7] = this->node(25);

          return;

          case 6:

          conn[0] = this->node(24);
          conn[1] = this->node(26);
          conn[2] = this->node(23);
          conn[3] = this->node(15);
          conn[4] = this->node(19);
          conn[5] = this->node(25);
          conn[6] = this->node(18);
          conn[7] = this->node(7);

          return;

          case 7:

          conn[0] = this->node(26);
          conn[1] = this->node(22);
          conn[2] = this->node(14);
          conn[3] = this->node(23);
          conn[4] = this->node(25);
          conn[5] = this->node(17);
          conn[6] = this->node(6);
          conn[7] = this->node(18);

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





unsigned int Hex27::n_second_order_adjacent_vertices (const unsigned int n) const
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

  libmesh_error_msg("We'll never get here!");
  return libMesh::invalid_uint;
}



unsigned short int Hex27::second_order_adjacent_vertex (const unsigned int n,
                                                        const unsigned int v) const
{
  libmesh_assert_greater_equal (n, this->n_vertices());
  libmesh_assert_less (n, this->n_nodes());

  switch (n)
    {
      /*
       * these are all nodes that are unique to Hex27,
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
       * Hex20 and Hex27.  Therefore use the
       * matrix stored in cell_hex.C
       */
    default:
      {
        libmesh_assert_less (v, 2);
        return _second_order_adjacent_vertices[n-this->n_vertices()][v];
      }
    }
}



const unsigned short int Hex27::_remaining_second_order_adjacent_vertices[6][4] =
  {
    { 0,  1,  2,  3}, // vertices adjacent to node 20   face nodes
    { 0,  1,  4,  5}, // vertices adjacent to node 21
    { 1,  2,  5,  6}, // vertices adjacent to node 22
    { 2,  3,  6,  7}, // vertices adjacent to node 23
    { 0,  3,  4,  7}, // vertices adjacent to node 24
    { 4,  5,  6,  7}, // vertices adjacent to node 25
  };



std::pair<unsigned short int, unsigned short int>
Hex27::second_order_child_vertex (const unsigned int n) const
{
  libmesh_assert_greater_equal (n, this->n_vertices());
  libmesh_assert_less (n, this->n_nodes());
  /*
   * the _second_order_vertex_child_* vectors are
   * stored in cell_hex.C, since they are identical
   * for Hex20 and Hex27 (for the first 12 higher-order nodes)
   */
  return std::pair<unsigned short int, unsigned short int>
    (_second_order_vertex_child_number[n],
     _second_order_vertex_child_index[n]);
}







#ifdef LIBMESH_ENABLE_AMR

const float Hex27::_embedding_matrix[8][27][27] =
  {
    // embedding matrix for child 0
    {
      //  0     1     2     3     4     5     6     7    8  9  10  11  12  13  14  15  16  17  18  19  20   21  22  23  24  25  26
      {    1.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 0
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 1
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 2
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 3
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 4
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 5
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000 }, // 6
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000 }, // 7
      {   0.375000,  -0.125000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 8
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.375000,    0.00000,  -0.125000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 9
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,  -0.125000,    0.00000,   0.375000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 10
      {   0.375000,    0.00000,    0.00000,  -0.125000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 11
      {   0.375000,    0.00000,    0.00000,    0.00000,  -0.125000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 12
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.375000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,  -0.125000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 13
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.375000,    0.00000,    0.00000,    0.00000,    0.00000,  -0.125000,   0.750000 }, // 14
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.375000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,  -0.125000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000 }, // 15
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.375000,  -0.125000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 16
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.375000,    0.00000,  -0.125000,    0.00000,    0.00000,   0.750000 }, // 17
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,  -0.125000,    0.00000,   0.375000,    0.00000,   0.750000 }, // 18
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.375000,    0.00000,    0.00000,  -0.125000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000 }, // 19
      {   0.140625, -0.0468750,  0.0156250, -0.0468750,    0.00000,    0.00000,    0.00000,    0.00000,   0.281250, -0.0937500, -0.0937500,   0.281250,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.562500,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 20
      {   0.140625, -0.0468750,    0.00000,    0.00000, -0.0468750,  0.0156250,    0.00000,    0.00000,   0.281250,    0.00000,    0.00000,    0.00000,   0.281250, -0.0937500,    0.00000,    0.00000, -0.0937500,    0.00000,    0.00000,    0.00000,    0.00000,   0.562500,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 21
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.140625,    0.00000, -0.0468750,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000, -0.0468750,    0.00000,  0.0156250,    0.00000,   0.281250,   0.281250,    0.00000, -0.0937500,    0.00000, -0.0937500,   0.562500 }, // 22
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000, -0.0468750,    0.00000,   0.140625,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,  0.0156250,    0.00000, -0.0468750,   0.281250,    0.00000, -0.0937500,    0.00000,   0.281250, -0.0937500,   0.562500 }, // 23
      {   0.140625,    0.00000,    0.00000, -0.0468750, -0.0468750,    0.00000,    0.00000,  0.0156250,    0.00000,    0.00000,    0.00000,   0.281250,   0.281250,    0.00000,    0.00000, -0.0937500,    0.00000,    0.00000,    0.00000, -0.0937500,    0.00000,    0.00000,    0.00000,    0.00000,   0.562500,    0.00000,    0.00000 }, // 24
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.140625, -0.0468750,  0.0156250, -0.0468750,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.281250, -0.0937500, -0.0937500,   0.281250,    0.00000,   0.562500 }, // 25
      { 0.052734375, -0.017578125, 0.005859375, -0.017578125, -0.017578125, 0.005859375, -0.001953125, 0.005859375, 0.10546875, -0.03515625, -0.03515625, 0.10546875, 0.10546875, -0.03515625, 0.01171875, -0.03515625, -0.03515625, 0.01171875, 0.01171875, -0.03515625, 0.2109375, 0.2109375, -0.0703125, -0.0703125, 0.2109375, -0.0703125, 0.421875 } // 26
    },

    // embedding matrix for child 1
    {
      //  0     1     2     3     4     5     6     7    8  9  11  12  13  14  15  16  17  18  19  20   21  22  23  24  25  26
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 0
      {    0.00000,    1.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 1
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 2
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 3
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 4
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 5
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 6
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000 }, // 7
      {  -0.125000,   0.375000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 8
      {    0.00000,   0.375000,  -0.125000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 9
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.375000,    0.00000,  -0.125000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 10
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.375000,    0.00000,  -0.125000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 11
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.375000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,  -0.125000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 12
      {    0.00000,   0.375000,    0.00000,    0.00000,    0.00000,  -0.125000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 13
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.375000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,  -0.125000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 14
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.375000,    0.00000,    0.00000,    0.00000,    0.00000,  -0.125000,   0.750000 }, // 15
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,  -0.125000,   0.375000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 16
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.375000,  -0.125000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 17
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.375000,    0.00000,  -0.125000,    0.00000,   0.750000 }, // 18
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.375000,    0.00000,  -0.125000,    0.00000,    0.00000,   0.750000 }, // 19
      { -0.0468750,   0.140625, -0.0468750,  0.0156250,    0.00000,    0.00000,    0.00000,    0.00000,   0.281250,   0.281250, -0.0937500, -0.0937500,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.562500,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 20
      { -0.0468750,   0.140625,    0.00000,    0.00000,  0.0156250, -0.0468750,    0.00000,    0.00000,   0.281250,    0.00000,    0.00000,    0.00000, -0.0937500,   0.281250,    0.00000,    0.00000, -0.0937500,    0.00000,    0.00000,    0.00000,    0.00000,   0.562500,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 21
      {    0.00000,   0.140625, -0.0468750,    0.00000,    0.00000, -0.0468750,  0.0156250,    0.00000,    0.00000,   0.281250,    0.00000,    0.00000,    0.00000,   0.281250, -0.0937500,    0.00000,    0.00000, -0.0937500,    0.00000,    0.00000,    0.00000,    0.00000,   0.562500,    0.00000,    0.00000,    0.00000,    0.00000 }, // 22
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.140625,    0.00000, -0.0468750,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000, -0.0468750,    0.00000,  0.0156250,   0.281250,    0.00000,   0.281250,    0.00000, -0.0937500, -0.0937500,   0.562500 }, // 23
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.140625,    0.00000, -0.0468750,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000, -0.0468750,    0.00000,  0.0156250,    0.00000,   0.281250,   0.281250,    0.00000, -0.0937500,    0.00000, -0.0937500,   0.562500 }, // 24
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000, -0.0468750,   0.140625, -0.0468750,  0.0156250,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.281250,   0.281250, -0.0937500, -0.0937500,    0.00000,   0.562500 }, // 25
      { -0.017578125, 0.052734375, -0.017578125, 0.005859375, 0.005859375, -0.017578125, 0.005859375, -0.001953125, 0.10546875, 0.10546875, -0.03515625, -0.03515625, -0.03515625, 0.10546875, -0.03515625, 0.01171875, -0.03515625, -0.03515625, 0.01171875, 0.01171875, 0.2109375, 0.2109375, 0.2109375, -0.0703125, -0.0703125, -0.0703125, 0.421875 } // 26

    },

    // embedding matrix for child 2
    {
      //  0     1     2     3     4     5     6     7    8  9  11  12  13  14  15  16  17  18  19  20   21  22  23  24  25  26
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 0
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 1
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 2
      {    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 3
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000 }, // 4
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000 }, // 5
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000 }, // 6
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 7
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,  -0.125000,    0.00000,   0.375000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 8
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,  -0.125000,    0.00000,   0.375000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 9
      {    0.00000,    0.00000,  -0.125000,   0.375000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 10
      {  -0.125000,    0.00000,    0.00000,   0.375000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 11
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.375000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,  -0.125000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000 }, // 12
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.375000,    0.00000,    0.00000,    0.00000,    0.00000,  -0.125000,   0.750000 }, // 13
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.375000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,  -0.125000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000 }, // 14
      {    0.00000,    0.00000,    0.00000,   0.375000,    0.00000,    0.00000,    0.00000,  -0.125000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 15
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,  -0.125000,    0.00000,   0.375000,    0.00000,   0.750000 }, // 16
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,  -0.125000,    0.00000,   0.375000,    0.00000,    0.00000,   0.750000 }, // 17
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,  -0.125000,   0.375000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000 }, // 18
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,  -0.125000,    0.00000,    0.00000,   0.375000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000 }, // 19
      { -0.0468750,  0.0156250, -0.0468750,   0.140625,    0.00000,    0.00000,    0.00000,    0.00000, -0.0937500, -0.0937500,   0.281250,   0.281250,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.562500,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 20
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000, -0.0468750,    0.00000,   0.140625,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,  0.0156250,    0.00000, -0.0468750,   0.281250,    0.00000, -0.0937500,    0.00000,   0.281250, -0.0937500,   0.562500 }, // 21
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000, -0.0468750,    0.00000,   0.140625,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,  0.0156250,    0.00000, -0.0468750,    0.00000,   0.281250, -0.0937500,    0.00000,   0.281250,    0.00000, -0.0937500,   0.562500 }, // 22
      {    0.00000,    0.00000, -0.0468750,   0.140625,    0.00000,    0.00000,  0.0156250, -0.0468750,    0.00000,    0.00000,   0.281250,    0.00000,    0.00000,    0.00000, -0.0937500,   0.281250,    0.00000,    0.00000, -0.0937500,    0.00000,    0.00000,    0.00000,    0.00000,   0.562500,    0.00000,    0.00000,    0.00000 }, // 23
      { -0.0468750,    0.00000,    0.00000,   0.140625,  0.0156250,    0.00000,    0.00000, -0.0468750,    0.00000,    0.00000,    0.00000,   0.281250, -0.0937500,    0.00000,    0.00000,   0.281250,    0.00000,    0.00000,    0.00000, -0.0937500,    0.00000,    0.00000,    0.00000,    0.00000,   0.562500,    0.00000,    0.00000 }, // 24
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000, -0.0468750,  0.0156250, -0.0468750,   0.140625,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000, -0.0937500, -0.0937500,   0.281250,   0.281250,    0.00000,   0.562500 }, // 25
      { -0.017578125, 0.005859375, -0.017578125, 0.052734375, 0.005859375, -0.001953125, 0.005859375, -0.017578125, -0.03515625, -0.03515625, 0.10546875, 0.10546875, -0.03515625, 0.01171875, -0.03515625, 0.10546875, 0.01171875, 0.01171875, -0.03515625, -0.03515625, 0.2109375, -0.0703125, -0.0703125, 0.2109375, 0.2109375, -0.0703125, 0.421875 } // 26
    },

    // embedding matrix for child 3
    {
      //  0     1     2     3     4     5     6     7    8  9  11  12  13  14  15  16  17  18  19  20   21  22  23  24  25  26
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 0
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 1
      {    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 2
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 3
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000 }, // 4
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 5
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 6
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000 }, // 7
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.375000,    0.00000,  -0.125000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 8
      {    0.00000,  -0.125000,   0.375000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 9
      {    0.00000,    0.00000,   0.375000,  -0.125000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 10
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,  -0.125000,    0.00000,   0.375000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 11
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.375000,    0.00000,    0.00000,    0.00000,    0.00000,  -0.125000,   0.750000 }, // 12
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.375000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,  -0.125000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 13
      {    0.00000,    0.00000,   0.375000,    0.00000,    0.00000,    0.00000,  -0.125000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 14
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.375000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,  -0.125000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000 }, // 15
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.375000,    0.00000,  -0.125000,    0.00000,   0.750000 }, // 16
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,  -0.125000,   0.375000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 17
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.375000,  -0.125000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000 }, // 18
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,  -0.125000,    0.00000,   0.375000,    0.00000,    0.00000,   0.750000 }, // 19
      {  0.0156250, -0.0468750,   0.140625, -0.0468750,    0.00000,    0.00000,    0.00000,    0.00000, -0.0937500,   0.281250,   0.281250, -0.0937500,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.562500,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 20
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.140625,    0.00000, -0.0468750,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000, -0.0468750,    0.00000,  0.0156250,   0.281250,    0.00000,   0.281250,    0.00000, -0.0937500, -0.0937500,   0.562500 }, // 21
      {    0.00000, -0.0468750,   0.140625,    0.00000,    0.00000,  0.0156250, -0.0468750,    0.00000,    0.00000,   0.281250,    0.00000,    0.00000,    0.00000, -0.0937500,   0.281250,    0.00000,    0.00000, -0.0937500,    0.00000,    0.00000,    0.00000,    0.00000,   0.562500,    0.00000,    0.00000,    0.00000,    0.00000 }, // 22
      {    0.00000,    0.00000,   0.140625, -0.0468750,    0.00000,    0.00000, -0.0468750,  0.0156250,    0.00000,    0.00000,   0.281250,    0.00000,    0.00000,    0.00000,   0.281250, -0.0937500,    0.00000,    0.00000, -0.0937500,    0.00000,    0.00000,    0.00000,    0.00000,   0.562500,    0.00000,    0.00000,    0.00000 }, // 23
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000, -0.0468750,    0.00000,   0.140625,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,  0.0156250,    0.00000, -0.0468750,    0.00000,   0.281250, -0.0937500,    0.00000,   0.281250,    0.00000, -0.0937500,   0.562500 }, // 24
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,  0.0156250, -0.0468750,   0.140625, -0.0468750,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000, -0.0937500,   0.281250,   0.281250, -0.0937500,    0.00000,   0.562500 }, // 25
      { 0.005859375, -0.017578125, 0.052734375, -0.017578125, -0.001953125, 0.005859375, -0.017578125, 0.005859375, -0.03515625, 0.10546875, 0.10546875, -0.03515625, 0.01171875, -0.03515625, 0.10546875, -0.03515625, 0.01171875, -0.03515625, -0.03515625, 0.01171875, 0.2109375, -0.0703125, 0.2109375, 0.2109375, -0.0703125, -0.0703125, 0.421875 } // 26
    },

    // embedding matrix for child 4
    {
      //  0     1     2     3     4     5     6     7    8  9  11  12  13  14  15  16  17  18  19  20   21  22  23  24  25  26
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 0
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 1
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000 }, // 2
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000 }, // 3
      {    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 4
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 5
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000 }, // 6
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 7
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.375000,  -0.125000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 8
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.375000,    0.00000,  -0.125000,    0.00000,    0.00000,   0.750000 }, // 9
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,  -0.125000,    0.00000,   0.375000,    0.00000,   0.750000 }, // 10
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.375000,    0.00000,    0.00000,  -0.125000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000 }, // 11
      {  -0.125000,    0.00000,    0.00000,    0.00000,   0.375000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 12
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,  -0.125000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.375000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 13
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,  -0.125000,    0.00000,    0.00000,    0.00000,    0.00000,   0.375000,   0.750000 }, // 14
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,  -0.125000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.375000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000 }, // 15
      {    0.00000,    0.00000,    0.00000,    0.00000,   0.375000,  -0.125000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 16
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.375000,    0.00000,  -0.125000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000 }, // 17
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,  -0.125000,    0.00000,   0.375000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000 }, // 18
      {    0.00000,    0.00000,    0.00000,    0.00000,   0.375000,    0.00000,    0.00000,  -0.125000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 19
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.140625, -0.0468750,  0.0156250, -0.0468750,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.281250, -0.0937500, -0.0937500,   0.281250,    0.00000,   0.562500 }, // 20
      { -0.0468750,  0.0156250,    0.00000,    0.00000,   0.140625, -0.0468750,    0.00000,    0.00000, -0.0937500,    0.00000,    0.00000,    0.00000,   0.281250, -0.0937500,    0.00000,    0.00000,   0.281250,    0.00000,    0.00000,    0.00000,    0.00000,   0.562500,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 21
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000, -0.0468750,    0.00000,  0.0156250,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.140625,    0.00000, -0.0468750,    0.00000, -0.0937500,   0.281250,    0.00000, -0.0937500,    0.00000,   0.281250,   0.562500 }, // 22
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,  0.0156250,    0.00000, -0.0468750,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000, -0.0468750,    0.00000,   0.140625, -0.0937500,    0.00000, -0.0937500,    0.00000,   0.281250,   0.281250,   0.562500 }, // 23
      { -0.0468750,    0.00000,    0.00000,  0.0156250,   0.140625,    0.00000,    0.00000, -0.0468750,    0.00000,    0.00000,    0.00000, -0.0937500,   0.281250,    0.00000,    0.00000, -0.0937500,    0.00000,    0.00000,    0.00000,   0.281250,    0.00000,    0.00000,    0.00000,    0.00000,   0.562500,    0.00000,    0.00000 }, // 24
      {    0.00000,    0.00000,    0.00000,    0.00000,   0.140625, -0.0468750,  0.0156250, -0.0468750,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.281250, -0.0937500, -0.0937500,   0.281250,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.562500,    0.00000 }, // 25
      { -0.017578125, 0.005859375, -0.001953125, 0.005859375, 0.052734375, -0.017578125, 0.005859375, -0.017578125, -0.03515625, 0.01171875, 0.01171875, -0.03515625, 0.10546875, -0.03515625, 0.01171875, -0.03515625, 0.10546875, -0.03515625, -0.03515625, 0.10546875, -0.0703125, 0.2109375, -0.0703125, -0.0703125, 0.2109375, 0.2109375, 0.421875 } // 26
    },

    // embedding matrix for child 5
    {
      //  0     1     2     3     4     5     6     7    8  9  11  12  13  14  15  16  17  18  19  20   21  22  23  24  25  26
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 0
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 1
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 2
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000 }, // 3
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 4
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 5
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 6
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000 }, // 7
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,  -0.125000,   0.375000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 8
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.375000,  -0.125000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 9
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.375000,    0.00000,  -0.125000,    0.00000,   0.750000 }, // 10
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.375000,    0.00000,  -0.125000,    0.00000,    0.00000,   0.750000 }, // 11
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,  -0.125000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.375000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 12
      {    0.00000,  -0.125000,    0.00000,    0.00000,    0.00000,   0.375000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 13
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,  -0.125000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.375000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 14
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,  -0.125000,    0.00000,    0.00000,    0.00000,    0.00000,   0.375000,   0.750000 }, // 15
      {    0.00000,    0.00000,    0.00000,    0.00000,  -0.125000,   0.375000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 16
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.375000,  -0.125000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 17
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.375000,    0.00000,  -0.125000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000 }, // 18
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.375000,    0.00000,  -0.125000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000 }, // 19
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000, -0.0468750,   0.140625, -0.0468750,  0.0156250,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.281250,   0.281250, -0.0937500, -0.0937500,    0.00000,   0.562500 }, // 20
      {  0.0156250, -0.0468750,    0.00000,    0.00000, -0.0468750,   0.140625,    0.00000,    0.00000, -0.0937500,    0.00000,    0.00000,    0.00000, -0.0937500,   0.281250,    0.00000,    0.00000,   0.281250,    0.00000,    0.00000,    0.00000,    0.00000,   0.562500,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 21
      {    0.00000, -0.0468750,  0.0156250,    0.00000,    0.00000,   0.140625, -0.0468750,    0.00000,    0.00000, -0.0937500,    0.00000,    0.00000,    0.00000,   0.281250, -0.0937500,    0.00000,    0.00000,   0.281250,    0.00000,    0.00000,    0.00000,    0.00000,   0.562500,    0.00000,    0.00000,    0.00000,    0.00000 }, // 22
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000, -0.0468750,    0.00000,  0.0156250,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.140625,    0.00000, -0.0468750, -0.0937500,    0.00000,   0.281250,    0.00000, -0.0937500,   0.281250,   0.562500 }, // 23
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000, -0.0468750,    0.00000,  0.0156250,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.140625,    0.00000, -0.0468750,    0.00000, -0.0937500,   0.281250,    0.00000, -0.0937500,    0.00000,   0.281250,   0.562500 }, // 24
      {    0.00000,    0.00000,    0.00000,    0.00000, -0.0468750,   0.140625, -0.0468750,  0.0156250,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.281250,   0.281250, -0.0937500, -0.0937500,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.562500,    0.00000 }, // 25
      { 0.005859375, -0.017578125, 0.005859375, -0.001953125, -0.017578125, 0.052734375, -0.017578125, 0.005859375, -0.03515625, -0.03515625, 0.01171875, 0.01171875, -0.03515625, 0.10546875, -0.03515625, 0.01171875, 0.10546875, 0.10546875, -0.03515625, -0.03515625, -0.0703125, 0.2109375, 0.2109375, -0.0703125, -0.0703125, 0.2109375, 0.421875 } // 26
    },

    // embedding matrix for child 6
    {
      //  0     1     2     3     4     5     6     7    8  9  11  12  13  14  15  16  17  18  19  20   21  22  23  24  25  26
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000 }, // 0
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000 }, // 1
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000 }, // 2
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 3
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 4
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000 }, // 5
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 6
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 7
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,  -0.125000,    0.00000,   0.375000,    0.00000,   0.750000 }, // 8
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,  -0.125000,    0.00000,   0.375000,    0.00000,    0.00000,   0.750000 }, // 9
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,  -0.125000,   0.375000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000 }, // 10
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,  -0.125000,    0.00000,    0.00000,   0.375000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000 }, // 11
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,  -0.125000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.375000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000 }, // 12
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,  -0.125000,    0.00000,    0.00000,    0.00000,    0.00000,   0.375000,   0.750000 }, // 13
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,  -0.125000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.375000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000 }, // 14
      {    0.00000,    0.00000,    0.00000,  -0.125000,    0.00000,    0.00000,    0.00000,   0.375000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 15
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,  -0.125000,    0.00000,   0.375000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000 }, // 16
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,  -0.125000,    0.00000,   0.375000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000 }, // 17
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,  -0.125000,   0.375000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 18
      {    0.00000,    0.00000,    0.00000,    0.00000,  -0.125000,    0.00000,    0.00000,   0.375000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 19
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000, -0.0468750,  0.0156250, -0.0468750,   0.140625,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000, -0.0937500, -0.0937500,   0.281250,   0.281250,    0.00000,   0.562500 }, // 20
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,  0.0156250,    0.00000, -0.0468750,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000, -0.0468750,    0.00000,   0.140625, -0.0937500,    0.00000, -0.0937500,    0.00000,   0.281250,   0.281250,   0.562500 }, // 21
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,  0.0156250,    0.00000, -0.0468750,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000, -0.0468750,    0.00000,   0.140625,    0.00000, -0.0937500, -0.0937500,    0.00000,   0.281250,    0.00000,   0.281250,   0.562500 }, // 22
      {    0.00000,    0.00000,  0.0156250, -0.0468750,    0.00000,    0.00000, -0.0468750,   0.140625,    0.00000,    0.00000, -0.0937500,    0.00000,    0.00000,    0.00000, -0.0937500,   0.281250,    0.00000,    0.00000,   0.281250,    0.00000,    0.00000,    0.00000,    0.00000,   0.562500,    0.00000,    0.00000,    0.00000 }, // 23
      {  0.0156250,    0.00000,    0.00000, -0.0468750, -0.0468750,    0.00000,    0.00000,   0.140625,    0.00000,    0.00000,    0.00000, -0.0937500, -0.0937500,    0.00000,    0.00000,   0.281250,    0.00000,    0.00000,    0.00000,   0.281250,    0.00000,    0.00000,    0.00000,    0.00000,   0.562500,    0.00000,    0.00000 }, // 24
      {    0.00000,    0.00000,    0.00000,    0.00000, -0.0468750,  0.0156250, -0.0468750,   0.140625,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000, -0.0937500, -0.0937500,   0.281250,   0.281250,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.562500,    0.00000 }, // 25
      { 0.005859375, -0.001953125, 0.005859375, -0.017578125, -0.017578125, 0.005859375, -0.017578125, 0.052734375, 0.01171875, 0.01171875, -0.03515625, -0.03515625, -0.03515625, 0.01171875, -0.03515625, 0.10546875, -0.03515625, -0.03515625, 0.10546875, 0.10546875, -0.0703125, -0.0703125, -0.0703125, 0.2109375, 0.2109375, 0.2109375, 0.421875 } // 26
    },

    // embedding matrix for child 7
    {
      //  0     1     2     3     4     5     6     7    8  9  11  12  13  14  15  16  17  18  19  20   21  22  23  24  25  26
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000 }, // 0
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 1
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 2
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000 }, // 3
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000 }, // 4
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 5
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 6
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 7
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.375000,    0.00000,  -0.125000,    0.00000,   0.750000 }, // 8
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,  -0.125000,   0.375000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 9
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.375000,  -0.125000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000 }, // 10
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,  -0.125000,    0.00000,   0.375000,    0.00000,    0.00000,   0.750000 }, // 11
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,  -0.125000,    0.00000,    0.00000,    0.00000,    0.00000,   0.375000,   0.750000 }, // 12
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,  -0.125000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.375000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 13
      {    0.00000,    0.00000,  -0.125000,    0.00000,    0.00000,    0.00000,   0.375000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 14
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,  -0.125000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.375000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000 }, // 15
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.375000,    0.00000,  -0.125000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000 }, // 16
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,  -0.125000,   0.375000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 17
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.375000,  -0.125000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 18
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,  -0.125000,    0.00000,   0.375000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000 }, // 19
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,  0.0156250, -0.0468750,   0.140625, -0.0468750,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000, -0.0937500,   0.281250,   0.281250, -0.0937500,    0.00000,   0.562500 }, // 20
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000, -0.0468750,    0.00000,  0.0156250,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.140625,    0.00000, -0.0468750, -0.0937500,    0.00000,   0.281250,    0.00000, -0.0937500,   0.281250,   0.562500 }, // 21
      {    0.00000,  0.0156250, -0.0468750,    0.00000,    0.00000, -0.0468750,   0.140625,    0.00000,    0.00000, -0.0937500,    0.00000,    0.00000,    0.00000, -0.0937500,   0.281250,    0.00000,    0.00000,   0.281250,    0.00000,    0.00000,    0.00000,    0.00000,   0.562500,    0.00000,    0.00000,    0.00000,    0.00000 }, // 22
      {    0.00000,    0.00000, -0.0468750,  0.0156250,    0.00000,    0.00000,   0.140625, -0.0468750,    0.00000,    0.00000, -0.0937500,    0.00000,    0.00000,    0.00000,   0.281250, -0.0937500,    0.00000,    0.00000,   0.281250,    0.00000,    0.00000,    0.00000,    0.00000,   0.562500,    0.00000,    0.00000,    0.00000 }, // 23
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,  0.0156250,    0.00000, -0.0468750,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000, -0.0468750,    0.00000,   0.140625,    0.00000, -0.0937500, -0.0937500,    0.00000,   0.281250,    0.00000,   0.281250,   0.562500 }, // 24
      {    0.00000,    0.00000,    0.00000,    0.00000,  0.0156250, -0.0468750,   0.140625, -0.0468750,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000, -0.0937500,   0.281250,   0.281250, -0.0937500,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.562500,    0.00000 }, // 25
      { -0.001953125, 0.005859375, -0.017578125, 0.005859375, 0.005859375, -0.017578125, 0.052734375, -0.017578125, 0.01171875, -0.03515625, -0.03515625, 0.01171875, 0.01171875, -0.03515625, 0.10546875, -0.03515625, -0.03515625, 0.10546875, 0.10546875, -0.03515625, -0.0703125, -0.0703125, 0.2109375, 0.2109375, -0.0703125, 0.2109375, 0.421875 } // 26
    }
  };

#endif

} // namespace libMesh
