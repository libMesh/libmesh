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
#include "libmesh/cell_prism15.h"
#include "libmesh/edge_edge3.h"
#include "libmesh/face_quad8.h"
#include "libmesh/face_tri6.h"

namespace libMesh
{



// ------------------------------------------------------------
// Prism15 class static member initializations
const unsigned int Prism15::side_nodes_map[5][8] =
  {
    {0, 2, 1,  8,  7,  6, 99, 99}, // Side 0
    {0, 1, 4,  3,  6, 10, 12,  9}, // Side 1
    {1, 2, 5,  4,  7, 11, 13, 10}, // Side 2
    {2, 0, 3,  5,  8,  9, 14, 11}, // Side 3
    {3, 4, 5, 12, 13, 14, 99, 99}  // Side 4
  };

const unsigned int Prism15::edge_nodes_map[9][3] =
  {
    {0, 1, 6},  // Side 0
    {1, 2, 7},  // Side 1
    {0, 2, 8},  // Side 2
    {0, 3, 9},  // Side 3
    {1, 4, 10}, // Side 4
    {2, 5, 11}, // Side 5
    {3, 4, 12}, // Side 6
    {4, 5, 13}, // Side 7
    {3, 5, 14}  // Side 8
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
  for (unsigned int i = 0; i != 8; ++i)
    if (side_nodes_map[s][i] == n)
      return true;
  return false;
}

bool Prism15::is_node_on_edge(const unsigned int n,
                              const unsigned int e) const
{
  libmesh_assert_less (e, n_edges());
  for (unsigned int i = 0; i != 3; ++i)
    if (edge_nodes_map[e][i] == n)
      return true;
  return false;
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



UniquePtr<Elem> Prism15::build_side (const unsigned int i,
                                     bool proxy) const
{
  libmesh_assert_less (i, this->n_sides());

  if (proxy)
    {
      switch (i)
        {
        case 0:  // the triangular face at z=-1
        case 4:
          return UniquePtr<Elem>(new Side<Tri6,Prism15>(this,i));

        case 1:
        case 2:
        case 3:
          return UniquePtr<Elem>(new Side<Quad8,Prism15>(this,i));

        default:
          libmesh_error_msg("Invalid side i = " << i);
        }
    }

  else
    {
      // Create NULL pointer to be initialized, returned later.
      Elem* face = NULL;

      switch (i)
        {
        case 0:  // the triangular face at z=-1
          {
            face = new Tri6;

            face->set_node(0) = this->get_node(0);
            face->set_node(1) = this->get_node(2);
            face->set_node(2) = this->get_node(1);
            face->set_node(3) = this->get_node(8);
            face->set_node(4) = this->get_node(7);
            face->set_node(5) = this->get_node(6);

            break;
          }
        case 1:  // the quad face at y=0
          {
            face = new Quad8;

            face->set_node(0) = this->get_node(0);
            face->set_node(1) = this->get_node(1);
            face->set_node(2) = this->get_node(4);
            face->set_node(3) = this->get_node(3);
            face->set_node(4) = this->get_node(6);
            face->set_node(5) = this->get_node(10);
            face->set_node(6) = this->get_node(12);
            face->set_node(7) = this->get_node(9);

            break;
          }
        case 2:  // the other quad face
          {
            face = new Quad8;

            face->set_node(0) = this->get_node(1);
            face->set_node(1) = this->get_node(2);
            face->set_node(2) = this->get_node(5);
            face->set_node(3) = this->get_node(4);
            face->set_node(4) = this->get_node(7);
            face->set_node(5) = this->get_node(11);
            face->set_node(6) = this->get_node(13);
            face->set_node(7) = this->get_node(10);

            break;
          }
        case 3: // the quad face at x=0
          {
            face = new Quad8;

            face->set_node(0) = this->get_node(2);
            face->set_node(1) = this->get_node(0);
            face->set_node(2) = this->get_node(3);
            face->set_node(3) = this->get_node(5);
            face->set_node(4) = this->get_node(8);
            face->set_node(5) = this->get_node(9);
            face->set_node(6) = this->get_node(14);
            face->set_node(7) = this->get_node(11);

            break;
          }
        case 4: // the triangular face at z=1
          {
            face = new Tri6;

            face->set_node(0) = this->get_node(3);
            face->set_node(1) = this->get_node(4);
            face->set_node(2) = this->get_node(5);
            face->set_node(3) = this->get_node(12);
            face->set_node(4) = this->get_node(13);
            face->set_node(5) = this->get_node(14);

            break;
          }
        default:
          libmesh_error_msg("Invalid side i = " << i);
        }

      face->subdomain_id() = this->subdomain_id();
      return UniquePtr<Elem>(face);
    }

  libmesh_error_msg("We'll never get here!");
  return UniquePtr<Elem>();
}


UniquePtr<Elem> Prism15::build_edge (const unsigned int i) const
{
  libmesh_assert_less (i, this->n_edges());

  return UniquePtr<Elem>(new SideEdge<Edge3,Prism15>(this,i));
}


void Prism15::connectivity(const unsigned int libmesh_dbg_var(sc),
                           const IOPackage iop,
                           std::vector<dof_id_type>& conn) const
{
  libmesh_assert(_nodes);
  libmesh_assert_less (sc, this->n_sub_elem());
  libmesh_assert_not_equal_to (iop, INVALID_IO_PACKAGE);

  switch (iop)
    {
    case TECPLOT:
      {
        conn.resize(8);
        conn[0] = this->node(0)+1;
        conn[1] = this->node(1)+1;
        conn[2] = this->node(2)+1;
        conn[3] = this->node(2)+1;
        conn[4] = this->node(3)+1;
        conn[5] = this->node(4)+1;
        conn[6] = this->node(5)+1;
        conn[7] = this->node(5)+1;
        return;
      }

    case VTK:
      {
        /*
          conn.resize(6);
          conn[0] = this->node(0);
          conn[1] = this->node(2);
          conn[2] = this->node(1);
          conn[3] = this->node(3);
          conn[4] = this->node(5);
          conn[5] = this->node(4);
        */

        // VTK's VTK_QUADRATIC_WEDGE first 9 nodes match, then their
        // middle and top layers of mid-edge nodes are reversed from
        // LibMesh's.
        conn.resize(15);
        for (unsigned i=0; i<9; ++i)
          conn[i] = this->node(i);

        // top "ring" of mid-edge nodes
        conn[9]  = this->node(12);
        conn[10] = this->node(13);
        conn[11] = this->node(14);

        // middle "ring" of mid-edge nodes
        conn[12] = this->node(9);
        conn[13] = this->node(10);
        conn[14] = this->node(11);


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

} // namespace libMesh
