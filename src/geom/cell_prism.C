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
#include "libmesh/cell_prism.h"
#include "libmesh/cell_prism6.h"
#include "libmesh/face_quad4.h"
#include "libmesh/face_tri3.h"

namespace libMesh
{


// ------------------------------------------------------------
// Prism class member functions
dof_id_type Prism::key (const unsigned int s) const
{
  libmesh_assert_less (s, this->n_sides());

  switch (s)
    {
    case 0:  // the triangular face at z=0

      return
        this->compute_key (this->node(0),
                           this->node(2),
                           this->node(1));

    case 1:  // the quad face at y=0

      return
        this->compute_key (this->node(0),
                           this->node(1),
                           this->node(4),
                           this->node(3));

    case 2:  // the other quad face

      return
        this->compute_key (this->node(1),
                           this->node(2),
                           this->node(5),
                           this->node(4));

    case 3: // the quad face at x=0

      return
        this->compute_key (this->node(2),
                           this->node(0),
                           this->node(3),
                           this->node(5));
    case 4: // the triangular face at z=1

      return
        this->compute_key (this->node(3),
                           this->node(4),
                           this->node(5));

    default:
      libmesh_error_msg("Invalid side " << s);
    }

  libmesh_error_msg("We'll never get here!");
  return 0;
}



AutoPtr<Elem> Prism::side (const unsigned int i) const
{
  libmesh_assert_less (i, this->n_sides());

  Elem* face = NULL;

  switch (i)
    {
    case 0:  // the triangular face at z=0
      {
        face = new Tri3;

        face->set_node(0) = this->get_node(0);
        face->set_node(1) = this->get_node(2);
        face->set_node(2) = this->get_node(1);

        break;
      }
    case 1:  // the quad face at y=0
      {
        face = new Quad4;

        face->set_node(0) = this->get_node(0);
        face->set_node(1) = this->get_node(1);
        face->set_node(2) = this->get_node(4);
        face->set_node(3) = this->get_node(3);

        break;
      }
    case 2:  // the other quad face
      {
        face = new Quad4;

        face->set_node(0) = this->get_node(1);
        face->set_node(1) = this->get_node(2);
        face->set_node(2) = this->get_node(5);
        face->set_node(3) = this->get_node(4);

        break;
      }
    case 3: // the quad face at x=0
      {
        face = new Quad4;

        face->set_node(0) = this->get_node(2);
        face->set_node(1) = this->get_node(0);
        face->set_node(2) = this->get_node(3);
        face->set_node(3) = this->get_node(5);

        break;
      }
    case 4: // the triangular face at z=1
      {
        face = new Tri3;

        face->set_node(0) = this->get_node(3);
        face->set_node(1) = this->get_node(4);
        face->set_node(2) = this->get_node(5);

        break;
      }
    default:
      libmesh_error_msg("Invalid side i = " << i);
    }

  return AutoPtr<Elem>(face);
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

  return (is_node_on_side(Prism6::edge_nodes_map[e][0],s) &&
          is_node_on_side(Prism6::edge_nodes_map[e][1],s));
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
