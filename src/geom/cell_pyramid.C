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
#include "libmesh/cell_pyramid.h"
#include "libmesh/cell_pyramid5.h"
#include "libmesh/face_tri3.h"
#include "libmesh/face_quad4.h"

namespace libMesh
{


// ------------------------------------------------------------
// Pyramid class member functions
dof_id_type Pyramid::key (const unsigned int s) const
{
  libmesh_assert_less (s, this->n_sides());

  switch (s)
    {
    case 0:  // triangular face 1

      return
        this->compute_key (this->node(0),
                           this->node(1),
                           this->node(4));

    case 1:  // triangular face 2

      return
        this->compute_key (this->node(1),
                           this->node(2),
                           this->node(4));

    case 2:  // triangular face 3

      return
        this->compute_key (this->node(2),
                           this->node(3),
                           this->node(4));

    case 3:  // triangular face 4

      return
        this->compute_key (this->node(3),
                           this->node(0),
                           this->node(4));

    case 4:  // the quad face at z=0

      return
        this->compute_key (this->node(0),
                           this->node(3),
                           this->node(2),
                           this->node(1));

    default:
      libmesh_error_msg("Invalid side s = " << s);
    }

  libmesh_error_msg("We'll never get here!");
  return 0;
}



UniquePtr<Elem> Pyramid::side (const unsigned int i) const
{
  libmesh_assert_less (i, this->n_sides());

  // To be returned wrapped in an UniquePtr
  Elem* face = NULL;

  switch (i)
    {
    case 0:  // triangular face 1
      {
        face = new Tri3;

        face->set_node(0) = this->get_node(0);
        face->set_node(1) = this->get_node(1);
        face->set_node(2) = this->get_node(4);

        break;
      }
    case 1:  // triangular face 2
      {
        face = new Tri3;

        face->set_node(0) = this->get_node(1);
        face->set_node(1) = this->get_node(2);
        face->set_node(2) = this->get_node(4);

        break;
      }
    case 2:  // triangular face 3
      {
        face = new Tri3;

        face->set_node(0) = this->get_node(2);
        face->set_node(1) = this->get_node(3);
        face->set_node(2) = this->get_node(4);

        break;
      }
    case 3:  // triangular face 4
      {
        face = new Tri3;

        face->set_node(0) = this->get_node(3);
        face->set_node(1) = this->get_node(0);
        face->set_node(2) = this->get_node(4);

        break;
      }
    case 4:  // the quad face at z=0
      {
        face = new Quad4;

        face->set_node(0) = this->get_node(0);
        face->set_node(1) = this->get_node(3);
        face->set_node(2) = this->get_node(2);
        face->set_node(3) = this->get_node(1);

        break;
      }
    default:
      libmesh_error_msg("Invalid side i = " << i);
    }

  return UniquePtr<Elem>(face);
}



bool Pyramid::is_child_on_side(const unsigned int c,
                               const unsigned int s) const
{
  libmesh_assert_less (c, this->n_children());
  libmesh_assert_less (s, this->n_sides());

  for (unsigned int i = 0; i != 4; ++i)
    if (Pyramid5::side_nodes_map[s][i] == c)
      return true;
  return false;
}



bool Pyramid::is_edge_on_side(const unsigned int e,
                              const unsigned int s) const
{
  libmesh_assert_less (e, this->n_edges());
  libmesh_assert_less (s, this->n_sides());

  return (is_node_on_side(Pyramid5::edge_nodes_map[e][0],s) &&
          is_node_on_side(Pyramid5::edge_nodes_map[e][1],s));
}



} // namespace libMesh
