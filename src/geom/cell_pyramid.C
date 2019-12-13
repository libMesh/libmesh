// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
// Pyramid class static member initializations


// We need to require C++11...
const Real Pyramid::_master_points[14][3] =
  {
    {-1, -1, 0},
    {1, -1, 0},
    {1, 1, 0},
    {-1, 1, 0},
    {0, 0, 1},
    {0, -1, 0},
    {1, 0, 0},
    {0, 1, 0},
    {-1, 0, 0},
    {0, -0.5, 0.5},
    {0.5, 0, 0.5},
    {0, 0.5, 0.5},
    {-0.5, 0, 0.5}
  };




// ------------------------------------------------------------
// Pyramid class member functions
dof_id_type Pyramid::key (const unsigned int s) const
{
  libmesh_assert_less (s, this->n_sides());

  switch (s)
    {
    case 0: // triangular face 1
    case 1: // triangular face 2
    case 2: // triangular face 3
    case 3: // triangular face 4
      return this->compute_key (this->node_id(Pyramid5::side_nodes_map[s][0]),
                                this->node_id(Pyramid5::side_nodes_map[s][1]),
                                this->node_id(Pyramid5::side_nodes_map[s][2]));

    case 4:  // the quad face at z=0
      return this->compute_key (this->node_id(Pyramid5::side_nodes_map[s][0]),
                                this->node_id(Pyramid5::side_nodes_map[s][1]),
                                this->node_id(Pyramid5::side_nodes_map[s][2]),
                                this->node_id(Pyramid5::side_nodes_map[s][3]));

    default:
      libmesh_error_msg("Invalid side s = " << s);
    }
}



unsigned int Pyramid::which_node_am_i(unsigned int side,
                                      unsigned int side_node) const
{
  libmesh_assert_less (side, this->n_sides());

  // Never more than 4 nodes per side.
  libmesh_assert_less(side_node, 4);

  // Some sides have 3 nodes.
  libmesh_assert(side == 4 || side_node < 3);

  return Pyramid5::side_nodes_map[side][side_node];
}



std::unique_ptr<Elem> Pyramid::side_ptr (const unsigned int i)
{
  libmesh_assert_less (i, this->n_sides());

  // Return value
  std::unique_ptr<Elem> face;

  // Set up the type of element
  switch (i)
    {
    case 0: // triangular face 1
    case 1: // triangular face 2
    case 2: // triangular face 3
    case 3: // triangular face 4
      {
        face = libmesh_make_unique<Tri3>();
        break;
      }
    case 4:  // the quad face at z=0
      {
        face = libmesh_make_unique<Quad4>();
        break;
      }
    default:
      libmesh_error_msg("Invalid side i = " << i);
    }

  // Set the nodes
  for (auto n : face->node_index_range())
    face->set_node(n) = this->node_ptr(Pyramid5::side_nodes_map[i][n]);

  return face;
}



void Pyramid::side_ptr (std::unique_ptr<Elem> & side,
                        const unsigned int i)
{
  libmesh_assert_less (i, this->n_sides());

  switch (i)
    {
    case 0: // triangular face 1
    case 1: // triangular face 2
    case 2: // triangular face 3
    case 3: // triangular face 4
      {
        if (!side.get() || side->type() != TRI3)
          {
            side = this->side_ptr(i);
            return;
          }
        break;
      }

    case 4:  // the quad face at z=0
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
    side->set_node(n) = this->node_ptr(Pyramid5::side_nodes_map[i][n]);
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
