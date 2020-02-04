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

#ifndef LIBMESH_CELL_PRISM_IMPL_H
#define LIBMESH_CELL_PRISM_IMPL_H

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
template <typename RealType>
dof_id_type PrismTempl<RealType>::key (const unsigned int s) const
{
  libmesh_assert_less (s, this->n_sides());

  switch (s)
    {
    case 0: // the triangular face at z=0
    case 4: // the triangular face at z=1
      return this->compute_key (this->node_id(Prism6::side_nodes_map[s][0]),
                                this->node_id(Prism6::side_nodes_map[s][1]),
                                this->node_id(Prism6::side_nodes_map[s][2]));

    case 1: // the quad face at y=0
    case 2: // the other quad face
    case 3: // the quad face at x=0
      return this->compute_key (this->node_id(Prism6::side_nodes_map[s][0]),
                                this->node_id(Prism6::side_nodes_map[s][1]),
                                this->node_id(Prism6::side_nodes_map[s][2]),
                                this->node_id(Prism6::side_nodes_map[s][3]));

    default:
      libmesh_error_msg("Invalid side " << s);
    }
}



template <typename RealType>
unsigned int PrismTempl<RealType>::which_node_am_i(unsigned int side,
                                    unsigned int side_node) const
{
  libmesh_assert_less (side, this->n_sides());

  // Never more than 4 nodes per side.
  libmesh_assert_less(side_node, 4);

  // Some sides have 3 nodes.
  libmesh_assert(!(side==0 || side==4) || side_node < 3);

  return Prism6::side_nodes_map[side][side_node];
}



template <typename RealType>
std::unique_ptr<ElemTempl<RealType>> PrismTempl<RealType>::side_ptr (const unsigned int i)
{
  libmesh_assert_less (i, this->n_sides());

  std::unique_ptr<Elem> face;

  // Set up the type of element
  switch (i)
    {
    case 0: // the triangular face at z=0
    case 4: // the triangular face at z=1
      {
        face = libmesh_make_unique<Tri3>();
        break;
      }
    case 1: // the quad face at y=0
    case 2: // the other quad face
    case 3: // the quad face at x=0
      {
        face = libmesh_make_unique<Quad4>();
        break;
      }
    default:
      libmesh_error_msg("Invalid side i = " << i);
    }

  // Set the nodes
  for (auto n : face->node_index_range())
    face->set_node(n) = this->node_ptr(Prism6::side_nodes_map[i][n]);

  return face;
}



template <typename RealType>
void PrismTempl<RealType>::side_ptr (std::unique_ptr<Elem> & side,
                      const unsigned int i)
{
  libmesh_assert_less (i, this->n_sides());

  switch (i)
    {
      // the base face
    case 0: // the triangular face at z=0
    case 4: // the triangular face at z=1
      {
        if (!side.get() || side->type() != TRI3)
          {
            side = this->side_ptr(i);
            return;
          }
        break;
      }

    case 1: // the quad face at y=0
    case 2: // the other quad face
    case 3: // the quad face at x=0
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
    side->set_node(n) = this->node_ptr(Prism6::side_nodes_map[i][n]);
}



template <typename RealType>
bool PrismTempl<RealType>::is_child_on_side(const unsigned int c,
                             const unsigned int s) const
{
  libmesh_assert_less (c, this->n_children());
  libmesh_assert_less (s, this->n_sides());

  for (unsigned int i = 0; i != 4; ++i)
    if (Prism6::side_elems_map[s][i] == c)
      return true;
  return false;
}



template <typename RealType>
bool PrismTempl<RealType>::is_edge_on_side(const unsigned int e,
                            const unsigned int s) const
{
  libmesh_assert_less (e, this->n_edges());
  libmesh_assert_less (s, this->n_sides());

  return (this->is_node_on_side(Prism6::edge_nodes_map[e][0],s) &&
          this->is_node_on_side(Prism6::edge_nodes_map[e][1],s));
}

} // namespace libMesh

#endif // LIBMESH_CELL_PRISM_IMPL_H
