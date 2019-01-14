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



// Local includes
#include "libmesh/edge.h"
#include "libmesh/node_elem.h"

namespace libMesh
{

unsigned int Edge::which_node_am_i(unsigned int side,
                                   unsigned int /*side_node*/) const
{
  libmesh_assert_less (side, this->n_sides());
  return side;
}



std::unique_ptr<Elem> Edge::side_ptr (const unsigned int i)
{
  libmesh_assert_less (i, 2);
  std::unique_ptr<Elem> nodeelem = libmesh_make_unique<NodeElem>(this);
  nodeelem->set_node(0) = this->node_ptr(i);
  return nodeelem;
}


void Edge::side_ptr (std::unique_ptr<Elem> & side,
                           const unsigned int i)
{
  libmesh_assert_less (i, this->n_sides());

  if (!side.get() || side->type() != NODEELEM)
    side = this->build_side_ptr(i, false);
  else
    {
      side->subdomain_id() = this->subdomain_id();

      side->set_node(0) = this->node_ptr(i);
    }
}




std::unique_ptr<Elem> Edge::build_side_ptr (const unsigned int i, bool)
{
  libmesh_assert_less (i, 2);
  std::unique_ptr<Elem> nodeelem = libmesh_make_unique<NodeElem>(this);
  nodeelem->set_node(0) = this->node_ptr(i);
  return nodeelem;
}


void Edge::build_side_ptr (std::unique_ptr<Elem> & side,
                           const unsigned int i)
{
  this->side_ptr(side, i);
}




bool Edge::is_child_on_side(const unsigned int c,
                            const unsigned int s) const
{
  libmesh_assert_less (c, this->n_children());
  libmesh_assert_less (s, this->n_sides());

  return (c == s);
}



unsigned int Edge::opposite_side(const unsigned int side_in) const
{
  libmesh_assert_less (side_in, 2);
  return 1 - side_in;
}



unsigned int Edge::opposite_node(const unsigned int node_in,
                                 const unsigned int libmesh_dbg_var(side_in)) const
{
  libmesh_assert_less (node_in, 2);
  libmesh_assert_less (side_in, this->n_sides());
  libmesh_assert(this->is_node_on_side(node_in, side_in));

  return 1 - node_in;
}

std::vector<unsigned>
Edge::nodes_on_side(const unsigned int s) const
{
  libmesh_assert_less(s, 2);
  return {s};
}

} // namespace libMesh
