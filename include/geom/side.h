// The libMesh Finite Element Library.
// Copyright (C) 2002-2016 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_SIDE_H
#define LIBMESH_SIDE_H

// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/elem.h"

// C++ includes
#include <cstddef>

namespace libMesh
{

// Forward declarations
class Point;
class Node;


/**
 * This defines the \p Side class.  A \p Side is basically a proxy
 * (or stand-in replacement) class for an element's side.  It acts
 * like a standard \p Elem, but allocates no additional memory for
 * storing connectivity.  Instead, its nodes are mapped directly from
 * the parent element (the element for which the side is created).
 * Similarly, you cannot access the neighbors of a side since it
 * does not store any.
 *
 * \author  Benjamin S. Kirk
 */
template <class SideType, class ParentType>
class Side : public SideType
{
public:

  /**
   * Constructor.  Creates a side from an element.
   */
  Side (const Elem * parent_in,
        const unsigned int side_in) :
    SideType(const_cast<Elem *>(parent_in)),
    _side_number(side_in)
  {
    libmesh_assert(parent_in);
    // may not be true when building infinite element sides
    // libmesh_assert_less (_side_number, this->parent()->n_sides());
    libmesh_assert_equal_to ((this->dim()+1), this->parent()->dim());

    for (unsigned int n=0; n != this->n_nodes(); ++n)
      this->_nodes[n] = this->parent()->get_node
        (ParentType::side_nodes_map[_side_number][n]);
  }

  /**
   * Setting a side node changes the node on the parent
   */
  virtual Node * & set_node (const unsigned int i) libmesh_override
  {
    libmesh_assert_less (i, this->n_nodes());
    return this->parent()->set_node (ParentType::side_nodes_map[_side_number][i]);
  }

  /**
   * Sides effectively do not have sides
   */
  virtual unsigned int n_sides () const libmesh_override
  { return 0; }

  virtual bool is_child_on_side(const unsigned int,
                                const unsigned int) const libmesh_override
  { libmesh_not_implemented(); return false; }


private:

  /**
   * The side on the parent element
   */
  const unsigned int _side_number;
};



/**
 * This defines the \p SideEdge class.  Like \p Side, \p SideEdge is basically
 * a proxy (or stand-in replacement) class, this time for an element's edge.
 * It acts like a standard \p Elem, but allocates no additional memory for
 * storing connectivity.  Instead, its nodes are mapped directly from the
 * parent element (the element for which the side is created).  Similarly, you
 * cannot access the neighbors of a side since it does not store any.
 *
 * \author  Roy H. Stogner
 */
template <class EdgeType, class ParentType>
class SideEdge : public EdgeType
{
public:

  /**
   * Constructor.  Creates a side from an element.
   */
  SideEdge (const Elem * my_parent,
            const unsigned int my_edge) :
    EdgeType(const_cast<Elem *>(my_parent)),
    _edge_number(my_edge)
  {
    libmesh_assert(my_parent);
    libmesh_assert_less (_edge_number, this->parent()->n_edges());
    libmesh_assert_equal_to (this->dim(), 1);

    for (unsigned int n=0; n != this->n_nodes(); ++n)
      this->_nodes[n] = this->parent()->get_node
        (ParentType::edge_nodes_map[_edge_number][n]);
  }

  /**
   * Setting an edge node changes the node on the parent
   */
  virtual Node * & set_node (const unsigned int i) libmesh_override
  {
    libmesh_assert_less (i, this->n_nodes());
    return this->parent()->set_node (ParentType::edge_nodes_map[_edge_number][i]);
  }

  /**
   * @returns 0. Edges effectively do not have sides, so
   * don't even ask!
   */
  virtual unsigned int n_sides () const libmesh_override { return 0; }

private:

  /**
   * The side on the parent element
   */
  const unsigned int _edge_number;
};


} // namespace libMesh

#endif // LIBMESH_SIDE_H
