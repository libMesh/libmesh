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



#ifndef LIBMESH_FACE_TRI3_SD_H
#define LIBMESH_FACE_TRI3_SD_H


// Local includes
#include "libmesh/face_tri3.h"

// C++ includes

namespace libMesh
{

/**
 * The \p Tri3 is an element in 2D composed of 3 nodes.
 */

// ------------------------------------------------------------
// Tri3SD class definition
class Tri3SD : public Tri3
{
public:

  /**
   * Constructor without parent specification.
   */
  Tri3SD() :
    Tri3(), _subdiv_updated(false), _is_ghost(false) {}

  /**
   * Constructor with parent specification.
   */
  Tri3SD(Elem *p);

  /**
   * @returns \p TRI3SD
   */
  ElemType type () const { return TRI3SD; }
  
  /**
   * @returns true iff the element map is definitely affine within
   * numerical tolerances
   */
  virtual bool has_affine_map () const { return false; }

  /**
   * @returns true iff the Lagrange shape functions on this element
   * are linear
   */
  virtual bool is_linear () const { return false; }

  /**
   * @returns FOURTH
   */
  Order default_order() const { return FOURTH; }

  /**
   * Prepares the element for use by reordering the nodes such that
   * the irregular node (valence != 6), if there is one, is the first.
   * The nodes are ordered once in advance for efficiency.
   */
  void prepare_subdiv_properties();

  /**
   * @returns \p true iff the subdivision element is ready for use,
   * i.e. the nodes have been reordered.
   */
  bool is_subdiv_updated() const { return _subdiv_updated; }

  /**
   * @returns a pointer to the node whose ordered id is \p node_id.
   */
  Node* get_ordered_node(unsigned int node_id) const;

  /**
   * @returns the number of nodes connected to the ordered node
   * whose id is \p node_id.
   */
  unsigned int get_ordered_valence(unsigned int node_id) const;

  /**
   * @returns the order number of the node whose unordered id is
   * \p node_id. This is the inverse of an \p _ordered_nodes lookup.
   */
  unsigned int local_node_number(unsigned int node_id) const;

  /**
   * @returns \p true iff the element is a ghost element.
   */
  bool is_ghost() const { return _is_ghost; }

  /**
   * Sets the boolean flag identifying ghost elements.
   */
  void set_ghost(bool ghosted) { _is_ghost = ghosted; }
  
private:

  /**
   * A list containing the ordered nodes such that the irregular
   * node (valence != 6), if there is one, is the first.
   */
  Node* _ordered_nodes[3];

  /**
   * \p true iff the subdivision element is ready for use,
   * i.e. the nodes have been reordered.
   */
  bool _subdiv_updated;

  /**
   * \p true iff the element is a ghost element
   * (e.g. for boundary conditions).
   */
  bool _is_ghost;
};


} // namespace libMesh

#endif // LIBMESH_FACE_TRI3_SD_H
