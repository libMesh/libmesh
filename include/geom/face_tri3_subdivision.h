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



#ifndef LIBMESH_FACE_TRI3_SUBDIVISION_H
#define LIBMESH_FACE_TRI3_SUBDIVISION_H


// Local includes
#include "libmesh/face_tri3.h"

namespace libMesh
{

/**
 * The \p Tri3Subdivision element is a three-noded subdivision surface
 * shell element used in mechanics calculations.  See also,
 * miscellaneous_ex11.
 *
 * \author Roman Vetter
 * \author Norbert Stoop
 * \date 2014
 * \brief A surface shell element used in mechanics calculations.
 */
class Tri3Subdivision libmesh_final : public Tri3
{
public:

  /**
   * Constructor without parent specification.
   */
  Tri3Subdivision() :
    Tri3(), _subdivision_updated(false), _is_ghost(false) {}

  /**
   * Constructor with parent specification.
   */
  Tri3Subdivision(Elem * p);

  /**
   * \returns \p TRI3SUBDIVISION.
   */
  virtual ElemType type () const libmesh_override { return TRI3SUBDIVISION; }

  /**
   * \returns \p true if the element map is definitely affine within
   * numerical tolerances.
   */
  virtual bool has_affine_map () const libmesh_override { return false; }

  /**
   * \returns \p true if the Lagrange shape functions on this element
   * are linear.
   */
  virtual bool is_linear () const libmesh_override { return false; }

  /**
   * \returns FOURTH.
   */
  virtual Order default_order() const libmesh_override { return FOURTH; }

  /**
   * Prepares the element for use by reordering the nodes such that
   * the irregular node (valence != 6), if there is one, is the first.
   * The nodes are ordered once in advance for efficiency.
   */
  void prepare_subdivision_properties();

  /**
   * \returns \p true if the subdivision element is ready for use,
   * i.e. the nodes have been reordered.
   */
  bool is_subdivision_updated() const { return _subdivision_updated; }

  /**
   * \returns A pointer to the node whose ordered id is \p node_id.
   */
  Node * get_ordered_node(unsigned int node_id) const;

  /**
   * \returns The number of nodes connected to the ordered node
   * whose id is \p node_id.
   */
  unsigned int get_ordered_valence(unsigned int node_id) const;

  /**
   * \returns The order number of the node whose unordered id is
   * \p node_id. This is the inverse of an \p _ordered_nodes lookup.
   */
  unsigned int local_node_number(unsigned int node_id) const;

  /**
   * \returns \p true if the element is a ghost element.
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
  Node * _ordered_nodes[3];

  /**
   * \p true if the subdivision element is ready for use,
   * i.e. the nodes have been reordered.
   */
  bool _subdivision_updated;

  /**
   * \p true if the element is a ghost element
   * (e.g. for boundary conditions).
   */
  bool _is_ghost;
};


} // namespace libMesh

#endif // LIBMESH_FACE_TRI3_SUBDIVISION_H
