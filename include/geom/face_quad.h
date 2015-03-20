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



#ifndef LIBMESH_FACE_QUAD_H
#define LIBMESH_FACE_QUAD_H


// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/face.h"

// C++ includes

namespace libMesh
{


// Forward declarations


/**
 * The \p QUAD is an element in 2D composed of 4 sides.
 * It looks like this:
 * \verbatim
 *
 *         -----------
 *        |           |
 *        |           |
 *        |           |
 *        |           |
 *        |           |
 *         -----------
 *
 * \endverbatim
 */
class Quad : public Face
{
public:

  /**
   * Default quadrilateral element, takes number of nodes and
   * parent. Derived classes implement 'true' elements.
   */
  Quad (const unsigned int nn, Elem* p, Node** nodelinkdata) :
    Face(nn, Quad::n_sides(), p, _elemlinks_data, nodelinkdata) {}

  /**
   * @returns 4.  All quad-derivatives are guaranteed to have at
   * least 4 nodes.
   */
  unsigned int n_nodes() const { return 4; }

  /**
   * @returns 4
   */
  unsigned int n_sides() const { return 4; }

  /**
   * @returns 4.  All quadrilaterals have 4 vertices.
   */
  unsigned int n_vertices() const { return 4; }

  /**
   * @returns 4.  All quadrilaterals have 4 edges.
   */
  unsigned int n_edges() const { return 4; }

  /**
   * @returns 4
   */
  unsigned int n_children() const { return 4; }

  /*
   * @returns true iff the specified child is on the
   * specified side
   */
  virtual bool is_child_on_side(const unsigned int c,
                                const unsigned int s) const;

  /**
   * @returns the side number opposite to \p s (for a tensor product
   * element), or throws an error otherwise.
   */
  virtual unsigned int opposite_side(const unsigned int s) const;

  /**
   * @returns the local node number for the node opposite to node n
   * on side \p opposite_side(s) (for a tensor product element), or
   * throws an error otherwise.
   */
  virtual unsigned int opposite_node(const unsigned int n,
                                     const unsigned int s) const;

  /**
   * @returns an id associated with the \p s side of this element.
   * The id is not necessariy unique, but should be close.  This is
   * particularly useful in the \p MeshBase::find_neighbors() routine.
   */
  dof_id_type key (const unsigned int s) const;

  /**
   * @returns a primitive (2-noded) edge for
   * edge i.
   */
  UniquePtr<Elem> side (const unsigned int i) const;

  /**
   * Based on the quality metric q specified by the user,
   * returns a quantitative assessment of element quality.
   */
  Real quality (const ElemQuality q) const;

  /**
   * Returns the suggested quality bounds for
   * the hex based on quality measure q.  These are
   * the values suggested by the CUBIT User's Manual.
   */
  std::pair<Real, Real> qual_bounds (const ElemQuality q) const;

protected:

  /**
   * Data for links to parent/neighbor/interior_parent elements.
   */
  Elem* _elemlinks_data[5+(LIBMESH_DIM>2)];

  /**
   * Matrix that tells which vertices define the location
   * of mid-side (or second-order) nodes.  Since most
   * second-order nodes are identical for \p Quad8 and \p Quad9,
   * we keep this matrix here in \p Quad
   */
  static const unsigned short int _second_order_adjacent_vertices[4][2];

  /**
   * Vector that names a child sharing each second order node.
   */
  static const unsigned short int _second_order_vertex_child_number[9];

  /**
   * Vector that names the child vertex index for each second order node.
   */
  static const unsigned short int _second_order_vertex_child_index[9];
};


// ------------------------------------------------------------
// Quad class member functions


} // namespace libMesh

#endif // LIBMESH_FACE_QUAD_H
