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



#ifndef LIBMESH_FACE_QUAD8_H
#define LIBMESH_FACE_QUAD8_H


// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/face_quad.h"

// C++ includes
#include <cstddef>

namespace libMesh
{



// Forward declarations
class Mesh;



/**
 * The \p QUAD8 is an element in 2D composed of 8 nodes.
 * It is numbered like this:
 * \verbatim
 *        3     6     2
 * QUAD8: o-----o-----o
 *        |           |
 *        |           |
 *      7 o           o 5
 *        |           |
 *        |           |
 *        o-----o-----o
 *        0     4     1
 * \endverbatim
 */
class Quad8 : public Quad
{
public:

  /**
   * Constructor. By default this element has no parent.
   */
  explicit
  Quad8 (Elem* p=NULL) :
    Quad(Quad8::n_nodes(), p, _nodelinks_data) {}

  /**
   * @returns \p QUAD8
   */
  ElemType type () const { return QUAD8; }

  /**
   * @returns 8
   */
  unsigned int n_nodes() const { return 8; }

  /**
   * @returns 5
   */
  unsigned int n_sub_elem() const { return 5; }

  /**
   * @returns true iff the specified (local) node number is a vertex.
   */
  virtual bool is_vertex(const unsigned int i) const;

  /**
   * @returns true iff the specified (local) node number is an edge.
   */
  virtual bool is_edge(const unsigned int i) const;

  /**
   * @returns true iff the specified (local) node number is a face.
   */
  virtual bool is_face(const unsigned int i) const;

  /*
   * @returns true iff the specified (local) node number is on the
   * specified side
   */
  virtual bool is_node_on_side(const unsigned int n,
                               const unsigned int s) const;

  /*
   * @returns true iff the specified (local) node number is on the
   * specified edge (== is_node_on_side in 2D)
   */
  virtual bool is_node_on_edge(const unsigned int n,
                               const unsigned int e) const
  { return this->is_node_on_side(n,e); }

  /*
   * @returns true iff the element map is definitely affine within
   * numerical tolerances
   */
  virtual bool has_affine_map () const;

  /**
   * @returns SECOND
   */
  Order default_order() const { return SECOND; }

  /**
   * @returns an id associated with the \p s side of this element.
   * The id is not necessariy unique, but should be close.  This is
   * particularly useful in the \p MeshBase::find_neighbors() routine.
   *
   * We reimplemenet this method here for the \p Quad8 since we can
   * use the center node of each edge to provide a perfect (unique)
   * key.
   */
  dof_id_type key (const unsigned int s) const;

  UniquePtr<Elem> build_side (const unsigned int i,
                              bool proxy) const;

  virtual void connectivity(const unsigned int sf,
                            const IOPackage iop,
                            std::vector<dof_id_type>& conn) const;

  /**
   * @returns 2 for all \p n
   */
  unsigned int n_second_order_adjacent_vertices (const unsigned int) const
  { return 2; }

  /**
   * @returns the element-local number of the  \f$ v^{th} \f$ vertex
   * that defines the \f$ n^{th} \f$ second-order node.
   * Note that \p n is counted as depicted above, \f$ 4 \le n < 8 \f$.
   */
  unsigned short int second_order_adjacent_vertex (const unsigned int n,
                                                   const unsigned int v) const;

  /**
   * @returns the child number \p c and element-local index \p v of the
   * \f$ n^{th} \f$ second-order node on the parent element.  Note that
   * the return values are always less \p this->n_children() and
   * \p this->child(c)->n_vertices(), while \p n has to be greater or equal
   * to \p * this->n_vertices().  For linear elements this returns 0,0.
   * On refined second order elements, the return value will satisfy
   * \p this->get_node(n)==this->child(c)->get_node(v)
   */
  virtual std::pair<unsigned short int, unsigned short int>
  second_order_child_vertex (const unsigned int n) const;

  /**
   * This maps the \f$ j^{th} \f$ node of the \f$ i^{th} \f$ side to
   * element node numbers.
   */
  static const unsigned int side_nodes_map[4][3];


protected:

  /**
   * Data for links to nodes
   */
  Node* _nodelinks_data[8];



#ifdef LIBMESH_ENABLE_AMR

  /**
   * Matrix used to create the elements children.
   */
  float embedding_matrix (const unsigned int i,
                          const unsigned int j,
                          const unsigned int k) const
  { return _embedding_matrix[i][j][k]; }

  /**
   * Matrix that computes new nodal locations/solution values
   * from current nodes/solution.
   */
  static const float _embedding_matrix[4][8][8];

#endif

};


} // namespace libMesh

#endif // LIBMESH_FACE_QUAD8_H
