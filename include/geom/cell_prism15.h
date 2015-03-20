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



#ifndef LIBMESH_CELL_PRISM15_H
#define LIBMESH_CELL_PRISM15_H

// Local includes
#include "libmesh/cell_prism.h"

// C++ includes
#include <cstddef>

namespace libMesh
{




/**
 * The \p Prism15 is an element in 3D composed of 15 nodes.
 * It is numbered like this:
 * \verbatim
 * PRISM15:
 *              5
 *              o
 *             /:\
 *            / : \
 *           /  :  \
 *          /   :   \
 *      14 o    :    o 13
 *        /     :     \
 *       /      :      \
 *      /       o 11    \
 *   3 /        :        \4
 *    o---------o---------o
 *    |         :12       |
 *    |         :         |
 *    |         :         |
 *    |         o         |
 *    |        .2.        |
 *    |       .   .       |
 *  9 o      .     .      o 10
 *    |     .       .     |
 *    |  8 o         o 7  |
 *    |   .           .   |
 *    |  .             .  |
 *    | .               . |
 *    |.                 .|
 *    o---------o---------o
 *    0         6         1
 *
 * \endverbatim
 */
class Prism15 : public Prism
{
public:

  /**
   * Constructor.  By default this element has no parent.
   */
  explicit
  Prism15  (Elem* p=NULL);

  /**
   * @returns \p PRISM15
   */
  ElemType     type () const   { return PRISM15; }

  /**
   * @returns 15
   */
  unsigned int n_nodes() const { return 15; }

  /**
   * @returns 1
   */
  unsigned int n_sub_elem() const { return 1; }

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
   * specified edge
   */
  virtual bool is_node_on_edge(const unsigned int n,
                               const unsigned int e) const;

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
   * Builds a \p QUAD8 or \p TRI6 built coincident with face i.
   * The \p UniquePtr<Elem> handles the memory aspect.
   */
  UniquePtr<Elem> build_side (const unsigned int i,
                              bool proxy) const;

  /**
   * Builds a \p EDGE3 or \p INFEDGE2 coincident with edge i.
   * The \p UniquePtr<Elem> handles the memory aspect.
   */
  UniquePtr<Elem> build_edge (const unsigned int i) const;

  virtual void connectivity(const unsigned int sc,
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
   * Note that \p n is counted as depicted above, \f$ 6 \le n < 15 \f$.
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
  static const unsigned int side_nodes_map[5][8];

  /**
   * This maps the \f$ j^{th} \f$ node of the \f$ i^{th} \f$ edge to
   * element node numbers.
   */
  static const unsigned int edge_nodes_map[9][3];


protected:

  /**
   * Data for links to nodes
   */
  Node* _nodelinks_data[15];



#ifdef LIBMESH_ENABLE_AMR

  /**
   * Matrix used to create the elements children.
   */
  float embedding_matrix (const unsigned int,
                          const unsigned int,
                          const unsigned int) const
  { libmesh_not_implemented(); return 0.; }

#endif


};



// ------------------------------------------------------------
// Prism15 class member functions
inline
Prism15::Prism15(Elem* p) :
  Prism(Prism15::n_nodes(), p, _nodelinks_data)
{
}


} // namespace libMesh

#endif // LIBMESH_CELL_PRISM15_H
