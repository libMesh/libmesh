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



#ifndef LIBMESH_FACE_INF_QUAD6_H
#define LIBMESH_FACE_INF_QUAD6_H


#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

// Local includes
#include "libmesh/face_inf_quad.h"

// C++ includes
#include <cstddef>

namespace libMesh
{



/**
 * The \p INFQUAD6 is an infinite element in 2D composed of 6 nodes.
 * It is numbered like this:
 * \verbatim
 *           2     5     3
 * INFQUAD6: o     o     o   closer to infinity
 *           |           |
 *           |           |
 *           |           |
 *           |           |
 *           |           |
 *           o-----o-----o   base side
 *           0     4     1
 * \endverbatim
 */
class InfQuad6 : public InfQuad
{
public:

  /**
   * Constructor.  By default this element has no parent.
   */
  explicit
  InfQuad6 (Elem * p=libmesh_nullptr):
    InfQuad(InfQuad6::n_nodes(), p, _nodelinks_data) {}

  /**
   * @returns 6
   */
  virtual unsigned int n_nodes() const libmesh_override { return 6; }

  /**
   * @returns \p INFQUAD6
   */
  virtual ElemType type () const libmesh_override { return INFQUAD6; }

  /**
   * @returns 2
   */
  virtual unsigned int n_sub_elem() const libmesh_override { return 2; }

  /**
   * @returns true iff the specified (local) node number is a vertex.
   */
  virtual bool is_vertex(const unsigned int i) const libmesh_override;

  /**
   * @returns true iff the specified (local) node number is an edge.
   */
  virtual bool is_edge(const unsigned int i) const libmesh_override;

  /**
   * @returns true iff the specified (local) node number is a face.
   */
  virtual bool is_face(const unsigned int i) const libmesh_override;

  /*
   * @returns true iff the specified (local) node number is on the
   * specified side
   */
  virtual bool is_node_on_side(const unsigned int n,
                               const unsigned int s) const libmesh_override;

  /*
   * @returns true iff the specified (local) node number is on the
   * specified edge (== is_node_on_side in 2D)
   */
  virtual bool is_node_on_edge(const unsigned int n,
                               const unsigned int e) const libmesh_override
  { return this->is_node_on_side(n,e); }

  /**
   * @returns \p SECOND
   */
  virtual Order default_order() const libmesh_override { return SECOND; }

  /**
   * Don't hide Elem::key() defined in the base class.
   */
  using Elem::key;

  /**
   * @returns an id associated with the \p s side of this element.
   * The id is not necessarily unique, but should be close.  This is
   * particularly useful in the \p MeshBase::find_neighbors() routine.
   *
   * We override this function to return a key for the Edge3 side
   * which is consistent which Edge3::key().
   */
  virtual dof_id_type key (const unsigned int s) const libmesh_override;

  /**
   * Creates and returns an \p Edge3 for the base (0) side, and an \p InfEdge2 for
   * the sides 1, 2.
   */
  virtual UniquePtr<Elem> build_side (const unsigned int i,
                                      bool proxy) const libmesh_override;

  virtual void connectivity(const unsigned int sf,
                            const IOPackage iop,
                            std::vector<dof_id_type> & conn) const libmesh_override;

  /**
   * @returns 2 for all \p n
   */
  virtual unsigned int n_second_order_adjacent_vertices (const unsigned int) const libmesh_override
  { return 2; }

  /**
   * @returns the element-local number of the  \f$ v^{th} \f$ vertex
   * that defines the \f$ n^{th} \f$ second-order node.
   * Note that \p n is counted as depicted above, \f$ 4 \le n < 6 \f$.
   */
  virtual unsigned short int second_order_adjacent_vertex (const unsigned int n,
                                                           const unsigned int v) const libmesh_override;

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
  second_order_child_vertex (const unsigned int n) const libmesh_override;

  /**
   * This maps the \f$ j^{th} \f$ node of the \f$ i^{th} \f$ side to
   * element node numbers.
   */
  static const unsigned int side_nodes_map[3][3];



protected:

  /**
   * Data for links to nodes
   */
  Node * _nodelinks_data[6];



#ifdef LIBMESH_ENABLE_AMR

  /**
   * Matrix used to create the elements children.
   */
  virtual float embedding_matrix (const unsigned int i,
                                  const unsigned int j,
                                  const unsigned int k) const libmesh_override
  { return _embedding_matrix[i][j][k]; }

  /**
   * Matrix that computes new nodal locations/solution values
   * from current nodes/solution.
   */
  static const float _embedding_matrix[2][6][6];

  LIBMESH_ENABLE_TOPOLOGY_CACHES;

#endif // LIBMESH_ENABLE_AMR

private:

  /**
   * Matrix that tells which vertices define the location
   * of mid-side (or second-order) nodes
   */
  static const unsigned short int _second_order_adjacent_vertices[2][2];
};


} // namespace libMesh

#endif // ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

#endif // LIBMESH_FACE_INF_QUAD6_H
