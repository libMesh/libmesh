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



#ifndef LIBMESH_FACE_QUAD8_H
#define LIBMESH_FACE_QUAD8_H

// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/face_quad.h"

namespace libMesh
{

/**
 * The \p QUAD8 is an element in 2D composed of 8 nodes.
 * It is numbered like this:
 * \verbatim
 *          3     6     2
 *   QUAD8: o-----o-----o
 *          |           |
 *          |           |
 *        7 o           o 5
 *          |           |
 *          |           |
 *          o-----o-----o
 *          0     4     1
 * \endverbatim
 *
 * \author Benjamin S. Kirk
 * \date 2002
 * \brief A 2D quadrilateral element with 8 nodes.
 */
class Quad8 : public Quad
{
public:

  /**
   * Constructor. By default this element has no parent.
   */
  explicit
  Quad8 (Elem * p=libmesh_nullptr) :
    Quad(Quad8::n_nodes(), p, _nodelinks_data) {}

  /**
   * @returns \p QUAD8.
   */
  virtual ElemType type () const libmesh_override { return QUAD8; }

  /**
   * @returns 8.
   */
  virtual unsigned int n_nodes() const libmesh_override { return 8; }

  /**
   * @returns 5.
   */
  virtual unsigned int n_sub_elem() const libmesh_override { return 5; }

  /**
   * @returns true if the specified (local) node number is a vertex.
   */
  virtual bool is_vertex(const unsigned int i) const libmesh_override;

  /**
   * @returns true if the specified (local) node number is an edge.
   */
  virtual bool is_edge(const unsigned int i) const libmesh_override;

  /**
   * @returns true if the specified (local) node number is a face.
   */
  virtual bool is_face(const unsigned int i) const libmesh_override;

  /**
   * @returns true if the specified (local) node number is on the
   * specified side.
   */
  virtual bool is_node_on_side(const unsigned int n,
                               const unsigned int s) const libmesh_override;

  /**
   * @returns true if the specified (local) node number is on the
   * specified edge (== is_node_on_side in 2D).
   */
  virtual bool is_node_on_edge(const unsigned int n,
                               const unsigned int e) const libmesh_override
  { return this->is_node_on_side(n,e); }

  /**
   * @returns true if the element map is definitely affine within
   * numerical tolerances.
   */
  virtual bool has_affine_map () const libmesh_override;

  /**
   * @returns SECOND.
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
   * We reimplemenet this method here for the \p Quad8 since we can
   * use the center node of each edge to provide a perfect (unique)
   * key.
   */
  virtual dof_id_type key (const unsigned int s) const libmesh_override;

  virtual UniquePtr<Elem> build_side_ptr (const unsigned int i,
                                          bool proxy) libmesh_override;

  virtual void connectivity(const unsigned int sf,
                            const IOPackage iop,
                            std::vector<dof_id_type> & conn) const libmesh_override;

  /**
   * @returns 2 for all \p n.
   */
  virtual unsigned int n_second_order_adjacent_vertices (const unsigned int) const libmesh_override
  { return 2; }

  /**
   * @returns the element-local number of the  \f$ v^{th} \f$ vertex
   * that defines the \f$ n^{th} \f$ second-order node.
   * Note that \p n is counted as depicted above, \f$ 4 \le n < 8 \f$.
   */
  virtual unsigned short int second_order_adjacent_vertex (const unsigned int n,
                                                           const unsigned int v) const libmesh_override;

  /**
   * @returns the child number \p c and element-local index \p v of the
   * \f$ n^{th} \f$ second-order node on the parent element.  See
   * elem.h for further details.
   */
  virtual std::pair<unsigned short int, unsigned short int>
  second_order_child_vertex (const unsigned int n) const libmesh_override;

  /**
   * This maps the \f$ j^{th} \f$ node of the \f$ i^{th} \f$ side to
   * element node numbers.
   */
  static const unsigned int side_nodes_map[4][3];

  /**
   * An optimized method for approximating the area of a
   * QUAD8 using quadrature.
   */
  virtual Real volume () const libmesh_override;

  /**
   * @return a bounding box (not necessarily the minimal bounding box)
   * containing the geometric element.
   */
  virtual BoundingBox loose_bounding_box () const libmesh_override;

protected:

  /**
   * Data for links to nodes.
   */
  Node * _nodelinks_data[8];



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
  static const float _embedding_matrix[4][8][8];

  LIBMESH_ENABLE_TOPOLOGY_CACHES;

#endif // LIBMESH_ENABLE_AMR

};


} // namespace libMesh

#endif // LIBMESH_FACE_QUAD8_H
