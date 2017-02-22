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



#ifndef LIBMESH_FACE_TRI6_H
#define LIBMESH_FACE_TRI6_H

// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/face_tri.h"

namespace libMesh
{

/**
 * The \p Tri6 is an element in 2D composed of 6 nodes.
 * It is numbered like this:
 * \verbatim
 *   TRI6:  2
 *          o
 *         / \
 *        /   \
 *     5 o     o 4
 *      /       \
 *     /         \
 *    o-----o-----o
 *    0     3     1
 * \endverbatim
 *
 * \author Benjamin S. Kirk
 * \date 2002
 * \brief A 2D triangular element with 6 nodes.
 */
class Tri6 : public Tri
{
public:

  /**
   * Constructor.  By default this element has no parent.
   */
  explicit
  Tri6 (Elem * p=libmesh_nullptr) :
    Tri(Tri6::n_nodes(), p, _nodelinks_data) {}

  /**
   * @returns \p TRI6.
   */
  virtual ElemType type () const libmesh_override { return TRI6; }

  /**
   * @returns 6.
   */
  virtual unsigned int n_nodes() const libmesh_override { return 6; }

  /**
   * @returns 4.
   */
  virtual unsigned int n_sub_elem() const libmesh_override { return 4; }

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
   * Note that \p n is counted as depicted above, \f$ 3 \le n < 6 \f$.
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
  static const unsigned int side_nodes_map[3][3];

  /**
   * An optimized method for approximating the area of a
   * TRI6 using quadrature.
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
  static const float _embedding_matrix[4][6][6];

  LIBMESH_ENABLE_TOPOLOGY_CACHES;

#endif // LIBMESH_ENABLE_AMR

private:

  /**
   * Matrix that tells which vertices define the location
   * of mid-side (or second-order) nodes
   */
  static const unsigned short int _second_order_adjacent_vertices[3][2];

  /**
   * Vector that names a child sharing each second order node.
   */
  static const unsigned short int _second_order_vertex_child_number[6];

  /**
   * Vector that names the child vertex index for each second order node.
   */
  static const unsigned short int _second_order_vertex_child_index[6];
};


} // namespace libMesh

#endif // LIBMESH_FACE_TRI6_H
