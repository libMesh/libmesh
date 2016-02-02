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



#ifndef LIBMESH_CELL_TET10_H
#define LIBMESH_CELL_TET10_H

// Local includes
#include "libmesh/cell_tet.h"

// C++ includes
#include <cstddef>

namespace libMesh
{




/**
 * The \p Tet10 is an element in 3D composed of 10 nodes.
 * It is numbered like this:
 * \verbatim
 *             3
 * TET10:      o
 *            /|\
 *           / | \
 *       7  /  |  \9
 *         o   |   o
 *        /    |8   \
 *       /     o     \
 *      /    6 |      \
 *   0 o.....o.|.......o 2
 *      \      |      /
 *       \     |     /
 *        \    |    /
 *       4 o   |   o 5
 *          \  |  /
 *           \ | /
 *            \|/
 *             o
 *             1
 *  \endverbatim
 */
class Tet10 libmesh_final : public Tet
{
public:

  /**
   * Constructor.  By default this element has no parent.
   */
  explicit
  Tet10 (Elem * p=libmesh_nullptr) :
    Tet(Tet10::n_nodes(), p, _nodelinks_data)
  {}

  /**
   * @returns \p TET10
   */
  virtual ElemType type () const libmesh_override { return TET10; }

  /**
   * @returns 10
   */
  virtual unsigned int n_nodes() const libmesh_override { return 10; }

  /**
   * @returns 8
   */
  virtual unsigned int n_sub_elem() const libmesh_override { return 8; }

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
   * specified edge
   */
  virtual bool is_node_on_edge(const unsigned int n,
                               const unsigned int e) const libmesh_override;

  /*
   * @returns true iff the specified child is on the
   * specified side
   */
  virtual bool is_child_on_side(const unsigned int c,
                                const unsigned int s) const libmesh_override;

  /*
   * @returns true iff the element map is definitely affine within
   * numerical tolerances
   */
  virtual bool has_affine_map () const libmesh_override;

  /**
   * @returns SECOND
   */
  virtual Order default_order() const libmesh_override { return SECOND; }

  /**
   * Builds a \p TRI6 built coincident with face i.
   * The \p UniquePtr<Elem> handles the memory aspect.
   */
  virtual UniquePtr<Elem> build_side (const unsigned int i,
                                      bool proxy) const libmesh_override;

  /**
   * Builds a \p EDGE3 built coincident with edge i.
   * The \p UniquePtr<Elem> handles the memory aspect.
   */
  virtual UniquePtr<Elem> build_edge (const unsigned int i) const libmesh_override;

  virtual void connectivity(const unsigned int sc,
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
   * Note that \p n is counted as depicted above, \f$ 4 \le n < 10 \f$.
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
  static const unsigned int side_nodes_map[4][6];

  /**
   * This maps the \f$ j^{th} \f$ node of the \f$ i^{th} \f$ edge to
   * element node numbers.
   */
  static const unsigned int edge_nodes_map[6][3];

  /**
   * A specialization for computing the volume of a Tet10.
   */
  virtual Real volume () const libmesh_override;

protected:

  /**
   * Data for links to nodes
   */
  Node * _nodelinks_data[10];



#ifdef LIBMESH_ENABLE_AMR

  /**
   * Matrix used to create the elements children.
   */
  virtual float embedding_matrix (const unsigned int i,
                                  const unsigned int j,
                                  const unsigned int k) const libmesh_override;

  /**
   * Matrix that computes new nodal locations/solution values
   * from current nodes/solution.
   */
  static const float _embedding_matrix[8][10][10];

  LIBMESH_ENABLE_TOPOLOGY_CACHES;

#endif // LIBMESH_ENABLE_AMR

private:

  /**
   * Matrix that tells which vertices define the location
   * of mid-side (or second-order) nodes
   */
  static const unsigned short int _second_order_adjacent_vertices[6][2];

  /**
   * Vector that names a child sharing each second order node.
   */
  static const unsigned short int _second_order_vertex_child_number[10];

  /**
   * Vector that names the child vertex index for each second order node.
   */
  static const unsigned short int _second_order_vertex_child_index[10];
};

} // namespace libMesh


#endif // LIBMESH_CELL_TET10_H
