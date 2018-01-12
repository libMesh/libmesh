// The libMesh Finite Element Library.
// Copyright (C) 2002-2018 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_CELL_INF_PRISM12_H
#define LIBMESH_CELL_INF_PRISM12_H

#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

// Local includes
#include "libmesh/cell_inf_prism.h"

namespace libMesh
{

/**
 * The \p InfPrism12 is an infinite element in 3D composed of 12 nodes.
 * It is numbered like this:
 * \verbatim
 *   INFPRISM12:
 *            5
 *            o
 *            :
 *            :
 *            :
 *     11 o   :   o 10
 *        :  2:   :
 *        :   o   :        closer to infinity
 *        :  . .  :
 *   3o   : . o9. :   o4
 *    |   :.  |  .:   |
 *    |   o   |   o   |
 *    |  . 8  |  7 .  |
 *    | .     |     . |
 *    |.      |      .|     base face
 *    o-------o-------o
 *    0       6       1
 * \endverbatim
 *
 * \author Daniel Dreyer
 * \date 2002
 * \brief A 3D infinite prismatic element with 12 nodes.
 */
class InfPrism12 libmesh_final : public InfPrism
{
public:

  /**
   * Constructor.  By default this element has no parent.
   */
  explicit
  InfPrism12 (Elem * p=libmesh_nullptr) :
    InfPrism(InfPrism12::n_nodes(), p, _nodelinks_data)
  {}

  /**
   * \returns 12.  The \p InfPrism12 has 12 nodes.
   */
  virtual unsigned int n_nodes() const libmesh_override { return 12; }

  /**
   * \returns \p INFPRISM12.
   */
  virtual ElemType type () const libmesh_override { return INFPRISM12; }

  /**
   * \returns 4.
   */
  virtual unsigned int n_sub_elem() const libmesh_override { return 4; }

  /**
   * \returns \p true if the specified (local) node number is a vertex.
   */
  virtual bool is_vertex(const unsigned int i) const libmesh_override;

  /**
   * \returns \p true if the specified (local) node number is an edge.
   */
  virtual bool is_edge(const unsigned int i) const libmesh_override;

  /**
   * \returns \p true if the specified (local) node number is a face.
   */
  virtual bool is_face(const unsigned int i) const libmesh_override;

  /**
   * \returns \p true if the specified (local) node number is on the
   * specified side.
   */
  virtual bool is_node_on_side(const unsigned int n,
                               const unsigned int s) const libmesh_override;

  /**
   * \returns \p true if the specified (local) node number is on the
   * specified edge.
   */
  virtual bool is_node_on_edge(const unsigned int n,
                               const unsigned int e) const libmesh_override;

  /**
   * \returns SECOND.
   */
  virtual Order default_order() const libmesh_override { return SECOND; }

  /**
   * \returns \p InfPrism12::side_nodes_map[side][side_node] after doing some range checking.
   */
  virtual unsigned int which_node_am_i(unsigned int side,
                                       unsigned int side_node) const libmesh_override;

  /**
   * \returns A \p TRI6 built coincident with face 0, or an \p
   * INFQUAD6 built coincident with faces 1 to 3.
   *
   * \note that the \p std::unique_ptr<Elem> takes care of freeing memory.
   */
  virtual std::unique_ptr<Elem> build_side_ptr (const unsigned int i,
                                                bool proxy) libmesh_override;

  /**
   * \returns An \p EDGE3 built coincident with edges 0 to 2, or an \p INFEDGE2
   * built coincident with edges 3 to 5.
   *
   * \note The \p std::unique_ptr<Elem> takes care of freeing memory.
   */
  virtual std::unique_ptr<Elem> build_edge_ptr (const unsigned int i) libmesh_override;

  virtual void connectivity(const unsigned int sc,
                            const IOPackage iop,
                            std::vector<dof_id_type> & conn) const libmesh_override;

  /**
   * \returns 2 for all \p n.
   */
  virtual unsigned int n_second_order_adjacent_vertices (const unsigned int) const libmesh_override
  { return 2; }

  /**
   * \returns The element-local number of the \f$ v^{th} \f$ vertex
   * that defines the \f$ n^{th} \f$ second-order node.
   *
   * \note \p n is counted as depicted above, \f$ 6 \le n < 12 \f$.
   */
  virtual unsigned short int second_order_adjacent_vertex (const unsigned int n,
                                                           const unsigned int v) const libmesh_override;

  /**
   * \returns The child number \p c and element-local index \p v of the
   * \f$ n^{th} \f$ second-order node on the parent element.  See
   * elem.h for further details.
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



protected:

  /**
   * Data for links to nodes.
   */
  Node * _nodelinks_data[12];


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
  static const float _embedding_matrix[4][12][12];

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
  static const unsigned short int _second_order_vertex_child_number[12];

  /**
   * Vector that names the child vertex index for each second order node.
   */
  static const unsigned short int _second_order_vertex_child_index[12];
};

} // namespace libMesh

#endif  // ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

#endif // LIBMESH_CELL_INF_PRISM12_H
