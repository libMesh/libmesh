// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

namespace libMesh
{

/**
 * The \p INFQUAD6 is an infinite element in 2D composed of 6 nodes.
 * It is numbered like this:
 * \verbatim
 *             2     5     3
 *   INFQUAD6: o     o     o   closer to infinity
 *             |           |
 *             |           |
 *             |           |
 *             |           |
 *             |           |
 *             o-----o-----o   base side
 *             0     4     1
 * \endverbatim
 *
 * \author Daniel Dreyer
 * \date 2002
 * \brief A 2D infinite quadrilateral element with 6 nodes.
 */
class InfQuad6 : public InfQuad
{
public:

  /**
   * Constructor.  By default this element has no parent.
   */
  explicit
  InfQuad6 (Elem * p=nullptr):
    InfQuad(InfQuad6::n_nodes(), p, _nodelinks_data) {}

  InfQuad6 (InfQuad6 &&) = delete;
  InfQuad6 (const InfQuad6 &) = delete;
  InfQuad6 & operator= (const InfQuad6 &) = delete;
  InfQuad6 & operator= (InfQuad6 &&) = delete;
  virtual ~InfQuad6() = default;

  /**
   * \returns 6.
   */
  virtual unsigned int n_nodes() const override { return num_nodes; }

  /**
   * \returns \p INFQUAD6.
   */
  virtual ElemType type () const override { return INFQUAD6; }

  /**
   * \returns 2.
   */
  virtual unsigned int n_sub_elem() const override { return 2; }

  /**
   * \returns \p true if the specified (local) node number is a vertex.
   */
  virtual bool is_vertex(const unsigned int i) const override;

  /**
   * \returns \p true if the specified (local) node number is an edge.
   */
  virtual bool is_edge(const unsigned int i) const override;

  /**
   * \returns \p true if the specified (local) node number is a face.
   */
  virtual bool is_face(const unsigned int i) const override;

  /**
   * \returns \p true if the specified (local) node number is on the
   * specified side.
   */
  virtual bool is_node_on_side(const unsigned int n,
                               const unsigned int s) const override;

  virtual std::vector<unsigned int> nodes_on_side(const unsigned int s) const override;

  /**
   * \returns \p true if the specified (local) node number is on the
   * specified edge (== is_node_on_side in 2D).
   */
  virtual bool is_node_on_edge(const unsigned int n,
                               const unsigned int e) const override
  { return this->is_node_on_side(n,e); }

  /**
   * \returns \p SECOND.
   */
  virtual Order default_order() const override;

  /**
   * Don't hide Elem::key() defined in the base class.
   */
  using Elem::key;

  /**
   * \returns An id associated with the \p s side of this element.
   * The id is not necessarily unique, but should be close.  This is
   * particularly useful in the \p MeshBase::find_neighbors() routine.
   *
   * We override this function to return a key for the Edge3 side
   * which is consistent which Edge3::key().
   */
  virtual dof_id_type key (const unsigned int s) const override;

  /**
   * \returns \p InfQuad6::side_nodes_map[side][side_node] after doing some range checking.
   */
  virtual unsigned int which_node_am_i(unsigned int side,
                                       unsigned int side_node) const override;

  /**
   * \returns An \p Edge3 for the base (0) side, and an \p InfEdge2 for
   * the sides 1, 2.
   */
  virtual std::unique_ptr<Elem> build_side_ptr (const unsigned int i,
                                                bool proxy=true) override;

  /**
   * Rebuilds an EDGE3 or INFEDGE2 coincident with face i.
   */
  virtual void build_side_ptr (std::unique_ptr<Elem> & elem,
                               const unsigned int i) override;

  virtual void connectivity(const unsigned int sf,
                            const IOPackage iop,
                            std::vector<dof_id_type> & conn) const override;

  /**
   * \returns 2 for all \p n.
   */
  virtual unsigned int n_second_order_adjacent_vertices (const unsigned int) const override
  { return 2; }

  /**
   * \returns The element-local number of the \f$ v^{th} \f$ vertex
   * that defines the \f$ n^{th} \f$ second-order node.
   *
   * \note \p n is counted as depicted above, \f$ 4 \le n < 6 \f$.
   */
  virtual unsigned short int second_order_adjacent_vertex (const unsigned int n,
                                                           const unsigned int v) const override;

  /**
   * \returns The child number \p c and element-local index \p v of
   * the \f$ n^{th} \f$ second-order node on the parent element.  See
   * elem.h for further details.
   */
  virtual std::pair<unsigned short int, unsigned short int>
  second_order_child_vertex (const unsigned int n) const override;

  /**
   * Geometric constants for InfQuad6.
   */
  static const int num_nodes = 6;
  static const int num_sides = 3;
  static const int num_children = 2;
  static const int nodes_per_side = 3;

  /**
   * This maps the \f$ j^{th} \f$ node of the \f$ i^{th} \f$ side to
   * element node numbers.
   */
  static const unsigned int side_nodes_map[num_sides][nodes_per_side];



protected:

  /**
   * Data for links to nodes.
   */
  Node * _nodelinks_data[num_nodes];



#ifdef LIBMESH_ENABLE_AMR

  /**
   * Matrix used to create the elements children.
   */
  virtual float embedding_matrix (const unsigned int i,
                                  const unsigned int j,
                                  const unsigned int k) const override
  { return _embedding_matrix[i][j][k]; }

  /**
   * Matrix that computes new nodal locations/solution values
   * from current nodes/solution.
   */
  static const float _embedding_matrix[num_children][num_nodes][num_nodes];

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
