// The libMesh Finite Element Library.
// Copyright (C) 2002-2020 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_CELL_INF_HEX18_H
#define LIBMESH_CELL_INF_HEX18_H

#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

// Local includes
#include "libmesh/cell_inf_hex.h"

namespace libMesh
{

/**
 * The \p InfHex18 is an infinite element in 3D composed of 18 nodes.
 * It is numbered like this:
 * \verbatim
 *   INFHEX18:   7              14             6
 *               o              o              o     closer to infinity
 *               :              :              |
 *               :              :              |
 *               :              :              |
 *         15    :        17    :        13    |
 *          o    :         o    :         o    |
 *          :    :              :         |    |
 *          :    :              :         |    |
 *          :    :              :         |    |
 *     4    :    :   12         :    5    |    |
 *     o    :    :    o         :    o    |    |
 *     |    :    :    |         :    |    |    |
 *     |    :    :    |         :    |    |    |
 *     |    :    :    |       10:    |    |    |
 *     |    :   3o....|.........o....|....|....o
 *     |    :   .     |              |    |   / 2
 *     |    :  .      |              |    |  /
 *     |    : .       |              |    | /
 *     |    :.        |              |    |/
 *     |  11o         |  16o         |    o           base face
 *     |   .          |              |   / 9
 *     |  .           |              |  /
 *     | .            |              | /
 *     |.             |              |/
 *     o--------------o--------------o
 *     0              8              1
 * \endverbatim
 *
 * \author Daniel Dreyer
 * \date 2002
 * \brief A 3D infinite hexahedral element with 18 nodes.
 */
class InfHex18 final : public InfHex
{
public:

  /**
   * Constructor.  By default this element has no parent.
   */
  explicit
  InfHex18 (Elem * p=nullptr) :
    InfHex(InfHex18::n_nodes(), p, _nodelinks_data)
  {}

  InfHex18 (InfHex18 &&) = delete;
  InfHex18 (const InfHex18 &) = delete;
  InfHex18 & operator= (const InfHex18 &) = delete;
  InfHex18 & operator= (InfHex18 &&) = delete;
  virtual ~InfHex18() = default;

  /**
   * \returns 18.  The \p InfHex18 has 18 nodes.
   */
  virtual unsigned int n_nodes() const override { return num_nodes; }

  /**
   * \returns \p INFHEX18.
   */
  virtual ElemType type () const override { return INFHEX18; }

  /**
   * \returns 4.
   */
  virtual unsigned int n_sub_elem() const override { return 4; }

  /**
   * \returns SECOND.
   */
  virtual Order default_order() const override;

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

  virtual std::vector<unsigned int> nodes_on_edge(const unsigned int e) const override;

  /**
   * \returns \p true if the specified (local) node number is on the
   * specified edge.
   */
  virtual bool is_node_on_edge(const unsigned int n,
                               const unsigned int e) const override;

  /**
   * \returns A \p QUAD9 built coincident with face 0, or an \p INFQUAD6
   * built coincident with faces 1 to 4.
   *
   * \note The \p std::unique_ptr<Elem> takes care of freeing memory.
   */
  virtual std::unique_ptr<Elem> build_side_ptr (const unsigned int i,
                                                bool proxy=true) override;

  /**
   * Rebuilds a \p QUAD9 built coincident with face 0, or an \p
   * INFQUAD6 built coincident with faces 1 to 4.
   */
  virtual void build_side_ptr (std::unique_ptr<Elem> & elem,
                               const unsigned int i) override;

  /**
   * \returns An \p EDGE3 built coincident with edges 0-3, or an \p INFEDGE2
   * built coincident with edges 4 to 11.
   *
   * \note that the \p std::unique_ptr<Elem> takes care of freeing memory.
   */
  virtual std::unique_ptr<Elem> build_edge_ptr (const unsigned int i) override;

  /**
   * Don't hide Elem::key() defined in the base class.
   */
  using Elem::key;

  /**
   * \returns An id associated with the \p s side of this element.
   * The id is not necessarily unique, but should be close.  This is
   * particularly useful in the \p MeshBase::find_neighbors() routine.
   *
   * We reimplement this method here for the \p InfHex18 since we can
   * use the center node of the bottom face to provide a perfect (unique)
   * key.
   */
  virtual dof_id_type key (const unsigned int s) const override;

  /**
   * \returns \p InfHex18::side_nodes_map[side][side_node] after doing some range checking.
   */
  virtual unsigned int local_side_node(unsigned int side,
                                       unsigned int side_node) const override;

  /**
   * \returns \p InfHex18::edge_nodes_map[edge][edge_node] after doing some range checking.
   */
  virtual unsigned int local_edge_node(unsigned int edge,
                                       unsigned int edge_node) const override;

  virtual void connectivity(const unsigned int sc,
                            const IOPackage iop,
                            std::vector<dof_id_type> & conn) const override;

  unsigned int vtk_element_type (const unsigned int) const
  { return 12; }

  /**
   * \returns 2 for all edge nodes, 4 for face nodes.
   */
  virtual unsigned int n_second_order_adjacent_vertices (const unsigned int) const override;

  /**
   * \returns The element-local number of the \f$ v^{th} \f$ vertex
   * that defines the \f$ n^{th} \f$ second-order node.
   *
   * \note \p n is counted as depicted above, \f$ 8 \le n < 18 \f$.
   */
  virtual unsigned short int second_order_adjacent_vertex (const unsigned int n,
                                                           const unsigned int v) const override;

  /**
   * \returns The child number \p c and element-local index \p v of the
   * \f$ n^{th} \f$ second-order node on the parent element.  See
   * elem.h for further details.
   */
  virtual std::pair<unsigned short int, unsigned short int>
  second_order_child_vertex (const unsigned int n) const override;

  /**
   * Geometric constants for InfHex18.
   */
  static const int num_nodes = 18;
  static const int num_sides = 5;
  static const int num_edges = 8;
  static const int num_children = 4;
  static const int nodes_per_side = 9;
  static const int nodes_per_edge = 3;

  /**
   * This maps the \f$ j^{th} \f$ node of the \f$ i^{th} \f$ side to
   * element node numbers.
   */
  static const unsigned int side_nodes_map[num_sides][nodes_per_side];

  /**
   * This maps the \f$ j^{th} \f$ node of the \f$ i^{th} \f$ edge to
   * element node numbers.
   */
  static const unsigned int edge_nodes_map[num_edges][nodes_per_edge];

  /**
   * This maps each edge to the sides that contain said edge.
   */
  static const unsigned int edge_sides_map[num_edges][2];

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

};

} // namespace libMesh

#endif  // ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS


#endif // LIBMESH_CELL_INF_HEX18_H
