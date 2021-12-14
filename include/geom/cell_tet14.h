// The libMesh Finite Element Library.
// Copyright (C) 2002-2021 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_CELL_TET14_H
#define LIBMESH_CELL_TET14_H

// Local includes
#include "libmesh/cell_tet.h"

namespace libMesh
{

/**
 * The \p Tet14 is an element in 3D composed of 14 nodes.  Its edges
 * and vertices are numbered as in Tet10, like this:
 *
 * \verbatim
 *               3
 *   TET14:      o
 *              /|\
 *             / | \
 *         7  /  |  \9
 *           o   |   o              zeta
 *          /    |8   \              ^
 *         /     o     \             |
 *        /    6 |      \            |
 *     0 o.....o.|.......o 2         o---> eta
 *        \      |      /             \
 *         \     |     /               \
 *          \    |    /                 xi (out of page)
 *         4 o   |   o 5
 *            \  |  /
 *             \ | /
 *              \|/
 *               o
 *               1
 * \endverbatim
 *
 * And it also includes four face nodes:
 * Node 10, on side 0, equidistant from 0/1/2 or 4/5/6
 * Node 11, on side 1, equidistant from 0/1/3 or 4/7/8
 * Node 12, on side 2, equidistant from 1/2/3 or 5/8/9
 * Node 13, on side 3, equidistant from 0/2/3 or 6/7/9
 *
 * (xi, eta, zeta): { 0 <= xi   <= 1
 *                  { 0 <= eta  <= 1
 *                  { 0 <= zeta <= 1
 *                  { xi + eta + zeta <= 1
 *
 * \author Roy H. Stogner
 * \date 2021
 * \brief A 3D tetrahedral element with 14 nodes.
 */
class Tet14 final : public Tet
{
public:

  /**
   * Constructor.  By default this element has no parent.
   */
  explicit
  Tet14 (Elem * p=nullptr) :
    Tet(Tet14::n_nodes(), p, _nodelinks_data)
  {}

  Tet14 (Tet14 &&) = delete;
  Tet14 (const Tet14 &) = delete;
  Tet14 & operator= (const Tet14 &) = delete;
  Tet14 & operator= (Tet14 &&) = delete;
  virtual ~Tet14() = default;

  /**
   * \returns \p TET14.
   */
  virtual ElemType type () const override { return TET14; }

  /**
   * \returns 10.
   */
  virtual unsigned int n_nodes() const override { return num_nodes; }

  /**
   * \returns 8.
   */
  virtual unsigned int n_sub_elem() const override { return 8; }

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
   * \returns \p true if the specified child is on the
   * specified side.
   */
  virtual bool is_child_on_side(const unsigned int c,
                                const unsigned int s) const override;

  /**
   * \returns \p true if the element map is definitely affine within
   * numerical tolerances.
   */
  virtual bool has_affine_map () const override;

  /**
   * \returns SECOND.
   */
  virtual Order default_order() const override;

  /**
   * \returns \p Tet14::side_nodes_map[side][side_node] after doing some range checking.
   */
  virtual unsigned int local_side_node(unsigned int side,
                                       unsigned int side_node) const override;

  /**
   * \returns \p Tet14::edge_nodes_map[edge][edge_node] after doing some range checking.
   */
  virtual unsigned int local_edge_node(unsigned int edge,
                                       unsigned int edge_node) const override;

  /**
   * Builds a \p TRI7 built coincident with face i.
   * The \p std::unique_ptr<Elem> handles the memory aspect.
   */
  virtual std::unique_ptr<Elem> build_side_ptr (const unsigned int i,
                                                bool proxy=false) override;

  /**
   * Rebuilds a TRI7 built coincident with face i.
   */
  virtual void build_side_ptr (std::unique_ptr<Elem> & elem,
                               const unsigned int i) override;

  /**
   * Builds a \p EDGE3 built coincident with edge i.
   * The \p std::unique_ptr<Elem> handles the memory aspect.
   */
  virtual std::unique_ptr<Elem> build_edge_ptr (const unsigned int i) override;

  /**
   * Rebuilds a \p EDGE3 coincident with edge i.
   */
  virtual void build_edge_ptr (std::unique_ptr<Elem> & edge, const unsigned int i) override;

  virtual void connectivity(const unsigned int sc,
                            const IOPackage iop,
                            std::vector<dof_id_type> & conn) const override;

  /**
   * \returns 2 for edge nodes and 3 for face nodes
   */
  virtual unsigned int n_second_order_adjacent_vertices (const unsigned int) const override;

  /**
   * \returns The element-local number of the \f$ v^{th} \f$ vertex
   * that defines the \f$ n^{th} \f$ second-order node.
   *
   * \note \p n is counted as depicted above, \f$ 4 \le n < 14 \f$.
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
   * Geometric constants for Tet14.
   */
  static const int num_nodes = 14;
  static const int num_sides = 4;
  static const int num_edges = 6;
  static const int num_children = 8;
  static const int nodes_per_side = 7;
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

  virtual void permute(unsigned int perm_num) override final;

  ElemType side_type (const unsigned int s) const override final;

#ifdef LIBMESH_ENABLE_AMR
  virtual
  const std::vector<std::pair<unsigned char, unsigned char>> &
  parent_bracketing_nodes(unsigned int c,
                          unsigned int n) const override
  { return _parent_bracketing_nodes[c][n]; }
#endif

protected:

  /**
   * Data for links to nodes.
   */
  Node * _nodelinks_data[num_nodes];



#ifdef LIBMESH_ENABLE_AMR

  /**
   * Matrix used to create the elements children.
   */
  virtual Real embedding_matrix (const unsigned int i,
                                 const unsigned int j,
                                 const unsigned int k) const override;

  /**
   * Matrix that computes new nodal locations/solution values
   * from current nodes/solution.
   */
  static const Real _embedding_matrix[num_children][num_nodes][num_nodes];

  /**
   * Pairs of nodes that bracket child nodes when doing mesh
   * refinement.
   */
  static const std::vector<std::pair<unsigned char, unsigned char>>
    _parent_bracketing_nodes[num_children][num_nodes];

  LIBMESH_ENABLE_TOPOLOGY_CACHES;

#endif // LIBMESH_ENABLE_AMR

private:

  /**
   * Matrix that tells which vertices define the location
   * of mid-side or mid-face nodes, indexed by node_num-4
   */
  static const unsigned short int _second_order_adjacent_vertices[10][3];

  /**
   * Vector that names a child sharing each second order node.
   */
  static const unsigned short int _second_order_vertex_child_number[14];

  /**
   * Vector that names the child vertex index for each second order node.
   */
  static const unsigned short int _second_order_vertex_child_index[14];
};

} // namespace libMesh


#endif // LIBMESH_CELL_TET14_H
