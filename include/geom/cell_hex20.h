// The libMesh Finite Element Library.
// Copyright (C) 2002-2026 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_CELL_HEX20_H
#define LIBMESH_CELL_HEX20_H

// Local includes
#include "libmesh/cell_hex.h"

namespace libMesh
{

/**
 * The \p Hex20 is an element in 3D composed of 20 nodes.  It is
 * numbered like this:
 *
 * \verbatim
 *   HEX20:      7              18             6
 *               o--------------o--------------o
 *              /:                            /|
 *             / :                           / |
 *            /  :                          /  |
 *         19/   :                       17/   |
 *          o    :                        o    |
 *         /     :                       /     |
 *        /    15o                      /    14o
 *       /       :                     /       |           zeta
 *     4/        :    16             5/        |            ^   eta (into page)
 *     o--------------o--------------o         |            | /
 *     |         :                   |         |            |/
 *     |         :                   |         |            o---> xi
 *     |         :              10   |         |
 *     |        3o..............o....|.........o
 *     |        .                    |        / 2
 *     |       .                   13|       /
 *  12 o      .                      o      /
 *     |     .                       |     /
 *     |  11o                        |    o
 *     |   .                         |   / 9
 *     |  .                          |  /
 *     | .                           | /
 *     |.                            |/
 *     o--------------o--------------o
 *     0              8              1
 * \endverbatim
 *
 * (xi, eta, zeta) in [-1,1]^3 are the reference element coordinates
 * associated with the given numbering.
 *
 * \author Benjamin S. Kirk
 * \date 2002
 * \brief A 3D hexahedral element with 20 nodes.
 */
class Hex20 final : public Hex
{
public:

  /**
   * Constructor.  By default this element has no parent.
   */
  explicit
  Hex20 (Elem * p=nullptr) :
    Hex(num_nodes, p, _nodelinks_data)
  {}

  Hex20 (Hex20 &&) = delete;
  Hex20 (const Hex20 &) = delete;
  Hex20 & operator= (const Hex20 &) = delete;
  Hex20 & operator= (Hex20 &&) = delete;
  virtual ~Hex20() = default;

  /**
   * \returns \p HEX20.
   */
  virtual ElemType type () const override { return HEX20; }

  /**
   * \returns 20.
   */
  virtual unsigned int n_nodes() const override { return num_nodes; }

  /**
   * \returns 1.
   */
  virtual unsigned int n_sub_elem() const override { return 1; }

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
   * \returns \p true if the element map is definitely affine within
   * numerical tolerances.
   */
  virtual bool has_affine_map () const override;

  /**
   * \returns SECOND.
   */
  virtual Order default_order() const override;

  /**
   * \returns \p Hex20::side_nodes_map[side][side_node] after doing some range checking.
   */
  virtual unsigned int local_side_node(unsigned int side,
                                       unsigned int side_node) const override;

  /**
   * \returns \p Hex20::edge_nodes_map[edge][edge_node] after doing some range checking.
   */
  virtual unsigned int local_edge_node(unsigned int edge,
                                       unsigned int edge_node) const override;

  /**
   * Builds a \p QUAD8 built coincident with face i.
   * The \p std::unique_ptr<Elem> handles the memory aspect.
   */
  virtual std::unique_ptr<Elem> build_side_ptr (const unsigned int i) override;

  /**
   * Rebuilds a \p QUAD8 built coincident with face i.
   */
  virtual void build_side_ptr (std::unique_ptr<Elem> & elem,
                               const unsigned int i) override;

  // Avoid hiding deprecated version with different signature
  using Elem::build_side_ptr;

  /**
   * Builds a \p EDGE3 built coincident with edge i.
   * The \p std::unique_ptr<Elem> handles the memory aspect.
   */
  virtual std::unique_ptr<Elem> build_edge_ptr (const unsigned int i) override;

  /**
   * Rebuilds a \p EDGE3 built coincident with edge i.
   */
  virtual void build_edge_ptr (std::unique_ptr<Elem> & edge, const unsigned int i) override;

  virtual void connectivity(const unsigned int sc,
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
   * \note \p n uses the numbering shown above, \f$ 8 \le n < 20 \f$.
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
   * Geometric constants for Hex20.
   */
  static const int num_nodes = 20;
  static const int nodes_per_side = 8;
  static const int nodes_per_edge = 3;

  /**
   * This maps the \f$ j^{th} \f$ node of the \f$ i^{th} \f$ side to
   * element node numbers.
   */
  static constexpr unsigned int side_nodes_map[num_sides][nodes_per_side] =
    {
      {0, 3, 2, 1, 11, 10,  9,  8}, // Side 0
      {0, 1, 5, 4,  8, 13, 16, 12}, // Side 1
      {1, 2, 6, 5,  9, 14, 17, 13}, // Side 2
      {2, 3, 7, 6, 10, 15, 18, 14}, // Side 3
      {3, 0, 4, 7, 11, 12, 19, 15}, // Side 4
      {4, 5, 6, 7, 16, 17, 18, 19}  // Side 5
    };

  /**
   * This maps the \f$ j^{th} \f$ node of the \f$ i^{th} \f$ edge to
   * element node numbers.
   */
  static constexpr unsigned int edge_nodes_map[num_edges][nodes_per_edge] =
    {
      {0, 1, 8},  // Edge 0
      {1, 2, 9},  // Edge 1
      {2, 3, 10}, // Edge 2
      {0, 3, 11}, // Edge 3
      {0, 4, 12}, // Edge 4
      {1, 5, 13}, // Edge 5
      {2, 6, 14}, // Edge 6
      {3, 7, 15}, // Edge 7
      {4, 5, 16}, // Edge 8
      {5, 6, 17}, // Edge 9
      {6, 7, 18}, // Edge 10
      {4, 7, 19}  // Edge 11
    };

  /**
   * A specialization for computing the volume of a Hex20.
   */
  virtual Real volume () const override;

  virtual void permute(unsigned int perm_num) override final;

  virtual void flip(BoundaryInfo *) override final;

  ElemType side_type (const unsigned int s) const override final;

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
                                 const unsigned int k) const override
  { return _embedding_matrix[i][j][k]; }

  /**
   * Matrix that computes new nodal locations/solution values
   * from current nodes/solution.
   */
  static const Real _embedding_matrix[num_children][num_nodes][num_nodes];

  LIBMESH_ENABLE_TOPOLOGY_CACHES;

#endif // LIBMESH_ENABLE_AMR

};

} // namespace libMesh

#endif // LIBMESH_CELL_HEX20_H
