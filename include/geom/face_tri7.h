// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_FACE_TRI7_H
#define LIBMESH_FACE_TRI7_H

// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/face_tri.h"

namespace libMesh
{

/**
 * The \p Tri7 is an element in 2D composed of 7 nodes.  It is
 * numbered like this:
 *
 * \verbatim
 *   TRI7:
 *    2
 *    o
 *    |\
 *    | \
 *    |  \            eta
 *    |   \            ^
 *    |    \           |
 *  5 o     o 4        |
 *    |   o  \         o---> xi
 *    |   6   \
 *    |        \
 *    o----o----o
 *    0    3    1
 * \endverbatim
 *
 * (xi, eta): { 0 <= xi  <= 1
 *            { 0 <= eta <= 1
 *            { xi + eta <= 1
 * are the reference element coordinates associated with the given
 * numbering.
 *
 * \author Roy H Stogner
 * \date 2021
 * \brief A 2D triangular element with 7 nodes.
 */
class Tri7 : public Tri
{
public:

  /**
   * Constructor.  By default this element has no parent.
   */
  explicit
  Tri7 (Elem * p=nullptr) :
    Tri(num_nodes, p, _nodelinks_data) {}

  Tri7 (Tri7 &&) = delete;
  Tri7 (const Tri7 &) = delete;
  Tri7 & operator= (const Tri7 &) = delete;
  Tri7 & operator= (Tri7 &&) = delete;
  virtual ~Tri7() = default;

  /**
   * \returns \p TRI7.
   */
  virtual ElemType type () const override { return TRI7; }

  /**
   * \returns 7.
   */
  virtual unsigned int n_nodes() const override { return num_nodes; }

  /**
   * \returns 4.
   */
  virtual unsigned int n_sub_elem() const override { return 4; }

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
   * specified edge (== is_node_on_side in 2D).
   */
  virtual bool is_node_on_edge(const unsigned int n,
                               const unsigned int e) const override
  { return this->is_node_on_side(n,e); }

  /**
   * \returns \p true if the element map is definitely affine within
   * numerical tolerances.
   */
  virtual bool has_affine_map () const override;

  /**
   * \returns THIRD.
   */
  virtual Order default_order() const override;

  /**
   * \returns SECOND - the sides are EDGE3
   */
  virtual Order default_side_order() const override;

  /**
   * Don't hide Elem::key() defined in the base class.
   */
  using Elem::key;

  /**
   * \returns An id associated with the \p s side of this element.
   * The id is not necessarily unique, but should be close.
   *
   * We reimplement this method here for the \p Tri7 since we can
   * use the center node of each edge to provide a perfect (unique)
   * key.
   */
  virtual dof_id_type key (const unsigned int s) const override;

  /**
   * \returns \p Tri7::side_nodes_map[side][side_node] after doing some range checking.
   */
  virtual unsigned int local_side_node(unsigned int side,
                                       unsigned int side_node) const override;

  virtual std::unique_ptr<Elem> build_side_ptr (const unsigned int i) override;

  /**
   * Rebuilds an EDGE2 coincident with face i.
   */
  virtual void build_side_ptr (std::unique_ptr<Elem> & elem,
                               const unsigned int i) override;

  // Avoid hiding deprecated version with different signature
  using Elem::build_side_ptr;

  virtual void connectivity(const unsigned int sf,
                            const IOPackage iop,
                            std::vector<dof_id_type> & conn) const override;

  /**
   * \returns 2 for edge nodes and 3 for the interior node.
   */
  virtual unsigned int n_second_order_adjacent_vertices (const unsigned int) const override;

  /**
   * \returns The element-local number of the  \f$ v^{th} \f$ vertex
   * that defines the \f$ n^{th} \f$ second-order node.
   *
   * \note \p n is counted as depicted above, \f$ 3 \le n < 7 \f$.
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
   * Geometric constants for Tri7.
   */
  static const int num_nodes = 7;
  static const int nodes_per_side = 3;

  /**
   * This maps the \f$ j^{th} \f$ node of the \f$ i^{th} \f$ side to
   * element node numbers.
   */
  static const unsigned int side_nodes_map[num_sides][nodes_per_side];

  /**
   * \returns A bounding box (not necessarily the minimal bounding box)
   * containing the geometric element.
   */
  virtual BoundingBox loose_bounding_box () const override;

  virtual void permute(unsigned int perm_num) override final;

  virtual void flip(BoundaryInfo *) override final;

#ifdef LIBMESH_ENABLE_AMR
  virtual
  const std::vector<std::pair<unsigned char, unsigned char>> &
  parent_bracketing_nodes(unsigned int c,
                          unsigned int n) const override
  { return _parent_bracketing_nodes[c][n]; }
#endif

  unsigned int center_node_on_side(const unsigned short side) const override final;

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
   * of mid-side (or second-order) nodes
   */
  static const unsigned short int _second_order_adjacent_vertices[num_sides][2];

  /**
   * Vector that names a child sharing each second order node.
   */
  static const unsigned short int _second_order_vertex_child_number[num_nodes];

  /**
   * Vector that names the child vertex index for each second order node.
   */
  static const unsigned short int _second_order_vertex_child_index[num_nodes];
};


} // namespace libMesh

#endif // LIBMESH_FACE_TRI6_H
