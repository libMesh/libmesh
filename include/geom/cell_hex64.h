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



#ifndef LIBMESH_CELL_HEX64_H
#define LIBMESH_CELL_HEX64_H

// Local includes
#include "libmesh/cell_hex.h"
#include "libmesh/cell_hex8.h"
#include "libmesh/fe_reference_element_traits.h"

namespace libMesh
{

/**
 * The \p Hex64 is an element in 3D composed of 64 nodes.  It is
 * numbered like this:
 *
 * \verbatim
 *   HEX64 (opaque view: vertex 3 and the back/left/bottom faces are hidden):
 *
 *              7 --------29--------28--------6
 *             /         /         /         /|
 *            /         /         /         / |
 *           30--------55--------54--------27 |
 *          /         /         /         /|  21
 *         /         /         /         / | /|          zeta
 *        31--------52--------53--------26 |/ |            ^   eta (into page)
 *       /         /         /         /|  42 |            |  /
 *      /         /         /         / | /|  20           | /
 *     4 --------24--------25--------5  |/ | /|            |/
 *     |         |         |         |  43 |/ |            o---> xi
 *     |         |         |         | /|  41 |
 *     |         |         |         |/ | /|  2
 *     17--------39--------38--------19 |/ | /
 *     |         |         |         |  40 |/
 *     |         |         |         | /|  11
 *     |         |         |         |/ | /
 *     16--------36--------37--------18 |/
 *     |         |         |         |  10
 *     |         |         |         | /
 *     |         |         |         |/
 *     0 --------8 --------9 --------1
 * \endverbatim
 *
 *
 * \verbatim
 *   HEX64 nodes by zeta slice (xi to the right, eta up the page).
 *   Slices 1 and 2 each contain 4 of the 8 interior nodes (56-63).
 *
 *   zeta = 0 (bottom face)            zeta = 1
 *   3 ----13----12----2               22----45----44----20
 *   |     |     |     |               |     |     |     |
 *   14----33----34----11              48----59----58----41
 *   |     |     |     |               |     |     |     |
 *   15----32----35----10              49----56----57----40
 *   |     |     |     |               |     |     |     |
 *   0 ----8 ----9 ----1               16----36----37----18
 *
 *   zeta = 2                          zeta = 3 (top face)
 *   23----46----47----21              7 ----29----28----6
 *   |     |     |     |               |     |     |     |
 *   51----63----62----42              30----55----54----27
 *   |     |     |     |               |     |     |     |
 *   50----60----61----43              31----52----53----26
 *   |     |     |     |               |     |     |     |
 *   17----39----38----19              4 ----24----25----5
 * \endverbatim
 *
 * (xi, eta, zeta) in [-1,1]^3 are the reference element coordinates
 * associated with the given numbering.
 *
 * \author Benjamin S. Kirk, Micah Gale
 * \date 2026
 * \brief A 3D hexahedral element with 64 nodes.
 */
class Hex64 final : public Hex
{
public:

  /**
   * Constructor.  By default this element has no parent.
   */
  explicit
  Hex64 (Elem * p=nullptr) :
    Hex(num_nodes, p, _nodelinks_data)
  {}

  Hex64 (Hex64 &&) = delete;
  Hex64 (const Hex64 &) = delete;
  Hex64 & operator= (const Hex64 &) = delete;
  Hex64 & operator= (Hex64 &&) = delete;
  virtual ~Hex64() = default;

  virtual Point master_point (const unsigned int i) const override final
  {
    libmesh_assert_less(i, this->n_nodes());
    return Point(_master_points[i][0],
                 _master_points[i][1],
                 _master_points[i][2]);
  }

  /**
   * \returns \p HEX64.
   */
  virtual ElemType type () const override { return HEX64; }

  /**
   * \returns 27.
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
   * \returns SECOND.
   */
  virtual Order default_order() const override;

  /**
   * Don't hide Elem::key() defined in the base class.
   */
  using Elem::key;

  /**
   * \returns An id associated with the \p s side of this element.
   * The id is not necessarily unique, but should be close.
   *
   * We reimplement this method here for the \p Hex64 since we can
   * use the center node of each face to provide a perfect (unique)
   * key.
   */
  virtual dof_id_type key (const unsigned int s) const override;

  /**
   * \returns \p Hex64::side_nodes_map[side][side_node] after doing some range checking.
   */
  virtual unsigned int local_side_node(unsigned int side,
                                       unsigned int side_node) const override;

  /**
   * \returns \p Hex64::edge_nodes_map[edge][edge_node] after doing some range checking.
   */
  virtual unsigned int local_edge_node(unsigned int edge,
                                       unsigned int edge_node) const override;

  /**
   * Builds a \p QUAD9 built coincident with face i.
   * The \p std::unique_ptr<Elem> handles the memory aspect.
   */
  virtual std::unique_ptr<Elem> build_side_ptr (const unsigned int i) override;

  /**
   * Rebuilds a \p QUAD9 built coincident with face i.
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


  // TODO is this third order?

  /**
   * \returns 2 for all edge nodes, 4 for all face nodes, and
   * 8 for the bubble node.
   */
  virtual unsigned int n_second_order_adjacent_vertices (const unsigned int) const override;

  /**
   * \returns The element-local number of the \f$ v^{th} \f$ vertex
   * that defines the \f$ n^{th} \f$ second-order node.
   *
   * \note \p n is counted as depicted above, \f$ 8 \le n < 27 \f$.
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
   * Geometric constants for Hex64.
   */
  static const int num_nodes = 64;
  static const int nodes_per_side = 16;
  static const int nodes_per_edge = 4;

  /**
   * These map the \f$ j^{th} \f$ node of the \f$ i^{th} \f$ edge or
   * side to element node numbers.  They are derived from the
   * first-order Hex8 tables; see fe_reference_element_traits.h.
   */
  static constexpr NodeMapTable<num_edges, nodes_per_edge>
  _edge_nodes = derived_edge_nodes<Hex8, num_edges, nodes_per_edge - 2>();
  static constexpr const unsigned int (&edge_nodes_map)[num_edges][nodes_per_edge] = _edge_nodes.values;

  static constexpr unsigned int _first_face_node = Hex8::num_nodes + (nodes_per_edge - 2) * num_edges;

  /**
   * Master element node locations
   */
  static const Real _master_points[64][3];

  static constexpr NodeMapTable<num_sides, nodes_per_side>
  _side_nodes = derived_side_nodes<Hex8, num_sides, nodes_per_side>
    (_edge_nodes.values,
      [](unsigned int s, unsigned int q)
      { return q < 4 ? _first_face_node + 4 * s + q : invalid_uint; });

  static constexpr const unsigned int (&side_nodes_map)[num_sides][nodes_per_side] = _side_nodes.values;


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
  virtual Real embedding_matrix (const unsigned int,
                                 const unsigned int,
                                 const unsigned int) const override
  { libmesh_not_implemented(); return 0; }


  LIBMESH_ENABLE_TOPOLOGY_CACHES;

#endif // LIBMESH_ENABLE_AMR


private:


  // TODO
  /**
   * Matrix that tells which vertices define the location
   * of mid-side (or second-order) nodes.  This matrix only
   * covers the nodes that are unique to \p Hex64, while the
   * second-order-nodes that are identical with \p Hex20 are covered
   * through the \p _second_order_adjacent_vertices matrix in
   * \p cell_hex.C.
   *
   * \note This matrix also does not cover the bubble node.  The
   * interpolation is trivial and would only blow up the size of this
   * matrix.
   */
  static const unsigned short int _remaining_second_order_adjacent_vertices[6][4];

};

} // namespace libMesh

#endif // LIBMESH_CELL_HEX64_H
