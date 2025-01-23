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



#ifndef LIBMESH_CELL_PYRAMID18_H
#define LIBMESH_CELL_PYRAMID18_H

// Local includes
#include "libmesh/cell_pyramid.h"

namespace libMesh
{

/**
 * The \p Pyramid18 is an element in 3D composed of 18 nodes, designed
 * to interface with a QUAD9 element on the base and a TRI7 element on
 * each of the triangular faces.  Cubit will generate hybrid meshes
 * with linear pyramids, but as of version 18 will not export
 * quadratic pyramids.  Paraview may support 13-node pyramids, but
 * does not render 18-node pyramids correctly.  So even if this
 * element works in libmesh, we are currently limited in what we can do
 * with it outside the library...
 *
 * The node numbering for the pyramid18 is given below:
 *
 * \verbatim
 *   PYRAMID18:
 *                       o 4
 *                     //|\
 *                    // | \
 *                   //  |  \
 *                  //   |   \
 *              12 o/    |    o 11
 *                //     |     \
 *               /o 9    o 10   \
 *              //       |       \          zeta
 *             //        |        \          ^   eta (into page)
 *          3 o/.......o.|........o 2        | /
 *           ./       7  |       /           |/
 *          ./           |      /            o---> xi
 *         ./            |     /
 *        ./             |    /
 *     8 o/       o      |   o 6
 *      ./        13     |  /
 *     ./                | /
 *    ./                 |/
 *    o--------o---------o
 *    0        5         1
 * \endverbatim
 *
 * And it also includes four triangle face nodes:
 * Node 14, centroid on side 0, arithmetic mean of 0/1/4 or 5/9/10
 * Node 15, centroid on side 1, arithmetic mean of 1/2/4 or 6/10/11
 * Node 16, centroid on side 2, arithmetic mean of 2/3/4 or 7/11/12
 * Node 17, centroid on side 3, arithmetic mean of 0/3/4 or 8/9/12
 *
 * (xi, eta, zeta): { zeta-1 <= xi   <= 1-zeta
 *                  { zeta-1 <= eta  <= 1-zeta
 *                  {      0 <= zeta <= 1
 * are the reference element coordinates associated with the given
 * numbering.
 *
 * \author Roy H. Stogner
 * \date 2022
 * \brief A 3D pyramid element with 18 nodes.
 */
class Pyramid18 final : public Pyramid
{
public:

  /**
   * Constructor.  By default this element has no parent.
   */
  explicit
  Pyramid18 (Elem * p=nullptr) :
    Pyramid(num_nodes, p, _nodelinks_data)
  {}

  Pyramid18 (Pyramid18 &&) = delete;
  Pyramid18 (const Pyramid18 &) = delete;
  Pyramid18 & operator= (const Pyramid18 &) = delete;
  Pyramid18 & operator= (Pyramid18 &&) = delete;
  virtual ~Pyramid18() = default;

  /**
   * \returns 18.
   */
  virtual unsigned int n_nodes() const override { return num_nodes; }

  /**
   * \returns \p PYRAMID18.
   */
  virtual ElemType type () const override { return PYRAMID18; }

  /**
   * FIXME: we don't yet have a refinement pattern for pyramids...
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
   * Don't hide Pyramid::key() defined in the base class.
   */
  using Pyramid::key;

  /**
   * \returns An id associated with the \p s side of this element.
   * The id is not necessarily unique, but should be close.
   *
   * We reimplement this method here for the \p Pyramid18 since we can
   * use the center node of the base face to provide a perfect (unique)
   * key.
   */
  virtual dof_id_type key (const unsigned int s) const override;

  /**
   * \returns \p Pyramid18::side_nodes_map[side][side_node] after doing some range checking.
   */
  virtual unsigned int local_side_node(unsigned int side,
                                       unsigned int side_node) const override;

  /**
   * \returns \p Pyramid18::edge_nodes_map[edge][edge_node] after doing some range checking.
   */
  virtual unsigned int local_edge_node(unsigned int edge,
                                       unsigned int edge_node) const override;

  /**
   * Builds a \p QUAD9 or \p TRI7 coincident with face i.
   * The \p std::unique_ptr<Elem> handles the memory aspect.
   */
  virtual std::unique_ptr<Elem> build_side_ptr (const unsigned int i,
                                                bool proxy=false) override;

  /**
   * Rebuilds a \p QUAD9 or \p TRI7 built coincident with face i.
   */
  virtual void build_side_ptr (std::unique_ptr<Elem> & elem,
                               const unsigned int i) override;

  /**
   * Builds a \p EDGE3 coincident with edge i.
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
   * \returns 2 for all edge nodes, 4 for quad face nodes, 3 for tri
   * face nodes.
   */
  virtual unsigned int n_second_order_adjacent_vertices (const unsigned int n) const override;

  /**
   * \returns The element-local number of the \f$ v^{th} \f$ vertex
   * that defines the \f$ n^{th} \f$ second-order node.
   */
  virtual unsigned short int second_order_adjacent_vertex (const unsigned int n,
                                                           const unsigned int v) const override;

  /**
   * Geometric constants for Pyramid18.
   */
  static const int num_nodes = 18;
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

  virtual void permute(unsigned int perm_num) override final;

  virtual void flip(BoundaryInfo *) override final;

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
  virtual Real embedding_matrix (const unsigned int,
                                 const unsigned int,
                                 const unsigned int) const override
  { libmesh_not_implemented(); return 0.; }

  LIBMESH_ENABLE_TOPOLOGY_CACHES;

#endif // LIBMESH_ENABLE_AMR

};

} // namespace libMesh


#endif // LIBMESH_CELL_PYRAMID18_H
