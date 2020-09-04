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



#ifndef LIBMESH_CELL_PYRAMID14_H
#define LIBMESH_CELL_PYRAMID14_H

// Local includes
#include "libmesh/cell_pyramid.h"

namespace libMesh
{

/**
 * The \p Pyramid14 is an element in 3D composed of 14 nodes, designed
 * to interface with a QUAD9 element on the base and a TRI6 element on
 * each of the triangular faces.  Cubit will generate hybrid meshes
 * with linear pyramids, but as of version 14 will not export
 * quadratic pyramids.  Paraview may support 13-node pyramids, but
 * does not render 14-node pyramids correctly.  So even if this
 * element works in libmesh, we are currently limited in what we can do
 * with it outside the library...
 *
 * The node numbering for the pyramid14 is given below:
 * \verbatim
 *   PYRAMID14:
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
 * (xi, eta, zeta) are the reference element coordinates associated with
 * the given numbering.
 *
 * \author John W. Peterson
 * \date 2013
 * \brief A 3D pyramid element with 14 nodes.
 */
class Pyramid14 final : public Pyramid
{
public:

  /**
   * Constructor.  By default this element has no parent.
   */
  explicit
  Pyramid14 (Elem * p=nullptr) :
    Pyramid(Pyramid14::n_nodes(), p, _nodelinks_data)
  {}

  Pyramid14 (Pyramid14 &&) = delete;
  Pyramid14 (const Pyramid14 &) = delete;
  Pyramid14 & operator= (const Pyramid14 &) = delete;
  Pyramid14 & operator= (Pyramid14 &&) = delete;
  virtual ~Pyramid14() = default;

  /**
   * \returns 14.
   */
  virtual unsigned int n_nodes() const override { return num_nodes; }

  /**
   * \returns \p PRYAMID14.
   */
  virtual ElemType type () const override { return PYRAMID14; }

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
   * The id is not necessarily unique, but should be close.  This is
   * particularly useful in the \p MeshBase::find_neighbors() routine.
   *
   * We reimplement this method here for the \p Pyramid14 since we can
   * use the center node of the base face to provide a perfect (unique)
   * key.
   */
  virtual dof_id_type key (const unsigned int s) const override;

  /**
   * \returns \p Pyramid14::side_nodes_map[side][side_node] after doing some range checking.
   */
  virtual unsigned int local_side_node(unsigned int side,
                                       unsigned int side_node) const override;

  /**
   * \returns \p Pyramid14::edge_nodes_map[edge][edge_node] after doing some range checking.
   */
  virtual unsigned int local_edge_node(unsigned int edge,
                                       unsigned int edge_node) const override;

  /**
   * Builds a \p QUAD9 or \p TRI6 coincident with face i.
   * The \p std::unique_ptr<Elem> handles the memory aspect.
   */
  virtual std::unique_ptr<Elem> build_side_ptr (const unsigned int i,
                                                bool proxy=true) override;

  /**
   * Rebuilds a \p QUAD9 or \p TRI6 built coincident with face i.
   */
  virtual void build_side_ptr (std::unique_ptr<Elem> & elem,
                               const unsigned int i) override;

  /**
   * Builds a \p EDGE3 coincident with edge i.
   * The \p std::unique_ptr<Elem> handles the memory aspect.
   */
  virtual std::unique_ptr<Elem> build_edge_ptr (const unsigned int i) override;

  virtual void connectivity(const unsigned int sc,
                            const IOPackage iop,
                            std::vector<dof_id_type> & conn) const override;

  /**
   * \returns 2 for all edge nodes and 4 for face nodes.
   */
  virtual unsigned int n_second_order_adjacent_vertices (const unsigned int n) const override;

  /**
   * \returns The element-local number of the \f$ v^{th} \f$ vertex
   * that defines the \f$ n^{th} \f$ second-order node.
   */
  virtual unsigned short int second_order_adjacent_vertex (const unsigned int n,
                                                           const unsigned int v) const override;

  /**
   * Geometric constants for Pyramid14.
   */
  static const int num_nodes = 14;
  static const int num_sides = 5;
  static const int num_edges = 8;
  static const int num_children = 0; // not implemented
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

  /**
   * Specialization for computing the volume of a Pyramid14.
   */
  virtual Real volume () const override;

protected:

  /**
   * Data for links to nodes.
   */
  Node * _nodelinks_data[num_nodes];



#ifdef LIBMESH_ENABLE_AMR

  /**
   * Matrix used to create the elements children.
   */
  virtual float embedding_matrix (const unsigned int,
                                  const unsigned int,
                                  const unsigned int) const override
  { libmesh_not_implemented(); return 0.; }

  LIBMESH_ENABLE_TOPOLOGY_CACHES;

#endif // LIBMESH_ENABLE_AMR

};

} // namespace libMesh


#endif // LIBMESH_CELL_PYRAMID14_H
