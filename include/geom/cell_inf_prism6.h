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



#ifndef LIBMESH_CELL_INF_PRISM6_H
#define LIBMESH_CELL_INF_PRISM6_H

#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

// Local includes
#include "libmesh/cell_inf_prism.h"

namespace libMesh
{

/**
 * The \p InfPrism6 is an infinite element in 3D composed of 6 nodes.
 * It is numbered like this:
 * \verbatim
 * INFPRISM6:
 *         5
 *         o
 *         :
 *         :         closer to infinity
 *         :
 *   3 o   :   o 4
 *     |   :   |
 *     | 2 o   |
 *     |  . .  |
 *     | .   . |
 *     |.     .|
 *     o-------o     base face
 *     0       1
 * \endverbatim
 *
 * \author Daniel Dreyer
 * \date 2002
 * \brief A 3D infinite prismatic element with 6 nodes.
 */
class InfPrism6 final : public InfPrism
{
public:

  /**
   * Constructor.  By default this element has no parent.
   */
  explicit
  InfPrism6 (Elem * p=nullptr) :
    InfPrism(InfPrism6::n_nodes(), p, _nodelinks_data)
  {}

  InfPrism6 (InfPrism6 &&) = delete;
  InfPrism6 (const InfPrism6 &) = delete;
  InfPrism6 & operator= (const InfPrism6 &) = delete;
  InfPrism6 & operator= (InfPrism6 &&) = delete;
  virtual ~InfPrism6() = default;

  /**
   * \returns 6.  The \p InfPrism6 has 6 nodes.
   */
  virtual unsigned int n_nodes() const override { return num_nodes; }

  /**
   * \returns \p INFPRISM6.
   */
  virtual ElemType type() const override { return INFPRISM6; }

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

  virtual std::vector<unsigned int> sides_on_edge(const unsigned int e) const override;

  /**
   * \returns \p true if the specified (local) node number is on the
   * specified edge.
   */
  virtual bool is_node_on_edge(const unsigned int n,
                               const unsigned int e) const override;

  /**
   * \returns FIRST.
   */
  virtual Order default_order() const override;

  /**
   * \returns A \p TRI3 built coincident with face 0, or an \p INFQUAD4
   * built coincident with faces 1 to 3.
   *
   * \note The \p std::unique_ptr<Elem> takes care of freeing memory.
   */
  virtual std::unique_ptr<Elem> build_side_ptr (const unsigned int i,
                                                bool proxy=true) override;

  /**
   * Rebuilds a \p TRI3 built coincident with face 0, or an \p INFQUAD4
   * built coincident with faces 1 to 3.
   */
  virtual void build_side_ptr (std::unique_ptr<Elem> & elem,
                               const unsigned int i) override;

  /**
   * \returns An \p EDGE2 built coincident with edges 0 to 2, an \p INFEDGE2
   * built coincident with edges 3 to 5.
   *
   * \note that the \p std::unique_ptr<Elem> takes care of freeing memory.
   */
  virtual std::unique_ptr<Elem> build_edge_ptr (const unsigned int i) override;

  virtual void connectivity(const unsigned int sc,
                            const IOPackage iop,
                            std::vector<dof_id_type> & conn) const override;

  /**
   * Geometric constants for InfPrism6.
   */
  static const int num_nodes = 6;
  static const int num_sides = 4;
  static const int num_edges = 6;
  static const int num_children = 4;
  static const int nodes_per_side = 4;
  static const int nodes_per_edge = 2;

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

#endif // LIBMESH_CELL_INF_PRISM6_H
