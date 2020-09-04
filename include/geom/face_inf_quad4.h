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


#ifndef LIBMESH_FACE_INF_QUAD4_H
#define LIBMESH_FACE_INF_QUAD4_H


#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

// Local includes
#include "libmesh/face_inf_quad.h"

namespace libMesh
{

/**
 * The \p INFQUAD4 is an infinite element in 2D composed of 4 nodes.
 * It is numbered like this:
 * \verbatim
 *              2           3
 *   INFQUAD4: o           o   closer to infinity
 *             |           |
 *             |           |
 *             |           |
 *             |           |
 *             |           |
 *             o-----------o   base side
 *             0           1
 * \endverbatim
 *
 * \author Daniel Dreyer
 * \date 2002
 * \brief A 2D infinite quadrilateral element with 4 nodes.
 */
class InfQuad4 : public InfQuad
{
public:

  /**
   * Constructor.  By default this element has no parent.
   */
  explicit
  InfQuad4 (Elem * p=nullptr) :
    InfQuad(InfQuad4::n_nodes(), p, _nodelinks_data) {}

  InfQuad4 (InfQuad4 &&) = delete;
  InfQuad4 (const InfQuad4 &) = delete;
  InfQuad4 & operator= (const InfQuad4 &) = delete;
  InfQuad4 & operator= (InfQuad4 &&) = delete;
  virtual ~InfQuad4() = default;

  /**
   * \returns 4.
   */
  virtual unsigned int n_nodes() const override { return num_nodes; }

  /**
   * \returns \p INFQUAD4.
   */
  virtual ElemType type () const override { return INFQUAD4; }

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
   * specified edge (== is_node_on_side in 2D).
   */
  virtual bool is_node_on_edge(const unsigned int n,
                               const unsigned int e) const override
  { return this->is_node_on_side(n,e); }

  /**
   * \returns \p FIRST.
   */
  virtual Order default_order() const override;

  /**
   * \returns An \p Edge2 for the base side, or an \p InfEdge2 for
   * the sides 1, 2.
   */
  virtual std::unique_ptr<Elem> build_side_ptr (const unsigned int i,
                                                bool proxy=true) override;

  /**
   * Rebuilds an EDGE2 or INFEDGE2 coincident with face i.
   */
  virtual void build_side_ptr (std::unique_ptr<Elem> & elem,
                               const unsigned int i) override;


  virtual void connectivity(const unsigned int sf,
                            const IOPackage iop,
                            std::vector<dof_id_type> & conn) const override;

  /**
   * \returns \p true when this element contains the point
   * \p p.  Customized for this \p InfQuad4, since knowledge
   * about the envelope can help avoiding slightly more
   * expensive computations.
   */
  virtual bool contains_point (const Point & p, Real tol=TOLERANCE) const override;

  /**
   * Geometric constants for InfQuad4.
   */
  static const int num_nodes = 4;
  static const int num_sides = 3;
  static const int num_children = 2;
  static const int nodes_per_side = 2;

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

};


} // namespace libMesh

#endif // ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

#endif // LIBMESH_FACE_INF_QUAD4_H
