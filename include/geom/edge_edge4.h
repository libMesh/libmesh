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



#ifndef LIBMESH_EDGE_EDGE4_H
#define LIBMESH_EDGE_EDGE4_H

// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/edge.h"

namespace libMesh
{

/**
 * The \p Edge4 is an element in 1D composed of 4 nodes.
 *
 * It is numbered like this:
 *
 * \verbatim
 *   EGDE4: o----o----o----o
 *          0    2    3    1
 * \endverbatim
 *
 * \author David Knezevic
 * \date 2005
 * \brief A 1D geometric element with 4 nodes.
 */
class Edge4 final : public Edge
{
public:

  /**
   * Constructor. By default this element has no parent.
   */
  explicit
  Edge4 (Elem * p=nullptr) :
    Edge(Edge4::n_nodes(), p, _nodelinks_data) {}

  Edge4 (Edge4 &&) = delete;
  Edge4 (const Edge4 &) = delete;
  Edge4 & operator= (const Edge4 &) = delete;
  Edge4 & operator= (Edge4 &&) = delete;
  virtual ~Edge4() = default;

  /**
   * \returns The \p Point associated with local \p Node \p i,
   * in master element rather than physical coordinates.
   */
  virtual Point master_point (const unsigned int i) const override
  {
    libmesh_assert_less(i, this->n_nodes());
    if (i < 2)
      return Point(2.0f*Real(i)-1.0f,0,0);
    return Point((Real(2)*Real(i)-5)/3,0,0);
  }

  /**
   * \returns 4.
   */
  virtual unsigned int n_nodes() const override { return num_nodes; }

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

  /**
   * \returns \p true if the specified (local) node number is on the
   * specified edge (always true in 1D).
   */
  virtual bool is_node_on_edge(const unsigned int n,
                               const unsigned int e) const override;

  /**
   * \returns \p true if the element map is definitely affine within
   * numerical tolerances.
   */
  virtual bool has_affine_map () const override;

  /**
   * \returns \p EDGE4.
   */
  virtual ElemType type() const override { return EDGE4; }

  /**
   * \returns THIRD.
   */
  virtual Order default_order() const override;

  virtual void connectivity(const unsigned int sc,
                            const IOPackage iop,
                            std::vector<dof_id_type> & conn) const override;

  /**
   * FIXME: This function could be generalized to work for Edges.
   */
  virtual unsigned int n_second_order_adjacent_vertices (const unsigned int) const override
  { libmesh_not_implemented(); return 0;  }

  /**
   * FIXME: This function could be generalized to work for Edges.
   */
  virtual unsigned short int second_order_adjacent_vertex (const unsigned int,
                                                           const unsigned int) const override
  { libmesh_not_implemented(); return 0; }

  /**
   * \returns A bounding box (not necessarily the minimal bounding box)
   * containing the edge.
   */
  virtual BoundingBox loose_bounding_box () const override;

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

  /**
   * \returns \p false.  This is a finite element.
   */
  virtual bool infinite () const override { return false; }

#endif

  /**
   * Don't hide Edge::key(side) defined in the base class.
   */
  using Edge::key;

  /**
   * \returns An id associated with the global node ids of this
   * element.  The id is not necessarily unique, but should be
   * close.
   */
  virtual dof_id_type key () const override;

  /**
   * An optimized method for approximating the length of an
   * EDGE4 using quadrature.
   */
  virtual Real volume () const override;

  /**
   * Geometric constants for Edge4.
   */
  static const int num_nodes = 4;
  static const int num_sides = 2;
  static const int num_edges = 0;
  static const int num_children = 2;
  static const int nodes_per_side = 1;
  static const int nodes_per_edge = invalid_int;

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


#endif // LIBMESH_EDGE_EDGE4_H
