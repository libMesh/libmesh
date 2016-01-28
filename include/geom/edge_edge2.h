// The libMesh Finite Element Library.
// Copyright (C) 2002-2016 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_EDGE_EDGE2_H
#define LIBMESH_EDGE_EDGE2_H

// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/edge.h"

// C++ includes
#include <cstddef>

namespace libMesh
{


/**
 * The \p Edge2 is an element in 1D composed of 2 nodes. It is numbered
 * like this:
 *
 * \verbatim
 *  EDGE2: o--------o
 *         0        1
 * \endverbatim
 */
class Edge2 : public Edge
{
public:

  /**
   * Constructor.  By default this element has no parent.
   */
  explicit
  Edge2 (Elem * p=libmesh_nullptr) :
    Edge(Edge2::n_nodes(), p, _nodelinks_data) {}

  /**
   * @returns the \p Point associated with local \p Node \p i,
   * in master element rather than physical coordinates.
   */
  virtual Point master_point (const unsigned int i) const libmesh_override
  {
    libmesh_assert_less(i, this->n_nodes());
    return Point(2.0f*i-1,0,0);
  }

  /**
   * @returns 1
   */
  virtual unsigned int n_sub_elem() const libmesh_override { return 1; }

  /**
   * @returns true iff the specified (local) node number is a vertex.
   */
  virtual bool is_vertex(const unsigned int i) const libmesh_override;

  /**
   * @returns true iff the specified (local) node number is an edge.
   */
  virtual bool is_edge(const unsigned int i) const libmesh_override;

  /**
   * @returns true iff the specified (local) node number is a face.
   */
  virtual bool is_face(const unsigned int i) const libmesh_override;

  /*
   * @returns true iff the specified (local) node number is on the
   * specified side
   */
  virtual bool is_node_on_side(const unsigned int n,
                               const unsigned int s) const libmesh_override;

  /*
   * @returns true iff the specified (local) node number is on the
   * specified edge (i.e. "returns true" in 1D)
   */
  virtual bool is_node_on_edge(const unsigned int n,
                               const unsigned int e) const libmesh_override;

  /*
   * @returns true iff the element map is definitely affine within
   * numerical tolerances
   */
  virtual bool has_affine_map () const libmesh_override { return true; }

  /**
   * @returns true iff the Lagrange shape functions on this element
   * are linear
   */
  virtual bool is_linear () const libmesh_override { return true; }

  /**
   * @returns \p EDGE2
   */
  virtual ElemType type() const libmesh_override { return EDGE2; }

  /**
   * @returns FIRST
   */
  virtual Order default_order() const libmesh_override { return FIRST; }

  virtual void connectivity(const unsigned int sc,
                            const IOPackage iop,
                            std::vector<dof_id_type> & conn) const libmesh_override;

  /**
   * An optimized method for computing the length of a 2-node edge.
   */
  virtual Real volume () const libmesh_override;

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

  /**
   * @returns \p false.  This is a finite element.
   */
  virtual bool infinite () const libmesh_override { return false; }

#endif

  /**
   * Don't hide Edge::key(side) defined in the base class.
   */
  using Edge::key;

  /**
   * @returns an id associated with the global node ids of this
   * element.  The id is not necessarily unique, but should be
   * close.
   */
  virtual dof_id_type key () const libmesh_override;

protected:

  /**
   * Data for links to nodes
   */
  Node * _nodelinks_data[2];



#ifdef LIBMESH_ENABLE_AMR

  /**
   * Matrix used to create the elements children.
   */
  virtual float embedding_matrix (const unsigned int i,
                                  const unsigned int j,
                                  const unsigned int k) const libmesh_override
  { return _embedding_matrix[i][j][k]; }

  /**
   * Matrix that computes new nodal locations/solution values
   * from current nodes/solution.
   */
  static const float _embedding_matrix[2][2][2];

  LIBMESH_ENABLE_TOPOLOGY_CACHES;

#endif // LIBMESH_ENABLE_AMR

};

} // namespace libMesh


#endif // LIBMESH_EDGE_EDGE2_H
