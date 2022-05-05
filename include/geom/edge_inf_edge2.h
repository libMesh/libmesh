// The libMesh Finite Element Library.
// Copyright (C) 2002-2022 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_EDGE_INF_EDGE2_H
#define LIBMESH_EDGE_INF_EDGE2_H

#include "libmesh/libmesh_common.h"

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

// Local includes
#include "libmesh/edge.h"

namespace libMesh
{

/**
 * The \p InfEdge2 is an infinite element in 1D composed of 2 nodes.
 * It is numbered like this:
 *
 * \verbatim
 *   INFEDGE2:
 *       o         closer to infinity
 *       | 1
 *       |
 *       |
 *       |
 *       o         base node
 *         0
 * \endverbatim
 *
 * \author Daniel Dreyer
 * \date 2003
 * \brief A 1D infinite element with 2 nodes.
 */
class InfEdge2 : public Edge
{
public:

  /**
   * Constructor.  By default this element has no parent.
   */
  explicit
  InfEdge2 (Elem * p=nullptr) :
    Edge(InfEdge2::n_nodes(), p, _nodelinks_data) {}

  InfEdge2 (InfEdge2 &&) = delete;
  InfEdge2 (const InfEdge2 &) = delete;
  InfEdge2 & operator= (const InfEdge2 &) = delete;
  InfEdge2 & operator= (InfEdge2 &&) = delete;
  virtual ~InfEdge2() = default;

  /**
   * \returns The \p Point associated with local \p Node \p i,
   * in master element rather than physical coordinates.
   */
  virtual Point master_point (const unsigned int i) const override
  {
    libmesh_assert_less(i, this->n_nodes());
    return Point(0,i,0);
  }

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

  /**
   * \returns \p true if the specified (local) node number is on the
   * specified edge (always true in 1D).
   */
  virtual bool is_node_on_edge(const unsigned int n,
                               const unsigned int e) const override;

  /**
   * \returns \p INFEDGE2.
   */
  virtual ElemType type() const override { return INFEDGE2; }

  /**
   * \returns FIRST.
   */
  virtual Order default_order() const override;

  virtual void connectivity(const unsigned int se,
                            const IOPackage iop,
                            std::vector<dof_id_type> & conn) const override;

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

  /**
   * \returns \p true.  This is an infinite element.
   */
  virtual bool infinite () const override { return true; }

  /**
   * \returns The origin of this infinite element.
   */
  virtual Point origin () const override;

  /**
   * \returns \p true if the specified (local) node number is a
   * "mid-edge" node on an infinite element edge.
   */
  virtual bool is_mid_infinite_edge_node(const unsigned int i) const
    override { return (i > 0); }

#endif


protected:

  /**
   * Data for links to nodes.
   */
  Node * _nodelinks_data[2];



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




// ------------------------------------------------------------
// InfEdge2 class member functions
inline
Point InfEdge2::origin () const
{
  return ( this->point(0)*2 - this->point(1) );
}


} // namespace libMesh

#endif

#endif // LIBMESH_EDGE_INF_EDGE2_H
