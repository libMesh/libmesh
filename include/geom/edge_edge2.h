// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/enum_io_package.h"
#include "libmesh/enum_order.h"

namespace libMesh
{

/**
 * The \p Edge2 is an element in 1D composed of 2 nodes. It is numbered
 * like this:
 *
 * \verbatim
 *   EDGE2: o--------o
 *          0        1
 * \endverbatim
 *
 * \author Benjamin S. Kirk
 * \date 2002
 * \brief A 1D geometric element with 2 nodes.
 */
template <typename RealType = Real>
class Edge2Templ : public EdgeTempl<RealType>
{
public:
  typedef Edge2Templ<RealType> Edge2;
  using typename EdgeTempl<RealType>::Edge;
  using typename EdgeTempl<RealType>::Elem;

  /**
   * Constructor.  By default this element has no parent.
   */
  explicit
  Edge2Templ (Elem * p=nullptr) :
    Edge(Edge2::n_nodes(), p, _nodelinks_data) {}

  Edge2Templ (Edge2 &&) = delete;
  Edge2Templ (const Edge2 &) = delete;
  Edge2 & operator= (const Edge2 &) = delete;
  Edge2 & operator= (Edge2 &&) = delete;
  virtual ~Edge2Templ() = default;

  /**
   * \returns The \p Point associated with local \p Node \p i,
   * in master element rather than physical coordinates.
   */
  virtual Point master_point (const unsigned int i) const override
  {
    libmesh_assert_less(i, this->n_nodes());
    return Point(2.0f*Real(i)-1.0f,0,0);
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
  virtual bool has_affine_map () const override { return true; }

  /**
   * \returns \p true if the Lagrange shape functions on this element
   * are linear.
   */
  virtual bool is_linear () const override { return true; }

  /**
   * \returns \p EDGE2.
   */
  virtual ElemType type() const override { return EDGE2; }

  /**
   * \returns FIRST.
   */
  virtual Order default_order() const override;

  virtual void connectivity(const unsigned int sc,
                            const IOPackage iop,
                            std::vector<dof_id_type> & conn) const override;

  /**
   * An optimized method for computing the length of a 2-node edge.
   */
  virtual RealType volume () const override;

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
   * Geometric constants for Edge2.
   */
  static const int num_nodes = 2;
  static const int num_children = 2;

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

template <typename RealType>
const int Edge2Templ<RealType>::num_nodes;

#ifdef LIBMESH_ENABLE_AMR

template <typename RealType>
const float Edge2Templ<RealType>::
_embedding_matrix[Edge2Templ<RealType>::num_children]
                 [Edge2Templ<RealType>::num_nodes]
                 [Edge2Templ<RealType>::num_nodes] =
  {
    // embedding matrix for child 0
    {
      // 0    1
      {1.0, 0.0}, // 0
      {0.5, 0.5}  // 1
    },

    // embedding matrix for child 1
    {
      // 0    1
      {0.5, 0.5}, // 0
      {0.0, 1.0}  // 1
    }
  };

#endif

template <typename RealType>
bool Edge2Templ<RealType>::is_vertex(const unsigned int) const
{
  return true;
}

template <typename RealType>
bool Edge2Templ<RealType>::is_edge(const unsigned int) const
{
  return false;
}

template <typename RealType>
bool Edge2Templ<RealType>::is_face(const unsigned int) const
{
  return false;
}

template <typename RealType>
bool Edge2Templ<RealType>::is_node_on_side(const unsigned int n,
                            const unsigned int s) const
{
  libmesh_assert_less (s, Edge2Templ<RealType>::num_nodes);
  return (s == n);
}

template <typename RealType>
bool Edge2Templ<RealType>::is_node_on_edge(const unsigned int,
                            const unsigned int libmesh_dbg_var(e)) const
{
  libmesh_assert_equal_to (e, 0);
  return true;
}



template <typename RealType>
Order Edge2Templ<RealType>::default_order() const
{
  return FIRST;
}



template <typename RealType>
void Edge2Templ<RealType>::connectivity(const unsigned int libmesh_dbg_var(sc),
                         const IOPackage iop,
                         std::vector<dof_id_type> & conn) const
{
  libmesh_assert_equal_to (sc, 0);
  libmesh_assert_less (sc, this->n_sub_elem());
  libmesh_assert_not_equal_to (iop, INVALID_IO_PACKAGE);

  // Create storage
  conn.resize(2);

  switch (iop)
    {
    case TECPLOT:
      {
        conn[0] = this->node_id(0)+1;
        conn[1] = this->node_id(1)+1;
        return;
      }

    case VTK:
      {
        conn[0] = this->node_id(0);
        conn[1] = this->node_id(1);
        return;
      }

    default:
      libmesh_error_msg("Unsupported IO package " << iop);
    }
}


template <typename RealType>
RealType Edge2Templ<RealType>::volume () const
{
  // OK, so this is probably overkill, since it is equivalent to
  // Elem::hmax() for the Edge2, but here it is nonetheless...
  return (this->point(1) - this->point(0)).norm();
}



template <typename RealType>
dof_id_type Edge2Templ<RealType>::key () const
{
  return this->compute_key(this->node_id(0),
                           this->node_id(1));
}


typedef Edge2Templ<Real> Edge2;

} // namespace libMesh


#endif // LIBMESH_EDGE_EDGE2_H
