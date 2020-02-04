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



#ifndef LIBMESH_NODE_ELEM_H
#define LIBMESH_NODE_ELEM_H

// Local includes
#include "libmesh/elem.h"
#include "libmesh/enum_order.h"

namespace libMesh
{

/**
 * The \p NodeElem is a point element, generally used as
 * a side of a 1D element.
 *
 * \author Roy H. Stogner
 * \date 2006
 * \brief A zero-dimensional geometric entity implementing the Elem interface.
 */
template <typename RealType = Real>
class NodeElemTempl : public ElemTempl<RealType>
{
public:
  typedef NodeElemTempl<RealType> NodeElem;
  typedef ElemTempl<RealType> Elem;
  typedef NodeTempl<RealType> Node;
  typedef PointTempl<RealType> Point;

  /**
   * Constructor.  By default this element has no parent.
   */
  explicit
  NodeElemTempl (Elem * p=nullptr) :
    Elem(NodeElem::n_nodes(), NodeElem::n_sides(), p, _elemlinks_data,
         _nodelinks_data)
  {
    // Make sure the interior parent isn't undefined
    if (LIBMESH_DIM > 0)
      this->set_interior_parent(nullptr);
  }

  NodeElemTempl (NodeElem &&) = delete;
  NodeElemTempl (const NodeElem &) = delete;
  NodeElem & operator= (const NodeElem &) = delete;
  NodeElem & operator= (NodeElem &&) = delete;
  virtual ~NodeElemTempl() = default;

  /**
   * \returns The \p Point associated with local \p Node \p i,
   * in master element rather than physical coordinates.
   */
  virtual Point master_point (const unsigned int libmesh_dbg_var(i)) const override
  {
    libmesh_assert_equal_to (i, 0);
    return Point(0,0,0);
  }

  /**
   * \returns 0, the dimensionality of the object.
   */
  virtual unsigned short dim () const override { return 0; }

  /**
   * \returns 1.
   */
  virtual unsigned int n_nodes() const override { return 1; }

  /**
   * \returns 0.
   */
  virtual unsigned int n_sides() const override { return 0; }

  /**
   * \returns 1.  Every NodeElem is a vertex
   */
  virtual unsigned int n_vertices() const override { return 1; }

  /**
   * \returns 0.
   */
  virtual unsigned int n_edges() const override { return 0; }

  /**
   * \returns 0.
   */
  virtual unsigned int n_faces() const override { return 0; }

  /**
   * \returns 1.
   */
  virtual unsigned int n_children() const override { return 1; }

  /**
   * Don't hide Elem::key() defined in the base class.
   */
  using Elem::key;

  /**
   * \returns An id associated with the \p s side of this element.
   * This should never be important for NodeElems.
   */
  virtual dof_id_type key (const unsigned int) const override
  { libmesh_error_msg("Calling NodeElem::key(side) does not make sense."); return 0; }

  /**
   * NodeElems don't have sides, so they can't have nodes on sides.
   */
  virtual unsigned int which_node_am_i(unsigned int /*side*/,
                                       unsigned int /*side_node*/) const override
  { libmesh_error_msg("Calling NodeElem::which_node_am_i() does not make sense."); return 0; }

  /**
   * The \p Elem::side_ptr() member makes no sense for nodes.
   */
  virtual std::unique_ptr<Elem> side_ptr (const unsigned int) override
  { libmesh_not_implemented(); return std::unique_ptr<Elem>(); }

  virtual void side_ptr (std::unique_ptr<Elem> &, const unsigned int) override
  { libmesh_not_implemented(); }

  /**
   * The \p Elem::build_side_ptr() member makes no sense for nodes.
   */
  virtual std::unique_ptr<Elem> build_side_ptr (const unsigned int, bool) override
  { libmesh_not_implemented(); return std::unique_ptr<Elem>(); }

  virtual void build_side_ptr (std::unique_ptr<Elem> &, const unsigned int) override
  { libmesh_not_implemented(); }

  /**
   * The \p Elem::build_edge_ptr() member makes no sense for nodes.
   */
  virtual std::unique_ptr<Elem> build_edge_ptr (const unsigned int) override
  { libmesh_not_implemented(); return std::unique_ptr<Elem>(); }

  /**
   * \returns 1.
   */
  virtual unsigned int n_sub_elem() const override { return 1; }

  /**
   * \returns \p true if the specified (local) node number is a vertex.
   */
  virtual bool is_vertex(const unsigned int) const override { return true; }

  /**
   * NodeElem objects don't have faces or sides.
   */
  virtual bool is_edge(const unsigned int) const override { return false; }
  virtual bool is_face(const unsigned int) const override { return false; }

  virtual bool is_child_on_side(const unsigned int,
                                const unsigned int) const override
  { libmesh_not_implemented(); return false; }

  virtual bool is_node_on_side(const unsigned int,
                               const unsigned int) const override
  { libmesh_not_implemented(); return false; }

  virtual std::vector<unsigned int> nodes_on_side(const unsigned int) const override
  {
    libmesh_not_implemented();
    return {0};
  }

  virtual bool is_node_on_edge(const unsigned int,
                               const unsigned int) const override
  { libmesh_not_implemented(); return false; }

  virtual bool is_edge_on_side(const unsigned int,
                               const unsigned int) const override
  { libmesh_not_implemented(); return false; }

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
   * \returns \p NODEELEM.
   */
  virtual ElemType type() const override { return NODEELEM; }

  /**
   * \returns FIRST.
   */
  virtual Order default_order() const override;

  virtual void connectivity(const unsigned int sc,
                            const IOPackage iop,
                            std::vector<dof_id_type> & conn) const override;


#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

  /**
   * \returns \p false.
   */
  virtual bool infinite () const override { return false; }

#endif


protected:

  /**
   * Data for links to parent/neighbor/interior_parent elements.
   */
  Elem * _elemlinks_data[1+(LIBMESH_DIM>0)];

  /**
   * Data for links to nodes.
   */
  Node * _nodelinks_data[1];


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
  static const float _embedding_matrix[1][1][1];

  /**
   * Matrix that allows children to inherit boundary conditions.
   */
  unsigned int side_children_matrix (const unsigned int,
                                     const unsigned int) const
  { libmesh_not_implemented(); return 0; }

  LIBMESH_ENABLE_TOPOLOGY_CACHES;

#endif // LIBMESH_ENABLE_AMR

};

template <typename RealType>
Order NodeElemTempl<RealType>::default_order() const
{
  return FIRST;
}



template <typename RealType>
void NodeElemTempl<RealType>::connectivity(const unsigned int,
                            const IOPackage,
                            std::vector<dof_id_type> &) const
{
  libmesh_not_implemented();
}

#ifdef LIBMESH_ENABLE_AMR

template <typename RealType>
const float NodeElemTempl<RealType>::_embedding_matrix[1][1][1] =
  {
    // embedding matrix for child 0
    {
      // 0
      {1.0}, // 0
    }
  };

#endif

typedef NodeElemTempl<Real> NodeElem;

} // namespace libMesh

#endif // LIBMESH_NODE_ELEM_H
