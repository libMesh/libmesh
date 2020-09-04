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



#ifndef LIBMESH_EDGE_H
#define LIBMESH_EDGE_H

// Local includes
#include "libmesh/elem.h"

namespace libMesh
{

/**
 * The \p Edge is an element in 1D. It can be thought of as a
 * line segment.
 *
 * \author Benjamin S. Kirk
 * \date 2002
 * \brief The base class for all 1D geometric element types.
 */
class Edge : public Elem
{
public:

  /**
   * Default line element, takes number of nodes and
   * parent. Derived classes implement 'true' elements.
   */
  Edge (const unsigned int nn,
        Elem * p,
        Node ** nodelinkdata) :
    Elem(nn, Edge::n_sides(), p, _elemlinks_data, nodelinkdata)
  {
    // Make sure the interior parent isn't undefined
    if (LIBMESH_DIM > 1)
      this->set_interior_parent(nullptr);
  }

  Edge (Edge &&) = delete;
  Edge (const Edge &) = delete;
  Edge & operator= (const Edge &) = delete;
  Edge & operator= (Edge &&) = delete;
  virtual ~Edge() = default;

  /**
   * \returns 1, the dimensionality of the object.
   */
  virtual unsigned short dim () const override final { return 1; }

  /**
   * \returns 2. Every edge is guaranteed to have at least 2 nodes.
   */
  virtual unsigned int n_nodes() const override { return 2; }

  /**
   * \returns 2.
   */
  virtual unsigned int n_sides() const override final { return 2; }

  /**
   * \returns 2.  Every edge has exactly two vertices.
   */
  virtual unsigned int n_vertices() const override final { return 2; }

  /**
   * \returns 0.  All 1D elements have no edges.
   */
  virtual unsigned int n_edges() const override final { return 0; }

  /**
   * \returns 0.  All 1D elements have no faces.
   */
  virtual unsigned int n_faces() const override final { return 0; }

  /**
   * \returns 2.
   */
  virtual unsigned int n_children() const override final { return 2; }

  /**
   * \returns \p true if the specified child is on the specified side.
   */
  virtual bool is_child_on_side(const unsigned int c,
                                const unsigned int s) const override final;

  /**
   * \returns \p true if the specified edge is on the specified side.
   */
  virtual bool is_edge_on_side(const unsigned int,
                               const unsigned int) const override final
  { return false; }

  /**
   * \returns The side number opposite to \p s (for a tensor product
   * element), or throws an error otherwise.
   */
  virtual unsigned int opposite_side(const unsigned int s) const override final;

  /**
   * \returns The local node number for the node opposite to node n
   * on side \p opposite_side(s) (for a tensor product element), or
   * throws an error otherwise.
   */
  virtual unsigned int opposite_node(const unsigned int n,
                                     const unsigned int s) const override final;

  /**
   * Don't hide Elem::key() defined in the base class.
   */
  using Elem::key;

  /**
   * \returns An id associated with the \p s side of this element.
   * The id is not necessarily unique, but should be close.  This is
   * particularly useful in the \p MeshBase::find_neighbors() routine.
   */
  virtual dof_id_type key (const unsigned int s) const override final
  { return this->compute_key(this->node_id(s)); }

  /**
   * \returns \p side after doing some range checking. \p side_node is ignored.
   */
  virtual unsigned int local_side_node(unsigned int side,
                                       unsigned int /*side_node*/) const override final;

  /**
   * Throws an error. Edge elems have n_edges() == 0, so it does not
   * make sense to call local_edge_node().
   */
  virtual unsigned int local_edge_node(unsigned int edge,
                                       unsigned int edge_node) const override final;

  /**
   * \returns A pointer to a NodeElem for the specified node.
   */
  virtual std::unique_ptr<Elem> side_ptr (const unsigned int i) override final;

  /**
   * Rebuilds a pointer to a NodeElem for the specified node.
   */
  virtual void side_ptr (std::unique_ptr<Elem> & side, const unsigned int i) override final;

  /**
   * \returns A pointer to a NodeElem for the specified node.
   */
  virtual std::unique_ptr<Elem> build_side_ptr (const unsigned int i,
                                                bool proxy=true) override final;

  /**
   * Rebuilds a NODEELEM for the specified node.
   */
  virtual void build_side_ptr (std::unique_ptr<Elem> & elem,
                               const unsigned int i) override final;

  /**
   * The \p Elem::build_edge_ptr() member makes no sense for edges.
   */
  virtual std::unique_ptr<Elem> build_edge_ptr (const unsigned int) override final
  { libmesh_not_implemented(); return std::unique_ptr<Elem>(); }

  virtual std::vector<unsigned int> nodes_on_side(const unsigned int s) const override;

  virtual std::vector<unsigned int> nodes_on_edge(const unsigned int e) const override;

  virtual std::vector<unsigned int> sides_on_edge(const unsigned int) const override final
  { return {}; }

protected:

  /**
   * Data for links to parent/neighbor/interior_parent elements.
   */
  Elem * _elemlinks_data[3+(LIBMESH_DIM>1)];

#ifdef LIBMESH_ENABLE_AMR

  /**
   * Matrix that allows children to inherit boundary conditions.
   */
  unsigned int side_children_matrix (const unsigned int,
                                     const unsigned int) const
  { libmesh_not_implemented(); return 0; }

#endif

};

} // namespace libMesh

#endif // LIBMESH_EDGE_H
