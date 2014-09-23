// The libMesh Finite Element Library.
// Copyright (C) 2002-2014 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

// C++ includes
#include <cstddef>

namespace libMesh
{


// Forward declarations


/**
 * The \p NodeElem is a point element, generally used as
 * a side of a 1D element.
 */
class NodeElem : public Elem
{
public:

  /**
   * Constructor.  By default this element has no parent.
   */
  explicit
  NodeElem (Elem* p=NULL) :
    Elem(NodeElem::n_nodes(), NodeElem::n_sides(), p, _elemlinks_data, _nodelinks_data) {}

  /**
   * @returns 0, the dimensionality of the object.
   */
  unsigned int dim () const { return 0; }

  /**
   * @returns 1.
   */
  unsigned int n_nodes() const { return 1; }

  /**
   * @returns 0
   */
  unsigned int n_sides() const { return 0; }

  /**
   * @returns 1.  Every NodeElem is a vertex
   */
  unsigned int n_vertices() const { return 1; }

  /**
   * @returns 0.
   */
  unsigned int n_edges() const { return 0; }

  /**
   * @returns 0.
   */
  unsigned int n_faces() const { return 0; }

  /**
   * @returns 1
   */
  unsigned int n_children() const { return 1; }

  /**
   * @returns an id associated with the \p s side of this element.
   * This should never be important for NodeElems
   */
  dof_id_type key (const unsigned int) const
  { return 0; }

  /**
   * The \p Elem::side() member makes no sense for nodes.
   */
  UniquePtr<Elem> side (const unsigned int) const
  { libmesh_not_implemented(); return UniquePtr<Elem>(); }

  /**
   * The \p Elem::build_side() member makes no sense for nodes.
   */
  UniquePtr<Elem> build_side (const unsigned int, bool) const
  { libmesh_not_implemented(); return UniquePtr<Elem>(); }

  /**
   * The \p Elem::build_edge() member makes no sense for nodes.
   */
  UniquePtr<Elem> build_edge (const unsigned int) const
  { libmesh_not_implemented(); return UniquePtr<Elem>(); }

  /**
   * @returns 1
   */
  unsigned int n_sub_elem() const { return 1; }

  /**
   * @returns true iff the specified (local) node number is a vertex.
   */
  virtual bool is_vertex(const unsigned int) const { return true; }

  /**
   * NodeElem objects don't have faces or sides
   */
  virtual bool is_edge(const unsigned int) const { return false; }

  virtual bool is_face(const unsigned int) const { return false; }

  virtual bool is_child_on_side(const unsigned int,
                                const unsigned int) const
  { libmesh_not_implemented(); return false; }

  virtual bool is_node_on_side(const unsigned int,
                               const unsigned int) const
  { libmesh_not_implemented(); return false; }

  virtual bool is_node_on_edge(const unsigned int,
                               const unsigned int) const
  { libmesh_not_implemented(); return false; }

  virtual bool is_edge_on_side(const unsigned int,
                               const unsigned int) const
  { libmesh_not_implemented(); return false; }

  /*
   * @returns true iff the element map is definitely affine within
   * numerical tolerances
   */
  virtual bool has_affine_map () const { return true; }

  /**
   * @returns true iff the Lagrange shape functions on this element
   * are linear
   */
  virtual bool is_linear () const { return true; }

  /**
   * @returns \p NODEELEM
   */
  ElemType type()  const { return NODEELEM; }

  /**
   * @returns FIRST
   */
  Order default_order() const { return FIRST; }

  virtual void connectivity(const unsigned int sc,
                            const IOPackage iop,
                            std::vector<dof_id_type>& conn) const;


#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

  /**
   * @returns \p false.
   */
  bool infinite () const { return false; }

#endif


protected:

  /**
   * Data for links to parent/neighbor/interior_parent elements.
   */
  Elem* _elemlinks_data[1+(LIBMESH_DIM>0)];

  /**
   * Data for links to nodes
   */
  Node* _nodelinks_data[1];


#ifdef LIBMESH_ENABLE_AMR

  /**
   * Matrix used to create the elements children.
   */
  float embedding_matrix (const unsigned int i,
                          const unsigned int j,
                          const unsigned int k) const
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

#endif

};





// ------------------------------------------------------------
// NodeElem class member functions

} // namespace libMesh

#endif // LIBMESH_NODE_ELEM_H
