// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2007  Benjamin S. Kirk, John W. Peterson
  
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



#ifndef __node_elem_h__
#define __node_elem_h__

// C++ includes

// Local includes
#include "elem.h"


// Forward declarations


/**
 * The \p NodeElem is a point element, generally used as
 * a side of a 1D element.
 */

// ------------------------------------------------------------
// NodeElem class definition
class NodeElem : public Elem
{
 public:

  /**
   * Constructor.  By default this element has no parent.
   */
  NodeElem (Elem* p=NULL) :
    Elem(NodeElem::n_nodes(), NodeElem::n_sides(), p) {}

  /**
   * Constructor.  Explicitly specifies the number of
   * nodes and neighbors for which storage will be allocated.
   */
  NodeElem (const unsigned int nn,
         const unsigned int ns,
         Elem* p) :
    Elem(nn, ns, p) { assert(nn == 1); assert (ns == 0); }

  /**
   * Default node element, takes number of nodes and 
   * parent. Derived classes implement 'true' elements.
   */
  NodeElem (const unsigned int nn,
	Elem* p) :
    Elem(nn, NodeElem::n_sides(), p) {}
   
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
  unsigned int key (const unsigned int) const
  { return 0; }
  
  /**
   * The \p Elem::side() member makes no sense for nodes.
   */
  AutoPtr<DofObject> side (const unsigned int) const
  { libmesh_error(); AutoPtr<DofObject> ap(NULL); return ap; }

  /**
   * The \p Elem::build_side() member makes no sense for nodes.
   */
  AutoPtr<Elem> build_side (const unsigned int, bool) const
  { libmesh_error(); AutoPtr<Elem> ap(NULL); return ap; }

  /**
   * The \p Elem::build_edge() member makes no sense for nodes.
   */
  AutoPtr<Elem> build_edge (const unsigned int) const
  { libmesh_error(); AutoPtr<Elem> ap(NULL); return ap; }

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
  { libmesh_error(); return false; }
  
  virtual bool is_node_on_side(const unsigned int,
			       const unsigned int) const
  { libmesh_error(); return false; }
  
  virtual bool is_node_on_edge(const unsigned int,
			       const unsigned int) const 
  { libmesh_error(); return false; }
  
  /*
   * @returns true iff the element map is definitely affine within
   * numerical tolerances
   */
  virtual bool has_affine_map () const { return true; }

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
			    std::vector<unsigned int>& conn) const;


#ifdef ENABLE_INFINITE_ELEMENTS

  /**
   * @returns \p false.
   */
  bool infinite () const { return false; }

#endif

  
protected:

  
#ifdef ENABLE_AMR

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
  { libmesh_error(); return 0; }

#endif
  
};





// ------------------------------------------------------------
// NodeElem class member functions

#endif
