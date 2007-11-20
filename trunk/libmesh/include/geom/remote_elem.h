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



#ifndef __remote_elem_h__
#define __remote_elem_h__

// Local includes
#include "elem.h"

// Forward declarations



/**
 * In parallel meshes where a ghost element has neighbors which do 
 * not exist on the local processor, the ghost element's neighbors
 * are set to point to the singleton RemoteElement instead.
 * Library code can then distinguish between such elements and
 * boundary elements (with NULL neighbors).
 *
 * \author Roy H. Stogner, 2007
 */

// ------------------------------------------------------------
// RemoteElem class definition
class RemoteElem : public Elem
{
 public:

  /**
   * Constructor.
   */ 
  RemoteElem () : Elem() {}

  /**
   * Destructor.
   */
  virtual ~RemoteElem() {}

  virtual const Point & point (const unsigned int i) const
  { error(); return Elem::point(i); }

  virtual Point & point (const unsigned int i)
  { error(); return Elem::point(i); }

  virtual unsigned int node (const unsigned int i) const
  { error(); return Elem::node(i); }

  virtual Node* get_node (const unsigned int i) const
  { error(); return Elem::get_node(i); }

  virtual Node* & set_node (const unsigned int i)
  { error(); return Elem::set_node(i); }
  
  virtual unsigned int key (const unsigned int) const
  { error(); return 0; }
  
  virtual void connectivity(const unsigned int,
			    const IOPackage,
			    std::vector<unsigned int>&) const
  { error(); }

  virtual ElemType type () const
  { return REMOTEELEM; }
  
  virtual unsigned int dim () const
  { error(); return 0; }
  
  virtual unsigned int n_nodes () const
  { error(); return 0; }

  virtual unsigned int n_sides () const 
  { error(); return 0; }

  virtual unsigned int n_vertices () const
  { error(); return 0; }

  virtual unsigned int n_edges () const
  { error(); return 0; }

  virtual unsigned int n_faces () const
  { error(); return 0; }
  
  virtual unsigned int n_children () const
  { error(); return 0; }

  virtual bool is_vertex(const unsigned int) const
  { error(); return false; }

  virtual bool is_edge(const unsigned int) const
  { error(); return false; }

  virtual bool is_face(const unsigned int) const
  { error(); return false; }
  
  virtual bool is_node_on_side(const unsigned int,
			       const unsigned int) const
  { error(); return false; }
  
  virtual bool is_node_on_edge(const unsigned int,
			       const unsigned int) const
  { error(); return false; }

  virtual unsigned int n_sub_elem () const
  { error(); return 0; }

  virtual AutoPtr<DofObject> side (const unsigned int) const
  { error(); return AutoPtr<DofObject>(NULL); }
  
  virtual AutoPtr<Elem> build_side (const unsigned int,
				    bool) const
  { error(); return AutoPtr<Elem>(NULL); }

  virtual AutoPtr<Elem> build_edge (const unsigned int) const
  { error(); return AutoPtr<Elem>(NULL); }

  virtual Order default_order () const
  { error(); return FIRST; }
  
#ifdef ENABLE_INFINITE_ELEMENTS

  virtual bool infinite () const
  { error(); return false; }

#endif

#ifdef ENABLE_AMR
  
  /**
   * Matrix that transforms the parents nodes into the children's
   * nodes
   */
  virtual float embedding_matrix (const unsigned int,
				  const unsigned int,
				  const unsigned int) const
  { error(); return 0.; }

#endif

};

// Singleton RemoteElem
extern const RemoteElem* remote_elem;

#endif // end #ifndef __remote_elem_h__
