// $Id: edge.h,v 1.6 2003-02-13 22:56:07 benkirk Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002  Benjamin S. Kirk, John W. Peterson
  
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



#ifndef __edge_h__
#define __edge_h__

// C++ includes

// Local includes
#include "elem.h"


// Forward declarations
class Mesh;



/**
 * The \p Edge is an element in 1D. It can be thought of as a
 * line segment.
 */

// ------------------------------------------------------------
// Edge class definition
class Edge : public Elem
{
 public:

  /**
   * Default line element, takes number of nodes and 
   * parent. Derived classes implement 'true' elements.
   */
  Edge (const unsigned int nn, Edge* p);

  /**
   * @returns 1, the dimensionality of the object.
   */
  unsigned int dim () const { return 1; }
  
  /**
   * @returns 2. Every edge is guaranteed to have at least 2 nodes.
   */
  unsigned int n_nodes() const { return 2; }

  /**
   * @returns 2
   */
  unsigned int n_sides() const { return 2; }

  /**
   * @returns 2.  Every edge has exactly two vertices.
   */
  unsigned int n_vertices() const { return 2; }
  
  /**
   * @returns 0.  All 1D elements have no edges.
   */  
  unsigned int n_edges() const { return 0; }
  
  /**
   * @returns 0.  All 1D elements have no faces.
   */  
  unsigned int n_faces() const { return 0; }
  
  /**
   * @returns 2
   */
  unsigned int n_children() const { return 2; }

  /**
   * The \p Elem::side() member makes no sense for edges.
   */
  AutoPtr<Elem> side (const unsigned int) const
  { error(); AutoPtr<Elem> ap(NULL); return ap; }

  /**
   * The \p Elem::build_side() member makes no sense for edges.
   */
  AutoPtr<Elem> build_side (const unsigned int) const
  { error(); AutoPtr<Elem> ap(NULL); return ap; }
  
 private:
  
};





// ------------------------------------------------------------
// Edge class member functions
inline
Edge::Edge(const unsigned int nn, Edge* p) :
  Elem(nn, Edge::n_sides(), p) 
{
}




#endif
