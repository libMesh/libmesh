// $Id: edge_edge2.h,v 1.6 2005-02-19 19:07:31 roystgnr Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2004  Benjamin S. Kirk, John W. Peterson
  
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



#ifndef __edge2_h__
#define __edge2_h__

// C++ includes


// Local includes
#include "libmesh_common.h"
#include "edge.h"


/**
 * The \p Edge2 is an element in 1D composed of 2 nodes. It is numbered
 * like this:
 *
   \verbatim
    EDGE2: o--------o
           0        1
   \endverbatim
 */

// ------------------------------------------------------------
// Edge class definition
class Edge2 : public Edge
{
 public:

  /**
   * Constructor.  By default this element has no parent.
   */
  Edge2 (const Elem* p=NULL) :
    Edge(Edge2::n_nodes(), p) {} 

  /**
   * Constructor.  Explicitly specifies the number of
   * nodes and neighbors for which storage will be allocated.
   */
  Edge2 (const unsigned int nn,
	 const unsigned int ns,
	 const Elem* p) :
    Edge(nn, p) { assert (ns == 0); }

  /**
   * @returns 1
   */
  unsigned int n_sub_elem() const { return 1; }
  
  /**
   * @returns true iff the specified (local) node number is a vertex.
   */
  virtual bool is_vertex(const unsigned int i) const;

  /**
   * @returns true iff the specified (local) node number is an edge.
   */
  virtual bool is_edge(const unsigned int i) const;

  /**
   * @returns true iff the specified (local) node number is a face.
   */
  virtual bool is_face(const unsigned int i) const;
  
  /**
   * @returns \p EDGE2
   */
  ElemType type()  const { return EDGE2; }
  
  /**
   * @returns FIRST
   */
  Order default_order() const { return FIRST; }

  virtual void connectivity(const unsigned int sc,
			    const IOPackage iop,
			    std::vector<unsigned int>& conn) const;


#ifdef ENABLE_INFINITE_ELEMENTS

  /**
   * @returns \p false.  This is a finite element. 
   */
  bool infinite () const { return false; }

#endif

  
protected:

  
#ifdef ENABLE_AMR

  /**
   * Matrix used to create the elements children.
   */
  float embedding_matrix (const unsigned int,
			 const unsigned int,
			 const unsigned int) const
  { error(); return 0.; }

#endif
};


#endif
