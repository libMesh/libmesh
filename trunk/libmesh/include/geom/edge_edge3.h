// $Id: edge_edge3.h,v 1.8 2005-02-25 19:15:41 roystgnr Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2005  Benjamin S. Kirk, John W. Peterson
  
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



#ifndef __edge3_h__
#define __edge3_h__


// Local includes
#include "libmesh_common.h"
#include "edge.h"



/**
 * The \p Edge3 is an element in 1D composed of 3 nodes. It is numbered
 * like this:
 *
 * \verbatim
 *  EGDE3: o----o----o
 *         0    2    1
 * \endverbatim
 */

// ------------------------------------------------------------
// Edge3 class definition
class Edge3 : public Edge
{
 public:

  /**
   * Constructor.  By default this element has no parent.
   */
  Edge3 (const Elem* p=NULL) :
    Edge(Edge3::n_nodes(), p) {}
  
  /**
   * Constructor.  Explicitly specifies the number of
   * nodes and neighbors for which storage will be allocated.
   */
  Edge3 (const unsigned int nn,
	 const unsigned int ns,
	 const Elem* p) :
    Edge(nn, p) { assert (ns == 0); }

  
  /**
   * @returns 3
   */
  unsigned int n_nodes() const { return 3; }

  /**
   * @returns 2
   */
  unsigned int n_sub_elem() const { return 2; }

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
  
  /*
   * @returns true iff the specified (local) node number is on the
   * specified side
   */
  virtual bool is_node_on_side(const unsigned int n,
			       const unsigned int s) const;
  
  /**
   * @returns \p EDGE3
   */
  ElemType type()  const { return EDGE3; }
  
  /**
   * @returns SECOND
   */
  Order default_order() const { return SECOND; }

  virtual void connectivity(const unsigned int sc,
			    const IOPackage iop,
			    std::vector<unsigned int>& conn) const;

  /**
   * @returns 2 for all \p n
   */
  unsigned int n_second_order_adjacent_vertices (const unsigned int) const
      { return 2; }

  /**
   * @returns the element-local number of the  \f$ v^{th} \f$ vertex
   * that defines the \f$ n^{th} \f$ second-order node.
   */
  unsigned short int second_order_adjacent_vertex (const unsigned int,
						   const unsigned int v) const
      { return static_cast<unsigned short int>(v); }


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
