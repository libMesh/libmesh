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



#ifndef __edge4_h__
#define __edge4_h__

// Local includes
#include "libmesh_common.h"
#include "edge.h"

/**
 * The \p Edge4 is an element in 1D composed of 4 nodes. It is numbered
 * like this:
 *
 * \verbatim
 *  EGDE4: o----o----o----o
 *         0    2    3    1
 * \endverbatim
 */

// ------------------------------------------------------------
// Edge4 class definition
class Edge4 : public Edge
{
  public:

  /**
   * Constructor. By default this element has no parent.
   */
  Edge4 (Elem* p=NULL) :
    Edge(Edge4::n_nodes(), p) {}

 /**
   * Constructor.  Explicitly specifies the number of
   * nodes and neighbors for which storage will be allocated.
   */
  Edge4(const unsigned int nn,
        const unsigned int ns,
	Elem* p) :
    Edge(nn, p) { assert (ns == 0); }

  /**
   * @returns 4
   */
  unsigned int n_nodes() const { return 4; }

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
  
  /*
   * @returns true iff the specified (local) node number is on the
   * specified edge (i.e. "returns true" in 1D)
   */
  virtual bool is_node_on_edge(const unsigned int n,
			       const unsigned int e) const;
  
  /*
   * @returns true iff the element map is definitely affine within
   * numerical tolerances
   */
  virtual bool has_affine_map () const;

  /**
   * @returns \p EDGE4
   */
  ElemType type()  const { return EDGE4; }
  
  /**
   * @returns THIRD
   */
  Order default_order() const { return THIRD; }

  virtual void connectivity(const unsigned int sc,
			    const IOPackage iop,
			    std::vector<unsigned int>& conn) const;

  /**
   */
  unsigned int n_second_order_adjacent_vertices (const unsigned int) const
      { libmesh_error(); return 0;  }

  /**
   */
  unsigned short int second_order_adjacent_vertex (const unsigned int,
						   const unsigned int) const
      { libmesh_error(); return 0; }


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
  float embedding_matrix (const unsigned int i,
			 const unsigned int j,
			 const unsigned int k) const
  { return _embedding_matrix[i][j][k]; }
  
  /**
   * Matrix that computes new nodal locations/solution values
   * from current nodes/solution.
   */
  static const float _embedding_matrix[2][4][4];


#endif
};


#endif
