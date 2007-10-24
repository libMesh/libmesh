// $Id: face_inf_quad6.h,v 1.12 2006-12-27 07:21:27 roystgnr Exp $

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



#ifndef __inf_quad6_h__
#define __inf_quad6_h__


#include "libmesh_config.h"
#ifdef ENABLE_INFINITE_ELEMENTS


// C++ includes


// Local includes
#include "face_inf_quad.h"



/**
 * The \p INFQUAD6 is an infinite element in 2D composed of 6 nodes.
 * It is numbered like this:
   \verbatim
             2     5     3
   INFQUAD6: o     o     o   closer to infinity
             |           |
             |           |
             |           |  
             |           |
             |           |
             o-----o-----o   base side
             0     4     1
   \endverbatim
 */

// ------------------------------------------------------------
// InfQuad6 class definition
class InfQuad6 : public InfQuad
{
public:

  /**
   * Constructor.  By default this element has no parent.
   */
  InfQuad6 (Elem* p=NULL):
    InfQuad(InfQuad6::n_nodes(), p) {}

  /**
   * Constructor.  Explicitly specifies the number of
   * nodes and neighbors for which storage will be allocated.
   */
  InfQuad6 (const unsigned int nn,
	    const unsigned int ns,
	    Elem* p) :
    InfQuad(nn, ns, p) {}
  
  /**
   * @returns 6
   */
  unsigned int n_nodes() const { return 6; }
  
  /**
   * @returns \p INFQUAD6
   */
  ElemType type () const { return INFQUAD6; }
  
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
   * specified edge (== is_node_on_side in 2D)
   */
  virtual bool is_node_on_edge(const unsigned int n,
			       const unsigned int e) const
  { return this->is_node_on_side(n,e); }
  
  
  /**
   * @returns \p SECOND
   */
  Order default_order() const { return SECOND; }
  
  /**
   * Creates and returns an \p Edge3 for the base (0) side, and an \p InfEdge2 for
   * the sides 1, 2.
   */
  AutoPtr<Elem> build_side (const unsigned int i) const;

  virtual void connectivity(const unsigned int sf,
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
   * Note that \p n is counted as depicted above, \f$ 4 \le n < 6 \f$.
   */
  unsigned short int second_order_adjacent_vertex (const unsigned int n,
						   const unsigned int v) const;
  
  /**
   * @returns the child number \p c and element-local index \p v of the
   * \f$ n^{th} \f$ second-order node on the parent element.  Note that
   * the return values are always less \p this->n_children() and 
   * \p this->child(c)->n_vertices(), while \p n has to be greater or equal
   * to \p * this->n_vertices().  For linear elements this returns 0,0.
   * On refined second order elements, the return value will satisfy
   * \p this->get_node(n)==this->child(c)->get_node(v)
   */
  virtual std::pair<unsigned short int, unsigned short int> 
	  second_order_child_vertex (const unsigned int n) const;

  /**
   * This maps the \f$ j^{th} \f$ node of the \f$ i^{th} \f$ side to
   * element node numbers.
   */
  static const unsigned int side_nodes_map[3][3];

  
  
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
  static const float _embedding_matrix[2][6][6];
  
#endif


private:
  
  /**
   * Matrix that tells which vertices define the location
   * of mid-side (or second-order) nodes
   */
  static const unsigned short int _second_order_adjacent_vertices[2][2];

};


#endif // ifdef ENABLE_INFINITE_ELEMENTS

#endif
