// $Id: face_quad4.h,v 1.6 2005-02-19 19:07:31 roystgnr Exp $

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



#ifndef __quad4_h__
#define __quad4_h__


// C++ includes


// Local includes
#include "libmesh_common.h"
#include "face_quad.h"



// Forward declarations



/**
 * The \p QUAD4 is an element in 2D composed of 4 nodes.
 * It is numbered like this:
 * \verbatim
 *        3           2
 * QUAD4: o-----------o
 *        |           |
 *        |           |
 *        |           | 
 *        |           |
 *        |           |
 *        o-----------o
 *        0           1
 * \endverbatim
 */

// ------------------------------------------------------------
// Quad4 class definition
class Quad4 : public Quad
{
public:

  /**
   * Constructor.  By default this element has no parent.
   */
  Quad4 (const Elem* p=NULL) :
    Quad(Quad::n_nodes(), p) {}

  /**
   * Constructor.  Explicitly specifies the number of
   * nodes and neighbors for which storage will be allocated.
   */
  Quad4 (const unsigned int nn,
	 const unsigned int ns,
	 const Elem* p) :
    Quad(nn, ns, p) {}
  
  /**
   * @returns \p QUAD4
   */
  ElemType type () const { return QUAD4; }
  
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
   * @returns FIRST
   */
  Order default_order() const { return FIRST; }
  
  AutoPtr<Elem> build_side (const unsigned int i) const;

  virtual void connectivity(const unsigned int sf,
			    const IOPackage iop,
			    std::vector<unsigned int>& conn) const;

  /**
   * This maps the \f$ j^{th} \f$ node of the \f$ i^{th} \f$ side to
   * element node numbers.
   */
  static const unsigned int side_nodes_map[4][2];
  
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
  static const float _embedding_matrix[4][4][4];
  
#endif
    
};


#endif
