// $Id: face_quad9.h,v 1.13 2003-08-07 19:25:31 ddreyer Exp $

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



#ifndef __quad9_h__
#define __quad9_h__


// C++ includes

// Local includes
#include "mesh_common.h"
#include "face_quad.h"



// Forward declarations



/**
 * The \p QUAD9 is an element in 2D composed of 9 nodes.
 * It is numbered like this:
 * \verbatim
 *        3     6     2
 * QUAD9: o-----o-----o
 *        |           |
 *        |     8     |
 *      7 o     o     o 5
 *        |           |
 *        |           |
 *        o-----o-----o
 *        0     4     1
 * \endverbatim
 */

// ------------------------------------------------------------
// Quad9 class definition
class Quad9 : public Quad
{
public:

  /**
   * Constructor.  By default this element has no parent.
   */
  Quad9  (const Elem* p=NULL);

  /**
   * @returns \p QUAD9
   */
  ElemType type () const { return QUAD9; }

  /**
   * @returns 9
   */
  unsigned int n_nodes() const { return 9; }
  
  /**
   * @returns 4
   */
  unsigned int n_sub_elem() const { return 4; }
  
  /**
   * @returns SECOND
   */
  Order default_order() const { return SECOND; }

  /**
   * @returns an id associated with the \p s side of this element.
   * The id is not necessariy unique, but should be close.  This is
   * particularly useful in the \p MeshBase::find_neighbors() routine.
   *
   * We reimplemenet this method here for the \p Quad9 since we can
   * use the center node of each edge to provide a perfect (unique)
   * key.
   */
  unsigned int key (const unsigned int s) const;
  
  AutoPtr<Elem> build_side (const unsigned int i) const;

  const std::vector<unsigned int> tecplot_connectivity(const unsigned int sf=0) const;
  
  
  void vtk_connectivity(const unsigned int sc,
			std::vector<unsigned int> *conn = NULL) const;

  unsigned int vtk_element_type (const unsigned int) const
  { return 9; }
  
  /**
   * @returns 2 for edge nodes and 4 for the face node.
   */
  unsigned int n_second_order_adjacent_vertices (const unsigned int n) const;

  /**
   * @returns the element-local number of the  \f$ v^{th} \f$ vertex
   * that defines the \f$ n^{th} \f$ second-order node.
   * Note that \p n is counted as depicted above, \f$ 4 \le n < 8 \f$.
   */
  unsigned int second_order_adjacent_vertex (const unsigned int n,
					     const unsigned int v) const;

  
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
  static const float _embedding_matrix[4][9][9];
  
#endif


private:
  
  /**
   * Matrix that tells which vertices define the location
   * of mid-side (or second-order) nodes
   */
  static const unsigned int _second_order_adjacent_vertices[5][4];
    
};


// ------------------------------------------------------------
// Quad9 class member functions
inline
Quad9::Quad9(const Elem* p) :
  Quad(Quad9::n_nodes(), p) 
{
}


#endif
