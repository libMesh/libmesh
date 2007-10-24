// $Id: face_tri3.h,v 1.1 2003-11-05 22:26:43 benkirk Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002-2003  Benjamin S. Kirk, John W. Peterson
  
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



#ifndef __tri3_h__
#define __tri3_h__


// C++ includes


// Local includes
#include "libmesh_common.h"
#include "face_tri.h"


// Forward declarations

/**
 * The \p Tri3 is an element in 2D composed of 3 nodes.
 * It is numbered like this:
 * \verbatim
 *   TRI3:  2      
 *          o      
 *         / \     
 *        /   \    
 *       /     \   
 *      /       \  
 *     /         \ 
 *    o-----------o
 *    0           1
 * \endverbatim
 */

// ------------------------------------------------------------
// Tri3 class definition
class Tri3 : public Tri
{
public:

  /**
   * Constructor.  By default this element has no parent.
   */
  Tri3  (const Elem* p=NULL);

  /**
   * @returns \p TRI3
   */
  ElemType type () const { return TRI3; }

  /**
   * @returns 1
   */
  unsigned int n_sub_elem() const { return 1; }
  
  /**
   * @returns FIRST
   */
  Order default_order() const { return FIRST; }

  AutoPtr<Elem> build_side (const unsigned int i) const;
  
  const std::vector<unsigned int> tecplot_connectivity(const unsigned int sf=0) const;
  
  
  void vtk_connectivity(const unsigned int sc,
			std::vector<unsigned int> *conn = NULL) const;

  unsigned int vtk_element_type (const unsigned int) const
  { return 5; }

  
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
  static const float _embedding_matrix[4][3][3];
  
#endif
    
};



// ------------------------------------------------------------
// Tri3 class member functions
inline
Tri3::Tri3(const Elem* p) :
  Tri(Tri3::n_nodes(), p) 
{
}

#endif
