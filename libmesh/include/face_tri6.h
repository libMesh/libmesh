// $Id: face_tri6.h,v 1.2 2003-01-20 16:31:22 jwpeterson Exp $

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



#ifndef __tri6_h__
#define __tri6_h__


// C++ includes

// Local includes
#include "mesh_common.h"
#include "face_tri.h"



// Forward declarations




/**
 * The \p Tri6 is an element in 2D composed of 6 nodes.
 * It is numbered like this:
 * \verbatim
 *  TRI6:  2      
 *         o      
 *        / \     
 *       /   \    
 *    5 o     o 4 
 *     /       \  
 *    /         \ 
 *   o-----o-----o
 *   0     3     1
 * \endverbatim
 */

// ------------------------------------------------------------
// Tri6 class definition
class Tri6 : public Tri
{
public:

  /**
   * Constructor.  By default this element has no parent.
   */
  Tri6  (Face* p=NULL);

  /**
   * @returns \p TRI6
   */
  ElemType type ()   const { return TRI6; };

  /**
   * @returns 6
   */
  unsigned int n_nodes() const { return 6; };
  
  /**
   * @returns 4
   */
  unsigned int n_sub_elem() const { return 4; };
  
  /**
   * @returns SECOND
   */
  Order default_order() const { return SECOND; };
  
  AutoPtr<Elem> build_side (const unsigned int i) const;

  const std::vector<unsigned int> tecplot_connectivity(const unsigned int sf=0) const;
  
    
  void vtk_connectivity(const unsigned int sc,
			std::vector<unsigned int> *conn = NULL) const;

  unsigned int vtk_element_type (const unsigned int) const
  { return 6; };
  
  
#ifdef ENABLE_AMR

  /**
   * Refine the element.
   */
  void refine(Mesh& mesh);
  
  /**
   * Refine the element.
   */
  void coarsen();

#endif
  
private:
  
  
#ifdef ENABLE_AMR
  
  /**
   * Matrix that computes new nodal locations/solution values
   * from current nodes/solution.
   */
  static const real embedding_matrix[4][6][6];
  
  /**
   * Matrix that tells which children share which of
   * my sides.
   */
  static const unsigned int side_children_matrix[3][2];
  
#endif
    
};




// ------------------------------------------------------------
// Tri6 class member functions
inline
Tri6::Tri6(Face* p) :
  Tri(Tri6::n_nodes(), p) 
{
};


#endif
