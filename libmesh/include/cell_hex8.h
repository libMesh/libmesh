// $Id: cell_hex8.h,v 1.7 2003-02-13 22:56:06 benkirk Exp $

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



#ifndef __cell_hex8_h__
#define __cell_hex8_h__

// C++ includes

// Local includes
#include "cell_hex.h"




/**
 * The \p Hex8 is an element in 3D composed of 8 nodes.
 * It is numbered like this:
   \verbatim
  HEX8: 7        6           
        o--------o           
       /|       /|           
      / |      / |           
   4 /  |   5 /  |           
    o--------o   |           
    |   o----|-- o 2         
    |  /3    |  /            
    | /      | /             
    |/       |/              
    o--------o               
    0        1               
                             
   \endverbatim
 */

// ------------------------------------------------------------
// Hex class definition
class Hex8 : public Hex
{
public:

  /**
   * Constructor.  By default this element has no parent.
   */
  Hex8  (Cell* p=NULL);
  
  /**
   * @returns \p HEX8
   */
  ElemType type () const { return HEX8; }

  /**
   * @returns 1
   */
  unsigned int n_sub_elem() const { return 1; }
  
  /**
   * @returns FIRST
   */
  Order default_order() const { return FIRST; }
  
  /**
   * Builds a QUAD4 built coincident with face i.  This
   * method allocates memory, so be sure to delete
   * the returned pointer when it is no longer needed.
   */
  AutoPtr<Elem> build_side (const unsigned int i) const;

  const std::vector<unsigned int> tecplot_connectivity(const unsigned int sc=0) const;
    
  void vtk_connectivity(const unsigned int sc,
			std::vector<unsigned int> *conn = NULL) const;
  
  unsigned int vtk_element_type (const unsigned int) const
  { return 12; }
  
#ifdef ENABLE_AMR

  /**
   * Refine the element.
   */
  void refine(Mesh& mesh);
  
  /**
   * Coarsen the element.
   */
  void coarsen();
  
#endif

  
private:
  
  
#ifdef ENABLE_AMR

  /**
   * Matrix that computes new nodal locations/solution values
   * from current nodes/solution.
   */
  static const Real embedding_matrix[8][8][8];
  
  /**
   * Matrix that tells which children share which of
   * my sides.
   */
  static const unsigned int side_children_matrix[6][4];
  
#endif

};



// ------------------------------------------------------------
// Hex8 class member functions
inline
Hex8::Hex8(Cell* p) :
  Hex(Hex8::n_nodes(), p) 
{
}



#endif
