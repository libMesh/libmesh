// $Id: cell_hex8.h,v 1.15 2003-08-18 14:44:51 ddreyer Exp $

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
       /:       /|           
      / :      / |           
   4 /  :   5 /  |           
    o--------o   |           
    |   o....|...o 2         
    |  .3    |  /            
    | .      | /             
    |.       |/              
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
  Hex8  (const Elem* p=NULL);
  
  /**
   * @returns \p HEX8
   */
  ElemType type () const { return HEX8; }
  
  /**
   * @returns 8
   */
  unsigned int n_nodes() const { return 8; }

  /**
   * @returns 1
   */
  unsigned int n_sub_elem() const { return 1; }
  
  /**
   * @returns FIRST
   */
  Order default_order() const { return FIRST; }
  
  /**
   * Builds a QUAD4 built coincident with face i.  
   * The \p AutoPtr<Elem> handles the memory aspect.
   */
  AutoPtr<Elem> build_side (const unsigned int i) const;

  const std::vector<unsigned int> tecplot_connectivity(const unsigned int sc=0) const;
    
  void vtk_connectivity(const unsigned int sc,
			std::vector<unsigned int> *conn = NULL) const;
  
  unsigned int vtk_element_type (const unsigned int) const
  { return 12; }
  

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
  static const float _embedding_matrix[8][8][8];
  
#endif

};



// ------------------------------------------------------------
// Hex8 class member functions
inline
Hex8::Hex8(const Elem* p) :
  Hex(Hex8::n_nodes(), p) 
{
}



#endif
