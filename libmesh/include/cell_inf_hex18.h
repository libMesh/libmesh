// $Id: cell_inf_hex18.h,v 1.13 2003-04-01 14:19:45 ddreyer Exp $

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



#ifndef __cell_inf_hex18_h__
#define __cell_inf_hex18_h__

// C++ includes

// Local includes
#include "mesh_config.h"
#ifdef ENABLE_INFINITE_ELEMENTS

#include "cell_inf_hex.h"




/**
 * The \p InfHex18 is an infinite element in 3D composed of 18 nodes.
 * It is numbered like this:
   \verbatim
   INFHEX18:   7              14             6     			      
               o              o              o     closer to infinity
               :              :              |     	
               :              :              |     	
               :              :              |     	 
         15    :        17    :        13    |         
          o    :         o    :         o    |     	
          :    :              :         |    |     	
          :    :              :         |    |     
          :    :              :         |    |     
     4    :    :   12         :    5    |    |     
     o    :    :    o         :    o    |    |         
     |    :    :    |         :    |    |    |     
     |    :    :    |         :    |    |    |     
     |    :    :    |       10:    |    |    |           
     |    :   3o....|.........o....|....|....o    
     |    :   .     |              |    |   / 2      
     |    :  .      |              |    |  /        
     |    : .       |              |    | /         
     |    :.        |              |    |/          
     |  11o         |  16o         |    o           base face
     |   .          |              |   / 9          
     |  .           |              |  /             
     | .            |              | /              
     |.             |              |/               
     o--------------o--------------o                
     0              8              1                
   \endverbatim
 */

// ------------------------------------------------------------
// InfHex18 class definition
class InfHex18 : public InfHex
{
public:

  /**
   * Constructor.  By default this element has no parent.
   */
  InfHex18  (const Elem* p=NULL);
    
  /**
   * @returns 18.  The \p InfHex18 has 18 nodes.
   */
  unsigned int n_nodes() const { return 18; }
  
  /**
   * @returns \p INFHEX18
   */
  ElemType     type ()   const { return INFHEX18; }
  
  /**
   * @returns 4
   */
  unsigned int n_sub_elem() const { return 4; }
  
  /**
   * @returns SECOND
   */
  Order default_order() const { return SECOND; }

  /**
   * Returns a \p QUAD9 built coincident with face 0, an \p INFQUAD6 
   * built coincident with faces 1 to 4. Note that the \p AutoPtr<Elem>
   * takes care of freeing memory.
   */
  AutoPtr<Elem> build_side (const unsigned int i) const;

  const std::vector<unsigned int> tecplot_connectivity(const unsigned int sc=0) const;

  void vtk_connectivity(const unsigned int,
			std::vector<unsigned int>*) const
  { error(); }
  
  unsigned int vtk_element_type (const unsigned int) const
  { return 12; }
  
  

protected:
  
  
#ifdef ENABLE_AMR

  /**
   * Matrix used to create the elements children.
   */
  Real embedding_matrix (const unsigned int i,
			 const unsigned int j,
			 const unsigned int k) const
  { return static_cast<Real>(_embedding_matrix[i][j][k]); }

  /**
   * Matrix that computes new nodal locations/solution values
   * from current nodes/solution.
   */
  static const float _embedding_matrix[4][18][18];
  
  
#endif

};



// ------------------------------------------------------------
// InfHex18 class member functions
inline
InfHex18::InfHex18(const Elem* p) :
  InfHex(InfHex18::n_nodes(), p) 
{
}



#endif  // ifdef ENABLE_INFINITE_ELEMENTS


#endif
