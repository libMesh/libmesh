// $Id: cell_hex27.h,v 1.14 2003-05-24 22:49:46 benkirk Exp $

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



#ifndef __cell_hex27_h__
#define __cell_hex27_h__

// C++ includes

// Local includes
#include "cell_hex.h"




/**
 * The \p Hex27 is an element in 3D composed of 27 nodes.
 * It is numbered like this:
   \verbatim
   HEX27:      7              18             6     			      
               o--------------o--------------o     			      
              /:             /              /|     			      
             / :            /              / |     			      
            /  :           /              /  |     			      
         19/   :        25/            17/   |     			      
          o--------------o--------------o    |     			      
         /     :        /              /|    |     			      
        /    15o       /    23o       / |  14o     
       /       :      /              /  |   /|     
     4/        :   16/             5/   |  / |     
     o--------------o--------------o    | /  |     			      
     |         :    |   26         |    |/   |     
     |  24o    :    |    o         |  22o    |     
     |         :    |       10     |   /|    |                                
     |        3o....|.........o....|../.|....o     
     |        .     |              | /  |   / 2      
     |       .    21|            13|/   |  /        
  12 o--------------o--------------o    | /         
     |     .        |              |    |/          
     |  11o         | 20o          |    o           
     |   .          |              |   / 9          
     |  .           |              |  /             
     | .            |              | /              
     |.             |              |/               
     o--------------o--------------o                
     0              8              1                
   \endverbatim
 */

// ------------------------------------------------------------
// Hex27 class definition
class Hex27 : public Hex
{
public:

  /**
   * Constructor.  By default this element has no parent.
   */
  Hex27  (const Elem* p=NULL);
  
  /**
   * @returns \p HEX27
   */
  ElemType     type ()   const { return HEX27; }

  /**
   * @returns 27
   */
  unsigned int n_nodes() const { return 27; }
  
  /**
   * @returns 8
   */
  unsigned int n_sub_elem() const { return 8; }
  
  /**
   * @returns SECOND
   */
  Order default_order() const { return SECOND; }

  /**
   * @returns an id associated with the \p s side of this element.
   * The id is not necessariy unique, but should be close.  This is
   * particularly useful in the \p MeshBase::find_neighbors() routine.
   *
   * We reimplemenet this method here for the \p Hex27 since we can
   * use the center node of each face to provide a perfect (unique)
   * key.
   */
  unsigned int key (const unsigned int s) const;
  
  /**
   * Builds a \p QUAD9 built coincident with face i.  
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
  static const float _embedding_matrix[8][27][27];
  
#endif

};



// ------------------------------------------------------------
// Hex27 class member functions
inline
Hex27::Hex27(const Elem* p) :
  Hex(Hex27::n_nodes(), p) 
{
}



#endif
