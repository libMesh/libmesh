// $Id: cell_hex27.h,v 1.1.1.1 2003-01-10 16:17:48 libmesh Exp $

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
              /|             /              /|     			      
             / |            /              / |     			      
            /  |           /              /  |     			      
         19/   |        25/            17/   |     			      
          o--------------o--------------o    |     			      
         /     |        /              /|    |     			      
        /    15o       /    23o       / |  14o     
       /       |      /              /  |   /|     
     4/        |   16/             5/   |  / |     
     o--------------o--------------o    | /  |     			      
     |         |    |   26         |    |/   |     
     |  24o    |    |    o         |  22o    |     
     |         |    |       10     |   /|    |                                
     |        3o----|---------o----|--/-|----o     
     |        /     |              | /  |   / 2      
     |       /    21|            13|/   |  /        
  12 o--------------o--------------o    | /         
     |     /        |              |    |/          
     |  11o         | 20o          |    o           
     |   /          |              |   / 9          
     |  /           |              |  /             
     | /            |              | /              
     |/             |              |/               
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
  Hex27  (Cell* p=NULL);
  
  /**
   * @returns \p HEX27
   */
  ElemType     type ()   const { return HEX27; };

  /**
   * @returns 27
   */
  unsigned int n_nodes() const { return 27; };
  
  /**
   * @returns 8
   */
  unsigned int n_sub_elem() const { return 8; };
  
  /**
   * @returns SECOND
   */
  Order default_order() const { return SECOND; };

  /**
   * Builds a QUAD9 built coincident with face i.  This
   * method allocates memory, so be sure to delete
   * the returned pointer when it is no longer needed.
   */
  std::auto_ptr<Elem> build_side (const unsigned int i) const;

  const std::vector<unsigned int> tecplot_connectivity(const unsigned int sc=0) const;
  
  
  void vtk_connectivity(const unsigned int sc,
			std::vector<unsigned int> *conn = NULL) const;
  
  unsigned int vtk_element_type (const unsigned int sc) const;
  
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
  static const real embedding_matrix[8][27][27];
  
  /**
   * Matrix that tells which children share which of
   * my sides.
   */
  static const unsigned int side_children_matrix[6][4];
  
#endif

};



// ------------------------------------------------------------
// Hex27 class member functions
inline
Hex27::Hex27(Cell* p) :
  Hex(Hex27::n_nodes(), p) 
{
};



#endif
