// $Id: cell_inf_hex8.h,v 1.8 2003-02-20 23:18:00 benkirk Exp $

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



#ifndef __cell_inf_hex8_h__
#define __cell_inf_hex8_h__

// C++ includes

// Local includes
#include "mesh_config.h"
#include "cell_hex.h"




#ifdef ENABLE_INFINITE_ELEMENTS

/**
 * The \p InfHex8 is an infinite element in 3D composed of 8 nodes.
 * It is numbered like this:
   \verbatim                
  INFHEX8: 7        6         z^  / y
           o--------o          | /
          /|       /|          |/ 
         / |      / |          +----> x
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
class InfHex8 : public Hex
{
public:

  /**
   * Constructor.  By default this element has no parent.
   */
  InfHex8  (Cell* p=NULL);
  
  /**
   * @returns \p INFHEX8
   */
  ElemType type () const { return INFHEX8; }

  /**
   * @returns 1
   */
  unsigned int n_sub_elem() const { return 1; }
  
  /**
   * @returns FIRST
   */
  Order default_order() const { return FIRST; }
  
  /**
   * Returns a QUAD4 built coincident with face 0, an INFQUAD4 
   * built coincident with faces 1 to 4.  Face 5 not supported. 
   * This method allocates memory, so be sure to delete
   * the returned pointer when it is no longer needed.
   */
  AutoPtr<Elem> build_side (const unsigned int i) const;

  const std::vector<unsigned int> tecplot_connectivity(const unsigned int sc=0) const;
  
  void vtk_connectivity(const unsigned int,
			std::vector<unsigned int>*) const
  { error(); }
  
  unsigned int vtk_element_type (const unsigned int) const
  { return 12; }
  
  void write_tecplot_connectivity(std::ostream &out) const;
  
#ifdef ENABLE_AMR

  /**
   * Refine the element.
   */
  void refine(Mesh& mesh);
  
#endif

  
private:
  
  
#ifdef ENABLE_AMR

  /**
   * Matrix that computes new nodal locations/solution values
   * from current nodes/solution.
   */
  static const float embedding_matrix[4][8][8];
  
  /**
   * Matrix that tells which children share which of
   * my sides. Note that infinite elements use different
   * storage scheme than conventional elements.
   */
  static const unsigned int side_children_matrix[6][5];
  
#endif

};



// ------------------------------------------------------------
// InfHex8 class member functions
inline
InfHex8::InfHex8(Cell* p) :
  Hex(InfHex8::n_nodes(), p) 
{
}



#endif

#endif
