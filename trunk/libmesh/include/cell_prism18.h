// $Id: cell_prism18.h,v 1.1 2003-02-25 18:34:41 benkirk Exp $

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



#ifndef __cell_prism18_h__
#define __cell_prism18_h__

// C++ includes

// Local includes
#include "cell_prism.h"




/**
 * The \p Prism18 is an element in 3D composed of 18 nodes.
 * It is numbered like this:
   \verbatim
   PRISM18:
            5
            o
           /|\
          /   \
         /  |  \
        /       \
    14 o    |    o 13
      /           \ 
     /      |      \
    /       o 11    \
 3 /                 \4
  o---------o---------o
  |         |12       |
  |                   |
  |    o    |    o    |
  |   17    o   16    |
  |        /2         |
  |           \       |
9 o      /  o         o 10
  |        15   \     |
  |  8 o         o 7  |
  |   /               |
  |                \  |
  | /                 |
  |                  \|
  o---------o---------o
  0         6         1

  
   \endverbatim
 */

// ------------------------------------------------------------
// Prism class definition
class Prism18 : public Prism
{
public:

  /**
   * Constructor.  By default this element has no parent.
   */
  Prism18  (Cell* p=NULL);
  
  /**
   * @returns \p PRISM18
   */
  ElemType     type () const   { return PRISM18; }

  /**
   * @returns 18
   */
  unsigned int n_nodes() const { return 18; }

  /**
   * @returns 2
   */
  unsigned int n_sub_elem() const { return 8; }
  
  /**
   * @returns SECOND
   */
  Order default_order() const { return SECOND; }
  
  /**
   * Builds a QUAD4 or TRI3 built coincident with face i.  This
   * method allocates memory, so be sure to delete
   * the returned pointer when it is no longer needed.
   */
  AutoPtr<Elem> build_side (const unsigned int i) const;

  const std::vector<unsigned int> tecplot_connectivity(const unsigned int sc=0) const;
  
  
  void vtk_connectivity(const unsigned int sc,
			std::vector<unsigned int> *conn = NULL) const;
  
  unsigned int vtk_element_type (const unsigned int) const
  { return 13; }
  
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
//  static const float embedding_matrix[8][18][18];
  
  /**
   * Matrix that tells which children share which of
   * my sides.
   */
  static const unsigned int side_children_matrix[5][4];
  
#endif
  
};



// ------------------------------------------------------------
// Prism18 class member functions
inline
Prism18::Prism18(Cell* p) :
  Prism(Prism18::n_nodes(), p) 
{
}



#endif // #define __cell_prism18_h__

