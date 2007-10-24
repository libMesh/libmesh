// $Id: cell_prism6.h,v 1.13 2003-05-22 21:18:02 benkirk Exp $

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



#ifndef __cell_prism6_h__
#define __cell_prism6_h__

// C++ includes

// Local includes
#include "cell_prism.h"




/**
 * The \p Prism6 is an element in 3D composed of 6 nodes.
 * It is numbered like this:
   \verbatim
   PRISM6:
           5
           o
          /:\
         / : \
        /  o  \
     3 o-------o 4
       | . 2 . |
       |.     .|
       o-------o
       0       1
   \endverbatim
 */

// ------------------------------------------------------------
// Prism class definition
class Prism6 : public Prism
{
public:

  /**
   * Constructor.  By default this element has no parent.
   */
  Prism6  (const Elem* p=NULL);
  
  /**
   * @returns \p PRISM6
   */
  ElemType     type () const   { return PRISM6; }

  /**
   * @returns 1
   */
  unsigned int n_sub_elem() const { return 1; }
  
  /**
   * @returns FIRST
   */
  Order default_order() const { return FIRST; }
  
  /**
   * Builds a \p QUAD4 or \p TRI3 built coincident with face i.  
   * The \p AutoPtr<Elem> handles the memory aspect.
   */
  AutoPtr<Elem> build_side (const unsigned int i) const;

  const std::vector<unsigned int> tecplot_connectivity(const unsigned int sc=0) const;
  
  
  void vtk_connectivity(const unsigned int sc,
			std::vector<unsigned int> *conn = NULL) const;
  
  unsigned int vtk_element_type (const unsigned int) const
  { return 13; }
   
  
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
  static const float _embedding_matrix[8][6][6];
  
#endif
  
};



// ------------------------------------------------------------
// Prism6 class member functions
inline
Prism6::Prism6(const Elem* p) :
  Prism(Prism6::n_nodes(), p) 
{
}



#endif
