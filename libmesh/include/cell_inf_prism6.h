// $Id: cell_inf_prism6.h,v 1.14 2003-05-16 00:30:08 benkirk Exp $

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



#ifndef __cell_inf_prism6_h__
#define __cell_inf_prism6_h__

// C++ includes

// Local includes
#include "mesh_config.h"
#ifdef ENABLE_INFINITE_ELEMENTS

#include "cell_inf_prism.h"




/**
 * The \p InfPrism6 is an infinite element in 3D composed of 6 nodes.
 * It is numbered like this:
   \verbatim
   INFPRISM6:
           5
           o
           : 
           :         closer to infinity
           :
     3 o   :   o 4
       |   :   |
       | 2 o   |
       |  . .  |  
       | .   . |
       |.     .|
       o-------o     base face
       0       1
   \endverbatim
 */

// ------------------------------------------------------------
// InfPrism6 class definition
class InfPrism6 : public InfPrism
{
public:

  /**
   * Constructor.  By default this element has no parent.
   */
  InfPrism6  (const Elem* p=NULL);
    
  /**
   * @returns 6.  The \p InfPrism6 has 6 nodes.
   */
  unsigned int n_nodes() const { return 6; }
  
  /**
   * @returns \p INFPRISM6
   */
  ElemType     type() const { return INFPRISM6; }

  /**
   * @returns 1
   */
  unsigned int n_sub_elem() const { return 1; }
  
  /**
   * @returns FIRST
   */
  Order        default_order() const { return FIRST; }
  
  /**
   * Returns a \p TRI3 built coincident with face 0, an \p INFQUAD4 
   * built coincident with faces 1 to 3.  Note that the \p AutoPtr<Elem>
   * takes care of freeing memory.  Here, this method is the same as 
   * the \p InfPrism::side(i).
   */
  AutoPtr<Elem> build_side (const unsigned int i) const 
  { AutoPtr<Elem> ap(this->side(i)); return ap; }

  const std::vector<unsigned int> tecplot_connectivity(const unsigned int sc=0) const;
  
  void vtk_connectivity(const unsigned int,
			std::vector<unsigned int>*) const
  { error(); }
  
  unsigned int vtk_element_type (const unsigned int) const
  { return 13; }

  /**
   * @returns \p true when this element contains the point
   * \p p.  Customized for infinite elements, since knowledge
   * about the envelope can be helpful.
   */
  bool contains_point (const Point& p) const;
  
  
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
  static const float _embedding_matrix[4][6][6];
  
#endif
  
};



// ------------------------------------------------------------
// InfPrism6 class member functions
inline
InfPrism6::InfPrism6(const Elem* p) :
  InfPrism(InfPrism6::n_nodes(), p) 
{
}



#endif  // ifdef ENABLE_INFINITE_ELEMENTS

#endif
