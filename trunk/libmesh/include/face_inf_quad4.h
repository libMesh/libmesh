// $Id: face_inf_quad4.h,v 1.3 2003-01-20 17:06:09 jwpeterson Exp $

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


#ifndef __inf_quad4_h__
#define __inf_quad4_h__




#include "mesh_common.h"

#ifdef ENABLE_INFINITE_ELEMENTS



// C++ includes


// Local includes
#include "face_quad.h"







/**
 * The \p INFQUAD4 is an infinite element in 2D composed of 4 nodes.
 * It is numbered like this:
 * \verbatim
 *           3           2
 * INFQUAD4: o-----------o
 *           |           |
 *           |           |
 *           |           |  
 *           |           |
 *           |           |
 *           o-----------o
 *           0           1
 * \endverbatim
 */
// ------------------------------------------------------------
// InfQuad4 class definition
class InfQuad4 : public Quad
{
public:

  /**
   * Constructor.  By default this element has no parent.
   */
  InfQuad4  (Face* p=NULL);
  
  /**
   * @returns \p INFQUAD4 
   */
  ElemType type () const { return INFQUAD4; };
  
  /**
   * @returns 4
   */
  unsigned int n_sides() const { return 4; };
  
  /**
   * @returns 2
   */
  unsigned int n_children() const { return 2; };
  
  /**
   * @returns 1
   */
  unsigned int n_sub_elem() const { return 1; };
  
  /**
   * @returns FIRST
   */
  Order default_order() const { return FIRST; };
  
  /**
   * @returns an EDGE2 for the base (0) side, an INFEDGE2 for
   * the sides 1, 3. Side 2 no supported.
   */
  AutoPtr<Elem> build_side (const unsigned int i) const;

  const std::vector<unsigned int> tecplot_connectivity(const unsigned int sf=0) const;
  
  void vtk_connectivity(const unsigned int,
			std::vector<unsigned int>*) const
  { error(); };
  
  unsigned int vtk_element_type (const unsigned int) const
  { return 9; };
  
  void write_tecplot_connectivity(std::ostream&) const
  { error(); };
    
  
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
  static const real embedding_matrix[2][4][4];
  
  /**
   * Matrix that tells which children share which of
   * my sides.
   */
  static const unsigned int side_children_matrix[4][3];
  
#endif
    
};









// ------------------------------------------------------------
// InfQuad4 class member functions
inline
InfQuad4::InfQuad4(Face* p) :
  Quad(InfQuad4::n_nodes(), p) 
{
};


#endif

#endif
