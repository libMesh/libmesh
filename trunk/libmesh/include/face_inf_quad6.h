// $Id: face_inf_quad6.h,v 1.7 2003-02-13 22:56:07 benkirk Exp $

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



#ifndef __inf_quad6_h__
#define __inf_quad6_h__


#include "mesh_common.h"

#ifdef ENABLE_INFINITE_ELEMENTS


// C++ includes


// Local includes
#include "face_quad.h"



// Forward declarations





/**
 * The \p INFQUAD6 is an infinite element in 2D composed of 6 nodes.
 * It is numbered like this:
 * \verbatim
 *           3     5     2
 * INFQUAD6: o-----o-----o
 *           |           |
 *           |           |
 *           |           |  
 *           |           |
 *           |           |
 *           o-----o-----o
 *           0     4     1
 * \endverbatim
 */

// ------------------------------------------------------------
// InfQuad6 class definition
class InfQuad6 : public Quad
{
public:

  /**
   * Constructor.  By default this element has no parent.
   */
  InfQuad6  (Face* p=NULL);
  
  /**
   * @returns \p INFQUAD6
   */
  ElemType type () const { return INFQUAD6; }

  /**
   * @returns 6
   */
  unsigned int n_nodes() const { return 6; }
  
  /**
   * @returns 4
   */
  unsigned int n_sides() const { return 4; }
  
  /**
   * @returns 2
   */
  unsigned int n_children() const { return 2; }
  
  /**
   * @returns 2
   */
  unsigned int n_sub_elem() const { return 2; }
  
  /**
   * @returns SECOND
   */
  Order default_order() const { return SECOND; }
  
  /**
   * @returns an EDGE3 for the base (0) side, an INFEDGE2 for
   * the sides 1, 3. Side 2 no supported.
   */
  AutoPtr<Elem> build_side (const unsigned int i) const;

  const std::vector<unsigned int> tecplot_connectivity(const unsigned int sf=0) const;
  
  void vtk_connectivity(const unsigned int,
			std::vector<unsigned int>*) const
  { error(); }
  
  unsigned int vtk_element_type (const unsigned int) const
  { return 9; }
  
  void write_tecplot_connectivity(std::ostream&) const
  { error(); }
  
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
  static const Real embedding_matrix[2][6][6];
  
  /**
   * Matrix that tells which children share which of
   * my sides.
   */
  static const unsigned int side_children_matrix[4][3];
  
#endif
    
};





// ------------------------------------------------------------
// InfQuad6 class member functions
inline
InfQuad6::InfQuad6(Face* p) :
  Quad(InfQuad6::n_nodes(), p) 
{
}

#endif

#endif
