// $Id: face_quad4.h,v 1.6 2003-02-03 03:51:49 ddreyer Exp $

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



#ifndef __quad4_h__
#define __quad4_h__


// C++ includes


// Local includes
#include "mesh_common.h"
#include "face_quad.h"



// Forward declarations



/**
 * The \p QUAD4 is an element in 2D composed of 4 nodes.
 * It is numbered like this:
 * \verbatim
 *        3           2
 * QUAD4: o-----------o
 *        |           |
 *        |           |
 *        |           | 
 *        |           |
 *        |           |
 *        o-----------o
 *        0           1
 * \endverbatim
 */

// ------------------------------------------------------------
// Quad4 class definition
class Quad4 : public Quad
{
public:

  /**
   * Constructor.  By default this element has no parent.
   */
  Quad4 (Face* p=NULL);
  
  /**
   * @returns \p QUAD4
   */
  ElemType type () const { return QUAD4; };
  
  /**
   * @returns 1
   */
  unsigned int n_sub_elem() const { return 1; };
  
  /**
   * @returns FIRST
   */
  Order default_order() const { return FIRST; };
  
  AutoPtr<Elem> build_side (const unsigned int i) const;

  const std::vector<unsigned int> tecplot_connectivity(const unsigned int sf=0) const;
  
  
  void vtk_connectivity(const unsigned int sc,
			std::vector<unsigned int> *conn = NULL) const;

  unsigned int vtk_element_type (const unsigned int) const
  { return 9; };
  
  
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
  static const Real embedding_matrix[4][4][4];
  
  /**
   * Matrix that tells which children share which of
   * my sides.
   */
  static const unsigned int side_children_matrix[4][2];
  
#endif
    
};




// ------------------------------------------------------------
// Quad4 class member functions
inline
Quad4::Quad4(Face* p) :
  Quad(Quad::n_nodes(), p) 
{
};

#endif
