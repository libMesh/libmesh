// $Id: cell_prism.h,v 1.11 2003-09-02 18:02:36 benkirk Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002-2003  Benjamin S. Kirk, John W. Peterson
  
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



#ifndef __cell_prism_h__
#define __cell_prism_h__

// C++ includes

// Local includes
#include "cell.h"




/**
 * The \p Prism is an element in 3D with 5 sides.
 */

// ------------------------------------------------------------
// Prism class definition
class Prism : public Cell
{
public:

  /**
   * Default prismatic element, takes number of nodes and 
   * parent. Derived classes implement 'true' elements.
   */
  Prism(const unsigned int nn, const Elem* p);
  
  /**
   * @returns 6.  All prism-derivatives are guaranteed to have at
   * least 6 nodes.
   */
  unsigned int n_nodes() const { return 6; }

  /**
   * @returns 5
   */
  unsigned int n_sides() const { return 5; }

  /**
   * @returns 6.  All prisms have 6 vertices.
   */
  unsigned int n_vertices() const { return 6; }

  /**
   * @returns 9.  All prisms have 9 edges.
   */
  unsigned int n_edges() const { return 9; }

  /**
   * @returns 5.  All prisms have 5 faces.
   */
  unsigned int n_faces() const { return 5; }
  
  /**
   * @returns 8
   */
  unsigned int n_children() const { return 8; }
  
  /**
   * @returns an id associated with the \p s side of this element.
   * The id is not necessariy unique, but should be close.  This is
   * particularly useful in the \p MeshBase::find_neighbors() routine.
   */
  unsigned int key (const unsigned int s) const;

  /**
   * @returns a primitive triangle or quad for 
   * face i.
   */
  AutoPtr<Elem> side (const unsigned int i) const;



protected:


  
#ifdef ENABLE_AMR
  
  /**
   * Matrix that allows children to inherit boundary conditions.
   */
  unsigned int side_children_matrix (const unsigned int i,
				     const unsigned int j) const
  { return _side_children_matrix[i][j]; }

#endif  

  
  /**
   * Matrix that tells which vertices define the location
   * of mid-side (or second-order) nodes.  for the first
   * 9 second-order nodes, this matrix is identical for
   * \p Prism15 and \p Prism18, therefore store it here.
   */
  static const unsigned short int _second_order_adjacent_vertices[9][2];


private:


  
#ifdef ENABLE_AMR
  
  /**
   * Matrix that tells which children share which of
   * my sides.
   */
  static const unsigned int _side_children_matrix[5][4];
  
#endif  


};



// ------------------------------------------------------------
// Prism class member functions
inline
Prism::Prism(const unsigned int nn, const Elem* p) :
  Cell(nn, Prism::n_sides(), p) 
{
}

#endif
