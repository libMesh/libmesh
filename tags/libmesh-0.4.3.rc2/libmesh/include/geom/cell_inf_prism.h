// $Id: cell_inf_prism.h,v 1.3 2004-10-25 21:49:24 benkirk Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2004  Benjamin S. Kirk, John W. Peterson
  
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



#ifndef __cell_inf_prism_h__
#define __cell_inf_prism_h__

// C++ includes

// Local includes
#include "libmesh_config.h"
#ifdef ENABLE_INFINITE_ELEMENTS

#include "cell_inf.h"




/**
 * The \p InfPrism is an element in 3D with 4 sides.
 * The \f$ 5^{th} \f$ side is theoretically located at infinity, 
 * and therefore not accounted for.
 * However, one could say that the \f$ 5^{th} \f$ side actually 
 * @e does exist in the mesh, since the outer nodes are located 
 * at a specific distance from the mesh origin (and therefore
 * define a side).  Still, this face is not to be used!
 */

// ------------------------------------------------------------
// InfPrism class definition
class InfPrism : public InfCell
{
public:

  /**
   * Default infinite prism element, takes number of nodes and 
   * parent. Derived classes implement 'true' elements.
   */
  InfPrism(const unsigned int nn, const Elem* p);

  /**
   * @returns 4 for the base \p s=0 and 2 for side faces. 
   */
  unsigned int n_children_per_side(const unsigned int s) const;

  /**
   * @returns 4.  Infinite elements have one side less
   * than their conventional counterparts, since one
   * side is supposed to be located at infinity.
   */
  unsigned int n_sides() const { return 4; }

  /**
   * @returns 6.  All infinite prisms (in our 
   * setting) have 6 vertices.
   */
  unsigned int n_vertices() const { return 6; }

  /**
   * @returns 6.  All infinite prismahedrals have 6 edges,
   * 3 lying in the base, and 3 perpendicular to the base.
   */
  unsigned int n_edges() const { return 6; }

  /**
   * @returns 4.  All prisms have 4 faces.
   */
  unsigned int n_faces() const { return 4; }
  
  /**
   * @returns 4
   */
  unsigned int n_children() const { return 4; }
  
  /**
   * @returns an id associated with the \p s side of this element.
   * The id is not necessariy unique, but should be close.  This is
   * particularly useful in the \p MeshBase::find_neighbors() routine.
   */
  unsigned int key (const unsigned int s) const;

  /**
   * @returns a primitive (3-noded) tri or (4-noded) infquad for 
   * face i.
   */
  AutoPtr<Elem> side (const unsigned int i) const;

};



// ------------------------------------------------------------
// InfPrism class member functions
inline
InfPrism::InfPrism(const unsigned int nn, const Elem* p) :
  InfCell(nn, InfPrism::n_sides(), p) 
{
}


inline
unsigned int InfPrism::n_children_per_side(const unsigned int s) const
{
  assert (s < this->n_sides());

  switch (s)
  {
    case 0:
      // every infinite prism has 4 children in the base side
      return 4;

    default:
      // on infinite faces (sides), only 2 children exist
      //
      // note that the face at infinity is already caught by the assertion
      return 2;
  }
}


#endif // ifdef ENABLE_INFINITE_ELEMENTS

#endif
