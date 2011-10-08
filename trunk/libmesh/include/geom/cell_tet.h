// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
  
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



#ifndef __cell_tet_h__
#define __cell_tet_h__

// C++ includes

// Local includes
#include "cell.h"

namespace libMesh
{




/**
 * The \p Tet is an element in 3D composed of 4 sides.
 */

// ------------------------------------------------------------
// Tet class definition
class Tet : public Cell
{
public:

  /**
   * Default tetrahedral element, takes number of nodes and 
   * parent. Derived classes implement 'true' elements.
   */
  Tet  (const unsigned int nn, Elem* p, Node** nodelinkdata);

  /**
   * @returns 4
   */
  unsigned int n_sides() const { return 4; }

  /**
   * @returns 4.  All tetrahedrals have 4 vertices.
   */
  unsigned int n_vertices() const { return 4; }

  /**
   * @returns 6.  All tetrahedrals have 6 edges.
   */
  unsigned int n_edges() const { return 6; }

  /**
   * @returns 4.  All tetrahedrals have 4 faces.
   */
  unsigned int n_faces() const { return 4; }

  /**
   * @returns 8
   */
  unsigned int n_children() const { return 8; } 
  
  /*
   * @returns true iff the specified child is on the
   * specified side
   */
  virtual bool is_child_on_side(const unsigned int c,
			        const unsigned int s) const;
  
  /**
   * @returns an id associated with the \p s side of this element.
   * The id is not necessariy unique, but should be close.  This is
   * particularly useful in the \p MeshBase::find_neighbors() routine.
   */
  unsigned int key (const unsigned int s) const;

  /**
   * @returns a primitive (3-noded) triangle for 
   * face i.
   */
  AutoPtr<Elem> side (const unsigned int i) const;
  
  /**
   * Based on the quality metric q specified by the user,
   * returns a quantitative assessment of element quality.
   */
  Real quality (const ElemQuality q) const;

  /**
   * Returns the suggested quality bounds for
   * the hex based on quality measure q.  These are
   * the values suggested by the CUBIT User's Manual.
   */
  std::pair<Real, Real> qual_bounds (const ElemQuality q) const;
  

protected:

  /**
   * Data for links to parent/neighbor/interior_parent elements.
   */
  Elem* _elemlinks_data[5+(LIBMESH_DIM>3)];
};



// ------------------------------------------------------------
// Tet class member functions
inline
Tet::Tet(const unsigned int nn, Elem* p, Node** nodelinkdata) :
  Cell(nn, Tet::n_sides(), p, _elemlinks_data, nodelinkdata) 
{
}

} // namespace libMesh

#endif
