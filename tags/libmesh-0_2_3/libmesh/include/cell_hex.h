// $Id: cell_hex.h,v 1.6 2003-02-03 03:51:48 ddreyer Exp $

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



#ifndef __cell_hex_h__
#define __cell_hex_h__

// C++ includes

// Local includes
#include "cell.h"




/**
 * The \p Hex is an element in 3D with 6 sides.
 */

// ------------------------------------------------------------
// Hex class definition
class Hex : public Cell
{
public:

  /**
   * Default brick element, takes number of nodes and 
   * parent. Derived classes implement 'true' elements.
   */
  Hex(const unsigned int nn, Cell* p);

  /**
   * @returns 8.  All hex-derivatives are guaranteed to have at
   * least 8 nodes.
   */
  unsigned int n_nodes() const { return 8; };

  /**
   * @returns 6
   */
  unsigned int n_sides() const { return 6; };

  /**
   * @returns 8.  All hexahedrals have 8 vertices.
   */
  unsigned int n_vertices() const { return 8; };

  /**
   * @returns 12.  All hexahedrals have 12 edges.
   */
  unsigned int n_edges() const { return 12; };

  /**
   * @returns 6.  All hexahedrals have 6 faces.
   */
  unsigned int n_faces() const { return 6; };
  
  /**
   * @returns 8
   */
  unsigned int n_children() const { return 8; }

  /**
   * @returns a primitive (4-noded) quad for 
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

private:

};



// ------------------------------------------------------------
// Hex class member functions
inline
Hex::Hex(const unsigned int nn, Cell* p) :
  Cell(nn, Hex::n_sides(), p) 
{
};



#endif
