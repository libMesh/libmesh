// $Id: cell_pyramid.h,v 1.1.1.1 2003-01-10 16:17:48 libmesh Exp $

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



#ifndef __cell_pyramid_h__
#define __cell_pyramid_h__

// C++ includes

// Local includes
#include "cell.h"




/**
 * The \p Pyramid is an element in 3D with 5 sides.
 */

// ------------------------------------------------------------
// Pyramid class definition
class Pyramid : public Cell
{
public:

  /**
   * Default pyramid, one quad face, four triangular faces,
   * takes number of nodes and parent. 
   * Derived classes implement 'true' elements.
   */
  Pyramid(const unsigned int nn, Cell* p);

  /**
   * @returns 5.  All pyramid-derivatives are guaranteed to have at
   * least 5 nodes.
   */
  unsigned int n_nodes() const { return 5; };

  /**
   * @returns 5
   */
  unsigned int n_sides() const { return 5; };

  /**
   * @returns 5.  All pyramids have 5 vertices.
   */
  unsigned int n_vertices() const { return 5; };

  /**
   * @returns 8.  All pyramids have 8 edges.
   */
  unsigned int n_edges() const { return 8; };

  /**
   * @returns 5.  All pyramids have 5 faces.
   */
  unsigned int n_faces() const { return 5; };
  
  /**
   * @returns 10
   */
  unsigned int n_children() const { return 10; };

  /**
   * @returns a primitive triangle or quad for 
   * face i.
   */
  Elem side (const unsigned int i) const;
  
private:

};



// ------------------------------------------------------------
// Pyramid class member functions
inline
Pyramid::Pyramid(const unsigned int nn, Cell* p) :
  Cell(nn, Pyramid::n_sides(), p) 
{
};



#endif
