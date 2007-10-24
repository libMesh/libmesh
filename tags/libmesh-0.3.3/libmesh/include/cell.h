// $Id: cell.h,v 1.8 2003-02-28 23:37:36 benkirk Exp $

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



#ifndef __cell_h__
#define __cell_h__

// C++ includes

// Local includes
#include "elem.h"


// Forward declarations



/**
 * The \p Cell is an abstract element type that lives in
 * three dimensions.  A face could be a tetrahedron, a hexahedron,
 * a pyramid, a prism, etc...
 */

// ------------------------------------------------------------
// Cell class definition
class Cell : public Elem
{
public:

  /**
   * Constructor.
   */
  Cell (const unsigned int nn,
	const unsigned int ns,
	const Elem* p) :
    Elem (nn, ns, p)
  {}

  /**
   * @returns 3, the dimensionality of the object.
   */
  unsigned int dim () const { return 3; }
  
  /**
   * @returns 4
   */
  unsigned int n_children_per_side(const unsigned int) const { return 4; }

};



#endif
