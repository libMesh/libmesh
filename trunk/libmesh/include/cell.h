// $Id: cell.h,v 1.1.1.1 2003-01-10 16:17:48 libmesh Exp $

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
#include "mesh_common.h"
#include "elem.h"


// Forward declarations
class Mesh;



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
	Cell* p) :
    Elem (nn, ns, p)
  {};

  /**
   * @returns 3, the dimensionality of the object.
   */
  unsigned int dim () const { return 3; };
};



#endif
