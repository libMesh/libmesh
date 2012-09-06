// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef __cell_inf_h__
#define __cell_inf_h__

// C++ includes

// Local includes
#include "libmesh_config.h"
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
#include "elem.h"

namespace libMesh
{




/**
 * The \p InfCell is an abstract element type that lives in
 * three dimensions.  An infinite cell could be an infinite hexahedron,
 * or an infinite prism.
 */

// ------------------------------------------------------------
// InfCell class definition
class InfCell : public Elem
{
public:

  /**
   * Constructor.
   */
  InfCell (const unsigned int nn,
	   const unsigned int ns,
	   Elem* p,
	   Elem** elemlinkdata,
           Node** nodelinkdata) :
    Elem (nn, ns, p, elemlinkdata, nodelinkdata)
  {}

  /**
   * @returns 3, the dimensionality of the object.
   */
  unsigned int dim () const { return 3; }

  /**
   * @returns \p true.  All classes derived from \p InfCell
   * are infinite elements.
   */
  bool infinite () const { return true; }

  /**
   * @returns the origin of this infinite element.
   */
  Point origin () const;

};


// ------------------------------------------------------------
// InfCell inline functions
inline
Point InfCell::origin () const
{
  return ( this->point(0)*2 - this->point( this->n_vertices()/2 ) );
}


} // namespace libMesh


#endif // ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

#endif
