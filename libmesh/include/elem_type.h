// $Id: elem_type.h,v 1.1.1.1 2003-01-10 16:17:48 libmesh Exp $

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



#ifndef __elem_type_h__
#define __elem_type_h__

// C++ includes
#include <string>

// Local includes
#include "mesh_common.h"




/**
 * The \p MeshEnums namespace is the namespace all \p enum definitions
 * should be put into.
 */

// ------------------------------------------------------------
// enum ElemType definition
namespace MeshEnums {
  
  /**
   * Defines an \p enum for geometric element types.
   */
  enum ElemType {EDGE2=0,
		 EDGE3,
		 EDGE4,
		 
		 TRI3,
		 TRI6,
		 
		 QUAD4,
		 QUAD8,
		 QUAD9,
		 
		 TET4,
		 TET10,
		 
		 HEX8,
		 HEX20,
		 HEX27,
		 
		 PRISM6,
		 PRISM18,
		 
		 PYRAMID5,
		 
#ifdef ENABLE_INFINITE_ELEMENTS
		 INFEDGE2,
		 
		 INFQUAD4,
		 INFQUAD6,
		 
		 INFHEX8,
		 INFHEX16,
		 INFHEX18,
		 
		 INFPRISM6,
		 INFPRISM12,
#endif
		 
		 INVALID_ELEM};
};

using namespace MeshEnums;



// A namespace for element type utility
// functions.  Similar to the one used
// in elem_quality.h
namespace ElementTypes
{
  /**
   * The number of element types that are
   * defined (INVALD_ELEM excluded).  
   * You might have to update this
   * if you add a new one!
   */
#ifdef ENABLE_INFINITE_ELEMENTS
  const unsigned int num_types = 24;
#else
  const unsigned int num_types = 16;  
#endif

  /**
   * Returns a standard string representation
   * of the basic name for element type t.
   * For example, a HEX27 has the basic name
   * of "Hexahedron".
   */
  std::string basic_name (const ElemType t);

  /**
   * Returns a standard string representation
   * for the specific name of element type t.
   * For example, HEX27 returns "Hex 27".
   */
  std::string name (const ElemType t);
}

#endif




