// $Id: enum_elem_type.h,v 1.1 2003-01-20 16:31:22 jwpeterson Exp $

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



#ifndef __enum_elem_type_h__
#define __enum_elem_type_h__

// C++ includes

// Local includes
#include "mesh_config.h"




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



#endif




