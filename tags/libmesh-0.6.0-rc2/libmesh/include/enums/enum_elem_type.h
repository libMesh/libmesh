// $Id: enum_elem_type.h,v 1.5 2006-10-04 22:26:53 roystgnr Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2005  Benjamin S. Kirk, John W. Peterson
  
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

/**
 * The \p libMeshEnums namespace is the namespace all \p enum definitions
 * should be put into.
 */

// ------------------------------------------------------------
// enum ElemType definition
namespace libMeshEnums {
  
  /**
   * Defines an \p enum for geometric element types.
   */
  enum ElemType {EDGE2=0,    // 0
		 EDGE3,      // 1
		 EDGE4,      // 2
		 
		 TRI3,       // 3
		 TRI6,       // 4
		 
		 QUAD4,      // 5
		 QUAD8,      // 6
		 QUAD9,      // 7
		 
		 TET4,       // 8
		 TET10,      // 9
		 
		 HEX8,       // 10
		 HEX20,      // 11
		 HEX27,      // 12
		 
		 PRISM6,     // 13
		 PRISM15,    // 14
		 PRISM18,    // 15
		 
		 PYRAMID5,   // 16
		 
		 INFEDGE2,   // 17
		 
		 INFQUAD4,   // 18
		 INFQUAD6,   // 19
		 
		 INFHEX8,    // 20
		 INFHEX16,   // 21
		 INFHEX18,   // 22
		 
		 INFPRISM6,  // 23
		 INFPRISM12, // 24

		 NODEELEM,   // 25
		 
		 INVALID_ELEM};
}

using namespace libMeshEnums;

#endif // #define __enum_elem_type_h__




