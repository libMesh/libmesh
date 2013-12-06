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



#ifndef LIBMESH_ENUM_ELEM_TYPE_H
#define LIBMESH_ENUM_ELEM_TYPE_H

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
                 PYRAMID14,  // 17

		 INFEDGE2,   // 18

		 INFQUAD4,   // 19
		 INFQUAD6,   // 20

		 INFHEX8,    // 21
		 INFHEX16,   // 22
		 INFHEX18,   // 23

		 INFPRISM6,  // 24
		 INFPRISM12, // 25

		 NODEELEM,   // 26

		 REMOTEELEM,   // 27

		 INVALID_ELEM};  // 28 - should always be last
}

using namespace libMeshEnums;

#endif // LIBMESH_ENUM_ELEM_TYPE_H
