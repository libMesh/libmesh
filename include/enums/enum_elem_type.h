// The libMesh Finite Element Library.
// Copyright (C) 2002-2018 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

namespace libMesh {

/**
 * Defines an \p enum for geometric element types.
 *
 * The fixed type, i.e. ": int", enumeration syntax used here allows
 * this enum to be forward declared as
 * enum ElemType : int;
 * reducing header file dependencies.
 */
enum ElemType : int {
               // 1D
               EDGE2 = 0,
               EDGE3 = 1,
               EDGE4 = 2,
               // 2D
               TRI3 = 3,
               TRI6 = 4,
               QUAD4 = 5,
               QUAD8 = 6,
               QUAD9 = 7,
               // 3D
               TET4 = 8,
               TET10 = 9,
               HEX8 = 10,
               HEX20 = 11,
               HEX27 = 12,
               PRISM6 = 13,
               PRISM15 = 14,
               PRISM18 = 15,
               PYRAMID5 = 16,
               PYRAMID13 = 17,
               PYRAMID14 = 18,
               // Infinite Elems
               INFEDGE2 = 19,
               INFQUAD4 = 20,
               INFQUAD6 = 21,
               INFHEX8 = 22,
               INFHEX16 = 23,
               INFHEX18 = 24,
               INFPRISM6 = 25,
               INFPRISM12 = 26,
               // 0D
               NODEELEM = 27,
               // Miscellaneous Elems
               REMOTEELEM = 28,
               TRI3SUBDIVISION = 29,
               // Shell Elems
               TRISHELL3 = 30,
               QUADSHELL4 = 31,
               QUADSHELL8 = 32,
               // Invalid
               INVALID_ELEM};   // should always be last
}

#endif
