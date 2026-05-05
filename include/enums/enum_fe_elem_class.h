// The libMesh Finite Element Library.
// Copyright (C) 2002-2026 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_ENUM_FE_ELEM_CLASS_H
#define LIBMESH_ENUM_FE_ELEM_CLASS_H

namespace libMesh {

/**
 * \enum libMesh::FEElemClass groups element types by topological class,
 * independent of polynomial order.
 *
 * e.g. QUAD4, QUAD8, QUAD9 all map to QUAD; TRI3, TRI6, TRI7 all map to TRI.
 * Used together with FEFamily and polynomial order to uniquely identify a
 * physics finite element space.
 *
 * The fixed type allows forward declaration as:
 * enum class FEElemClass : unsigned int;
 */
enum class FEElemClass : unsigned int
{
  EDGE    = 0,
  TRI     = 1,
  QUAD    = 2,
  TET     = 3,
  HEX     = 4,
  PRISM   = 5,
  PYRAMID = 6,
  N_CLASSES
};

} // namespace libMesh

#endif // LIBMESH_ENUM_FE_ELEM_CLASS_H
