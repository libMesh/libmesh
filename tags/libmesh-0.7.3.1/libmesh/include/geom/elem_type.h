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



#ifndef __elem_type_h__
#define __elem_type_h__

// C++ includes
#include <string>

// Local includes
#include "libmesh_common.h"
#include "enum_elem_type.h"

namespace libMesh
{



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
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
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

} // namespace libMesh

#endif




