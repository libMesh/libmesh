// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



// Local includes
#include "libmesh/quadrature_clough.h"
#include "libmesh/string_to_enum.h"

namespace libMesh
{

void QClough::init_3D(const ElemType, unsigned int)
{
#if LIBMESH_DIM == 3
  // QClough not supported on any 3D elements
  libmesh_error_msg("ERROR: Unsupported type: " << Utility::enum_to_string(_type));
#endif
}

} // namespace libMesh
