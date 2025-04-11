// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#include "libmesh/variable.h"
#include "libmesh/fe_interface.h"

namespace libMesh
{
unsigned int Variable::n_components() const
{
  if (_type.family == SCALAR)
    return _type.order.get_order();

  if (const auto fe_field_type = FEInterface::field_type(_type.family);
      fe_field_type == TYPE_VECTOR)
    libmesh_error_msg("Cannot determine the number of components. Please call the n_components "
                      "overload with a mesh argument");
  else
    return 1;
}

unsigned int Variable::n_components(const MeshBase & mesh) const
{
  if (_type.family == SCALAR)
    return _type.order.get_order();
  else
    return FEInterface::n_vec_dim(mesh, _type);
}
}
