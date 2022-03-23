// The libMesh Finite Element Library.
// Copyright (C) 2002-2022 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#ifndef LIBMESH_AUTO_PTR_H
#define LIBMESH_AUTO_PTR_H

// libMesh includes
#include "libmesh/libmesh_config.h"

// C++ includes
#include <memory>

// We still provide the UniquePtr alias for backwards compatibility,
// but it should now be considered deprecated and users should just
// use std::unique_ptr instead.
#ifdef LIBMESH_ENABLE_DEPRECATED

namespace libMesh
{
template<typename T>
using UniquePtr = std::unique_ptr<T>;
}

// Set up the libmesh_make_unique macro, for backwards compatibility.
// We require C++17 now; users should just use std::make_unique
// directly.

#define libmesh_make_unique std::make_unique

#endif // LIBMESH_ENABLE_DEPRECATED

#endif // LIBMESH_AUTO_PTR_H
