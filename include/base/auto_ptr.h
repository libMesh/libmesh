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

#ifndef LIBMESH_AUTO_PTR_H
#define LIBMESH_AUTO_PTR_H

// libMesh includes
#include "libmesh/libmesh_config.h"

// C++ includes
#include <memory>

// We still provide the UniquePtr alias for backwards compatibility,
// but it should now be considered deprecated and users should just
// use std::unique_ptr instead.
namespace libMesh
{
template<typename T>
using UniquePtr = std::unique_ptr<T>;
}

// Set up the libmesh_make_unique macro. We don't yet require C++14, so this is a C++11 based workaround.
#ifdef LIBMESH_HAVE_CXX14_MAKE_UNIQUE

#define libmesh_make_unique std::make_unique

#else
namespace libMesh
{
template<typename T, typename... Args>
std::unique_ptr<T> make_unique(Args &&... args)
{
  return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}
}
#define libmesh_make_unique make_unique

#endif

#endif // LIBMESH_AUTO_PTR_H
