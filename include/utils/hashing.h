// The libMesh Finite Element Library.
// Copyright (C) 2002-2020 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#ifndef LIBMESH_HASHING_H
#define LIBMESH_HASHING_H

#include "libmesh/hashword.h"

#include <functional>

namespace libMesh
{
namespace boostcopy
{

inline void hash_combine_impl(std::size_t & seed, std::size_t value)
{
  seed ^= value + 0x9e3779b9 + (seed<<6) + (seed>>2);
}

template <typename T>
inline void hash_combine(std::size_t & seed, const T & value)
{
  hash_combine_impl(seed, std::hash<T>{}(value));
}

}
}

#endif // LIBMESH_HASHING_H
