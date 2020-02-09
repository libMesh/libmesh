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

#ifndef LIBMESH_RAW_VALUE_H
#define LIBMESH_RAW_VALUE_H

#include "libmesh/libmesh_common.h"

namespace libMesh
{
#ifdef LIBMESH_HAVE_METAPHYSICL
#include "metaphysicl/raw_type.h"

using MetaPhysicL::raw_value;
using typename MetaPhysicL::RawType;

#else

// RawType strips the derivatives, shadow magic, etc. off of
// types.  Vector types remain vector types.
template <typename T>
struct RawType
{
  typedef T value_type;

  static value_type value(const T& a) { return a; }
};

template <typename T>
struct RawType<const T>
{
  typedef const typename RawType<T>::value_type value_type;

  static value_type value(const T& a) { return RawType<T>::value(a); }
};


// Make the user syntax slightly nicer
template <typename T>
inline
typename RawType<T>::value_type
raw_value(const T& a) { return RawType<T>::value(a); }

#endif
}

#endif // LIBMESH_RAW_VALUE_H
