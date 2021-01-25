// The libMesh Finite Element Library.
// Copyright (C) 2002-2021 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_LIBMESH_AUGMENT_STD_NAMESPACE_H
#define LIBMESH_LIBMESH_AUGMENT_STD_NAMESPACE_H

#include "libmesh_common.h"


// For some reason the real std::max, std::min
// don't handle mixed compatible types
#define LIBMESH_DEFINE_MAX_MIN(type1, type2) \
inline auto max(type1 a, type2 b) -> decltype(a>b?a:b) \
{ return (a>b?a:b); } \
inline auto min(type1 a, type2 b) -> decltype(a<b?a:b) \
{ return (a<b?a:b); }

#define LIBMESH_DEFINE_BOTH_MAX_MIN(type1, type2) \
LIBMESH_DEFINE_MAX_MIN(type1, type2) \
LIBMESH_DEFINE_MAX_MIN(type2, type1)

namespace std {
LIBMESH_DEFINE_BOTH_MAX_MIN(long double, double)
LIBMESH_DEFINE_BOTH_MAX_MIN(long double, float)
LIBMESH_DEFINE_BOTH_MAX_MIN(double, float)

#ifdef LIBMESH_DEFAULT_QUADRUPLE_PRECISION
LIBMESH_DEFINE_BOTH_MAX_MIN(libMesh::Real, long double)
LIBMESH_DEFINE_BOTH_MAX_MIN(libMesh::Real, double)
LIBMESH_DEFINE_BOTH_MAX_MIN(libMesh::Real, float)
#endif

// fix for std::abs() overload ambiguity
#if defined (__SUNPRO_CC) || defined(__PGI)
inline double abs(double a)
{ return ::fabs(a); }

#endif

// fix for std::pow() overload ambiguity
#if defined (__SUNPRO_CC)
inline double pow(double a, int b)
{ return std::pow(a, static_cast<double>(b)); }
#endif
}

#endif // #define LIBMESH_LIBMESH_AUGMENT_STD_NAMESPACE_H
