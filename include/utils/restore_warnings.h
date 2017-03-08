// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


// Note: no include guards!  We want to be able to #include this
// header multiple times.

#include "libmesh/libmesh_config.h"

#ifdef __clang__
#pragma clang diagnostic pop
#endif

#if defined(__GNUC__) && !defined(__INTEL_COMPILER) && !defined(__clang__)
// GCC > 4.1 supports diagnostic pragmas
#if (__GNUC__ > 4) || (__GNUC__ == 4 && __GNUC_MINOR__ > 1)
// TODO: use the gcc 4.6 push/pop when available
#if (__GNUC__ > 5)
#pragma GCC diagnostic warning "-Wplacement-new"
#pragma GCC diagnostic warning "-Wmisleading-indentation"
#endif // GCC > 5
#pragma GCC diagnostic warning "-Wdeprecated-declarations"
#pragma GCC diagnostic warning "-Wunused-parameter"
#pragma GCC diagnostic warning "-Wdeprecated"
#pragma GCC diagnostic warning "-Wpedantic"
#endif // GCC > 4.1
#endif // __GNUC__ && !__INTEL_COMPILER
