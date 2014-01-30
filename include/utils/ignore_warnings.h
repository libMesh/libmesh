// The libMesh Finite Element Library.
// Copyright (C) 2002-2014 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#ifndef NDEBUG
// TODO: icpc options
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-parameter"
#endif

#if defined(__GNUC__) && !defined(__INTEL_COMPILER) && !defined(__clang__)
#if __GNUC__ > 3 && __GNUC_MINOR__ > 1
// These two don't work?
#pragma GCC diagnostic ignored "-pedantic"
#pragma GCC diagnostic ignored "-Wdeprecated"
// But this is helpful with some MPI stacks
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif //__GNUC__ > 3 && __GNUC_MINOR__ > 1
#endif // __GNUC__ && !__INTEL_COMPILER
#endif // NDEBUG
