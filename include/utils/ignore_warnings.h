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


// Note: no include guards!  We want to be able to #include this
// header multiple times.

#include "libmesh/libmesh_config.h"

// TODO: icpc options
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-parameter"
#pragma clang diagnostic ignored "-Wextra-semi"
#pragma clang diagnostic ignored "-Wvariadic-macros"
#pragma clang diagnostic ignored "-Wc++11-extensions"
#pragma clang diagnostic ignored "-Wmacro-redefined"
#pragma clang diagnostic ignored "-Wnested-anon-types"
#pragma clang diagnostic ignored "-Wsign-compare"
#pragma clang diagnostic ignored "-Wunused-private-field"
#if (__clang_major__ > 3) || (__clang_major__ == 3 && __clang_minor__ > 5)
// This was introduced in 3.6
#pragma clang diagnostic ignored "-Winconsistent-missing-override"
#endif // clang > 3.5
#endif

#if defined(__GNUC__) && !defined(__INTEL_COMPILER) && !defined(__clang__)
// GCC > 4.5 supports diagnostic pragmas with push/pop
#if (__GNUC__ > 4) || (__GNUC__ == 4 && __GNUC_MINOR__ > 5)
#pragma GCC diagnostic push
// These two don't work?
#pragma GCC diagnostic ignored "-Wpedantic"
#pragma GCC diagnostic ignored "-Wdeprecated"
// But this is helpful with some MPI stacks
#pragma GCC diagnostic ignored "-Wunused-parameter"
// Ignore warnings from code that uses deprecated members of std, like std::auto_ptr.
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#if (__GNUC__ > 5)
// Ignore warnings from code that does "if (foo) bar();"
#pragma GCC diagnostic ignored "-Wmisleading-indentation"
// Ignore warnings from bad placement new use
#pragma GCC diagnostic ignored "-Wplacement-new"
#if (__GNUC__ > 6)
#pragma GCC diagnostic ignored "-Wint-in-bool-context"
#pragma GCC diagnostic ignored "-Wimplicit-fallthrough"
#endif // GCC > 6
#endif // GCC > 5
#endif // GCC > 4.5
#endif // __GNUC__ && !__INTEL_COMPILER
