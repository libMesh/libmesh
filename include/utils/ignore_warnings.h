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
#pragma clang diagnostic ignored "-Wnested-anon-types"
#pragma clang diagnostic ignored "-Wsign-compare"
#pragma clang diagnostic ignored "-Wunused-private-field"
#pragma clang diagnostic ignored "-Wextra"
#pragma clang diagnostic ignored "-Wredundant-decls"
#pragma clang diagnostic ignored "-Wcast-qual"
#pragma clang diagnostic ignored "-Wswitch-default"
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
#pragma clang diagnostic ignored "-Woverloaded-virtual"
// This isn't supported in 3.4.2 at least
#if (__clang_major__ > 3)
#pragma clang diagnostic ignored "-Wmacro-redefined"
#endif
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
// And these are for cppunit
#pragma GCC diagnostic ignored "-Woverloaded-virtual"
#pragma GCC diagnostic ignored "-Wundef"
// And these are for Trilinos
#pragma GCC diagnostic ignored "-Wextra"
#pragma GCC diagnostic ignored "-Wshadow"
// And these are for PETSc
#pragma GCC diagnostic ignored "-Wredundant-decls"
#pragma GCC diagnostic ignored "-Wcast-qual"
#pragma GCC diagnostic ignored "-Wswitch-default"
// And these are for Eigen
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wstack-protector"
// And this for VTK
#pragma GCC diagnostic ignored "-Wlogical-op"
// Ignore warnings from code that uses deprecated members of std, like std::auto_ptr.
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#if (__GNUC__ > 5)
// Ignore this for VTK
#pragma GCC diagnostic ignored "-Wfloat-conversion"
// Ignore warnings from code that does "if (foo) bar();"
#pragma GCC diagnostic ignored "-Wmisleading-indentation"
// Ignore warnings from bad placement new use
#pragma GCC diagnostic ignored "-Wplacement-new"
#if (__GNUC__ > 6)
#pragma GCC diagnostic ignored "-Wint-in-bool-context"
#pragma GCC diagnostic ignored "-Wimplicit-fallthrough"
#if (__GNUC__ > 7)
#pragma GCC diagnostic ignored "-Wparentheses"
#endif // GCC > 7
#endif // GCC > 6
#endif // GCC > 5
#endif // GCC > 4.5
#endif // __GNUC__ && !__INTEL_COMPILER
