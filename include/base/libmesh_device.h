// The libMesh Finite Element Library.
// Copyright (C) 2002-2026 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#ifndef LIBMESH_LIBMESH_DEVICE_H
#define LIBMESH_LIBMESH_DEVICE_H

// Defines LIBMESH_DEVICE_INLINE, mirroring MetaPhysicL's METAPHYSICL_INLINE
// pattern (metaphysicl_device.h / METAPHYSICL_KOKKOS_COMPILATION).
//
// When compiling a .K translation unit (ACSM_KOKKOS_COMPILATION is defined
// by kokkos.mk), this expands to KOKKOS_INLINE_FUNCTION so that annotated
// methods are callable from both host and device code.  In all other
// translation units it expands to plain `inline`.
#ifdef ACSM_KOKKOS_COMPILATION
#  include "libmesh/ignore_warnings.h"
#  include <Kokkos_Macros.hpp>
#  include <Kokkos_Abort.hpp>
#  include <Kokkos_MathematicalFunctions.hpp>
#  include "libmesh/restore_warnings.h"
#  define LIBMESH_DEVICE_INLINE KOKKOS_INLINE_FUNCTION

// Backend-neutral device-code detection for Kokkos .K translation units.
// This lets error/exception plumbing share a single predicate instead of
// hardcoding per-backend checks in multiple headers.
#  if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__) || defined(__SYCL_DEVICE_ONLY__)
#    define LIBMESH_IN_DEVICE_CODE 1
#  else
#    define LIBMESH_IN_DEVICE_CODE 0
#  endif

// Two-level stringification, so that a macro argument such as __LINE__
// is expanded before it is turned into a string literal.
#  define LIBMESH_STRINGIFY_IMPL(x) #x
#  define LIBMESH_STRINGIFY(x) LIBMESH_STRINGIFY_IMPL(x)

// Device-safe assert and error. Kokkos::abort() prints its message and
// terminates portably on every backend, so the whole message, location
// included, is assembled as a single string literal at compile time;
// this avoids a separate formatted printf call in device code.
// Defined here (not in libmesh_common.h) because Kokkos headers are
// only available in .K translation units.
#  ifndef NDEBUG
#    define LIBMESH_DEVICE_ASSERT(asserted)                             \
       do { if (!(asserted)) {                                          \
         ::Kokkos::abort("libMesh assert failed: " #asserted            \
                         ", file " __FILE__                             \
                         ", line " LIBMESH_STRINGIFY(__LINE__));        \
       } } while (0)
#  else
#    define LIBMESH_DEVICE_ASSERT(asserted) ((void) 0)
#  endif

// The message is stringified rather than evaluated: libMesh error
// messages are often stream expressions (e.g. "text" << value), which
// cannot be formatted in device code, and host-only inline functions
// using them are still compiled in the device pass of a .K file.
#  define LIBMESH_DEVICE_ERROR_MSG(msg)                                 \
     ::Kokkos::abort("libMesh error: " #msg                             \
                     ", file " __FILE__                                 \
                     ", line " LIBMESH_STRINGIFY(__LINE__))

#  define LIBMESH_DEVICE_ERROR_MSG_IF(cond, msg)                        \
     do { if (cond) { LIBMESH_DEVICE_ERROR_MSG(msg); } } while (0)

#else
#  define LIBMESH_DEVICE_INLINE inline
#  define LIBMESH_IN_DEVICE_CODE 0

namespace std {}
#endif

namespace libMesh
{
// Namespace providing the math functions (sqrt, abs, ...) for
// LIBMESH_DEVICE_INLINE functions: std:: functions are not callable from
// device code, while their Kokkos counterparts are callable from both host
// and device.  Use as `using portable_math::sqrt; sqrt(x);` so that
// argument-dependent lookup still finds overloads for non-arithmetic types
// such as MetaPhysicL::DualNumber.
#ifdef ACSM_KOKKOS_COMPILATION
namespace portable_math = ::Kokkos;
#else
namespace portable_math = std;
#endif
}

#endif // LIBMESH_LIBMESH_DEVICE_H
