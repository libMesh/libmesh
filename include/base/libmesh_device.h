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
// When compiling a .K translation unit (LIBMESH_KOKKOS_COMPILATION is defined
// by kokkos.mk), this expands to KOKKOS_INLINE_FUNCTION so that annotated
// methods are callable from both host and device code.  In all other
// translation units it expands to plain `inline`.
#ifdef LIBMESH_KOKKOS_COMPILATION
#  include <Kokkos_Macros.hpp>
#  include <Kokkos_Abort.hpp>
#  define LIBMESH_DEVICE_INLINE KOKKOS_INLINE_FUNCTION

// Backend-neutral device-code detection for Kokkos .K translation units.
// This lets error/exception plumbing share a single predicate instead of
// hardcoding per-backend checks in multiple headers.
#  if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__) || defined(__SYCL_DEVICE_ONLY__)
#    define LIBMESH_IN_DEVICE_CODE 1
#  else
#    define LIBMESH_IN_DEVICE_CODE 0
#  endif

// Device-safe assert: uses printf (supported on CUDA/HIP) and
// Kokkos::abort() for backend-portable device termination.
// Defined here (not in libmesh_common.h) because Kokkos headers
// are only available in .K translation units.
#  ifndef NDEBUG
#    define LIBMESH_DEVICE_ASSERT(asserted)                             \
       do { if (!(asserted)) {                                          \
         printf("libMesh assert failed: %s, file %s, line %d\n",       \
                #asserted, __FILE__, __LINE__);                         \
         ::Kokkos::abort("libmesh_assert failed");                     \
       } } while (0)
#  else
#    define LIBMESH_DEVICE_ASSERT(asserted) ((void) 0)
#  endif

#  define LIBMESH_DEVICE_ERROR_MSG(msg)                                 \
     do {                                                               \
       printf("libMesh error: %s, file %s, line %d\n",                 \
              msg, __FILE__, __LINE__);                                 \
       ::Kokkos::abort(msg);                                            \
     } while (0)

#  define LIBMESH_DEVICE_ERROR_MSG_IF(cond, msg)                        \
     do { if (cond) { LIBMESH_DEVICE_ERROR_MSG(msg); } } while (0)

#else
#  define LIBMESH_DEVICE_INLINE inline
#  define LIBMESH_IN_DEVICE_CODE 0
#  define LIBMESH_DEVICE_ERROR_MSG(msg) libmesh_error_msg(msg)
#  define LIBMESH_DEVICE_ERROR_MSG_IF(cond, msg) libmesh_error_msg_if(cond, msg)
#endif

#endif // LIBMESH_LIBMESH_DEVICE_H
