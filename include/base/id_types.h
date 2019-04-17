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



#ifndef LIBMESH_ID_TYPES_H
#define LIBMESH_ID_TYPES_H

#include <limits>
#include <stdint.h>

#include "libmesh/libmesh_config.h"

namespace libMesh
{

// A useful way to debug:
#if 0
class TestClass {
  //int _c;
  unsigned int _c;
public:
  TestClass() : _c(0) {}
  TestClass(unsigned int c) : _c(c) {}
  TestClass & operator=(unsigned int c) { _c = c; return *this; }
  bool operator<(const TestClass & l) const { return _c < l._c; }
  operator int() const { return _c; }
};
typedef TestClass subdomain_id_type;
#endif


// How many bytes do we need to specify subsets of the boundary?  By
// default we'll allow for tens of thousands of different boundary
// ids.
#if LIBMESH_BOUNDARY_ID_BYTES == 1
typedef int8_t boundary_id_type;
#elif LIBMESH_BOUNDARY_ID_BYTES == 4
typedef int32_t boundary_id_type;
#elif LIBMESH_BOUNDARY_ID_BYTES == 8
typedef int64_t boundary_id_type;
#else // LIBMESH_BOUNDARY_ID_BYTES = 2 (default)
typedef int16_t boundary_id_type;
#endif


// How many bytes do we need to specify DoFObjects?  By default we'll
// allow for a few billion (each) nodes & elements.
//
// We'll need a sign bit in some internal buffers, so let's also
// define a signed type that we can convert to without loss.
#if LIBMESH_DOF_ID_BYTES == 1
typedef uint8_t dof_id_type;
typedef  int8_t dof_id_signed_type;
#elif LIBMESH_DOF_ID_BYTES == 2
typedef uint16_t dof_id_type;
typedef  int16_t dof_id_signed_type;
#elif LIBMESH_DOF_ID_BYTES == 8
typedef uint64_t dof_id_type;
typedef  int64_t dof_id_signed_type;
#else // LIBMESH_DOF_ID_BYTES = 4 (default)
typedef uint32_t dof_id_type;
typedef  int32_t dof_id_signed_type;
#endif



// How many bytes do we need to specify DoFObjects without ever
// renumbering?  By default we'll allow for quintillions of
// one-time-use ids.
#if LIBMESH_UNIQUE_ID_BYTES == 1
typedef uint8_t unique_id_type;
#elif LIBMESH_UNIQUE_ID_BYTES == 2
typedef uint16_t unique_id_type;
#elif LIBMESH_UNIQUE_ID_BYTES == 4
typedef uint32_t unique_id_type;
#else // LIBMESH_UNIQUE_ID_BYTES == 8 (default)
typedef uint64_t unique_id_type;

#endif


// We may want to specialize this later, but for now we'll assume
// numeric vector indices are the same as dof indices
typedef dof_id_type numeric_index_type;


// Define processor id storage type.
#if LIBMESH_PROCESSOR_ID_BYTES == 1
typedef uint8_t processor_id_type;
#elif LIBMESH_PROCESSOR_ID_BYTES == 2
typedef uint16_t processor_id_type;
#elif LIBMESH_PROCESSOR_ID_BYTES == 4
typedef uint32_t processor_id_type; // default
#elif LIBMESH_PROCESSOR_ID_BYTES == 8
typedef uint64_t processor_id_type;
#else
// We should not get here: it's not currently possible to set any
// other size.
DIE A HORRIBLE DEATH HERE...
#endif


// How many bytes do we need to specify subsets of the interior?  By
// default we'll allow for tens of thousands of different subdomain
// ids.
#if LIBMESH_SUBDOMAIN_ID_BYTES == 1
typedef uint8_t subdomain_id_type;
#elif LIBMESH_SUBDOMAIN_ID_BYTES == 4
/**
 * \note \p subdomain_id_type should be a positive integer, but due to
 * a limitation in the exodusII API, we are forced to use a signed
 * integer here to represent subdomains.  This gives us 2^31 possible
 * unique blocks.
 */
typedef int32_t subdomain_id_type;
#elif LIBMESH_SUBDOMAIN_ID_BYTES == 8
/**
 * Based on the 4-byte comment warning above, this probably doesn't
 * work with exodusII *at all*...
 */
typedef int64_t subdomain_id_type;
#else // LIBMESH_SUBDOMAIN_ID_BYTES = 2 (default)
typedef uint16_t subdomain_id_type;
#endif


// For serialization purposes we often like to pack the different
// kinds of ids together; how large a data type do we need to hold an
// arbitrary id?
#if (LIBMESH_BOUNDARY_ID_BYTES > 4) || (LIBMESH_DOF_ID_BYTES > 4) ||    \
  (LIBMESH_UNIQUE_ID_BYTES > 4) || (LIBMESH_PROCESSOR_ID_BYTES > 4) ||  \
  (LIBMESH_SUBDOMAIN_ID_BYTES > 4)
typedef uint64_t largest_id_type;
#elif (LIBMESH_BOUNDARY_ID_BYTES > 2) || (LIBMESH_DOF_ID_BYTES > 2) ||  \
  (LIBMESH_UNIQUE_ID_BYTES > 2) || (LIBMESH_PROCESSOR_ID_BYTES > 2) ||  \
  (LIBMESH_SUBDOMAIN_ID_BYTES > 2)
typedef uint32_t largest_id_type;
#elif (LIBMESH_BOUNDARY_ID_BYTES > 1) || (LIBMESH_DOF_ID_BYTES > 1) ||  \
  (LIBMESH_UNIQUE_ID_BYTES > 1) || (LIBMESH_PROCESSOR_ID_BYTES > 1) ||  \
  (LIBMESH_SUBDOMAIN_ID_BYTES > 1)
typedef uint16_t largest_id_type;
#else
typedef uint8_t largest_id_type;
#endif

} // namespace libMesh

#endif // LIBMESH_ID_TYPES_H
