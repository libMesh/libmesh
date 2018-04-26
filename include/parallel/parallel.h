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


#ifndef LIBMESH_PARALLEL_H
#define LIBMESH_PARALLEL_H

// Parallel includes
#include "libmesh/attributes.h"
#include "libmesh/communicator.h"
#include "libmesh/data_type.h"
#include "libmesh/libmesh_call_mpi.h"
#include "libmesh/message_tag.h"
#include "libmesh/op_function.h"
#include "libmesh/packing.h"
#include "libmesh/parallel_only.h"
#include "libmesh/post_wait_copy_buffer.h"
#include "libmesh/post_wait_delete_buffer.h"
#include "libmesh/post_wait_unpack_buffer.h"
#include "libmesh/post_wait_work.h"
#include "libmesh/request.h"
#include "libmesh/status.h"
#include "libmesh/standard_type.h"

// libMesh Includes
#include "libmesh/libmesh_common.h" // libmesh_assert, cast_int
#include "libmesh/libmesh_logging.h"
#include "libmesh/auto_ptr.h" // deprecated

// C++ includes
#include <cstddef>
#include <climits>
#include <iterator>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <vector>
#include <memory>

namespace libMesh
{

/**
 * The Parallel namespace is for wrapper functions
 * for common general parallel synchronization tasks.
 *
 * For MPI 1.1 compatibility, temporary buffers are used
 * instead of MPI 2's MPI_IN_PLACE
 */
namespace Parallel
{

// FakeCommunicator for debugging inappropriate CommWorld uses
class FakeCommunicator
{
  operator Communicator & ()
  {
    libmesh_not_implemented();
    static Communicator temp;
    return temp;
  }
};

} // namespace Parallel

  /**
   * The default libMesh communicator.
   *
   * If this communicator is disabled, we also disable it as a default
   * argument to functions which accept a default communicator
   * argument.  This should expose implicit uses of the default
   * communicator as compile-time rather than run-time errors.
   *
   * The macro LIBMESH_CAN_DEFAULT_TO_COMMWORLD effects this
   * functionality; it is empty (and so leaves arguments with no
   * default value) if the default is disabled, and it sets something
   * equal to the default otherwise.
   */
#ifdef LIBMESH_DISABLE_COMMWORLD
extern Parallel::FakeCommunicator CommWorld;
#define LIBMESH_CAN_DEFAULT_TO_COMMWORLD
#else
extern Parallel::Communicator CommWorld;
#define LIBMESH_CAN_DEFAULT_TO_COMMWORLD = libMesh::CommWorld
#endif

} // namespace libMesh

// Define all the implementations separately; users might want to look
// through this file for APIs, and it's long enough already.

#include "libmesh/parallel_implementation.h"

#endif // LIBMESH_PARALLEL_H
