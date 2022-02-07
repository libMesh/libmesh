// The libMesh Finite Element Library.
// Copyright (C) 2002-2022 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

// libMesh Includes
#include "libmesh/auto_ptr.h" // for backwards compatibility, not internal use
#include "libmesh/libmesh_call_mpi.h"
#include "libmesh/libmesh_common.h"
#include "libmesh/parallel_only.h"

// TIMPI includes
#include "timpi/attributes.h"
#include "timpi/communicator.h"
#include "timpi/data_type.h"
#include "timpi/message_tag.h"
#include "timpi/op_function.h"
#include "timpi/packing.h"
#include "timpi/timpi_call_mpi.h"
#include "timpi/post_wait_copy_buffer.h"
#include "timpi/post_wait_delete_buffer.h"
#include "timpi/post_wait_dereference_shared_ptr.h"
#include "timpi/post_wait_dereference_tag.h"
#include "timpi/post_wait_free_buffer.h"
#include "timpi/post_wait_unpack_buffer.h"
#include "timpi/post_wait_work.h"
#include "timpi/request.h"
#include "timpi/status.h"
#include "timpi/standard_type.h"

// Define all the implementations separately; users might want to look
// through this file for APIs, and it's long enough already.

#include "timpi/parallel_implementation.h"

#endif // LIBMESH_PARALLEL_H
