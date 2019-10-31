// The PassMess Message-Passing Parallelism Library.
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


#ifndef PASSMESS_PARALLEL_H
#define PASSMESS_PARALLEL_H

// Parallel includes
#include "libmesh/attributes.h"
#include "libmesh/communicator.h"
#include "libmesh/data_type.h"
#include "libmesh/passmess_call_mpi.h"
#include "libmesh/message_tag.h"
#include "libmesh/op_function.h"
#include "libmesh/packing.h"
#include "libmesh/parallel_only.h"
#include "libmesh/passmess_call_mpi.h"
#include "libmesh/post_wait_copy_buffer.h"
#include "libmesh/post_wait_delete_buffer.h"
#include "libmesh/post_wait_dereference_shared_ptr.h"
#include "libmesh/post_wait_dereference_tag.h"
#include "libmesh/post_wait_free_buffer.h"
#include "libmesh/post_wait_unpack_buffer.h"
#include "libmesh/post_wait_work.h"
#include "libmesh/request.h"
#include "libmesh/status.h"
#include "libmesh/standard_type.h"

// Define all the implementations separately; users might want to look
// through this file for APIs, and it's long enough already.

#include "libmesh/parallel_implementation.h"

#endif // PASSMESS_PARALLEL_H
