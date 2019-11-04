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

// PassMess includes
#include "passmess/attributes.h"
#include "passmess/communicator.h"
#include "passmess/data_type.h"
#include "passmess/message_tag.h"
#include "passmess/op_function.h"
#include "passmess/packing.h"
#include "passmess/passmess_call_mpi.h"
#include "passmess/passmess_init.h"
#include "passmess/post_wait_copy_buffer.h"
#include "passmess/post_wait_delete_buffer.h"
#include "passmess/post_wait_dereference_shared_ptr.h"
#include "passmess/post_wait_dereference_tag.h"
#include "passmess/post_wait_free_buffer.h"
#include "passmess/post_wait_unpack_buffer.h"
#include "passmess/post_wait_work.h"
#include "passmess/request.h"
#include "passmess/status.h"
#include "passmess/standard_type.h"

// Define all the implementations separately; users might want to look
// through this file for APIs, and it's long enough already.

#include "passmess/parallel_implementation.h"

#endif // PASSMESS_PARALLEL_H
