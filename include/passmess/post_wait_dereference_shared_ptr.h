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


#ifndef LIBMESH_POST_WAIT_DEREFERENCE_SHARED_PTR_H
#define LIBMESH_POST_WAIT_DEREFERENCE_SHARED_PTR_H

// Parallel includes
#include "libmesh/message_tag.h"
#include "libmesh/post_wait_work.h"

// libMesh Includes
#include "libmesh/libmesh_common.h"

// System Includes
#include <memory>

namespace libMesh
{

namespace Parallel
{

// PostWaitWork specialization for holding a shared_ptr reference
template <typename T>
struct PostWaitDereferenceSharedPtr : public PostWaitWork {
  PostWaitDereferenceSharedPtr(const std::shared_ptr<T> & ptr) : _ptr(ptr) {}

  // Our real work is done by the destructor;
  virtual void run() {}

private:
  const std::shared_ptr<T> _ptr;
};

} // namespace Parallel

} // namespace libMesh

#endif // LIBMESH_POST_WAIT_DEREFERENCE_SHARED_PTR_H
