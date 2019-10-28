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


#ifndef LIBMESH_POST_WAIT_DEREFERENCE_TAG_H
#define LIBMESH_POST_WAIT_DEREFERENCE_TAG_H

// Parallel includes
#include "libmesh/message_tag.h"
#include "libmesh/post_wait_work.h"

// libMesh Includes
#include "libmesh/libmesh_common.h"

namespace libMesh
{

namespace Parallel
{

// PostWaitWork specialization for holding a MessageTag.  This
// prevents the MessageTag from being completely dereferenced and thus
// prevents a unique tag number from being reused until after the
// Request has been cleaned up.
struct PostWaitDereferenceTag : public PostWaitWork {
  PostWaitDereferenceTag(const MessageTag & tag) : _tag(tag) {}

  // Our real work is done by the MessageTag member's destructor;
  // there's no advantage to dereferencing our tag slightly earlier.
  virtual void run() {}

private:
  MessageTag _tag;
};

} // namespace Parallel

} // namespace libMesh

#endif // LIBMESH_POST_WAIT_DEREFERENCE_TAG_H
