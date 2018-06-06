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


#ifndef LIBMESH_POST_WAIT_UNPACK_BUFFER_H
#define LIBMESH_POST_WAIT_UNPACK_BUFFER_H

// Parallel includes
#include "libmesh/post_wait_work.h"
#include "libmesh/packing.h"

// libMesh Includes
#include "libmesh/libmesh_common.h"

namespace libMesh
{

namespace Parallel
{

// PostWaitWork specialization for unpacking received buffers.
template <typename Container, typename Context, typename OutputIter,
          typename T>
struct PostWaitUnpackBuffer : public PostWaitWork {
  PostWaitUnpackBuffer(const Container & buffer, Context * context, OutputIter out) :
    _buf(buffer), _context(context), _out(out) {}

  virtual void run() override {

    Parallel::unpack_range(_buf, _context, _out, (T*)libmesh_nullptr);
  }

private:
  const Container & _buf;
  Context * _context;
  OutputIter _out;
};

} // namespace Parallel

} // namespace libMesh

#endif // LIBMESH_POST_WAIT_UNPACK_BUFFER_H
