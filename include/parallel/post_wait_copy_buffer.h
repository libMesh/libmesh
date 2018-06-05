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


#ifndef LIBMESH_POST_WAIT_COPY_BUFFER_H
#define LIBMESH_POST_WAIT_COPY_BUFFER_H

// Parallel includes
#include "libmesh/post_wait_work.h"

// libMesh Includes
#include "libmesh/libmesh_common.h"

// C++ includes
#include <algorithm>

namespace libMesh
{

namespace Parallel
{

// PostWaitWork specialization for copying from temporary to
// output containers
template <typename Container, typename OutputIter>
struct PostWaitCopyBuffer : public PostWaitWork {
  PostWaitCopyBuffer(const Container & buffer, const OutputIter out)
    : _buf(buffer), _out(out) {}

  virtual void run() override { std::copy(_buf.begin(), _buf.end(), _out); }

private:
  const Container & _buf;
  OutputIter _out;
};

} // namespace Parallel

} // namespace libMesh

#endif // LIBMESH_POST_WAIT_COPY_BUFFER_H
