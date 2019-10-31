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


#ifndef PASSMESS_POST_WAIT_UNPACK_BUFFER_H
#define PASSMESS_POST_WAIT_UNPACK_BUFFER_H

// Parallel includes
#include "passmess/post_wait_work.h"
#include "passmess/packing.h"

namespace PassMess
{

// PostWaitWork specialization for unpacking received buffers.
template <typename Container, typename Context, typename OutputIter,
          typename T>
struct PostWaitUnpackBuffer : public PostWaitWork {
  PostWaitUnpackBuffer(const Container & buffer, Context * context, OutputIter out) :
    _buf(buffer), _context(context), _out(out) {}

  virtual void run() override {

    unpack_range(_buf, _context, _out, (T*)nullptr);
  }

private:
  const Container & _buf;
  Context * _context;
  OutputIter _out;
};

} // namespace PassMess

#endif // PASSMESS_POST_WAIT_UNPACK_BUFFER_H
