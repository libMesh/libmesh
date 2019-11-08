// The TIMPI Message-Passing Parallelism Library.
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


#ifndef TIMPI_POST_WAIT_WORK_H
#define TIMPI_POST_WAIT_WORK_H

namespace TIMPI
{
//-------------------------------------------------------------------
/**
 * An abstract base class that can be subclassed to allow other code
 * to perform work after a MPI_Wait succeeds.  This makes it possible
 * to automatically schedule deserialization, cleanup, or other
 * operations such that they automatically take place after an
 * asynchronous I/O operation has completed.
 *
 * See PostWaitCopyBuffer, PostWaitDeleteBuffer, and
 * PostWaitUnpackBuffer for examples of useful subclasses that are
 * automatically employed by the library.  See the
 * Request::add_post_wait_work method for documentation of how to use
 * these subclasses.  See Communicator method implementations for
 * examples of all this in use, including chaining of multiple
 * PostWaitWork operations.
 */
struct PostWaitWork {
  virtual ~PostWaitWork() {}

  // Do *something* after a communication wait has succeeded.
  virtual void run() = 0;
};

} // namespace TIMPI

#endif // TIMPI_POST_WAIT_WORK_H
