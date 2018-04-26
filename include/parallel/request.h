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


#ifndef LIBMESH_REQUEST_H
#define LIBMESH_REQUEST_H

// Parallel includes
#include "libmesh/post_wait_work.h"
#include "libmesh/status.h"

// libMesh Includes
#include "libmesh/libmesh_common.h"

// C++ includes
#include <memory>
#include <vector>
#include <utility>

namespace libMesh
{

namespace Parallel
{

#ifdef LIBMESH_HAVE_MPI

//-------------------------------------------------------------------
/**
 * Request object for non-blocking I/O
 */
typedef MPI_Request request;

#else

// This shouldn't actually be needed, but must be
// a unique type for function overloading to work
// properly.
struct request      { /* unsigned int r; */ };
#endif // LIBMESH_HAVE_MPI


//-------------------------------------------------------------------
/**
 * Encapsulates the MPI_Request
 */
class Request
{
public:
  Request ();

  Request (const request & r);

  Request (const Request & other);

  void cleanup();

  Request & operator = (const Request & other);

  Request & operator = (const request & r);

  ~Request ();

  request * get() { return &_request; }

  const request * get() const { return &_request; }

  Status wait ();

  bool test ();

  bool test (status & status);

  void add_prior_request(const Request & req);

  void add_post_wait_work(PostWaitWork * work);

private:
  request _request;

  // Breaking non-blocking sends into multiple requests can require chaining
  // multiple requests into a single Request
  std::unique_ptr<Request> _prior_request;

  // post_wait_work->first is a vector of work to do after a wait
  // finishes; post_wait_work->second is a reference count so that
  // Request objects will behave roughly like a shared_ptr and be
  // usable in STL containers
  std::pair<std::vector <PostWaitWork * >, unsigned int> * post_wait_work;
};

/**
 * Wait for a non-blocking send or receive to finish
 */
inline Status wait (Request & r) { return r.wait(); }

/**
 * Wait for a non-blocking send or receive to finish
 */
inline void wait (std::vector<Request> & r)
{ for (std::size_t i=0; i<r.size(); i++) r[i].wait(); }


} // namespace Parallel

} // namespace libMesh

#endif // LIBMESH_REQUEST_H
