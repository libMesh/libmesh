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


#ifndef LIBMESH_STATUS_H
#define LIBMESH_STATUS_H

// Parallel includes
#include "libmesh/data_type.h"

// libMesh Includes
#include "libmesh/libmesh_common.h"

namespace libMesh
{

namespace Parallel
{

#ifdef LIBMESH_HAVE_MPI

//-------------------------------------------------------------------

/**
 * Status object for querying messages
 */
typedef MPI_Status status;

#else

// This shouldn't actually be needed, but must be
// unique types for function overloading to work
// properly.
struct status       { /* unsigned int s; */ };

#endif // LIBMESH_HAVE_MPI



//-------------------------------------------------------------------
/**
 * Encapsulates the MPI_Status struct.  Allows the source and size
 * of the message to be determined.
 */
class Status
{
public:
  Status ();

  explicit Status (const data_type & type);

  explicit Status (const status & status);

  Status (const status    & status,
          const data_type & type);

  Status (const Status & status);

  Status (const Status    & status,
          const data_type & type);

  status * get() { return &_status; }

  status const * get() const { return &_status; }

  int source () const;

  int tag () const;

  data_type & datatype () { return _datatype; }

  const data_type & datatype () const { return _datatype; }

  unsigned int size (const data_type & type) const;

  unsigned int size () const;

private:

  status    _status;
  data_type _datatype;
};

} // namespace Parallel

} // namespace libMesh

#endif // LIBMESH_STATUS_H
