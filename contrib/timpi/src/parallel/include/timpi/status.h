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


#ifndef TIMPI_STATUS_H
#define TIMPI_STATUS_H

// TIMPI includes
#include "timpi/data_type.h"
#include "timpi/timpi_assert.h"
#include "timpi/timpi_config.h"

// C/C++ includes
#ifdef TIMPI_HAVE_MPI
#  include "timpi/ignore_warnings.h"
#  include "mpi.h"
#  include "timpi/restore_warnings.h"
#endif // TIMPI_HAVE_MPI

namespace TIMPI
{

#ifdef TIMPI_HAVE_MPI

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

#endif // TIMPI_HAVE_MPI



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

// ------------------------------------------------------------
// Status member functions

inline Status::Status () :
  _status(),
  _datatype()
{}

inline Status::Status (const data_type & type) :
  _status(),
  _datatype(type)
{}

inline Status::Status (const status & stat) :
  _status(stat),
  _datatype()
{}

inline Status::Status (const status & stat,
                       const data_type & type) :
  _status(stat),
  _datatype(type)
{}

inline Status::Status (const Status & stat) :
  _status(stat._status),
  _datatype(stat._datatype)
{}

inline Status::Status (const Status    & stat,
                       const data_type & type) :
  _status(stat._status),
  _datatype(type)
{}

inline int Status::source () const
{
#ifdef TIMPI_HAVE_MPI
  return _status.MPI_SOURCE;
#else
  return 0;
#endif
}

inline int Status::tag () const
{
#ifdef TIMPI_HAVE_MPI
  return _status.MPI_TAG;
#else
  timpi_not_implemented();
  return 0;
#endif
}

inline unsigned int Status::size (const data_type & type) const
{
  ignore(type); // We don't use this ifndef TIMPI_HAVE_MPI
  int msg_size = 1;
  timpi_call_mpi
    (MPI_Get_count (const_cast<MPI_Status*>(&_status), type,
                    &msg_size));

  timpi_assert_greater_equal (msg_size, 0);
  return msg_size;
}

inline unsigned int Status::size () const
{ return this->size (this->datatype()); }


} // namespace TIMPI

#endif // TIMPI_STATUS_H
