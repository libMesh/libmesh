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


#ifndef LIBMESH_DATA_TYPE_H
#define LIBMESH_DATA_TYPE_H

// Parallel includes
#include "libmesh/libmesh_call_mpi.h"

// libMesh includes
#include "libmesh/libmesh_common.h"

namespace libMesh
{

/**
 * The Parallel namespace is for wrapper functions
 * for common general parallel synchronization tasks.
 *
 * For MPI 1.1 compatibility, temporary buffers are used
 * instead of MPI 2's MPI_IN_PLACE
 */
namespace Parallel
{

#ifdef LIBMESH_HAVE_MPI
//-------------------------------------------------------------------
/**
 * Data types for communication
 */
typedef MPI_Datatype data_type;

#else

// These shouldn't actually be needed, but must be
// unique types for function overloading to work
// properly.
struct data_type    { /* unsigned int t; */ };

#endif // LIBMESH_HAVE_MPI



//-------------------------------------------------------------------
/**
 * Encapsulates the MPI_Datatype.
 */
class DataType
{
public:
  DataType () : _datatype() {}

  DataType (const DataType & other) :
    _datatype(other._datatype)
  {}

  DataType (const data_type & type) :
    _datatype(type)
  {}

  DataType (const DataType & other, unsigned int count)
  {
    // FIXME - if we nest an inner type here will we run into bug
    // https://github.com/libMesh/libmesh/issues/631 again?
    libmesh_call_mpi(MPI_Type_contiguous(count, other._datatype, &_datatype));
    libmesh_ignore(other, count); // ifndef LIBMESH_HAVE_MPI
    this->commit();
  }

  DataType & operator = (const DataType & other)
  { _datatype = other._datatype; return *this; }

  DataType & operator = (const data_type & type)
  { _datatype = type; return *this; }

  operator const data_type & () const
  { return _datatype; }

  operator data_type & ()
  { return _datatype; }

  //     operator data_type const * () const
  //     { return &_datatype; }

  //     operator data_type * ()
  //     { return &_datatype; }

  void commit ()
  {
    libmesh_call_mpi
      (MPI_Type_commit (&_datatype));
  }

  void free ()
  {
    libmesh_call_mpi
      (MPI_Type_free (&_datatype));
  }

protected:

  data_type _datatype;
};

} // namespace Parallel

} // namespace libMesh

#endif // LIBMESH_DATA_TYPE_H
