// The libMesh Finite Element Library.
// Copyright (C) 2002-2016 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


#ifndef LIBMESH_PARALLEL_HILBERT_H
#define LIBMESH_PARALLEL_HILBERT_H

// This class contains all the functionality for bin sorting
// Templated on the type of keys you will be sorting and the
// type of iterator you will be using.

#include "libmesh/libmesh_config.h"

#if defined(LIBMESH_HAVE_LIBHILBERT)

// Local includes
#include "hilbert.h"
#include "libmesh/parallel.h"

// C++ includes
#include <cstddef>

namespace libMesh {
namespace Parallel {

#ifdef LIBMESH_HAVE_MPI
// A StandardType<> specialization to return a derived MPI datatype
// to handle communication of HilbertIndices.  We use a singleton
// pattern here because a global variable would have tried to call
// MPI functions before MPI got initialized.
template <>
class StandardType<Hilbert::HilbertIndices> : public DataType
{
public:
  explicit
  StandardType(const Hilbert::HilbertIndices * =libmesh_nullptr) {
    // _static_type never gets freed, but it only gets committed once
    // so it's not a *huge* memory leak...
    static DataType _static_type;
    static bool _is_initialized = false;
    if (!_is_initialized)
      {
        _static_type = DataType(Parallel::StandardType<Hilbert::inttype>(), 3);
        _is_initialized = true;
      }
    _datatype = _static_type;
  }
};

#endif // LIBMESH_HAVE_MPI

#ifdef LIBMESH_ENABLE_UNIQUE_ID
  typedef
  std::pair<Hilbert::HilbertIndices, unique_id_type> DofObjectKey;
#else
  typedef
  Hilbert::HilbertIndices DofObjectKey;
#endif


} // namespace Parallel


} // namespace libMesh


namespace Hilbert {

// This has to be in the Hilbert namespace for Koenig lookup to work?
// g++ doesn't find it if it's in the global namespace.
// XCode didn't find it in the libMesh namespace.
#ifdef LIBMESH_ENABLE_UNIQUE_ID
inline
std::ostream&
operator <<
  (std::ostream& os,
   const libMesh::Parallel::DofObjectKey & hilbert_pair)
{
  os << '(' << hilbert_pair.first << ',' << hilbert_pair.second << ')' << std::endl;
  return os;
}
#endif

}


// Appropriate operator< definitions for std::pair let the same code handle
// both DofObjectKey types

inline
void dofobjectkey_max_op (libMesh::Parallel::DofObjectKey *in,
                          libMesh::Parallel::DofObjectKey *inout,
                          int *len, void *)
{
  // When (*in <= *inout), then inout already contains max(*in,*inout)
  // Otherwise we need to copy from in.
  for (int i=0; i<*len; i++, in++, inout++)
    if (*inout < *in)
      *inout = *in;
}

inline
void dofobjectkey_min_op (libMesh::Parallel::DofObjectKey *in,
                          libMesh::Parallel::DofObjectKey *inout,
                          int *len, void *)
{
  // When (*in >= *inout), then inout already contains min(*in,*inout)
  // Otherwise we need to copy from in.
  for (int i=0; i<*len; i++, in++, inout++)
    if (*in < *inout)
      *inout = *in;
}

#endif // LIBMESH_HAVE_LIBHILBERT && LIBMESH_HAVE_MPI

#endif // LIBMESH_PARALLEL_HILBERT_H
