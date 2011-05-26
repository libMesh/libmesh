// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
  
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


#ifndef __parallel_hilbert_h__
#define __parallel_hilbert_h__

// This class contains all the functionality for bin sorting
// Templated on the type of keys you will be sorting and the
// type of iterator you will be using.


#include "libmesh_config.h"

#if defined(LIBMESH_HAVE_LIBHILBERT) && defined(LIBMESH_HAVE_MPI)

#include "hilbert.h"
#include "parallel.h"

namespace libMesh {
namespace Parallel {
  // A StandardType<> specialization to return a derived MPI datatype
  // to handle communication of HilbertIndices.  We use a singleton
  // pattern here because a global variable would have tried to call
  // MPI functions before MPI got initialized.
  template <>
  class StandardType<Hilbert::HilbertIndices> : public DataType
  {
  public:
    inline StandardType(const Hilbert::HilbertIndices *example=NULL) {
      // Make stupid compiler think "example" is used
      libmesh_ignore(example);

      // _static_type never gets freed, but it only gets committed once
      // so it's not a *huge* memory leak...
      static DataType _static_type;
      static bool _is_initialized = false;
      if (!_is_initialized)
        {
          _static_type = DataType(Parallel::StandardType<unsigned int>(), 3);
          _is_initialized = true;
        }
      _datatype = _static_type;
    }
  };
} // namespace Parallel
} // namespace libMesh

#endif // LIBMESH_HAVE_LIBHILBERT && LIBMESH_HAVE_MPI

#endif // __parallel_hilbert_h__

