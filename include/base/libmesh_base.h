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



#ifndef LIBMESH_LIBMESH_BASE_H
#define LIBMESH_LIBMESH_BASE_H

#include "libmesh/id_types.h"

namespace libMesh {

#ifndef LIBMESH_DISABLE_COMMWORLD
/**
 * @returns the number of processors used in the current simulation.
 */
processor_id_type n_processors();

/**
 * @returns the index of the local processor.
 */
processor_id_type processor_id();
#endif

/**
 * @returns the number of processors libMesh was initialized with.
 */
processor_id_type global_n_processors();

/**
 * @returns the index of the local processor with respect to the
 * original MPI pool libMesh was initialized with.
 */
processor_id_type global_processor_id();

/**
 * @returns the maximum number of threads used in the simulation.
 */
unsigned int n_threads();

/**
 * Namespaces don't provide private data,
 * so let's take the data we would like
 * private and put it in an obnoxious
 * namespace.  At least that way it is a
 * pain to use, thus discouraging errors.
 */
namespace libMeshPrivateData {
#ifdef LIBMESH_HAVE_MPI
/**
 * Total number of processors used.
 */
extern processor_id_type _n_processors;

/**
 * The local processor id.
 */
extern processor_id_type _processor_id;
#endif

/**
 * Total number of threads possible.
 */
extern int _n_threads;
}
}



// ------------------------------------------------------------
// libMesh inline member functions
#ifndef LIBMESH_DISABLE_COMMWORLD
inline
libMesh::processor_id_type libMesh::n_processors()
{
  return libMesh::global_n_processors();
}



inline
libMesh::processor_id_type libMesh::processor_id()
{
  return libMesh::global_processor_id();
}
#endif // LIBMESH_DISABLE_COMMWORLD


inline
libMesh::processor_id_type libMesh::global_n_processors()
{
#ifdef LIBMESH_HAVE_MPI
  return libMeshPrivateData::_n_processors;
#else
  return 1;
#endif
}

inline
libMesh::processor_id_type libMesh::global_processor_id()
{
#ifdef LIBMESH_HAVE_MPI
  return libMeshPrivateData::_processor_id;
#else
  return 0;
#endif
}


inline
unsigned int libMesh::n_threads()
{
  return static_cast<unsigned int>(libMeshPrivateData::_n_threads);
}


// We now put everything we can into a separate libMesh namespace;
// code which forward declares libMesh classes or which specializes
// libMesh templates may want to know whether it is compiling under
// such conditions, to be backward compatible with older libMesh
// versions:
#define LIBMESH_USE_SEPARATE_NAMESPACE 1


// Unless configured otherwise, we import all of namespace libMesh,
// for backwards compatibility with pre-namespaced codes.

#ifndef LIBMESH_REQUIRE_SEPARATE_NAMESPACE
using namespace libMesh;
#endif


#endif // LIBMESH_LIBMESH_BASE_H
