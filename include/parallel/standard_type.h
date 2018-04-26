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


#ifndef LIBMESH_STANDARD_TYPE_H
#define LIBMESH_STANDARD_TYPE_H

// Parallel includes
#include "libmesh/data_type.h"

// libMesh Includes
#include "libmesh/libmesh_common.h"

namespace libMesh
{

/**
 * The Parallel namespace is for wrapper functions
 * for common general parallel synchronization tasks.
 */
namespace Parallel
{

//-------------------------------------------------------------------

// Templated helper class to be used with static_assert.
template<typename T>
struct standardtype_dependent_false : std::false_type
{};

/**
 * Templated class to provide the appropriate MPI datatype
 * for use with built-in C types or simple C++ constructions.
 *
 * More complicated data types may need to provide a pointer-to-T so
 * that we can use MPI_Address without constructing a new T.
 */
template <typename T>
class StandardType : public DataType
{
  // Get a slightly better compiler diagnostic
  static_assert(standardtype_dependent_false<T>::value,
                "Only specializations of StandardType may be used, did you forget to include a header file (e.g. parallel_algebra.h)?");

  /*
   * The unspecialized class is useless, so we make its constructor
   * private to catch mistakes at compile-time rather than link-time.
   * Specializations should have a public constructor of the same
   * form.
   */
private:
  StandardType(const T * example = libmesh_nullptr);
};

} // namespace Parallel

} // namespace libMesh

#endif // LIBMESH_STANDARD_TYPE_H
