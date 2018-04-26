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


#ifndef LIBMESH_OP_FUNCTION_H
#define LIBMESH_OP_FUNCTION_H

// C++ includes
#include <type_traits>

namespace libMesh
{

namespace Parallel
{
//-------------------------------------------------------------------

// Templated helper class to be used with static_assert.
template<typename T>
struct opfunction_dependent_false : std::false_type
{};

/**
 * Templated class to provide the appropriate MPI reduction operations
 * for use with built-in C types or simple C++ constructions.
 *
 * More complicated data types may need to provide a pointer-to-T so
 * that we can use MPI_Address without constructing a new T.
 */
template <typename T>
class OpFunction
{
  // Get a slightly better compiler diagnostic if we have C++11
  static_assert(opfunction_dependent_false<T>::value,
                "Only specializations of OpFunction may be used, did you forget to include a header file (e.g. parallel_algebra.h)?");

  /*
   * The unspecialized class defines none of these functions;
   * specializations will need to define any functions that need to be
   * usable.
   *
   * Most specializations will just return MPI_MIN, etc, but we'll use
   * a whitelist rather than a default implementation, so that any
   * attempt to perform a reduction on an unspecialized type will be a
   * compile-time rather than a run-time failure.
   */
  // static MPI_Op max();
  // static MPI_Op min();
  // static MPI_Op sum();
  // static MPI_Op product();
  // static MPI_Op logical_and();
  // static MPI_Op bitwise_and();
  // static MPI_Op logical_or();
  // static MPI_Op bitwise_or();
  // static MPI_Op logical_xor();
  // static MPI_Op bitwise_xor();
  // static MPI_Op max_loc();
  // static MPI_Op min_loc();
};

} // namespace Parallel

} // namespace libMesh

#endif // LIBMESH_OP_FUNCTION_H
