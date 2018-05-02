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

// C++ includes
#include <numeric>
#include <utility>

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



// ------------------------------------------------------------
// Declare StandardType specializations for C++ built-in types

#ifdef LIBMESH_HAVE_MPI

#define LIBMESH_STANDARD_TYPE(cxxtype,mpitype)                          \
  template<>                                                            \
  class StandardType<cxxtype> : public DataType                         \
  {                                                                     \
  public:                                                               \
    explicit                                                            \
      StandardType(const cxxtype * = libmesh_nullptr) : DataType(mpitype) {} \
  }

#else

#define LIBMESH_STANDARD_TYPE(cxxtype,mpitype)                          \
  template<>                                                            \
  class StandardType<cxxtype> : public DataType                         \
  {                                                                     \
  public:                                                               \
    explicit                                                            \
      StandardType(const cxxtype * = libmesh_nullptr) : DataType() {}   \
  }

#endif

LIBMESH_STANDARD_TYPE(char,MPI_CHAR);
LIBMESH_STANDARD_TYPE(signed char,MPI_SIGNED_CHAR);
LIBMESH_STANDARD_TYPE(unsigned char,MPI_UNSIGNED_CHAR);
LIBMESH_STANDARD_TYPE(short int,MPI_SHORT);
LIBMESH_STANDARD_TYPE(unsigned short int,MPI_UNSIGNED_SHORT);
LIBMESH_STANDARD_TYPE(int,MPI_INT);
LIBMESH_STANDARD_TYPE(unsigned int,MPI_UNSIGNED);
LIBMESH_STANDARD_TYPE(long,MPI_LONG);
LIBMESH_STANDARD_TYPE(long long,MPI_LONG_LONG_INT);
LIBMESH_STANDARD_TYPE(unsigned long,MPI_UNSIGNED_LONG);
LIBMESH_STANDARD_TYPE(unsigned long long,MPI_UNSIGNED_LONG_LONG);
LIBMESH_STANDARD_TYPE(float,MPI_FLOAT);
LIBMESH_STANDARD_TYPE(double,MPI_DOUBLE);
LIBMESH_STANDARD_TYPE(long double,MPI_LONG_DOUBLE);

template<typename T1, typename T2>
class StandardType<std::pair<T1, T2>> : public DataType
{
public:
  explicit
  StandardType(const std::pair<T1, T2> * example = libmesh_nullptr) {
    // We need an example for MPI_Address to use
    static const std::pair<T1, T2> p;
    if (!example)
      example = &p;

#ifdef LIBMESH_HAVE_MPI

    // Get the sub-data-types, and make sure they live long enough
    // to construct the derived type
    StandardType<T1> d1(&example->first);
    StandardType<T2> d2(&example->second);

    MPI_Datatype types[] = { (data_type)d1, (data_type)d2 };
    int blocklengths[] = {1,1};
    MPI_Aint displs[2], start;

    libmesh_call_mpi
      (MPI_Get_address (const_cast<std::pair<T1,T2> *>(example),
                        &start));
    libmesh_call_mpi
      (MPI_Get_address (const_cast<T1*>(&example->first),
                        &displs[0]));
    libmesh_call_mpi
      (MPI_Get_address (const_cast<T2*>(&example->second),
                        &displs[1]));
    displs[0] -= start;
    displs[1] -= start;

    // create a prototype structure
    MPI_Datatype tmptype;
    libmesh_call_mpi
      (MPI_Type_create_struct (2, blocklengths, displs, types,
                               &tmptype));
    libmesh_call_mpi
      (MPI_Type_commit (&tmptype));

    // resize the structure type to account for padding, if any
    libmesh_call_mpi
      (MPI_Type_create_resized (tmptype, 0,
                                sizeof(std::pair<T1,T2>),
                                &_datatype));
    libmesh_call_mpi
      (MPI_Type_free (&tmptype));

    this->commit();

#endif // LIBMESH_HAVE_MPI

  }

  StandardType(const StandardType<std::pair<T1, T2>> & t)
  {
#ifdef LIBMESH_HAVE_MPI
    libmesh_call_mpi
      (MPI_Type_dup (t._datatype, &_datatype));
#endif
  }

  ~StandardType() { this->free(); }
};

template<typename T>
class StandardType<std::complex<T>> : public DataType
{
public:
  explicit
  StandardType(const std::complex<T> * /*example*/ = libmesh_nullptr) :
    DataType(StandardType<T>(libmesh_nullptr), 2) {}

  ~StandardType() { this->free(); }
};



} // namespace Parallel

} // namespace libMesh

#endif // LIBMESH_STANDARD_TYPE_H
