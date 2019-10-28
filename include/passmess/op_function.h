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


#ifndef PASSMESS_OP_FUNCTION_H
#define PASSMESS_OP_FUNCTION_H

#include "libmesh_config.h"

#ifdef LIBMESH_HAVE_MPI
#  include "mpi.h"
#endif // LIBMESH_HAVE_MPI

// C++ includes
#include <functional>
#include <type_traits>



namespace libMesh
{

namespace Parallel
{
#ifdef LIBMESH_DEFAULT_QUADRUPLE_PRECISION
# ifdef LIBMESH_HAVE_MPI
# define LIBMESH_MPI_BINARY(funcname) \
inline void \
libmesh_mpi_##funcname(void * a, void * b, int * len, MPI_Datatype *) \
{ \
  const int size = *len; \
 \
  Real *in = static_cast<Real*>(a); \
  Real *inout = static_cast<Real*>(b); \
  for (int i=0; i != size; ++i) \
    inout[i] = std::funcname(in[i],inout[i]); \
}

# define LIBMESH_MPI_LOCATOR(funcname) \
inline void \
libmesh_mpi_##funcname##_location(void * a, void * b, int * len, MPI_Datatype *) \
{ \
  const int size = *len; \
 \
  typedef std::pair<Real, int> dtype; \
 \
  dtype *in = static_cast<dtype*>(a); \
  dtype *inout = static_cast<dtype*>(b); \
  for (int i=0; i != size; ++i) \
    { \
      Real old_inout = inout[i].first; \
      inout[i].first = std::funcname(in[i].first,inout[i].first); \
      if (old_inout != inout[i].first) \
        inout[i].second = in[i].second; \
    } \
}


# define LIBMESH_MPI_BINARY_FUNCTOR(funcname) \
inline void \
libmesh_mpi_##funcname(void * a, void * b, int * len, MPI_Datatype *) \
{ \
  const int size = *len; \
 \
  Real *in = static_cast<Real*>(a); \
  Real *inout = static_cast<Real*>(b); \
  for (int i=0; i != size; ++i) \
    inout[i] = std::funcname<Real>()(in[i],inout[i]); \
}


LIBMESH_MPI_BINARY(max)
LIBMESH_MPI_BINARY(min)
LIBMESH_MPI_LOCATOR(max)
LIBMESH_MPI_LOCATOR(min)
LIBMESH_MPI_BINARY_FUNCTOR(plus)
LIBMESH_MPI_BINARY_FUNCTOR(multiplies)

# endif // LIBMESH_HAVE_MPI
#endif // LIBMESH_DEFAULT_QUADRUPLE_PRECISION


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



// ------------------------------------------------------------
// Declare OpFunction specializations for C++ built-in types

#ifdef LIBMESH_HAVE_MPI

#define LIBMESH_PARALLEL_INTEGER_OPS(cxxtype)           \
  template<>                                            \
  class OpFunction<cxxtype>                             \
  {                                                     \
  public:                                               \
    static MPI_Op max()          { return MPI_MAX; }    \
    static MPI_Op min()          { return MPI_MIN; }    \
    static MPI_Op sum()          { return MPI_SUM; }    \
    static MPI_Op product()      { return MPI_PROD; }   \
    static MPI_Op logical_and()  { return MPI_LAND; }   \
    static MPI_Op bitwise_and()  { return MPI_BAND; }   \
    static MPI_Op logical_or()   { return MPI_LOR; }    \
    static MPI_Op bitwise_or()   { return MPI_BOR; }    \
    static MPI_Op logical_xor()  { return MPI_LXOR; }   \
    static MPI_Op bitwise_xor()  { return MPI_BXOR; }   \
    static MPI_Op max_location() { return MPI_MAXLOC; } \
    static MPI_Op min_location() { return MPI_MINLOC; } \
  }

#define LIBMESH_PARALLEL_FLOAT_OPS(cxxtype)             \
  template<>                                            \
  class OpFunction<cxxtype>                             \
  {                                                     \
  public:                                               \
    static MPI_Op max()          { return MPI_MAX; }    \
    static MPI_Op min()          { return MPI_MIN; }    \
    static MPI_Op sum()          { return MPI_SUM; }    \
    static MPI_Op product()      { return MPI_PROD; }   \
    static MPI_Op max_location() { return MPI_MAXLOC; } \
    static MPI_Op min_location() { return MPI_MINLOC; } \
  }

#else

#define LIBMESH_PARALLEL_INTEGER_OPS(cxxtype)   \
  template<>                                    \
  class OpFunction<cxxtype>                     \
  {                                             \
  }

#define LIBMESH_PARALLEL_FLOAT_OPS(cxxtype)     \
  template<>                                    \
  class OpFunction<cxxtype>                     \
  {                                             \
  }

#endif

LIBMESH_PARALLEL_INTEGER_OPS(char);
LIBMESH_PARALLEL_INTEGER_OPS(signed char);
LIBMESH_PARALLEL_INTEGER_OPS(unsigned char);
LIBMESH_PARALLEL_INTEGER_OPS(short int);
LIBMESH_PARALLEL_INTEGER_OPS(unsigned short int);
LIBMESH_PARALLEL_INTEGER_OPS(int);
LIBMESH_PARALLEL_INTEGER_OPS(unsigned int);
LIBMESH_PARALLEL_INTEGER_OPS(long);
LIBMESH_PARALLEL_INTEGER_OPS(long long);
LIBMESH_PARALLEL_INTEGER_OPS(unsigned long);
LIBMESH_PARALLEL_INTEGER_OPS(unsigned long long);                                \

LIBMESH_PARALLEL_FLOAT_OPS(float);
LIBMESH_PARALLEL_FLOAT_OPS(double);
LIBMESH_PARALLEL_FLOAT_OPS(long double);

#ifdef LIBMESH_DEFAULT_QUADRUPLE_PRECISION
# ifdef LIBMESH_HAVE_MPI

// FIXME - we'll still need to free these to keep valgrind happy
#define LIBMESH_MPI_OPFUNCTION(mpiname, funcname) \
  static MPI_Op mpiname() { \
    static MPI_Op LIBMESH_MPI_##mpiname = MPI_OP_NULL; \
    if (LIBMESH_MPI_##mpiname == MPI_OP_NULL) \
      { \
        libmesh_call_mpi \
          (MPI_Op_create(libmesh_mpi_##funcname, true, &LIBMESH_MPI_##mpiname)); \
      } \
    return LIBMESH_MPI_##mpiname;  \
  }

  template<>
  class OpFunction<Real>
  {
  public:
    LIBMESH_MPI_OPFUNCTION(max, max)
    LIBMESH_MPI_OPFUNCTION(min, min)
    LIBMESH_MPI_OPFUNCTION(sum, plus)
    LIBMESH_MPI_OPFUNCTION(product, multiplies)

    LIBMESH_MPI_OPFUNCTION(max_location, max_location)
    LIBMESH_MPI_OPFUNCTION(min_location, min_location)
  };

# else
  LIBMESH_PARALLEL_FLOAT_OPS(Real);
# endif
#endif

} // namespace Parallel

} // namespace libMesh

#endif // PASSMESS_OP_FUNCTION_H
