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


#ifndef TIMPI_OP_FUNCTION_H
#define TIMPI_OP_FUNCTION_H

#include "timpi/timpi_config.h"

#ifdef TIMPI_HAVE_MPI
#  include "timpi/ignore_warnings.h"
#  include "mpi.h"
#  include "timpi/restore_warnings.h"
#endif // TIMPI_HAVE_MPI

// C++ includes
#include <functional>
#include <type_traits>



namespace TIMPI
{
#ifdef TIMPI_DEFAULT_QUADRUPLE_PRECISION
# ifdef TIMPI_HAVE_MPI
# define TIMPI_MPI_BINARY(funcname) \
inline void \
timpi_mpi_##funcname(void * a, void * b, int * len, MPI_Datatype *) \
{ \
  const int size = *len; \
 \
  Real *in = static_cast<Real*>(a); \
  Real *inout = static_cast<Real*>(b); \
  for (int i=0; i != size; ++i) \
    inout[i] = std::funcname(in[i],inout[i]); \
}

# define TIMPI_MPI_LOCATOR(funcname) \
inline void \
timpi_mpi_##funcname##_location(void * a, void * b, int * len, MPI_Datatype *) \
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


# define TIMPI_MPI_BINARY_FUNCTOR(funcname) \
inline void \
timpi_mpi_##funcname(void * a, void * b, int * len, MPI_Datatype *) \
{ \
  const int size = *len; \
 \
  Real *in = static_cast<Real*>(a); \
  Real *inout = static_cast<Real*>(b); \
  for (int i=0; i != size; ++i) \
    inout[i] = std::funcname<Real>()(in[i],inout[i]); \
}


TIMPI_MPI_BINARY(max)
TIMPI_MPI_BINARY(min)
TIMPI_MPI_LOCATOR(max)
TIMPI_MPI_LOCATOR(min)
TIMPI_MPI_BINARY_FUNCTOR(plus)
TIMPI_MPI_BINARY_FUNCTOR(multiplies)

# endif // TIMPI_HAVE_MPI
#endif // TIMPI_DEFAULT_QUADRUPLE_PRECISION


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

#ifdef TIMPI_HAVE_MPI

#define TIMPI_PARALLEL_INTEGER_OPS(cxxtype)          \
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

#define TIMPI_PARALLEL_FLOAT_OPS(cxxtype)            \
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

#define TIMPI_PARALLEL_INTEGER_OPS(cxxtype)  \
  template<>                                    \
  class OpFunction<cxxtype>                     \
  {                                             \
  }

#define TIMPI_PARALLEL_FLOAT_OPS(cxxtype)    \
  template<>                                    \
  class OpFunction<cxxtype>                     \
  {                                             \
  }

#endif

TIMPI_PARALLEL_INTEGER_OPS(char);
TIMPI_PARALLEL_INTEGER_OPS(signed char);
TIMPI_PARALLEL_INTEGER_OPS(unsigned char);
TIMPI_PARALLEL_INTEGER_OPS(short int);
TIMPI_PARALLEL_INTEGER_OPS(unsigned short int);
TIMPI_PARALLEL_INTEGER_OPS(int);
TIMPI_PARALLEL_INTEGER_OPS(unsigned int);
TIMPI_PARALLEL_INTEGER_OPS(long);
TIMPI_PARALLEL_INTEGER_OPS(long long);
TIMPI_PARALLEL_INTEGER_OPS(unsigned long);
TIMPI_PARALLEL_INTEGER_OPS(unsigned long long);

TIMPI_PARALLEL_FLOAT_OPS(float);
TIMPI_PARALLEL_FLOAT_OPS(double);
TIMPI_PARALLEL_FLOAT_OPS(long double);

#ifdef TIMPI_DEFAULT_QUADRUPLE_PRECISION
# ifdef TIMPI_HAVE_MPI

// FIXME - we'll still need to free these to keep valgrind happy
#define TIMPI_MPI_OPFUNCTION(mpiname, funcname) \
  static MPI_Op mpiname() { \
    static MPI_Op TIMPI_MPI_##mpiname = MPI_OP_NULL; \
    if (TIMPI_MPI_##mpiname == MPI_OP_NULL) \
      { \
        timpi_call_mpi \
          (MPI_Op_create(timpi_mpi_##funcname, true, &TIMPI_MPI_##mpiname)); \
      } \
    return TIMPI_MPI_##mpiname;  \
  }

  template<>
  class OpFunction<Real>
  {
  public:
    TIMPI_MPI_OPFUNCTION(max, max)
    TIMPI_MPI_OPFUNCTION(min, min)
    TIMPI_MPI_OPFUNCTION(sum, plus)
    TIMPI_MPI_OPFUNCTION(product, multiplies)

    TIMPI_MPI_OPFUNCTION(max_location, max_location)
    TIMPI_MPI_OPFUNCTION(min_location, min_location)
  };

# else
  TIMPI_PARALLEL_FLOAT_OPS(Real);
# endif
#endif // TIMPI_DEFAULT_QUADRUPLE_PRECISION

} // namespace TIMPI

#endif // TIMPI_OP_FUNCTION_H
