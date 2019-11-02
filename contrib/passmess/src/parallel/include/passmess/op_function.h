// The PassMess Message-Passing Parallelism Library.
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

#include "passmess_config.h"

#ifdef PASSMESS_HAVE_MPI
#  include "passmess/ignore_warnings.h"
#  include "mpi.h"
#  include "passmess/restore_warnings.h"
#endif // PASSMESS_HAVE_MPI

// C++ includes
#include <functional>
#include <type_traits>



namespace PassMess
{
#ifdef PASSMESS_DEFAULT_QUADRUPLE_PRECISION
# ifdef PASSMESS_HAVE_MPI
# define PASSMESS_MPI_BINARY(funcname) \
inline void \
passmess_mpi_##funcname(void * a, void * b, int * len, MPI_Datatype *) \
{ \
  const int size = *len; \
 \
  Real *in = static_cast<Real*>(a); \
  Real *inout = static_cast<Real*>(b); \
  for (int i=0; i != size; ++i) \
    inout[i] = std::funcname(in[i],inout[i]); \
}

# define PASSMESS_MPI_LOCATOR(funcname) \
inline void \
passmess_mpi_##funcname##_location(void * a, void * b, int * len, MPI_Datatype *) \
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


# define PASSMESS_MPI_BINARY_FUNCTOR(funcname) \
inline void \
passmess_mpi_##funcname(void * a, void * b, int * len, MPI_Datatype *) \
{ \
  const int size = *len; \
 \
  Real *in = static_cast<Real*>(a); \
  Real *inout = static_cast<Real*>(b); \
  for (int i=0; i != size; ++i) \
    inout[i] = std::funcname<Real>()(in[i],inout[i]); \
}


PASSMESS_MPI_BINARY(max)
PASSMESS_MPI_BINARY(min)
PASSMESS_MPI_LOCATOR(max)
PASSMESS_MPI_LOCATOR(min)
PASSMESS_MPI_BINARY_FUNCTOR(plus)
PASSMESS_MPI_BINARY_FUNCTOR(multiplies)

# endif // PASSMESS_HAVE_MPI
#endif // PASSMESS_DEFAULT_QUADRUPLE_PRECISION


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

#ifdef PASSMESS_HAVE_MPI

#define PASSMESS_PARALLEL_INTEGER_OPS(cxxtype)          \
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

#define PASSMESS_PARALLEL_FLOAT_OPS(cxxtype)            \
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

#define PASSMESS_PARALLEL_INTEGER_OPS(cxxtype)  \
  template<>                                    \
  class OpFunction<cxxtype>                     \
  {                                             \
  }

#define PASSMESS_PARALLEL_FLOAT_OPS(cxxtype)    \
  template<>                                    \
  class OpFunction<cxxtype>                     \
  {                                             \
  }

#endif

PASSMESS_PARALLEL_INTEGER_OPS(char);
PASSMESS_PARALLEL_INTEGER_OPS(signed char);
PASSMESS_PARALLEL_INTEGER_OPS(unsigned char);
PASSMESS_PARALLEL_INTEGER_OPS(short int);
PASSMESS_PARALLEL_INTEGER_OPS(unsigned short int);
PASSMESS_PARALLEL_INTEGER_OPS(int);
PASSMESS_PARALLEL_INTEGER_OPS(unsigned int);
PASSMESS_PARALLEL_INTEGER_OPS(long);
PASSMESS_PARALLEL_INTEGER_OPS(long long);
PASSMESS_PARALLEL_INTEGER_OPS(unsigned long);
PASSMESS_PARALLEL_INTEGER_OPS(unsigned long long);

PASSMESS_PARALLEL_FLOAT_OPS(float);
PASSMESS_PARALLEL_FLOAT_OPS(double);
PASSMESS_PARALLEL_FLOAT_OPS(long double);

#ifdef PASSMESS_DEFAULT_QUADRUPLE_PRECISION
# ifdef PASSMESS_HAVE_MPI

// FIXME - we'll still need to free these to keep valgrind happy
#define PASSMESS_MPI_OPFUNCTION(mpiname, funcname) \
  static MPI_Op mpiname() { \
    static MPI_Op PASSMESS_MPI_##mpiname = MPI_OP_NULL; \
    if (PASSMESS_MPI_##mpiname == MPI_OP_NULL) \
      { \
        passmess_call_mpi \
          (MPI_Op_create(passmess_mpi_##funcname, true, &PASSMESS_MPI_##mpiname)); \
      } \
    return PASSMESS_MPI_##mpiname;  \
  }

  template<>
  class OpFunction<Real>
  {
  public:
    PASSMESS_MPI_OPFUNCTION(max, max)
    PASSMESS_MPI_OPFUNCTION(min, min)
    PASSMESS_MPI_OPFUNCTION(sum, plus)
    PASSMESS_MPI_OPFUNCTION(product, multiplies)

    PASSMESS_MPI_OPFUNCTION(max_location, max_location)
    PASSMESS_MPI_OPFUNCTION(min_location, min_location)
  };

# else
  PASSMESS_PARALLEL_FLOAT_OPS(Real);
# endif
#endif // PASSMESS_DEFAULT_QUADRUPLE_PRECISION

} // namespace PassMess

#endif // PASSMESS_OP_FUNCTION_H
