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


#ifndef TIMPI_STANDARD_TYPE_H
#define TIMPI_STANDARD_TYPE_H

// TIMPI includes
#include "timpi/data_type.h"
#include "timpi/timpi_config.h"

// C/C++ includes
#ifdef TIMPI_HAVE_MPI
#  include "timpi/ignore_warnings.h"
#  include "mpi.h"
#  include "timpi/restore_warnings.h"
#endif // TIMPI_HAVE_MPI

#include <array>
#include <complex>
#include <memory>
#include <numeric>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace TIMPI
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
  StandardType(const T * example = nullptr);
};



// ------------------------------------------------------------
// Declare StandardType specializations for C++ built-in types

#ifdef TIMPI_HAVE_MPI

#define TIMPI_STANDARD_TYPE(cxxtype,mpitype)                      \
  template<>                                                         \
  class StandardType<cxxtype> : public DataType                      \
  {                                                                  \
  public:                                                            \
    explicit                                                         \
      StandardType(const cxxtype * = nullptr) : DataType(mpitype) {} \
  }

#else

#define TIMPI_STANDARD_TYPE(cxxtype,mpitype)               \
  template<>                                                  \
  class StandardType<cxxtype> : public DataType               \
  {                                                           \
  public:                                                     \
    explicit                                                  \
      StandardType(const cxxtype * = nullptr) : DataType() {} \
  }

#endif

TIMPI_STANDARD_TYPE(char,MPI_CHAR);
TIMPI_STANDARD_TYPE(signed char,MPI_SIGNED_CHAR);
TIMPI_STANDARD_TYPE(unsigned char,MPI_UNSIGNED_CHAR);
TIMPI_STANDARD_TYPE(short int,MPI_SHORT);
TIMPI_STANDARD_TYPE(unsigned short int,MPI_UNSIGNED_SHORT);
TIMPI_STANDARD_TYPE(int,MPI_INT);
TIMPI_STANDARD_TYPE(unsigned int,MPI_UNSIGNED);
TIMPI_STANDARD_TYPE(long,MPI_LONG);
TIMPI_STANDARD_TYPE(long long,MPI_LONG_LONG_INT);
TIMPI_STANDARD_TYPE(unsigned long,MPI_UNSIGNED_LONG);
TIMPI_STANDARD_TYPE(unsigned long long,MPI_UNSIGNED_LONG_LONG);
TIMPI_STANDARD_TYPE(float,MPI_FLOAT);
TIMPI_STANDARD_TYPE(double,MPI_DOUBLE);
TIMPI_STANDARD_TYPE(long double,MPI_LONG_DOUBLE);

// Quad and float128 types aren't standard C++, so only work with them
// if configure and PETSc encapsulated the non-standard issues.
#ifdef TIMPI_HAVE_MPI
# ifdef TIMPI_DEFAULT_QUADRUPLE_PRECISION
  template<>
  class StandardType<Real> : public DataType
  {
  public:
    explicit
      StandardType(const Real * = nullptr) {
        timpi_call_mpi(MPI_Type_contiguous(2, MPI_DOUBLE, &_datatype));
        this->commit();
      }

    StandardType(const StandardType<Real> & t)
      : DataType()
    {
      timpi_call_mpi (MPI_Type_dup (t._datatype, &_datatype));
    }

    ~StandardType() {
      this->free();
    }
  };
# endif
#else
# ifdef TIMPI_DEFAULT_QUADRUPLE_PRECISION
  TIMPI_STANDARD_TYPE(Real,);
# endif
#endif

template<typename T1, typename T2>
class StandardType<std::pair<T1, T2>> : public DataType
{
public:
  explicit
  StandardType(const std::pair<T1, T2> * example = nullptr) {
    // We need an example for MPI_Address to use
    static const std::pair<T1, T2> p;
    if (!example)
      example = &p;

#ifdef TIMPI_HAVE_MPI

    // Get the sub-data-types, and make sure they live long enough
    // to construct the derived type
    StandardType<T1> d1(&example->first);
    StandardType<T2> d2(&example->second);

    MPI_Datatype types[] = { (data_type)d1, (data_type)d2 };
    int blocklengths[] = {1,1};
    MPI_Aint displs[2], start;

    timpi_call_mpi
      (MPI_Get_address (const_cast<std::pair<T1,T2> *>(example),
                        &start));
    timpi_call_mpi
      (MPI_Get_address (const_cast<T1*>(&example->first),
                        &displs[0]));
    timpi_call_mpi
      (MPI_Get_address (const_cast<T2*>(&example->second),
                        &displs[1]));
    displs[0] -= start;
    displs[1] -= start;

    // create a prototype structure
    MPI_Datatype tmptype;
    timpi_call_mpi
      (MPI_Type_create_struct (2, blocklengths, displs, types,
                               &tmptype));
    timpi_call_mpi
      (MPI_Type_commit (&tmptype));

    // resize the structure type to account for padding, if any
    timpi_call_mpi
      (MPI_Type_create_resized (tmptype, 0,
                                sizeof(std::pair<T1,T2>),
                                &_datatype));
    timpi_call_mpi
      (MPI_Type_free (&tmptype));

    this->commit();

#endif // TIMPI_HAVE_MPI

  }

  StandardType(const StandardType<std::pair<T1, T2>> & timpi_mpi_var(t))
  {
    timpi_call_mpi
      (MPI_Type_dup (t._datatype, &_datatype));
  }

  ~StandardType() { this->free(); }
};


// Helper functions for creating type/displacement arrays for tuples
//
// These are classes since we can't partially specialize functions
template<std::size_t n_minus_i>
struct BuildStandardTypeVector
{
  template<typename... Types>
  static void build(std::vector<std::unique_ptr<DataType>> & out_vec,
                    const std::tuple<Types...> & example);
};

template <>
struct BuildStandardTypeVector<0>
{
  template<typename... Types>
  static void build(std::vector<std::unique_ptr<DataType>> & /*out_vec*/,
                    const std::tuple<Types...> & /*example*/) {}
};

template<std::size_t n_minus_i>
template<typename... Types>
void BuildStandardTypeVector<n_minus_i>::build
  (std::vector<std::unique_ptr<DataType>> & out_vec,
   const std::tuple<Types...> & example)
{
  typedef typename
    std::tuple_element<sizeof...(Types)-n_minus_i, std::tuple<Types...>>::type
    ith_type;

  out_vec.emplace_back
    (new StandardType<ith_type>
     (&std::get<sizeof...(Types)-n_minus_i>(example)));

  BuildStandardTypeVector<n_minus_i-1>::build(out_vec, example);
}


template<std::size_t n_minus_i>
struct FillDisplacementArray
{
  template <typename OutArray, class... Types>
  static void fill(OutArray & out,
                   const std::tuple<Types...> & example);
};

template<>
struct FillDisplacementArray<0>
{
  template <typename OutArray, typename... Types>
  static void fill(OutArray & /*out*/,
                   const std::tuple<Types...> & /*example*/) {}
};


template<std::size_t n_minus_i>
template<typename OutArray, typename... Types>
void FillDisplacementArray<n_minus_i>::fill
  (OutArray & out_vec,
   const std::tuple<Types...> & example)
{
  timpi_call_mpi
    (MPI_Get_address
      (&std::get<sizeof...(Types)-n_minus_i>(example),
       &out_vec[sizeof...(Types)-n_minus_i]));

  FillDisplacementArray<n_minus_i-1>::fill(out_vec, example);
}


template<typename... Types>
class StandardType<std::tuple<Types...>> : public DataType
{
public:
  explicit
  StandardType(const std::tuple<Types...> * example = nullptr) {
    // We need an example for MPI_Address to use
    static const std::tuple<Types...> t;
    if (!example)
      example = &t;

#ifdef TIMPI_HAVE_MPI
    MPI_Aint start;

    timpi_call_mpi
      (MPI_Get_address (example, &start));

    const std::size_t tuplesize = sizeof...(Types);

    std::vector<std::unique_ptr<DataType>> subtypes;
    BuildStandardTypeVector<sizeof...(Types)>::build(subtypes, *example);

    std::array<MPI_Aint, sizeof...(Types)> displs;
    FillDisplacementArray<sizeof...(Types)>::fill(displs, *example);

    std::array<MPI_Datatype, sizeof...(Types)> types;
    std::array<int, sizeof...(Types)> blocklengths;

    for (std::size_t i = 0; i != tuplesize; ++i)
    {
      displs[i] -= start;
      types[i] = (data_type)(*subtypes[i]);
      blocklengths[i] = 1;
    }

    // create a prototype structure
    MPI_Datatype tmptype;
    timpi_call_mpi
      (MPI_Type_create_struct (tuplesize, blocklengths.data(), displs.data(), types.data(),
                               &tmptype));
    timpi_call_mpi
      (MPI_Type_commit (&tmptype));

    // resize the structure type to account for padding, if any
    timpi_call_mpi
      (MPI_Type_create_resized (tmptype, 0,
                                sizeof(std::tuple<Types...>),
                                &_datatype));
    timpi_call_mpi
      (MPI_Type_free (&tmptype));

    this->commit();

#endif // TIMPI_HAVE_MPI

  }

  StandardType(const StandardType<std::tuple<Types...>> & timpi_mpi_var(t))
  {
    timpi_call_mpi
      (MPI_Type_dup (t._datatype, &_datatype));
  }

  ~StandardType() { this->free(); }
};


template<typename T>
class StandardType<std::complex<T>> : public DataType
{
public:
  explicit
  StandardType(const std::complex<T> * /*example*/ = nullptr) :
    DataType(StandardType<T>(nullptr), 2) {}

  ~StandardType() { this->free(); }
};

} // namespace TIMPI

#endif // TIMPI_STANDARD_TYPE_H
