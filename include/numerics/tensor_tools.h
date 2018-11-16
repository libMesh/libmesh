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



#ifndef LIBMESH_TENSOR_TOOLS_H
#define LIBMESH_TENSOR_TOOLS_H

// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/compare_types.h"
#ifdef LIBMESH_HAVE_METAPHYSICL
#include "metaphysicl/dualnumber.h"
#include "metaphysicl/dynamicsparsenumberbase.h"
#endif

namespace libMesh
{
// Forward declarations
template <typename T> class TypeVector;
template <typename T> class VectorValue;
template <typename T> class TypeTensor;
template <typename T> class TensorValue;
template <unsigned int N, typename T> class TypeNTensor;

namespace TensorTools
{
// Any tensor-rank-independent code will need to include
// tensor_tools.h, so we define a product/dot-product here, starting
// with the generic case to apply to scalars.
// Vector specializations will follow.

template <typename T, typename T2>
inline
typename boostcopy::enable_if_c<ScalarTraits<T>::value && ScalarTraits<T2>::value,
                                typename CompareTypes<T, T2>::supertype>::type
inner_product(const T & a, const T2& b)
{ return a * b; }

template <typename T, typename T2>
inline
typename CompareTypes<T, T2>::supertype
inner_product(const TypeVector<T> & a, const TypeVector<T2> & b)
{ return a * b; }

template <typename T, typename T2>
inline
typename CompareTypes<T, T2>::supertype
inner_product(const TypeTensor<T> & a, const TypeTensor<T2> & b)
{ return a.contract(b); }

template <unsigned int N, typename T, typename T2>
inline
typename CompareTypes<T, T2>::supertype
inner_product(const TypeNTensor<N,T> & a, const TypeNTensor<N,T2> & b)
{ return a.contract(b); }

template<typename T>
inline
T norm_sq(T a) { return a*a; }

template<typename T>
inline
T norm_sq(std::complex<T> a) { return std::norm(a); }

#ifdef LIBMESH_HAVE_METAPHYSICL
template<typename T, typename D>
inline
auto norm_sq(const MetaPhysicL::DualNumber<T, D> & a) -> decltype(std::norm(a))
{ return std::norm(a); }
#endif

template <typename T>
inline
auto norm_sq(const TypeVector<T> & a) -> decltype(std::norm(T()))
{return a.norm_sq();}

template <typename T>
inline
auto norm_sq(const VectorValue<T> & a) -> decltype(std::norm(T()))
{return a.norm_sq();}

// Any tensor-rank-independent code will need to include
// tensor_tools.h, so we define rank-increasing and real-to-number type
// conversion functions here, starting with the generic case to apply
// to scalars.
// Tensor(and higher?) specializations will go in the tensor
// header(s).
template <typename T>
struct IncrementRank
{
  typedef VectorValue<T> type;
};

template <typename T>
struct IncrementRank<VectorValue<T>>
{
  typedef TensorValue<T> type;
};


template <typename T>
struct IncrementRank<TypeVector<T>>
{
  typedef TensorValue<T> type;
};

template <typename T>
struct IncrementRank<TypeTensor<T>>
{
  typedef TypeNTensor<3,T> type;
};


template <typename T>
struct IncrementRank<TensorValue<T>>
{
  typedef TypeNTensor<3,T> type;
};

template <unsigned int N, typename T>
struct IncrementRank<TypeNTensor<N,T>>
{
  typedef TypeNTensor<N+1,T> type;
};


// Also need rank-decreasing case
template <typename T>
struct DecrementRank
{
  // The default case is typically an error, but for simpler
  // templated code we need it to be compatible with Number
  // operations...
  typedef T type;
};

template <typename T>
struct DecrementRank<VectorValue<T>>
{
  typedef T type;
};

template <typename T>
struct DecrementRank<TypeVector<T>>
{
  typedef T type;
};

template <typename T>
struct DecrementRank<TensorValue<T>>
{
  typedef VectorValue<T> type;
};

template <typename T>
struct DecrementRank<TypeTensor<T>>
{
  typedef VectorValue<T> type;
};

template <unsigned int N, typename T>
struct DecrementRank<TypeNTensor<N,T>>
{
  typedef TypeNTensor<N-1,T> type;
};

// Handle the complex-valued case
template <typename T>
struct MakeNumber
{
#ifdef LIBMESH_USE_COMPLEX_NUMBERS
  typedef std::complex<T> type;
#else
  typedef T type;
#endif
};

template <typename T>
struct MakeNumber<std::complex<T>>
{
  // Compile-time error: we shouldn't need to make numbers out of
  // numbers
  //typedef std::complex<T> type;
};


template <typename T>
struct MakeNumber<TypeVector<T>>
{
  typedef TypeVector<typename MakeNumber<T>::type > type;
};

template <typename T>
struct MakeNumber<VectorValue<T>>
{
  typedef VectorValue<typename MakeNumber<T>::type > type;
};

template <typename T>
struct MakeNumber<TypeTensor<T>>
{
  typedef TypeTensor<typename MakeNumber<T>::type> type;
};

template <typename T>
struct MakeNumber<TensorValue<T>>
{
  typedef TypeTensor<typename MakeNumber<T>::type> type;
};

template <unsigned int N, typename T>
struct MakeNumber<TypeNTensor<N,T>>
{
#ifdef LIBMESH_USE_COMPLEX_NUMBERS
  typedef TypeNTensor<N,std::complex<T>> type;
#else
  typedef TypeNTensor<N,T> type;
#endif
};

// A utility for determining real-valued (e.g. shape function)
// types from corresponding complex-valued types
template <typename T>
struct MakeReal
{
  typedef T type;
};

template <typename T>
struct MakeReal<std::complex<T>>
{
  typedef T type;
};

template <typename T>
struct MakeReal<TypeVector<T>>
{
  typedef TypeVector<typename MakeReal<T>::type > type;
};

template <typename T>
struct MakeReal<VectorValue<T>>
{
  typedef VectorValue<typename MakeReal<T>::type > type;
};

template <typename T>
struct MakeReal<TypeTensor<T>>
{
  typedef TypeTensor<typename MakeReal<T>::type> type;
};

template <typename T>
struct MakeReal<TensorValue<T>>
{
  typedef TypeTensor<typename MakeReal<T>::type> type;
};

template <unsigned int N, typename T>
struct MakeReal<TypeNTensor<N,T>>
{
  typedef TypeNTensor<N,typename MakeReal<T>::type> type;
};

// Needed for ExactSolution to compile
Number curl_from_grad( const VectorValue<Number> & );

//! Computes the curl of a vector given the gradient of that vector
VectorValue<Number> curl_from_grad( const TensorValue<Number> & grad );

/*! Place holder needed for ExactSolution to compile. Will compute the
  curl of a tensor given the gradient of that tensor. */
TensorValue<Number> curl_from_grad( const TypeNTensor<3,Number> & grad );

//! Dummy. Divergence of a scalar not defined, but is needed for ExactSolution to compile
Number div_from_grad( const VectorValue<Number> & grad );

//! Computes the divergence of a vector given the gradient of that vector
Number div_from_grad( const TensorValue<Number> & grad );

/*! Place holder needed for ExactSolution to compile. Will compute the
  divergence of a tensor given the gradient of that tensor. */
VectorValue<Number> div_from_grad( const TypeNTensor<3,Number> & grad );

}//namespace TensorTools

}//namespace libMesh

#endif // LIBMESH_TENSOR_TOOLS_H
