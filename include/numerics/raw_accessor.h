// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_RAW_ACCESSOR_H
#define LIBMESH_RAW_ACCESSOR_H

// Local includes
#include "libmesh/libmesh_common.h"

#include "libmesh/tensor_value.h"
#include "libmesh/vector_value.h"
#include "libmesh/type_n_tensor.h"

namespace libMesh
{

/**
 * What underlying data type would we need to access in each field?
 */
template <typename FieldType>
struct RawFieldType {};

template <>
struct RawFieldType<Number>
{
  typedef Number type;
};

template <>
struct RawFieldType<Gradient>
{
  typedef Number type;
};

template <>
struct RawFieldType<Tensor>
{
  typedef Number type;
};

template<>
struct RawFieldType<TypeNTensor<3, Number> >
{
  typedef Number type;
};

#ifdef LIBMESH_USE_COMPLEX_NUMBERS
template <>
struct RawFieldType<Real>
{
  typedef Real type;
};

template <>
struct RawFieldType<RealGradient>
{
  typedef Real type;
};

template <>
struct RawFieldType<RealTensor>
{
  typedef Real type;
};

template<>
struct RawFieldType<TypeNTensor<3, Real> >
{
  typedef Real type;
};
#endif

/**
 * This class provides single index access to FieldType (i.e. Number, Gradient, Tensor, etc.).
 */
template <typename FieldType>
class RawAccessor
{
public:

  RawAccessor(FieldType & data, const unsigned int dim)
    : _data(data),
      _dim(dim)
  {}

  ~RawAccessor(){}

  typename RawFieldType<FieldType>::type & operator()( unsigned int i );
  const typename RawFieldType<FieldType>::type & operator()( unsigned int i ) const;

private:
  RawAccessor();

  FieldType & _data;
  const unsigned int _dim;
};

// Specialize for specific cases
template<>
inline
Number & RawAccessor<Number>::operator()( unsigned int libmesh_dbg_var(i) )
{
  libmesh_assert_equal_to (i, 0);
  return this->_data;
}

template<>
inline
Number & RawAccessor<Gradient>::operator()( unsigned int i )
{
  libmesh_assert_less (i, this->_dim);
  return this->_data(i);
}

template<>
inline
Number & RawAccessor<Tensor>::operator()( unsigned int k )
{
  libmesh_assert_less (k, this->_dim*this->_dim);

  // For tensors, each row is filled first, i.e. for 2-D
  // [ 0 1; 2 3]
  // Thus, k(i,j) = j + i*dim
  unsigned int ii = k/_dim;
  unsigned int jj = k - ii*_dim;

  return this->_data(ii,jj);
}

/**
 * Stub implementations for stub TypeNTensor object
 */
template <unsigned int N, typename ScalarType>
class RawAccessor<TypeNTensor<N, ScalarType> >
{
public:

  typedef TypeNTensor<N, ScalarType> FieldType;

  RawAccessor(FieldType & data, const unsigned int dim)
    : _data(data),
      _dim(dim)
  {}

  ~RawAccessor(){}

  typename RawFieldType<FieldType>::type & operator()( unsigned int /*i*/ )
  { return dummy; }

  const typename RawFieldType<FieldType>::type & operator()( unsigned int /*i*/ ) const
  { return dummy; }

private:
  RawAccessor();

  ScalarType dummy;

  FieldType & _data;
  const unsigned int _dim;
};

#ifdef LIBMESH_USE_COMPLEX_NUMBERS
template<>
inline
Real & RawAccessor<Real>::operator()( unsigned int libmesh_dbg_var(i) )
{
  libmesh_assert_equal_to (i, 0);
  return this->_data;
}

template<>
inline
Real & RawAccessor<RealGradient>::operator()( unsigned int i )
{
  libmesh_assert_less (i, this->_dim);
  return this->_data(i);
}

template<>
inline
Real & RawAccessor<RealTensor>::operator()( unsigned int k )
{
  libmesh_assert_less (k, this->_dim*this->_dim);

  // For tensors, each row is filled first, i.e. for 2-D
  // [ 0 1; 2 3]
  // Thus, k(i,j) = i + j*dim
  unsigned int jj = k/_dim;
  unsigned int ii = k - jj*_dim;

  return this->_data(ii,jj);
}

#endif

}
#endif // LIBMESH_RAW_ACCESSOR_H
