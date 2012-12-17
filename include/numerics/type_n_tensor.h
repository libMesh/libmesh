// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_TYPE_N_TENSOR_H
#define LIBMESH_TYPE_N_TENSOR_H

// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/type_vector.h"

// C++ includes
#include <cstdlib> // *must* precede <cmath> for proper std:abs() on PGI, Sun Studio CC
#include <cmath>

namespace libMesh
{


/**
 * This class will eventually define a rank-N tensor in \p LIBMESH_DIM
 * dimensional space of type T.
 * Right now it defines a shim to allow for rank-independent code to
 * compile (but not give correct results) in the case of vector-valued
 * elements and second derivatives.
 *
 * \author Roy Stogner, 2012.
 */

template <unsigned int N, typename T>
class TypeNTensor
{
public:
  TypeNTensor () {}

  TypeNTensor (const T&) {}

  TypeNTensor (const TypeVector<T>&) {}

  TypeNTensor (const TypeTensor<T>&) {}

  operator TypeVector<T> () const { libmesh_error(); return 0; }
  operator VectorValue<T> () const { libmesh_error(); return 0; }

  operator TypeTensor<T> () const { libmesh_error(); return 0; }
  operator TensorValue<T> () const { libmesh_error(); return 0; }

  /**
   * Destructor.
   */
  ~TypeNTensor() {}

  /**
   * Return a proxy for the \f$ i^{th} \f$ slice of the tensor.
   */
  const TypeNTensor<N-1,T> slice (const unsigned int /*i*/) const
    { return TypeNTensor<N-1,T>(); }

  /**
   * Return a writeable proxy for the \f$ i^{th} \f$ slice of the tensor.
   */
  TypeNTensor<N-1,T> slice (const unsigned int /*i*/)
    { return TypeNTensor<N-1,T>(); }

  /**
   * Add two tensors.
   */
  template<typename T2>
  TypeNTensor<N,typename CompareTypes<T, T2>::supertype>
  operator + (const TypeNTensor<N,T2> &) const
    { return TypeNTensor<N,typename CompareTypes<T,T2>::supertype>(); }

  /**
   * Add to this tensor.
   */
  template<typename T2>
  const TypeNTensor<N,T> & operator += (const TypeNTensor<N,T2> &/*rhs*/)
    { return *this; }

  /**
   * Subtract two tensors.
   */
  template<typename T2>
  TypeNTensor<N,typename CompareTypes<T, T2>::supertype>
  operator - (const TypeNTensor<N,T2> &) const
    { return TypeNTensor<N,typename CompareTypes<T,T2>::supertype>(); }

  /**
   * Subtract from this tensor.
   */
  template<typename T2>
  const TypeNTensor<N,T> & operator -= (const TypeNTensor<N,T2> &)
    { return *this; }

  /**
   * Return the opposite of a tensor
   */
  TypeNTensor<N,T> operator - () const
    { return *this; }

  /**
   * Multiply a tensor by a number, i.e. scale.
   */
  template <typename Scalar>
  typename boostcopy::enable_if_c<
    ScalarTraits<Scalar>::value,
    TypeNTensor<N,typename CompareTypes<T, Scalar>::supertype> >::type
  operator * (const Scalar) const
    { return TypeNTensor<N,typename CompareTypes<T, Scalar>::supertype>(); }

  /**
   * Multiply this tensor by a number, i.e. scale.
   */
  template <typename Scalar>
  const TypeNTensor<N,T> & operator *= (const Scalar) { return *this; }

  /**
   * Divide a tensor by a number, i.e. scale.
   */
  template <typename Scalar>
  typename boostcopy::enable_if_c<
    ScalarTraits<Scalar>::value,
    TypeNTensor<N,typename CompareTypes<T, Scalar>::supertype> >::type
  operator / (const Scalar) const { return *this; }

  /**
   * Divide this tensor by a number, i.e. scale.
   */
  const TypeNTensor<N,T> & operator /= (const T) { return *this; }

  /**
   * Multiply 2 tensors together, i.e. dyadic product
   * sum_ij Aij*Bij.
   * The tensors may be of different types.
   */
  template <typename T2>
  typename CompareTypes<T,T2>::supertype
  contract (const TypeNTensor<N,T2> &) const { return 0; }

  /**
   * Returns the Frobenius norm of the tensor squared, i.e.  sum of the
   * element magnitudes squared.
   */
  Real size_sq() const { return 0.;}

  /**
   * @returns \p true if two tensors are equal valued.
   */
  bool operator == (const TypeNTensor<N,T>& /*rhs*/) const
    { return true; }

  /**
   * @returns \p true if this tensor is "less"
   * than another.  Useful for sorting.
   */
  bool operator < (const TypeNTensor<N,T>& /*rhs*/) const
    { return false; }

  /**
   * @returns \p true if this tensor is "greater"
   * than another.
   */
  bool operator > (const TypeNTensor<N,T>& /*rhs*/) const
    { return false; }

  /**
   * Formatted print, by default to \p libMesh::out.
   */
  void print(std::ostream& /*os = libMesh::out*/) const {}

  /**
   * Formatted print as above but allows you to do
   * std::cout << t << std::endl;
   */
  friend std::ostream& operator << (std::ostream& os,
                                    const TypeNTensor<N,T>& t)
  {
    t.print(os);
    return os;
  }
};


} // namespace libMesh

#endif // LIBMESH_TYPE_N_TENSOR_H
