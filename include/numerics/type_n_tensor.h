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



#ifndef LIBMESH_TYPE_N_TENSOR_H
#define LIBMESH_TYPE_N_TENSOR_H

// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/type_vector.h"
#include "libmesh/tuple_of.h"

// C++ includes
#include <cstdlib> // *must* precede <cmath> for proper std:abs() on PGI, Sun Studio CC
#include <cmath>
#include <vector>

namespace libMesh
{

/**
 * This class will eventually define a rank-N tensor in \p LIBMESH_DIM
 * dimensional space of type T.
 *
 * Right now it defines a shim to allow for rank-independent code to
 * compile (but not give correct results) in the case of vector-valued
 * elements and second derivatives.
 *
 * \author Roy Stogner
 * \date 2012
 */
template <unsigned int N, typename T>
class TypeNTensor
{
public:
  /**
   * Helper typedef for generic index programming
   */
  typedef tuple_of<N, unsigned int> index_type;

  TypeNTensor () : _coords(std::vector<T>(int_pow(LIBMESH_DIM, N))) {}

  TypeNTensor (const T &) : _coords(std::vector<T>(int_pow(LIBMESH_DIM, N))) {}

  TypeNTensor (const TypeVector<T> &) : _coords(std::vector<T>(int_pow(LIBMESH_DIM, N))) {}

  TypeNTensor (const TypeTensor<T> &) : _coords(std::vector<T>(int_pow(LIBMESH_DIM, N))) {}

  operator TypeVector<T> () const { libmesh_not_implemented(); return 0; }
  operator VectorValue<T> () const { libmesh_not_implemented(); return 0; }

  operator TypeTensor<T> () const { libmesh_not_implemented(); return 0; }
  operator TensorValue<T> () const { libmesh_not_implemented(); return 0; }

  /**
   * Destructor.
   */
  ~TypeNTensor() {}

  /**
   * \returns A proxy for the \f$ i^{th} \f$ slice of the tensor.
   */
  const TypeNTensor<N-1,T> slice (const unsigned int /*i*/) const
  { return TypeNTensor<N-1,T>(); }

  /**
   * \returns A writable proxy for the \f$ i^{th} \f$ slice of the tensor.
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
   * \returns The negative of a tensor.
   */
  TypeNTensor<N,T> operator - () const
  { return *this; }

  /**
   * Multiply every entry of a tensor by a number.
   */
  template <typename Scalar>
  typename boostcopy::enable_if_c<
    ScalarTraits<Scalar>::value,
    TypeNTensor<N,typename CompareTypes<T, Scalar>::supertype>>::type
  operator * (const Scalar) const
  { return TypeNTensor<N,typename CompareTypes<T, Scalar>::supertype>(); }

  /**
   * Multiply every entry of this tensor by a number.
   */
  template <typename Scalar>
  const TypeNTensor<N,T> & operator *= (const Scalar) { return *this; }

  /**
   * Divide every entry of a tensor by a number.
   */
  template <typename Scalar>
  typename boostcopy::enable_if_c<
    ScalarTraits<Scalar>::value,
    TypeNTensor<N,typename CompareTypes<T, Scalar>::supertype>>::type
  operator / (const Scalar) const { return *this; }

  /**
   * Divide every entry of this tensor by a number.
   */
  const TypeNTensor<N,T> & operator /= (const T) { return *this; }

  /**
   * Multiply 2 tensors together to return a scalar, i.e.
   * \f$ \sum_{ij} A_{ij} B_{ij} \f$
   * The tensors may contain different numeric types.
   * Also known as the "double inner product" or "double dot product"
   * of tensors.
   *
   * \returns The scalar-valued result, this tensor is unchanged.
   */
  template <typename T2>
  typename CompareTypes<T,T2>::supertype
  contract (const TypeNTensor<N,T2> &) const { return 0; }

  /**
   * \returns The Frobenius norm of the tensor squared, i.e. the sum of the
   * entry magnitudes squared.
   *
   * \deprecated Use the norm_sq() function instead.
   */
#ifdef LIBMESH_ENABLE_DEPRECATED
  auto size_sq() const -> decltype(std::norm(T())) { libmesh_deprecated(); return 0.;}
#endif

  /**
   * \returns The Frobenius norm of the tensor squared, i.e. the sum of the
   * entry magnitudes squared.
   */
  auto norm_sq() const -> decltype(std::norm(T())) { return 0.;}

  /**
   * \returns \p true if two tensors are equal, \p false otherwise.
   */
  bool operator == (const TypeNTensor<N,T> & /*rhs*/) const
  { return true; }

  /**
   * \returns \p true if this tensor is "less" than another.
   *
   * Useful for sorting.
   */
  bool operator < (const TypeNTensor<N,T> & /*rhs*/) const
  { return false; }

  /**
   * \returns \p true if this tensor is "greater" than another.
   */
  bool operator > (const TypeNTensor<N,T> & /*rhs*/) const
  { return false; }

  /**
   * Do a formatted print of this tensor to a stream which defaults to
   * \p libMesh::out.
   */
  void print(std::ostream & /*os = libMesh::out*/) const {}

  /**
   * Does a formatted print (as above) but supports the syntax:
   * \code
   * std::cout << t << std::endl;
   * \endcode
   */
  friend std::ostream & operator << (std::ostream & os,
                                     const TypeNTensor<N,T> & t)
  {
    t.print(os);
    return os;
  }

  /**
   * Add a scaled type N tensor to this type N tensor without creating a temporary.
   */
  template<typename T2>
  void add_scaled (const TypeNTensor<N, T2> &, const T &);

  /**
   * The coordinates of the \p TypeNTensor
   */
  std::vector<T> _coords;

private:
  static constexpr int int_pow(int b, int e)
  {
    return (e == 0) ? 1 : b * int_pow(b, e - 1);
  }
};


template<unsigned int N, typename T>
template<typename T2>
inline
void TypeNTensor<N, T>::add_scaled (const TypeNTensor<N, T2> & p, const T & factor)
{
  unsigned int size = int_pow(LIBMESH_DIM, N);
  for (unsigned int i = 0; i < size ; i++)
    _coords[i] += factor*p._coords[i];
}

template <unsigned int N, typename T, typename Scalar>
typename boostcopy::enable_if_c<
  ScalarTraits<Scalar>::value,
  TypeNTensor<N,typename CompareTypes<Scalar, T>::supertype>>::type
operator * (const Scalar &, const TypeNTensor<N, T> &)
{
  libmesh_not_implemented();
  return TypeNTensor<N,typename CompareTypes<Scalar, T>::supertype>();
}

template <unsigned int N, typename T, typename Scalar>
typename boostcopy::enable_if_c<
  ScalarTraits<Scalar>::value,
  TypeNTensor<N,typename CompareTypes<Scalar, T>::supertype>>::type
operator / (const Scalar &, const TypeNTensor<N, T> &)
{
  libmesh_not_implemented();
  return TypeNTensor<N,typename CompareTypes<Scalar, T>::supertype>();
}


} // namespace libMesh


#endif // LIBMESH_TYPE_N_TENSOR_H
