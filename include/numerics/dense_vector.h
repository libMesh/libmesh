// The libMesh Finite Element Library.
// Copyright (C) 2002-2020 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_DENSE_VECTOR_H
#define LIBMESH_DENSE_VECTOR_H

// Local Includes
#include "libmesh/libmesh_common.h"
#include "libmesh/compare_types.h"
#include "libmesh/dense_vector_base.h"
#include "libmesh/tensor_tools.h"

#ifdef LIBMESH_HAVE_EIGEN
#include "libmesh/ignore_warnings.h"
#include <Eigen/Core>
#include "libmesh/restore_warnings.h"
#endif

#ifdef LIBMESH_HAVE_METAPHYSICL
#include "metaphysicl/raw_type.h"
#endif

// C++ includes
#include <vector>

namespace libMesh
{

/**
 * Defines a dense vector for use in Finite Element-type computations.
 * This class is to basically compliment the \p DenseMatrix class.  It
 * has additional capabilities over the \p std::vector that make it
 * useful for finite elements, particularly for systems of equations.
 * All overridden virtual functions are documented in dense_vector_base.h.
 *
 * \author Benjamin S. Kirk
 * \date 2003
 */
template<typename T>
class DenseVector : public DenseVectorBase<T>
{
public:

  /**
   * Constructor.  Creates a dense vector of dimension \p n.
   */
  explicit
  DenseVector(const unsigned int n=0);

  /**
   * Constructor.  Creates a dense vector of dimension \p n where all
   * entries have value \p val.
   */
  explicit
  DenseVector(const unsigned int n,
              const T & val);

  /**
   * Copy-constructor.
   */
  template <typename T2>
  DenseVector (const DenseVector<T2> & other_vector);

  /**
   * Copy-constructor, from a \p std::vector.
   */
  template <typename T2>
  DenseVector (const std::vector<T2> & other_vector);

  /**
   * The 5 special functions can be defaulted for this class, as it
   * does not manage any memory itself.
   */
  DenseVector (DenseVector &&) = default;
  DenseVector (const DenseVector &) = default;
  DenseVector & operator= (const DenseVector &) = default;
  DenseVector & operator= (DenseVector &&) = default;
  virtual ~DenseVector() = default;

  virtual unsigned int size() const override
  {
    return cast_int<unsigned int>(_val.size());
  }

  virtual bool empty() const override
  { return _val.empty(); }

  virtual void zero() override;

  /**
   * \returns Entry \p i of the vector as a const reference.
   */
  const T & operator() (const unsigned int i) const;

  /**
   * \returns Entry \p i of the vector as a writable reference.
   */
  T & operator() (const unsigned int i);

  virtual T el(const unsigned int i) const override
  { return (*this)(i); }

  virtual T & el(const unsigned int i) override
  { return (*this)(i); }

  /**
   * Assignment operator.
   *
   * \returns A reference to *this.
   */
  template <typename T2>
  DenseVector<T> & operator = (const DenseVector<T2> & other_vector);

  /**
   * STL-like swap method
   */
  void swap(DenseVector<T> & other_vector);

  /**
   * Resize the vector. Sets all elements to 0.
   */
  void resize (const unsigned int n);

  /**
   * Append additional entries to (resizing, but unchanging) the
   * vector.
   */
  template <typename T2>
  void append (const DenseVector<T2> & other_vector);

  /**
   * Multiplies every element in the vector by \p factor.
   */
  void scale (const T factor);

  /**
   * Multiplies every element in the vector by \p factor.
   *
   * \returns A reference to *this.
   */
  DenseVector<T> & operator*= (const T factor);

  /**
   * Adds \p factor times \p vec to this vector.
   * This should only work if T += T2 * T3 is valid C++ and
   * if T2 is scalar.  Return type is void
   *
   * \returns A reference to *this.
   */
  template <typename T2, typename T3>
  typename boostcopy::enable_if_c<
    ScalarTraits<T2>::value, void >::type
  add (const T2 factor,
       const DenseVector<T3> & vec);

  /**
   * \returns The dot product of *this with \p vec.
   *
   * In the complex-valued case, uses the complex conjugate of \p vec.
   */
  template <typename T2>
  typename CompareTypes<T, T2>::supertype dot (const DenseVector<T2> & vec) const;

  /**
   * \returns The dot product of *this with \p vec.
   *
   * In the complex-valued case, does not use the complex conjugate of
   * \p vec.
   */
  template <typename T2>
  typename CompareTypes<T, T2>::supertype indefinite_dot (const DenseVector<T2> & vec) const;

  /**
   * \returns \p true if \p vec is exactly equal to this vector, false otherwise.
   */
  template <typename T2>
  bool operator== (const DenseVector<T2> & vec) const;

  /**
   * \returns \p true if \p vec is not exactly equal to this vector, false otherwise.
   */
  template <typename T2>
  bool operator!= (const DenseVector<T2> & vec) const;

  /**
   * Adds \p vec to this vector.
   *
   * \returns A reference to *this.
   */
  template <typename T2>
  DenseVector<T> & operator+= (const DenseVector<T2> & vec);

  /**
   * Subtracts \p vec from this vector.
   *
   * \returns A reference to *this.
   */
  template <typename T2>
  DenseVector<T> & operator-= (const DenseVector<T2> & vec);

  /**
   * \returns The minimum entry of the vector, or the minimum real
   * part in the case of complex numbers.
   */
  Real min () const;

  /**
   * \returns The maximum entry of the vector, or the maximum real
   * part in the case of complex numbers.
   */
  Real max () const;

  /**
   * \returns The \f$l_1\f$-norm of the vector, i.e. the sum of the
   * absolute values of the entries.
   */
  Real l1_norm () const;

  /**
   * \returns The \f$l_2\f$-norm of the vector, i.e. the square root
   * of the sum of the squares of the entries.
   */
  Real l2_norm () const;

  /**
   * \returns The \f$l_\infty\f$-norm of the vector, i.e. the maximum
   * absolute value of the entries.
   */
  Real linfty_norm () const;

  /**
   * Puts the principal subvector of size \p sub_n (i.e. first sub_n
   * entries) into \p dest.
   */
  void get_principal_subvector (unsigned int sub_n, DenseVector<T> & dest) const;

  /**
   * \returns A reference to the underlying data storage vector.
   *
   * This should be used with caution (i.e. one should not change the
   * size of the vector, etc.) but is useful for interoperating with
   * low level BLAS routines which expect a simple array.
   */
  std::vector<T> & get_values() { return _val; }

  /**
   * \returns A constant reference to the underlying data storage vector.
   */
  const std::vector<T> & get_values() const { return _val; }

private:

  /**
   * The actual data values, stored as a 1D array.
   */
  std::vector<T> _val;
};



// ------------------------------------------------------------
// DenseVector member functions
template<typename T>
inline
DenseVector<T>::DenseVector(const unsigned int n) :
  _val (n, T{})
{
}


template<typename T>
inline
DenseVector<T>::DenseVector(const unsigned int n,
                            const T & val) :
  _val (n, val)
{
}



template<typename T>
template<typename T2>
inline
DenseVector<T>::DenseVector (const DenseVector<T2> & other_vector) :
  DenseVectorBase<T>()
{
  const std::vector<T2> & other_vals = other_vector.get_values();

  _val.clear();

  const int N = cast_int<int>(other_vals.size());
  _val.reserve(N);

  for (int i=0; i<N; i++)
    _val.push_back(other_vals[i]);
}



template<typename T>
template<typename T2>
inline
DenseVector<T>::DenseVector (const std::vector<T2> & other_vector) :
  _val(other_vector)
{
}





template<typename T>
template<typename T2>
inline
DenseVector<T> & DenseVector<T>::operator = (const DenseVector<T2> & other_vector)
{
  const std::vector<T2> & other_vals = other_vector.get_values();

  _val.clear();

  const int N = cast_int<int>(other_vals.size());
  _val.reserve(N);

  for (int i=0; i<N; i++)
    _val.push_back(other_vals[i]);

  return *this;
}



template<typename T>
inline
void DenseVector<T>::swap(DenseVector<T> & other_vector)
{
  _val.swap(other_vector._val);
}



template<typename T>
inline
void DenseVector<T>::resize(const unsigned int n)
{
  _val.resize(n);

  zero();
}



template<typename T>
template<typename T2>
inline
void DenseVector<T>::append (const DenseVector<T2> & other_vector)
{
  const std::vector<T2> & other_vals = other_vector.get_values();

  _val.reserve(this->size() + other_vals.size());
  _val.insert(_val.end(), other_vals.begin(), other_vals.end());
}



template<typename T>
inline
void DenseVector<T>::zero()
{
  std::fill (_val.begin(),
             _val.end(),
             T{});
}



template<typename T>
inline
const T & DenseVector<T>::operator () (const unsigned int i) const
{
  libmesh_assert_less (i, _val.size());

  return _val[i];
}



template<typename T>
inline
T & DenseVector<T>::operator () (const unsigned int i)
{
  libmesh_assert_less (i, _val.size());

  return _val[i];
}



template<typename T>
inline
void DenseVector<T>::scale (const T factor)
{
  const int N = cast_int<int>(_val.size());
  for (int i=0; i<N; i++)
    _val[i] *= factor;
}



template<typename T>
inline
DenseVector<T> & DenseVector<T>::operator*= (const T factor)
{
  this->scale(factor);
  return *this;
}



template<typename T>
template<typename T2, typename T3>
inline
typename boostcopy::enable_if_c<
  ScalarTraits<T2>::value, void >::type
DenseVector<T>::add (const T2 factor,
                     const DenseVector<T3> & vec)
{
  libmesh_assert_equal_to (this->size(), vec.size());

  const int N = cast_int<int>(_val.size());
  for (int i=0; i<N; i++)
    (*this)(i) += static_cast<T>(factor)*vec(i);
}



template<typename T>
template<typename T2>
inline
typename CompareTypes<T, T2>::supertype DenseVector<T>::dot (const DenseVector<T2> & vec) const
{
  if (!_val.size())
    return 0.;

  libmesh_assert_equal_to (this->size(), vec.size());

#ifdef LIBMESH_HAVE_EIGEN
  // We reverse the order of the arguments to dot() here since
  // the convention in Eigen is to take the complex conjugate of the
  // *first* argument, while ours is to take the complex conjugate of
  // the second.
  return Eigen::Map<const typename Eigen::Matrix<T2, Eigen::Dynamic, 1>>(vec.get_values().data(), vec.size())
    .dot(Eigen::Map<const typename Eigen::Matrix<T, Eigen::Dynamic, 1>>(_val.data(), _val.size()));
#else
  typename CompareTypes<T, T2>::supertype val = 0.;

  const int N = cast_int<int>(_val.size());
  // The following pragma tells clang's vectorizer that it is safe to
  // reorder floating point operations for this loop.
#ifdef __clang__
#pragma clang loop vectorize(enable)
#endif
  for (int i=0; i<N; i++)
    val += (*this)(i)*libmesh_conj(vec(i));

  return val;
#endif
}

template<typename T>
template<typename T2>
inline
typename CompareTypes<T, T2>::supertype DenseVector<T>::indefinite_dot (const DenseVector<T2> & vec) const
{
  libmesh_assert_equal_to (this->size(), vec.size());

  typename CompareTypes<T, T2>::supertype val = 0.;

  const int N = cast_int<int>(_val.size());
  for (int i=0; i<N; i++)
    val += (*this)(i)*(vec(i));

  return val;
}

template<typename T>
template<typename T2>
inline
bool DenseVector<T>::operator== (const DenseVector<T2> & vec) const
{
  libmesh_assert_equal_to (this->size(), vec.size());

  const int N = cast_int<int>(_val.size());
  for (int i=0; i<N; i++)
    if ((*this)(i) != vec(i))
      return false;

  return true;
}



template<typename T>
template<typename T2>
inline
bool DenseVector<T>::operator!= (const DenseVector<T2> & vec) const
{
  libmesh_assert_equal_to (this->size(), vec.size());

  const int N = cast_int<int>(_val.size());
  for (int i=0; i<N; i++)
    if ((*this)(i) != vec(i))
      return true;

  return false;
}



template<typename T>
template<typename T2>
inline
DenseVector<T> & DenseVector<T>::operator+= (const DenseVector<T2> & vec)
{
  libmesh_assert_equal_to (this->size(), vec.size());

  const int N = cast_int<int>(_val.size());
  for (int i=0; i<N; i++)
    (*this)(i) += vec(i);

  return *this;
}



template<typename T>
template<typename T2>
inline
DenseVector<T> & DenseVector<T>::operator-= (const DenseVector<T2> & vec)
{
  libmesh_assert_equal_to (this->size(), vec.size());

  const int N = cast_int<int>(_val.size());
  for (int i=0; i<N; i++)
    (*this)(i) -= vec(i);

  return *this;
}



template<typename T>
inline
Real DenseVector<T>::min () const
{
  libmesh_assert (this->size());
  Real my_min = libmesh_real((*this)(0));

  const int N = cast_int<int>(_val.size());
  for (int i=1; i!=N; i++)
    {
      Real current = libmesh_real((*this)(i));
      my_min = (my_min < current? my_min : current);
    }
  return my_min;
}



template<typename T>
inline
Real DenseVector<T>::max () const
{
  libmesh_assert (this->size());
  Real my_max = libmesh_real((*this)(0));

  const int N = cast_int<int>(_val.size());
  for (int i=1; i!=N; i++)
    {
      Real current = libmesh_real((*this)(i));
      my_max = (my_max > current? my_max : current);
    }
  return my_max;
}



template<typename T>
inline
Real DenseVector<T>::l1_norm () const
{
  if (!_val.size())
    return 0.;

#ifdef LIBMESH_HAVE_EIGEN
  return Eigen::Map<const typename Eigen::Matrix<T, Eigen::Dynamic, 1>>(_val.data(), _val.size()).template lpNorm<1>();
#else
  Real my_norm = 0.;
  const int N = cast_int<int>(_val.size());
  for (int i=0; i!=N; i++)
    my_norm += std::abs((*this)(i));

  return my_norm;
#endif
}



template<typename T>
inline
Real DenseVector<T>::l2_norm () const
{
  if (!_val.size())
    return 0.;

#ifdef LIBMESH_HAVE_EIGEN
  return Eigen::Map<const typename Eigen::Matrix<T, Eigen::Dynamic, 1>>(_val.data(), _val.size()).norm();
#else
  Real my_norm = 0.;
  const int N = cast_int<int>(_val.size());
  // The following pragma tells clang's vectorizer that it is safe to
  // reorder floating point operations for this loop.
#ifdef __clang__
#pragma clang loop vectorize(enable)
#endif
  for (int i=0; i!=N; i++)
    my_norm += TensorTools::norm_sq((*this)(i));

  return sqrt(my_norm);
#endif
}



template<typename T>
inline
Real DenseVector<T>::linfty_norm () const
{
  if (!_val.size())
    return 0.;

#ifdef LIBMESH_HAVE_EIGEN
  return Eigen::Map<const typename Eigen::Matrix<T, Eigen::Dynamic, 1>>(_val.data(), _val.size()).template lpNorm<Eigen::Infinity>();
#else
  Real my_norm = TensorTools::norm_sq((*this)(0));

  const int N = cast_int<int>(_val.size());
  for (int i=1; i!=N; i++)
    {
      Real current = TensorTools::norm_sq((*this)(i));
      my_norm = (my_norm > current? my_norm : current);
    }
  return sqrt(my_norm);
#endif
}



template<typename T>
inline
void DenseVector<T>::get_principal_subvector (unsigned int sub_n,
                                              DenseVector<T> & dest) const
{
  libmesh_assert_less_equal ( sub_n, this->size() );

  dest.resize(sub_n);
  const int N = cast_int<int>(sub_n);
  for (int i=0; i<N; i++)
    dest(i) = _val[i];
}

} // namespace libMesh

#ifdef LIBMESH_HAVE_METAPHYSICL
namespace MetaPhysicL
{
template <typename T>
struct RawType<libMesh::DenseVector<T>>
{
  typedef libMesh::DenseVector<typename RawType<T>::value_type> value_type;

  static value_type value (const libMesh::DenseVector<T> & in)
    {
      value_type ret(in.size());
      for (unsigned int i = 0; i < in.size(); ++i)
          ret(i) = raw_value(in(i));

      return ret;
    }
};
}
#endif

#endif // LIBMESH_DENSE_VECTOR_H
