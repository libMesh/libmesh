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



#ifndef LIBMESH_TYPE_TENSOR_H
#define LIBMESH_TYPE_TENSOR_H

// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/type_vector.h"

// C++ includes
#include <cstdlib> // *must* precede <cmath> for proper std:abs() on PGI, Sun Studio CC
#include <cmath>
#include <tuple>

namespace libMesh
{

// Forward declarations
template <typename T> class TypeTensorColumn;
template <typename T> class ConstTypeTensorColumn;
template <unsigned int N, typename T> class TypeNTensor;

/**
 * This class defines a tensor in \p LIBMESH_DIM dimensional space of
 * type T.  T may either be Real or Complex.
 *
 * \author Roy Stogner
 * \date 2004
 */
template <typename T>
class TypeTensor
{
  template <typename T2>
  friend class TypeTensor;

public:
  /**
   * Empty constructor. Gives the tensor 0 in \p LIBMESH_DIM
   * dimensions.
   */
  TypeTensor  ();

protected:
  /**
   * Constructor-from-T.  By default sets higher dimensional entries
   * to 0.  This is a poor constructor for 2D tensors - if the default
   * arguments are to be overridden it requires that the "xz = 0."
   * etc. arguments also be given explicitly.
   */
  explicit TypeTensor  (const T & xx,
                        const T & xy=0,
                        const T & xz=0,
                        const T & yx=0,
                        const T & yy=0,
                        const T & yz=0,
                        const T & zx=0,
                        const T & zy=0,
                        const T & zz=0);


  /**
   * Constructor-from-Scalar.
   */
  template <typename Scalar>
  explicit TypeTensor  (const Scalar & xx,
                        const Scalar & xy=0,
                        const Scalar & xz=0,
                        const Scalar & yx=0,
                        const Scalar & yy=0,
                        const Scalar & yz=0,
                        const Scalar & zx=0,
                        const Scalar & zy=0,
                        typename
                        boostcopy::enable_if_c<ScalarTraits<Scalar>::value,
                        const Scalar>::type & zz=0);

  /**
   * Constructor.  Assigns each vector to a different row of the
   * tensor.  We're in LIBMESH_DIM space dimensions and so LIBMESH_DIM
   * many vectors are needed.
   */
  template<typename T2>
  TypeTensor(const TypeVector<T2> & vx);

  template<typename T2>
  TypeTensor(const TypeVector<T2> & vx,
             const TypeVector<T2> & vy);

  template<typename T2>
  TypeTensor(const TypeVector<T2> & vx,
             const TypeVector<T2> & vy,
             const TypeVector<T2> & vz);

public:

  /**
   * Helper typedef for C++98 generic programming.
   */
  typedef T value_type;

  /**
   * Helper typedef for generic index programming
   */
  typedef std::tuple<unsigned int, unsigned int> index_type;

  /**
   * Copy-constructor.
   */
  template<typename T2>
  TypeTensor(const TypeTensor<T2> & p);

  /**
   * Destructor.
   */
  ~TypeTensor();

  /**
   * Assign to this tensor without creating a temporary.
   */
  template<typename T2>
  void assign (const TypeTensor<T2> &);

  /**
   * Assignment-from-scalar operator.  Used only to zero out vectors.
   *
   * \returns A reference to *this.
   */
  template <typename Scalar>
  typename boostcopy::enable_if_c<
    ScalarTraits<Scalar>::value,
    TypeTensor &>::type
  operator = (const Scalar & libmesh_dbg_var(p))
  { libmesh_assert_equal_to (p, Scalar(0)); this->zero(); return *this; }

  /**
   * \returns A const reference to the (i,j) entry of the tensor.
   */
  const T & operator () (const unsigned int i, const unsigned int j) const;

  /**
   * \returns A writable reference to the (i,j) entry of the tensor.
   */
  T & operator () (const unsigned int i, const unsigned int j);

  /**
   * \returns A proxy for the \f$ i^{th} \f$ column of the tensor.
   */
  ConstTypeTensorColumn<T> slice (const unsigned int i) const;

  /**
   * \returns A writable proxy for the \f$ i^{th} \f$ column of the tensor.
   */
  TypeTensorColumn<T> slice (const unsigned int i);

  /**
   * \returns A copy of one row of the tensor as a TypeVector.
   */
  TypeVector<T> row(const unsigned int r) const;

  /**
   * Add another tensor to this tensor.
   *
   * \returns A copy of the result, this tensor is unchanged.
   */
  template<typename T2>
  TypeTensor<typename CompareTypes<T, T2>::supertype>
  operator + (const TypeTensor<T2> &) const;

  /**
   * Add to this tensor.
   *
   * \returns A reference to *this.
   */
  template<typename T2>
  const TypeTensor<T> & operator += (const TypeTensor<T2> &);

  /**
   * Add to this tensor without creating a temporary.
   */
  template<typename T2>
  void add (const TypeTensor<T2> &);

  /**
   * Add a scaled tensor to this tensor without creating a temporary.
   */
  template <typename T2>
  void add_scaled (const TypeTensor<T2> &, const T &);

  /**
   * Subtract a tensor from this tensor.
   *
   * \returns A copy of the result, this tensor is unchanged.
   */
  template<typename T2>
  TypeTensor<typename CompareTypes<T, T2>::supertype>
  operator - (const TypeTensor<T2> &) const;

  /**
   * Subtract from this tensor.
   *
   * \returns A reference to *this.
   */
  template<typename T2>
  const TypeTensor<T> & operator -= (const TypeTensor<T2> &);

  /**
   * Subtract from this tensor without creating a temporary.
   */
  template<typename T2>
  void subtract (const TypeTensor<T2> &);

  /**
   * Subtract a scaled value from this tensor without creating a
   * temporary.
   */
  template <typename T2>
  void subtract_scaled (const TypeTensor<T2> &, const T &);

  /**
   * \returns The negative of this tensor in a separate copy.
   */
  TypeTensor<T> operator - () const;

  /**
   * Multiply this tensor by a scalar value.
   *
   * \returns A copy of the result, this tensor is unchanged.
   */
  template <typename Scalar>
  auto
  operator * (const Scalar & scalar) const -> typename boostcopy::enable_if_c<
    ScalarTraits<Scalar>::value,
    TypeTensor<decltype(T() * scalar)>>::type;

  /**
   * Multiply this tensor by a scalar value in place.
   *
   * \returns A reference to *this.
   */
  template <typename Scalar, typename boostcopy::enable_if_c<
                               ScalarTraits<Scalar>::value, int>::type = 0>
  const TypeTensor<T> & operator *= (const Scalar & factor)
    {
      for (unsigned int i=0; i<LIBMESH_DIM*LIBMESH_DIM; i++)
        _coords[i] *= factor;

      return *this;
    }

  /**
   * Divide each entry of this tensor by a scalar value.
   *
   * \returns A copy of the result, this tensor is unchanged.
   */
  template <typename Scalar>
  typename boostcopy::enable_if_c<
    ScalarTraits<Scalar>::value,
    TypeTensor<typename CompareTypes<T, Scalar>::supertype>>::type
  operator / (const Scalar &) const;

  /**
   * Divide each entry of this tensor by a scalar value.
   *
   * \returns A reference to *this.
   */
  const TypeTensor<T> & operator /= (const T &);

  /**
   * Multiply 2 tensors together, i.e. matrix-matrix product.
   * The tensors may contain different numeric types.
   *
   * \returns A copy of the result, this tensor is unchanged.
   */
  template <typename T2>
  TypeTensor<T> operator * (const TypeTensor<T2> &) const;

  /**
   * Multiply this tensor by a tensor value in place
   *
   * \returns A reference to *this
   */
  template <typename T2>
  const TypeTensor<T> & operator *= (const TypeTensor<T2> &);

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
  contract (const TypeTensor<T2> &) const;

  /**
   * Multiply this tensor by a vector, i.e. matrix-vector product.
   * The tensor and vector may contain different numeric types.
   *
   * \returns A copy of the result vector, this tensor is unchanged.
   */
  template <typename T2>
  TypeVector<typename CompareTypes<T,T2>::supertype>
  operator * (const TypeVector<T2> &) const;

  /**
   * \returns The transpose of this tensor (with complex numbers not conjugated).
   */
  TypeTensor<T> transpose() const;

  /**
   * \returns The inverse of this tensor as an independent object.
   */
  TypeTensor<T> inverse() const;

  /**
   * Solve the 2x2 or 3x3 system of equations A*x = b for x
   * by directly inverting A and multiplying it by b.
   *
   * \returns The solution in the \p x vector.
   */
  void solve(const TypeVector<T> & b, TypeVector<T> & x) const;

  /**
   * \returns The Frobenius norm of the tensor, i.e. the square-root of
   * the sum of the elements squared.
   *
   * \deprecated Use the norm() function instead.
   */
#ifdef LIBMESH_ENABLE_DEPRECATED
  auto size() const -> decltype(std::norm(T()));
#endif

  /**
   * \returns The Frobenius norm of the tensor, i.e. the square-root of
   * the sum of the elements squared.
   */
  auto norm() const -> decltype(std::norm(T()));

  /**
   * \returns The Frobenius norm of the tensor squared, i.e. sum of the
   * element magnitudes squared.
   *
   * \deprecated Use the norm_sq() function instead.
   */
#ifdef LIBMESH_ENABLE_DEPRECATED
  auto size_sq() const -> decltype(std::norm(T()));
#endif

  /**
   * \returns The Frobenius norm of the tensor squared, i.e. sum of the
   * element magnitudes squared.
   */
  auto norm_sq() const -> decltype(std::norm(T()));

  /**
   * \returns The determinant of the tensor.
   *
   * Because these are 3x3 tensors at most, we don't do an LU
   * decomposition like DenseMatrix does.
   */
  T det() const;

  /**
   * \returns The trace of the tensor.
   */
  T tr() const;

  /**
   * Set all entries of the tensor to 0.
   */
  void zero();

  /**
   * \returns \p true if two tensors are equal, \p false otherwise.
   */
  bool operator == (const TypeTensor<T> & rhs) const;

  /**
   * \returns \p true if this tensor is "less" than another.  Useful
   * for sorting.
   */
  bool operator < (const TypeTensor<T> & rhs) const;

  /**
   * \returns \p true if this tensor is "greater" than another.
   */
  bool operator > (const TypeTensor<T> & rhs) const;

  /**
   * Formatted print to a stream which defaults to \p libMesh::out.
   */
  void print(std::ostream & os = libMesh::out) const;

  /**
   * Formatted print as above but supports the syntax:
   *
   * \code
   * std::cout << t << std::endl;
   * \endcode
   */
  friend std::ostream & operator << (std::ostream & os, const TypeTensor<T> & t)
  {
    t.print(os);
    return os;
  }

  /**
   * Unformatted print to the stream \p out.  Simply prints the elements
   * of the tensor separated by spaces and newlines.
   */
  void write_unformatted (std::ostream & out, const bool newline = true) const;

protected:

  /**
   * The coordinates of the \p TypeTensor
   */
  T _coords[LIBMESH_DIM*LIBMESH_DIM];
};


template <typename T>
class TypeTensorColumn
{
public:

  TypeTensorColumn(TypeTensor<T> & tensor,
                   unsigned int j) :
    _tensor(&tensor), _j(j) {}

  /**
   * \returns A writable reference to the \f$ i,this \f$ element of the
   * tensor.
   */
  T & operator () (const unsigned int i)
  { return (*_tensor)(i,_j); }

  T & slice (const unsigned int i)
  { return (*_tensor)(i,_j); }

  /**
   * Assign values to this column of the tensor.
   */
  TypeTensorColumn<T> & operator = (const TypeVector<T> & rhs)
  {
    for (unsigned int i=0; i != LIBMESH_DIM; ++i)
      (*this)(i) = rhs(i);
    return *this;
  }

private:
  TypeTensor<T> * _tensor;
  const unsigned int _j;
};


template <typename T>
class ConstTypeTensorColumn
{
public:

  ConstTypeTensorColumn(const TypeTensor<T> & tensor,
                        unsigned int j) :
    _tensor(&tensor), _j(j) {}

  /**
   * Return the \f$ i,this \f$ element of the tensor.
   */
  const T & operator () (const unsigned int i) const
  { return (*_tensor)(i,_j); }

  const T & slice (const unsigned int i) const
  { return (*_tensor)(i,_j); }

private:
  const TypeTensor<T> * _tensor;
  const unsigned int _j;
};

//------------------------------------------------------
// Inline functions
template <typename T>
inline
TypeTensor<T>::TypeTensor ()
{
  _coords[0] = 0;

#if LIBMESH_DIM > 1
  _coords[1] = 0;
  _coords[2] = 0;
  _coords[3] = 0;
#endif

#if LIBMESH_DIM > 2
  _coords[4] = 0;
  _coords[5] = 0;
  _coords[6] = 0;
  _coords[7] = 0;
  _coords[8] = 0;
#endif
}



template <typename T>
inline
TypeTensor<T>::TypeTensor (const T & xx,
                           const T & xy,
                           const T & xz,
                           const T & yx,
                           const T & yy,
                           const T & yz,
                           const T & zx,
                           const T & zy,
                           const T & zz)
{
  _coords[0] = xx;

#if LIBMESH_DIM == 2
  _coords[1] = xy;
  _coords[2] = yx;
  _coords[3] = yy;
#elif LIBMESH_DIM == 1
  libmesh_assert_equal_to (xy, 0);
  libmesh_assert_equal_to (yx, 0);
  libmesh_assert_equal_to (yy, 0);
#endif

#if LIBMESH_DIM == 3
  _coords[1] = xy;
  _coords[2] = xz;
  _coords[3] = yx;
  _coords[4] = yy;
  _coords[5] = yz;
  _coords[6] = zx;
  _coords[7] = zy;
  _coords[8] = zz;
#else
  libmesh_assert_equal_to (xz, 0);
  libmesh_assert_equal_to (yz, 0);
  libmesh_assert_equal_to (zx, 0);
  libmesh_assert_equal_to (zy, 0);
  libmesh_assert_equal_to (zz, 0);
#endif
}


template <typename T>
template <typename Scalar>
inline
TypeTensor<T>::TypeTensor (const Scalar & xx,
                           const Scalar & xy,
                           const Scalar & xz,
                           const Scalar & yx,
                           const Scalar & yy,
                           const Scalar & yz,
                           const Scalar & zx,
                           const Scalar & zy,
                           typename
                           boostcopy::enable_if_c<ScalarTraits<Scalar>::value,
                           const Scalar>::type & zz)
{
  _coords[0] = xx;

#if LIBMESH_DIM == 2
  _coords[1] = xy;
  _coords[2] = yx;
  _coords[3] = yy;
#endif

#if LIBMESH_DIM == 3
  _coords[1] = xy;
  _coords[2] = xz;
  _coords[3] = yx;
  _coords[4] = yy;
  _coords[5] = yz;
  _coords[6] = zx;
  _coords[7] = zy;
  _coords[8] = zz;
#endif
}



template <typename T>
template<typename T2>
inline
TypeTensor<T>::TypeTensor (const TypeTensor <T2> & p)
{
  // copy the nodes from vector p to me
  for (unsigned int i=0; i<LIBMESH_DIM*LIBMESH_DIM; i++)
    _coords[i] = p._coords[i];
}


template <typename T>
template <typename T2>
TypeTensor<T>::TypeTensor(const TypeVector<T2> & vx)
{
  libmesh_assert_equal_to (LIBMESH_DIM, 1);
  _coords[0] = vx(0);
}

template <typename T>
template <typename T2>
TypeTensor<T>::TypeTensor(const TypeVector<T2> & vx,
                          const TypeVector<T2> & vy)
{
  libmesh_assert_equal_to (LIBMESH_DIM, 2);
  _coords[0] = vx(0);
  _coords[1] = vx(1);
  _coords[2] = vy(0);
  _coords[3] = vy(1);
}

template <typename T>
template <typename T2>
TypeTensor<T>::TypeTensor(const TypeVector<T2> & vx,
                          const TypeVector<T2> & vy,
                          const TypeVector<T2> & vz)
{
  libmesh_assert_equal_to (LIBMESH_DIM, 3);
  _coords[0] = vx(0);
  _coords[1] = vx(1);
  _coords[2] = vx(2);
  _coords[3] = vy(0);
  _coords[4] = vy(1);
  _coords[5] = vy(2);
  _coords[6] = vz(0);
  _coords[7] = vz(1);
  _coords[8] = vz(2);
}




template <typename T>
inline
TypeTensor<T>::~TypeTensor ()
{
}



template <typename T>
template<typename T2>
inline
void TypeTensor<T>::assign (const TypeTensor<T2> & p)
{
  for (unsigned int i=0; i<LIBMESH_DIM*LIBMESH_DIM; i++)
    _coords[i] = p._coords[i];
}



template <typename T>
inline
const T & TypeTensor<T>::operator () (const unsigned int i,
                                      const unsigned int j) const
{
  libmesh_assert_less (i, 3);
  libmesh_assert_less (j, 3);

#if LIBMESH_DIM < 3
  const static T my_zero = 0;
  if (i >= LIBMESH_DIM || j >= LIBMESH_DIM)
    return my_zero;
#endif

  return _coords[i*LIBMESH_DIM+j];
}



template <typename T>
inline
T & TypeTensor<T>::operator () (const unsigned int i,
                                const unsigned int j)
{
#if LIBMESH_DIM < 3

  if (i >= LIBMESH_DIM || j >= LIBMESH_DIM)
    libmesh_error_msg("ERROR:  You are assigning to a tensor component that is out of range for the compiled LIBMESH_DIM!");

#endif

  libmesh_assert_less (i, LIBMESH_DIM);
  libmesh_assert_less (j, LIBMESH_DIM);

  return _coords[i*LIBMESH_DIM+j];
}


template <typename T>
inline
ConstTypeTensorColumn<T>
TypeTensor<T>::slice (const unsigned int i) const
{
  libmesh_assert_less (i, LIBMESH_DIM);
  return ConstTypeTensorColumn<T>(*this, i);
}


template <typename T>
inline
TypeTensorColumn<T>
TypeTensor<T>::slice (const unsigned int i)
{
  libmesh_assert_less (i, LIBMESH_DIM);
  return TypeTensorColumn<T>(*this, i);
}


template <typename T>
inline
TypeVector<T>
TypeTensor<T>::row(const unsigned int r) const
{
  TypeVector<T> return_vector;

  for (unsigned int j=0; j<LIBMESH_DIM; j++)
    return_vector._coords[j] = _coords[r*LIBMESH_DIM + j];

  return return_vector;
}

template <typename T>
template<typename T2>
inline
TypeTensor<typename CompareTypes<T, T2>::supertype>
TypeTensor<T>::operator + (const TypeTensor<T2> & p) const
{

#if LIBMESH_DIM == 1
  return TypeTensor(_coords[0] + p._coords[0]);
#endif

#if LIBMESH_DIM == 2
  return TypeTensor(_coords[0] + p._coords[0],
                    _coords[1] + p._coords[1],
                    0.,
                    _coords[2] + p._coords[2],
                    _coords[3] + p._coords[3]);
#endif

#if LIBMESH_DIM == 3
  return TypeTensor(_coords[0] + p._coords[0],
                    _coords[1] + p._coords[1],
                    _coords[2] + p._coords[2],
                    _coords[3] + p._coords[3],
                    _coords[4] + p._coords[4],
                    _coords[5] + p._coords[5],
                    _coords[6] + p._coords[6],
                    _coords[7] + p._coords[7],
                    _coords[8] + p._coords[8]);
#endif

}



template <typename T>
template<typename T2>
inline
const TypeTensor<T> & TypeTensor<T>::operator += (const TypeTensor<T2> & p)
{
  this->add (p);

  return *this;
}



template <typename T>
template<typename T2>
inline
void TypeTensor<T>::add (const TypeTensor<T2> & p)
{
  for (unsigned int i=0; i<LIBMESH_DIM*LIBMESH_DIM; i++)
    _coords[i] += p._coords[i];
}



template <typename T>
template <typename T2>
inline
void TypeTensor<T>::add_scaled (const TypeTensor<T2> & p, const T & factor)
{
  for (unsigned int i=0; i<LIBMESH_DIM*LIBMESH_DIM; i++)
    _coords[i] += factor*p._coords[i];

}



template <typename T>
template<typename T2>
inline
TypeTensor<typename CompareTypes<T, T2>::supertype>
TypeTensor<T>::operator - (const TypeTensor<T2> & p) const
{

#if LIBMESH_DIM == 1
  return TypeTensor(_coords[0] - p._coords[0]);
#endif

#if LIBMESH_DIM == 2
  return TypeTensor(_coords[0] - p._coords[0],
                    _coords[1] - p._coords[1],
                    0.,
                    _coords[2] - p._coords[2],
                    _coords[3] - p._coords[3]);
#endif

#if LIBMESH_DIM == 3
  return TypeTensor(_coords[0] - p._coords[0],
                    _coords[1] - p._coords[1],
                    _coords[2] - p._coords[2],
                    _coords[3] - p._coords[3],
                    _coords[4] - p._coords[4],
                    _coords[5] - p._coords[5],
                    _coords[6] - p._coords[6],
                    _coords[7] - p._coords[7],
                    _coords[8] - p._coords[8]);
#endif

}



template <typename T>
template <typename T2>
inline
const TypeTensor<T> & TypeTensor<T>::operator -= (const TypeTensor<T2> & p)
{
  this->subtract (p);

  return *this;
}



template <typename T>
template <typename T2>
inline
void TypeTensor<T>::subtract (const TypeTensor<T2> & p)
{
  for (unsigned int i=0; i<LIBMESH_DIM*LIBMESH_DIM; i++)
    _coords[i] -= p._coords[i];
}



template <typename T>
template <typename T2>
inline
void TypeTensor<T>::subtract_scaled (const TypeTensor<T2> & p, const T & factor)
{
  for (unsigned int i=0; i<LIBMESH_DIM*LIBMESH_DIM; i++)
    _coords[i] -= factor*p._coords[i];
}



template <typename T>
inline
TypeTensor<T> TypeTensor<T>::operator - () const
{

#if LIBMESH_DIM == 1
  return TypeTensor(-_coords[0]);
#endif

#if LIBMESH_DIM == 2
  return TypeTensor(-_coords[0],
                    -_coords[1],
                    -_coords[2],
                    -_coords[3]);
#endif

#if LIBMESH_DIM == 3
  return TypeTensor(-_coords[0],
                    -_coords[1],
                    -_coords[2],
                    -_coords[3],
                    -_coords[4],
                    -_coords[5],
                    -_coords[6],
                    -_coords[7],
                    -_coords[8]);
#endif

}



template <typename T>
template <typename Scalar>
inline
auto
TypeTensor<T>::operator * (const Scalar & factor) const -> typename boostcopy::enable_if_c<
  ScalarTraits<Scalar>::value,
  TypeTensor<decltype(T() * factor)>>::type
{
  typedef decltype((*this)(0, 0) * factor) TS;


#if LIBMESH_DIM == 1
  return TypeTensor<TS>(_coords[0]*factor);
#endif

#if LIBMESH_DIM == 2
  return TypeTensor<TS>(_coords[0]*factor,
                        _coords[1]*factor,
                        _coords[2]*factor,
                        _coords[3]*factor);
#endif

#if LIBMESH_DIM == 3
  return TypeTensor<TS>(_coords[0]*factor,
                        _coords[1]*factor,
                        _coords[2]*factor,
                        _coords[3]*factor,
                        _coords[4]*factor,
                        _coords[5]*factor,
                        _coords[6]*factor,
                        _coords[7]*factor,
                        _coords[8]*factor);
#endif
}



template <typename T, typename Scalar>
inline
typename boostcopy::enable_if_c<
  ScalarTraits<Scalar>::value,
  TypeTensor<typename CompareTypes<T, Scalar>::supertype>>::type
operator * (const Scalar & factor,
            const TypeTensor<T> & t)
{
  return t * factor;
}

template <typename T>
template <typename Scalar>
inline
typename boostcopy::enable_if_c<
  ScalarTraits<Scalar>::value,
  TypeTensor<typename CompareTypes<T, Scalar>::supertype>>::type
TypeTensor<T>::operator / (const Scalar & factor) const
{
  libmesh_assert_not_equal_to (factor, static_cast<T>(0.));

  typedef typename CompareTypes<T, Scalar>::supertype TS;

#if LIBMESH_DIM == 1
  return TypeTensor<TS>(_coords[0]/factor);
#endif

#if LIBMESH_DIM == 2
  return TypeTensor<TS>(_coords[0]/factor,
                        _coords[1]/factor,
                        _coords[2]/factor,
                        _coords[3]/factor);
#endif

#if LIBMESH_DIM == 3
  return TypeTensor<TS>(_coords[0]/factor,
                        _coords[1]/factor,
                        _coords[2]/factor,
                        _coords[3]/factor,
                        _coords[4]/factor,
                        _coords[5]/factor,
                        _coords[6]/factor,
                        _coords[7]/factor,
                        _coords[8]/factor);
#endif

}



template <typename T>
inline
TypeTensor<T> TypeTensor<T>::transpose() const
{
#if LIBMESH_DIM == 1
  return TypeTensor(_coords[0]);
#endif

#if LIBMESH_DIM == 2
  return TypeTensor(_coords[0],
                    _coords[2],
                    _coords[1],
                    _coords[3]);
#endif

#if LIBMESH_DIM == 3
  return TypeTensor(_coords[0],
                    _coords[3],
                    _coords[6],
                    _coords[1],
                    _coords[4],
                    _coords[7],
                    _coords[2],
                    _coords[5],
                    _coords[8]);
#endif
}



template <typename T>
inline
TypeTensor<T> TypeTensor<T>::inverse() const
{
#if LIBMESH_DIM == 1
  if (_coords[0] == static_cast<T>(0.))
    libmesh_convergence_failure();
  return TypeTensor(1. / _coords[0]);
#endif

#if LIBMESH_DIM == 2
  // Get convenient reference to this.
  const TypeTensor<T> & A = *this;

  // Use temporary variables, avoid multiple accesses
  T a = A(0,0), b = A(0,1),
    c = A(1,0), d = A(1,1);

  // Make sure det = ad - bc is not zero
  T my_det = a*d - b*c;

  if (my_det == static_cast<T>(0.))
    libmesh_convergence_failure();

  return TypeTensor(d/my_det, -b/my_det, -c/my_det, a/my_det);
#endif

#if LIBMESH_DIM == 3
  // Get convenient reference to this.
  const TypeTensor<T> & A = *this;

  T a11 = A(0,0), a12 = A(0,1), a13 = A(0,2),
    a21 = A(1,0), a22 = A(1,1), a23 = A(1,2),
    a31 = A(2,0), a32 = A(2,1), a33 = A(2,2);

  T my_det = a11*(a33*a22-a32*a23) - a21*(a33*a12-a32*a13) + a31*(a23*a12-a22*a13);

  if (my_det == static_cast<T>(0.))
    libmesh_convergence_failure();

  // Inline comment characters are for lining up columns.
  return TypeTensor(/**/  (a33*a22-a32*a23)/my_det, -(a33*a12-a32*a13)/my_det,  (a23*a12-a22*a13)/my_det,
                    /**/ -(a33*a21-a31*a23)/my_det,  (a33*a11-a31*a13)/my_det, -(a23*a11-a21*a13)/my_det,
                    /**/  (a32*a21-a31*a22)/my_det, -(a32*a11-a31*a12)/my_det,  (a22*a11-a21*a12)/my_det);
#endif
}



template <typename T>
inline
void TypeTensor<T>::solve(const TypeVector<T> & b, TypeVector<T> & x) const
{
#if LIBMESH_DIM == 1
  if (_coords[0] == static_cast<T>(0.))
    libmesh_convergence_failure();
  x(0) = b(0) / _coords[0];
#endif

#if LIBMESH_DIM == 2
  T my_det = _coords[0]*_coords[3] - _coords[1]*_coords[2];

  if (my_det == static_cast<T>(0.))
    libmesh_convergence_failure();

  T my_det_inv = 1./my_det;

  x(0) = my_det_inv*( _coords[3]*b(0) - _coords[1]*b(1));
  x(1) = my_det_inv*(-_coords[2]*b(0) + _coords[0]*b(1));
#endif

#if LIBMESH_DIM == 3
  T my_det =
    //           a11*(a33       *a22        - a32       *a23)
    /**/  _coords[0]*(_coords[8]*_coords[4] - _coords[7]*_coords[5])
    //          -a21*(a33       *a12        - a32       *a13)
    /**/ -_coords[3]*(_coords[8]*_coords[1] - _coords[7]*_coords[2]) +
    //          +a31*(a23       *a12        - a22       *a13)
    /**/ +_coords[6]*(_coords[5]*_coords[1] - _coords[4]*_coords[2]);

  if (my_det == static_cast<T>(0.))
    libmesh_convergence_failure();

  T my_det_inv = 1./my_det;
  x(0) = my_det_inv*((_coords[8]*_coords[4] - _coords[7]*_coords[5])*b(0) -
                     (_coords[8]*_coords[1] - _coords[7]*_coords[2])*b(1) +
                     (_coords[5]*_coords[1] - _coords[4]*_coords[2])*b(2));

  x(1) = my_det_inv*((_coords[6]*_coords[5] - _coords[8]*_coords[3])*b(0) +
                     (_coords[8]*_coords[0] - _coords[6]*_coords[2])*b(1) -
                     (_coords[5]*_coords[0] - _coords[3]*_coords[2])*b(2));

  x(2) = my_det_inv*((_coords[7]*_coords[3] - _coords[6]*_coords[4])*b(0) -
                     (_coords[7]*_coords[0] - _coords[6]*_coords[1])*b(1) +
                     (_coords[4]*_coords[0] - _coords[3]*_coords[1])*b(2));
#endif
}



template <typename T>
inline
const TypeTensor<T> & TypeTensor<T>::operator /= (const T & factor)
{
  libmesh_assert_not_equal_to (factor, static_cast<T>(0.));

  for (unsigned int i=0; i<LIBMESH_DIM*LIBMESH_DIM; i++)
    _coords[i] /= factor;

  return *this;
}




template <typename T>
template <typename T2>
inline
TypeVector<typename CompareTypes<T,T2>::supertype>
TypeTensor<T>::operator * (const TypeVector<T2> & p) const
{
  TypeVector<typename CompareTypes<T,T2>::supertype> returnval;
  for (unsigned int i=0; i<LIBMESH_DIM; i++)
    for (unsigned int j=0; j<LIBMESH_DIM; j++)
      returnval(i) += (*this)(i,j)*p(j);

  return returnval;
}



template <typename T>
template <typename T2>
inline
TypeTensor<T> TypeTensor<T>::operator * (const TypeTensor<T2> & p) const
{
  TypeTensor<T> returnval;
  for (unsigned int i=0; i<LIBMESH_DIM; i++)
    for (unsigned int j=0; j<LIBMESH_DIM; j++)
      for (unsigned int k=0; k<LIBMESH_DIM; k++)
        returnval(i,j) += (*this)(i,k)*p(k,j);

  return returnval;
}

template <typename T>
template <typename T2>
inline
const TypeTensor<T> & TypeTensor<T>::operator *= (const TypeTensor<T2> & p)
{
  TypeTensor<T> temp;
  for (unsigned int i=0; i<LIBMESH_DIM; i++)
    for (unsigned int j=0; j<LIBMESH_DIM; j++)
      for (unsigned int k=0; k<LIBMESH_DIM; k++)
        temp(i,j) += (*this)(i,k)*p(k,j);

  this->assign(temp);
  return *this;
}


/**
 * Multiply 2 tensors together, i.e. sum Aij*Bij.
 * The tensors may be of different types.
 */
template <typename T>
template <typename T2>
inline
typename CompareTypes<T,T2>::supertype
TypeTensor<T>::contract (const TypeTensor<T2> & t) const
{
  typename CompareTypes<T,T2>::supertype sum = 0.;
  for (unsigned int i=0; i<LIBMESH_DIM*LIBMESH_DIM; i++)
    sum += _coords[i]*t._coords[i];
  return sum;
}



#ifdef LIBMESH_ENABLE_DEPRECATED
template <typename T>
inline
auto TypeTensor<T>::size() const -> decltype(std::norm(T()))
{
  libmesh_deprecated();
  return this->norm();
}
#endif



template <typename T>
inline
auto TypeTensor<T>::norm() const -> decltype(std::norm(T()))
{
  return std::sqrt(this->norm_sq());
}



template <typename T>
inline
T TypeTensor<T>::det() const
{
#if LIBMESH_DIM == 1
  return _coords[0];
#endif

#if LIBMESH_DIM == 2
  return (_coords[0] * _coords[3]
          - _coords[1] * _coords[2]);
#endif

#if LIBMESH_DIM == 3
  return (_coords[0] * _coords[4] * _coords[8]
          + _coords[1] * _coords[5] * _coords[6]
          + _coords[2] * _coords[3] * _coords[7]
          - _coords[0] * _coords[5] * _coords[7]
          - _coords[1] * _coords[3] * _coords[8]
          - _coords[2] * _coords[4] * _coords[6]);
#endif
}

template <typename T>
inline
T TypeTensor<T>::tr() const
{
#if LIBMESH_DIM == 1
  return _coords[0];
#endif

#if LIBMESH_DIM == 2
  return _coords[0] + _coords[3];
#endif

#if LIBMESH_DIM == 3
  return _coords[0] + _coords[4] + _coords[8];
#endif
}

template <typename T>
inline
void TypeTensor<T>::zero()
{
  for (unsigned int i=0; i<LIBMESH_DIM*LIBMESH_DIM; i++)
    _coords[i] = 0.;
}



#ifdef LIBMESH_ENABLE_DEPRECATED
template <typename T>
inline
auto TypeTensor<T>::size_sq () const -> decltype(std::norm(T()))
{
  libmesh_deprecated();
  return this->norm_sq();
}
#endif



template <typename T>
inline
auto TypeTensor<T>::norm_sq () const -> decltype(std::norm(T()))
{
  Real sum = 0.;
  for (unsigned int i=0; i<LIBMESH_DIM*LIBMESH_DIM; i++)
    sum += TensorTools::norm_sq(_coords[i]);
  return sum;
}



template <typename T>
inline
bool TypeTensor<T>::operator == (const TypeTensor<T> & rhs) const
{
#if LIBMESH_DIM == 1
  return (std::abs(_coords[0] - rhs._coords[0])
          < TOLERANCE);
#endif

#if LIBMESH_DIM == 2
  return ((std::abs(_coords[0] - rhs._coords[0]) +
           std::abs(_coords[1] - rhs._coords[1]) +
           std::abs(_coords[2] - rhs._coords[2]) +
           std::abs(_coords[3] - rhs._coords[3]))
          < 4.*TOLERANCE);
#endif

#if LIBMESH_DIM == 3
  return ((std::abs(_coords[0] - rhs._coords[0]) +
           std::abs(_coords[1] - rhs._coords[1]) +
           std::abs(_coords[2] - rhs._coords[2]) +
           std::abs(_coords[3] - rhs._coords[3]) +
           std::abs(_coords[4] - rhs._coords[4]) +
           std::abs(_coords[5] - rhs._coords[5]) +
           std::abs(_coords[6] - rhs._coords[6]) +
           std::abs(_coords[7] - rhs._coords[7]) +
           std::abs(_coords[8] - rhs._coords[8]))
          < 9.*TOLERANCE);
#endif

}

} // namespace libMesh

#endif // LIBMESH_TYPE_TENSOR_H
