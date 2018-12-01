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



#ifndef LIBMESH_TENSOR_VALUE_H
#define LIBMESH_TENSOR_VALUE_H

// Local includes
#include "libmesh/type_tensor.h"

// C++ includes

namespace libMesh
{

/**
 * This class defines a tensor in LIBMESH_DIM dimensional Real or Complex
 * space.  The typedef RealTensorValue always defines a real-valued tensor,
 * and NumberTensorValue defines a real or complex-valued tensor depending
 * on how the library was configured.
 *
 * \author Roy H. Stogner
 * \date 2004
 */
template <typename T>
class TensorValue : public TypeTensor<T>
{
public:

  /**
   * Empty constructor.
   * Gives the tensor 0 in \p LIBMESH_DIM dimensional T space.
   */
  TensorValue  ();

  /**
   * Constructor-from-T.  By default sets higher dimensional
   * entries to 0.
   */
  explicit TensorValue  (const T & xx,
                         const T & xy=0,
                         const T & xz=0,
                         const T & yx=0,
                         const T & yy=0,
                         const T & yz=0,
                         const T & zx=0,
                         const T & zy=0,
                         const T & zz=0);

  /**
   * Constructor-from-scalars.  By default sets higher dimensional
   * entries to 0.
   */
  template <typename Scalar>
  explicit TensorValue  (const Scalar & xx,
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
   * Constructor.  Takes 1 row vector for LIBMESH_DIM=1
   */
  template <typename T2>
  TensorValue (const TypeVector<T2> & vx);

  /**
   * Constructor.  Takes 2 row vectors for LIBMESH_DIM=2
   */
  template <typename T2>
  TensorValue (const TypeVector<T2> & vx,
               const TypeVector<T2> & vy);

  /**
   * Constructor.  Takes 3 row vectors for LIBMESH_DIM=3
   */
  template <typename T2>
  TensorValue (const TypeVector<T2> & vx,
               const TypeVector<T2> & vy,
               const TypeVector<T2> & vz);


  /**
   * Copy-constructor.
   */
  template <typename T2>
  TensorValue (const TensorValue<T2> & p);

  /**
   * Copy-constructor.
   */
  template <typename T2>
  TensorValue (const TypeTensor<T2> & p);


#ifdef LIBMESH_USE_COMPLEX_NUMBERS
  /**
   * Constructor that takes two \p TypeTensor<Real>
   * representing the real and imaginary part as
   * arguments.
   */
  TensorValue (const TypeTensor<Real> & p_re,
               const TypeTensor<Real> & p_im);
#endif


  /**
   * Assignment-from-scalar operator.  Used only to zero out tensors.
   */
  template <typename Scalar>
  typename boostcopy::enable_if_c<
    ScalarTraits<Scalar>::value,
    TensorValue &>::type
  operator = (const Scalar & libmesh_dbg_var(p) )
  { libmesh_assert_equal_to (p, Scalar(0)); this->zero(); return *this; }
};



/**
 * Useful typedefs to allow transparent switching
 * between Real and Complex data types.
 */
typedef TensorValue<Real>   RealTensorValue;
typedef TensorValue<Number> NumberTensorValue;
typedef RealTensorValue     RealTensor;
typedef NumberTensorValue   Tensor;



//------------------------------------------------------
// Inline functions
template <typename T>
inline
TensorValue<T>::TensorValue () :
  TypeTensor<T> ()
{
}



template <typename T>
inline
TensorValue<T>::TensorValue (const T & xx,
                             const T & xy,
                             const T & xz,
                             const T & yx,
                             const T & yy,
                             const T & yz,
                             const T & zx,
                             const T & zy,
                             const T & zz) :
  TypeTensor<T> (xx,xy,xz,yx,yy,yz,zx,zy,zz)
{
}


template <typename T>
template <typename Scalar>
inline
TensorValue<T>::TensorValue (const Scalar & xx,
                             const Scalar & xy,
                             const Scalar & xz,
                             const Scalar & yx,
                             const Scalar & yy,
                             const Scalar & yz,
                             const Scalar & zx,
                             const Scalar & zy,
                             typename
                             boostcopy::enable_if_c<ScalarTraits<Scalar>::value,
                             const Scalar>::type & zz) :
  TypeTensor<T> (xx,xy,xz,yx,yy,yz,zx,zy,zz)
{
}



template <typename T>
template <typename T2>
inline
TensorValue<T>::TensorValue (const TensorValue<T2> & p) :
  TypeTensor<T> (p)
{
}



template <typename T>
template <typename T2>
inline
TensorValue<T>::TensorValue (const TypeVector<T2> & vx) :
  TypeTensor<T> (vx)
{
}



template <typename T>
template <typename T2>
inline
TensorValue<T>::TensorValue (const TypeVector<T2> & vx,
                             const TypeVector<T2> & vy) :
  TypeTensor<T> (vx, vy)
{
}



template <typename T>
template <typename T2>
inline
TensorValue<T>::TensorValue (const TypeVector<T2> & vx,
                             const TypeVector<T2> & vy,
                             const TypeVector<T2> & vz) :
  TypeTensor<T> (vx, vy, vz)
{
}



template <typename T>
template <typename T2>
inline
TensorValue<T>::TensorValue (const TypeTensor<T2> & p) :
  TypeTensor<T> (p)
{
}


#ifdef LIBMESH_USE_COMPLEX_NUMBERS
template <typename T>
inline
TensorValue<T>::TensorValue (const TypeTensor<Real> & p_re,
                             const TypeTensor<Real> & p_im) :
  TypeTensor<T> (Complex (p_re(0,0), p_im(0,0)),
                 Complex (p_re(0,1), p_im(0,1)),
                 Complex (p_re(0,2), p_im(0,2)),
                 Complex (p_re(1,0), p_im(1,0)),
                 Complex (p_re(1,1), p_im(1,1)),
                 Complex (p_re(1,2), p_im(1,2)),
                 Complex (p_re(2,0), p_im(2,0)),
                 Complex (p_re(2,1), p_im(2,1)),
                 Complex (p_re(2,2), p_im(2,2)))
{
}
#endif


} // namespace libMesh

#endif // LIBMESH_TENSOR_VALUE_H
