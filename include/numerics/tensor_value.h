// The libMesh Finite Element Library.
// Copyright (C) 2002-2022 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/libmesh.h" // for pi

#ifdef LIBMESH_HAVE_METAPHYSICL
#include "metaphysicl/raw_type.h"
#endif

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
  typedef T value_type;

  template <typename T2>
  struct rebind
  {
    typedef TensorValue<T2> other;
  };

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

  /**
   * Generate the rotation matrix associated with the provided Euler angles.  We follow the
   * convention described at http://mathworld.wolfram.com/EulerAngles.html (equations 6-14 give the
   * entries of the composite transformation matrix). The rotations are performed sequentially about
   * the z, x', and z'' axes, in that order. A positive angle yields a counter-clockwise rotation
   * about the axis in question. Note that angles should be provided in degrees
   * @param phi rotation angle around original z-axis
   * @param theta rotation angle around "current" x-axis (post phi), e.g. x'
   * @param psi rotation angle around "current" z-axis (post phi and theta), e.g. z''
   * @return The associated rotation matrix
   */
  static TensorValue<Real> rotation_matrix(Real phi, Real theta, Real psi);

  /**
   * Invert the rotation that would occur if the same angles were provided to \p rotation_matrix,
   * e.g. return to the original starting point
   */
  static TensorValue<Real> inverse_rotation_matrix(Real phi, Real theta, Real psi);
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

template <typename T>
TensorValue<Real>
TensorValue<T>::rotation_matrix(const Real phi, const Real theta, const Real psi)
{
#if LIBMESH_DIM == 3
  const Real p = -phi / 180. * pi;
  const Real t = -theta / 180. * pi;
  const Real s = -psi / 180. * pi;
  const Real sp = std::sin(p), cp = std::cos(p);
  const Real st = std::sin(t), ct = std::cos(t);
  const Real ss = std::sin(s), cs = std::cos(s);

  return TensorValue<Real>(cp * cs - sp * ct * ss,
                           sp * cs + cp * ct * ss,
                           st * ss,
                           -cp * ss - sp * ct * cs,
                           -sp * ss + cp * ct * cs,
                           st * cs,
                           sp * st,
                           -cp * st,
                           ct);
#else
  libmesh_ignore(phi, theta, psi);
  libmesh_error_msg(
      "TensorValue<T>::rotation_matrix() requires libMesh to be compiled with LIBMESH_DIM==3");
  // We'll never get here
  return TensorValue<Real>();
#endif
}

template <typename T>
TensorValue<Real>
TensorValue<T>::inverse_rotation_matrix(const Real phi, const Real theta, const Real psi)
{
#if LIBMESH_DIM == 3
  using namespace std;

  const Real p = -phi / 180. * pi;
  const Real t = -theta / 180. * pi;
  const Real s = -psi / 180. * pi;

  // These inverse matrices are constructed from equations 3, 4, and 5 at
  // https://mathworld.wolfram.com/EulerAngles.html and combined using the knowledge that
  // - We want to apply the matrices in the reverse order, e.g. D^{-1}C^{-1}B^{-1}(BCD)
  // - We want to apply the negative of the original angle
  // - cos(-foo) = cos(foo) and sin(-foo) = -sin(foo)
  TensorValue<Real> Binv(cos(s), -sin(s), 0, sin(s), cos(s), 0, 0, 0, 1);
  TensorValue<Real> Cinv(1, 0, 0, 0, cos(t), -sin(t), 0, sin(t), cos(t));
  TensorValue<Real> Dinv(cos(p), -sin(p), 0, sin(p), cos(p), 0, 0, 0, 1);

  return Dinv * Cinv * Binv;

#else
  libmesh_ignore(phi, theta, psi);
  libmesh_error_msg(
      "TensorValue<T>::rotation_matrix() requires libMesh to be compiled with LIBMESH_DIM==3");
  // We'll never get here
  return TensorValue<Real>();
#endif
}

} // namespace libMesh

#ifdef LIBMESH_HAVE_METAPHYSICL
namespace MetaPhysicL
{
template <typename T>
struct RawType<libMesh::TensorValue<T>>
{
  typedef libMesh::TensorValue<typename RawType<T>::value_type> value_type;

  static value_type value (const libMesh::TensorValue<T> & in)
    {
      value_type ret;
      for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
        for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
          ret(i,j) = raw_value(in(i,j));

      return ret;
    }
};
}
#endif

#endif // LIBMESH_TENSOR_VALUE_H
