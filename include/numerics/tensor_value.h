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
   * Generate the intrinsic rotation matrix associated with the provided Euler angles. An instrinsic
   * rotation leaves bodies in the domain fixed while rotating the coordinate axes. We follow the
   * convention described at http://mathworld.wolfram.com/EulerAngles.html (equations 6-14 give the
   * entries of the composite transformation matrix). The rotations are performed sequentially about
   * the z, x', and z'' axes, in that order. A positive angle for a given step in the rotation
   * sequences gives the appearance of rotating an entity in the domain counter-clockwise around the
   * rotation axis, although in fact it is the coordinate axes themselves that are rotating. In
   * order to give the appearance of a body rotating counter-clockwise, we actually rotate the
   * coordinate axes by the \emph negative of the angle passed into the method. All angles should be
   * provided in degrees
   * @param phi The \emph negative of the angle we will rotate the coordinate axes around the
   * original z-axis
   * @param theta The \emph negative of the angle we will rotate the coordinate axes around the
   * "current" x-axis (post phi), e.g. x'
   * @param psi The \emph negative of the angle we will rotate the coordinate axes around the
   * "current" z-axis (post phi and theta), e.g. z''
   * @return The associated rotation matrix
   */
  static TensorValue<Real> intrinsic_rotation_matrix(Real phi, Real theta, Real psi);

  /**
   * Invert the rotation that would occur if the same angles were provided to \p
   * intrinsic_rotation_matrix, e.g. return to the original starting point. All angles should be
   * provided in degrees
   */
  static TensorValue<Real> inverse_intrinsic_rotation_matrix(Real phi, Real theta, Real psi);

  /**
   * Generate the extrinsic rotation matrix associated with the provided Euler angles. An extrinsic
   * rotation rotates the bodies in the domain and leaves the coordinate axes fixed. We follow the
   * convention described at https://en.wikipedia.org/wiki/Euler_angles#Rotation_matrix and we use
   * the matrix described by the 'Proper Euler Angles' column and Z1 X2 Z3 row, which indicates that
   * the rotations are performed sequentially about the z, x, and z axes, in that order. A
   * positive angle yields a counter-clockwise rotation about the axis in question. Note that angles
   * should be provided in degrees
   * @param angle1_deg rotation angle around z-axis
   * @param angle2_deg rotation angle around x-axis (post angle1)
   * @param angle3_deg rotation angle around z-axis (post angle1 and angle2)
   * @return The associated rotation matrix
   */
  static TensorValue<Real>
  extrinsic_rotation_matrix(Real angle1_deg, Real angle2_deg, Real angle3_deg);

  /**
   * Invert the rotation that would occur if the same angles were provided to \p
   * extrinsic_rotation_matrix, e.g. return to the original starting point
   */
  static TensorValue<Real>
  inverse_extrinsic_rotation_matrix(Real angle1_deg, Real angle2_deg, Real angle3_deg);
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
TensorValue<T>::intrinsic_rotation_matrix(const Real phi, const Real theta, const Real psi)
{
#if LIBMESH_DIM == 3
  // We apply a negative sign here or else we don't get the appearance of
  // counter-clockwise/right-hand-rule rotation of the bodies with respect to the coordinate axes
  // (but as explained in the method doxygen we are *actually* rotating the coordinate axes while
  // leaving the bodies fixed)
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
  libmesh_error_msg("TensorValue<T>::intrinsic_rotation_matrix() requires libMesh to be compiled "
                    "with LIBMESH_DIM==3");
  // We'll never get here
  return TensorValue<Real>();
#endif
}

template <typename T>
TensorValue<Real>
TensorValue<T>::inverse_intrinsic_rotation_matrix(const Real phi, const Real theta, const Real psi)
{
  // The inverse of a rotation matrix is just the transpose
  return TensorValue<T>::intrinsic_rotation_matrix(phi, theta, psi).transpose();
}

template <typename T>
TensorValue<Real>
TensorValue<T>::extrinsic_rotation_matrix(const Real angle1_deg,
                                          const Real angle2_deg,
                                          const Real angle3_deg)
{
#if LIBMESH_DIM == 3
  const auto angle1 = angle1_deg / 180. * pi;
  const auto angle2 = angle2_deg / 180. * pi;
  const auto angle3 = angle3_deg / 180. * pi;
  const auto s1 = std::sin(angle1), c1 = std::cos(angle1);
  const auto s2 = std::sin(angle2), c2 = std::cos(angle2);
  const auto s3 = std::sin(angle3), c3 = std::cos(angle3);

  return TensorValue<Real>(c1 * c3 - c2 * s1 * s3,
                           -c1 * s3 - c2 * c3 * s1,
                           s1 * s2,
                           c3 * s1 + c1 * c2 * s3,
                           c1 * c2 * c3 - s1 * s3,
                           -c1 * s2,
                           s2 * s3,
                           c3 * s2,
                           c2);
#else
  libmesh_ignore(angle1_deg, angle2_deg, angle3_deg);
  libmesh_error_msg("TensorValue<T>::extrinsic_rotation_matrix() requires libMesh to be compiled "
                    "with LIBMESH_DIM==3");
  // We'll never get here
  return TensorValue<Real>();
#endif
}

template <typename T>
TensorValue<Real>
TensorValue<T>::inverse_extrinsic_rotation_matrix(const Real angle1_deg,
                                                  const Real angle2_deg,
                                                  const Real angle3_deg)
{
  // The inverse of a rotation matrix is just the transpose
  return TensorValue<T>::extrinsic_rotation_matrix(angle1_deg, angle2_deg, angle3_deg).transpose();
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
