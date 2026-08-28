// The libMesh Finite Element Library.
// Copyright (C) 2002-2026 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#ifndef LIBMESH_KOKKOS_VECTOR_OPS_H
#define LIBMESH_KOKKOS_VECTOR_OPS_H

// libMesh Kokkos generic vector operations.
//
// These free functions build vector algebra on top of the primitive
// access/materialization layer in kokkos_linalg_base.h. They are written
// against vector-like inputs so both libMesh owning types and storage-backed
// refs can participate in the same math.
//
// The header is organized as:
//
// - construction/materialization helpers (copy_vector),
// - detail:: componentwise kernels shared by everything below,
// - reductions, geometry helpers, and libMesh-like wrappers, and
// - the out-of-line vector_ref definitions and the operator overloads.
//
// The operator overloads are constrained (via is_vector_ref_v) to
// require at least one storage-backed ref operand: expressions between
// two owning vectors keep using the operators those types already
// define, while vector-only expressions involving a ref resolve here
// (mixed vector/tensor products live in kokkos_tensor_ops.h).
//
// Functions taking an optional ResultVector template parameter return
// the input's semantic type (e.g. TypeVector) by default, but a caller
// can request a different owning type explicitly, e.g.
// vector_cross<Point>(a, b).

#include "libmesh/kokkos_linalg_base.h"

#include "libmesh/tensor_tools.h"

#include <cmath>

namespace libMesh::Kokkos
{

// Construction and materialization

/**
 * \returns An owning copy of the vector-like \p v: by default a value
 * of \p v 's semantic type, or of \p ResultVector if that parameter is
 * given explicitly.  This is how a kernel snapshots a storage-backed
 * ref's current value for register-resident work.
 */
template <typename ResultVector = void, typename VectorLike>
LIBMESH_DEVICE_INLINE
auto copy_vector(const VectorLike & v)
  -> std::conditional_t<std::is_void<ResultVector>::value,
                        vector_semantic_type_t<VectorLike>,
                        ResultVector>
{
  using output_type = std::conditional_t<std::is_void<ResultVector>::value,
                                         vector_semantic_type_t<VectorLike>,
                                         ResultVector>;
  return materialize_vector<output_type>(v);
}

namespace detail
{

// These helpers are shared by the public functions and ref operators so
// Kokkos-backed refs use direct component access without extra materialization.

/**
 * Componentwise copy: left = right.
 */
template <typename LeftVector, typename RightVector>
LIBMESH_DEVICE_INLINE
void assign_vector_components(LeftVector & left, const RightVector & right)
{
  for (unsigned int component = 0; component < LIBMESH_DIM; ++component)
    left(component) = right(component);
}

/**
 * Sets every component of \p v to \p value.
 */
template <typename VectorLike, typename Scalar>
LIBMESH_DEVICE_INLINE
void fill_vector_components(VectorLike & v, const Scalar & value)
{
  for (unsigned int component = 0; component < LIBMESH_DIM; ++component)
    v(component) = value;
}

/**
 * Componentwise axpy-style update: left += factor * right.
 */
template <typename LeftVector, typename RightVector, typename Scalar>
LIBMESH_DEVICE_INLINE
void update_vector_components(LeftVector & left, const RightVector & right, const Scalar & factor)
{
  for (unsigned int component = 0; component < LIBMESH_DIM; ++component)
    left(component) = left(component) + factor * right(component);
}

/**
 * Componentwise transform: out(i) = op(in(i)).  \p out and \p in may
 * alias, as they do in the compound assignment operators below.
 */
template <typename OutputVector, typename InputVector, typename TransformOp>
LIBMESH_DEVICE_INLINE
void transform_vector_components(OutputVector & out, const InputVector & in, const TransformOp & op)
{
  for (unsigned int component = 0; component < LIBMESH_DIM; ++component)
    out(component) = op(in(component));
}

/**
 * \returns A new owning \p ResultVector with components op(v(i)).
 */
template <typename ResultVector, typename VectorLike, typename TransformOp>
LIBMESH_DEVICE_INLINE
ResultVector transformed_vector(const VectorLike & v, const TransformOp & op)
{
  ResultVector out;
  out.zero();
  transform_vector_components(out, v, op);
  return out;
}

/**
 * \returns \p true if \p left and \p right agree exactly in every
 * component, \p false otherwise.
 */
template <typename LeftVector, typename RightVector>
LIBMESH_DEVICE_INLINE
bool vector_equal_impl(const LeftVector & left, const RightVector & right)
{
  for (unsigned int component = 0; component < LIBMESH_DIM; ++component)
    if (left(component) != right(component))
      return false;

  return true;
}

/**
 * Componentwise transform functors, shared by the vector helpers here
 * and the tensor helpers in kokkos_tensor_ops.h.  scale_value and
 * divide_value hold their scalar by reference, so they must not
 * outlive it; both headers only ever build and consume them within a
 * single expression.
 */
template <typename ValueType>
struct negate_value
{
  LIBMESH_DEVICE_INLINE
  auto operator()(const ValueType & value) const
  {
    return -value;
  }
};

template <typename Scalar>
struct scale_value
{
  const Scalar & alpha;

  LIBMESH_DEVICE_INLINE
  auto operator()(const Scalar & value) const -> decltype(value * alpha)
  {
    return value * alpha;
  }

  template <typename ValueType>
  LIBMESH_DEVICE_INLINE
  auto operator()(const ValueType & value) const -> decltype(value * alpha)
  {
    return value * alpha;
  }
};

template <typename Scalar>
struct divide_value
{
  const Scalar & alpha;

  LIBMESH_DEVICE_INLINE
  auto operator()(const Scalar & value) const -> decltype(value / alpha)
  {
    return value / alpha;
  }

  template <typename ValueType>
  LIBMESH_DEVICE_INLINE
  auto operator()(const ValueType & value) const -> decltype(value / alpha)
  {
    return value / alpha;
  }
};

} // namespace detail

// Reductions and normalization

/**
 * \returns The dot product of \p left and \p right, in whatever type
 * their componentwise products sum to.
 */
template <typename LeftVector, typename RightVector>
LIBMESH_DEVICE_INLINE
auto vector_dot(const LeftVector & left, const RightVector & right)
{
  static_assert(is_vector_like_v<LeftVector>, "vector_dot() requires a vector-like left input");
  static_assert(is_vector_like_v<RightVector>, "vector_dot() requires a vector-like right input");

  using sum_type =
    detail::remove_cvref_t<decltype(left(0) * right(0))>;

  sum_type sum = sum_type(0);
  for (unsigned int component = 0; component < LIBMESH_DIM; ++component)
    sum += left(component) * right(component);

  return sum;
}


/**
 * \returns An owning unit vector in the direction of \p v.  Asserts
 * (in debug modes) that \p v is not the zero vector.
 */
template <typename ResultVector = void, typename VectorLike>
LIBMESH_DEVICE_INLINE
auto vector_unit(const VectorLike & v)
  -> std::conditional_t<std::is_void<ResultVector>::value,
                        vector_semantic_type_t<VectorLike>,
                        ResultVector>
{
  const auto length = v.norm();
  libmesh_assert_not_equal_to(length, static_cast<Real>(0.));
  using output_type = std::conditional_t<std::is_void<ResultVector>::value,
                                         vector_semantic_type_t<VectorLike>,
                                         ResultVector>;
  return detail::transformed_vector<output_type>(v, detail::divide_value<decltype(length)>{length});
}

// Geometry

/**
 * \returns The cross product of \p left and \p right, as an owning
 * vector.
 *
 * \note The result is only nonzero when LIBMESH_DIM == 3; in lower
 * dimensions the zero vector is returned, without the debug-mode
 * assertion TypeVector::cross() performs.
 */
template <typename ResultVector = void, typename LeftVector, typename RightVector>
LIBMESH_DEVICE_INLINE
auto vector_cross(const LeftVector & left, const RightVector & right)
  -> std::conditional_t<std::is_void<ResultVector>::value,
                        vector_semantic_type_t<LeftVector>,
                        ResultVector>
{
  using output_type = std::conditional_t<std::is_void<ResultVector>::value,
                                         vector_semantic_type_t<LeftVector>,
                                         ResultVector>;
  output_type out;
  out.zero();

#if LIBMESH_DIM == 3
  out(0) = left(1) * right(2) - left(2) * right(1);
  out(1) = -left(0) * right(2) + left(2) * right(0);
  out(2) = left(0) * right(1) - left(1) * right(0);
#else
  libmesh_ignore(left);
  libmesh_ignore(right);
#endif

  return out;
}

/**
 * \returns The scalar triple product left . (middle x right), i.e. the
 * signed volume of the parallelepiped the three vectors span (zero
 * unless LIBMESH_DIM == 3).
 */
template <typename LeftVector, typename MiddleVector, typename RightVector>
LIBMESH_DEVICE_INLINE
auto vector_triple_product(const LeftVector & left,
                           const MiddleVector & middle,
                           const RightVector & right)
{
#if LIBMESH_DIM == 3
  return left(0) * (middle(1) * right(2) - middle(2) * right(1)) -
         left(1) * (middle(0) * right(2) - middle(2) * right(0)) +
         left(2) * (middle(0) * right(1) - middle(1) * right(0));
#else
  libmesh_ignore(left, middle, right);
  using value_type =
    detail::remove_cvref_t<decltype(left(0) * middle(0))>;
  return value_type(0);
#endif
}

/**
 * \returns |left x right|^2, computed without the extra temporary a
 * cross().norm_sq() chain would create, mirroring cross_norm_sq() in
 * type_vector.h.  In 2D only the out-of-plane component contributes.
 */
template <typename LeftVector, typename RightVector>
LIBMESH_DEVICE_INLINE
auto vector_cross_norm_sq(const LeftVector & left, const RightVector & right)
{
  const auto z = left(0) * right(1) - left(1) * right(0);

#if LIBMESH_DIM == 3
  const auto x = left(1) * right(2) - left(2) * right(1);
  const auto y = left(0) * right(2) - left(2) * right(0);
  return x * x + y * y + z * z;
#else
  return z * z;
#endif
}

/**
 * \returns The solid angle subtended by a tetrahedral vertex, as
 * defined by the edge vectors \p v01, \p v02, \p v03 - positive if
 * the vectors obey the right-hand rule, negative for a left-hand
 * orientation.  Computed by Van Oosterom and Strackee's formula,
 * mirroring the solid_angle() free function in type_vector.h.
 */
template <typename VectorA, typename VectorB, typename VectorC>
LIBMESH_DEVICE_INLINE
auto vector_solid_angle(const VectorA & v01, const VectorB & v02, const VectorC & v03)
{
  using std::atan;

  const auto norm01 = v01.norm();
  const auto norm02 = v02.norm();
  const auto norm03 = v03.norm();
  const auto tan_half_angle =
    vector_triple_product(v01, v02, v03) /
    (vector_dot(v01, v02) * norm03 +
     vector_dot(v01, v03) * norm02 +
     vector_dot(v02, v03) * norm01 +
     norm01 * norm02 * norm03);

  return Real(2) * atan(tan_half_angle);
}

// libMesh-like convenience wrappers

/**
 * \returns The dot product of \p left and \p right; a synonym for
 * vector_dot() matching the TypeVector::contract() name.
 */
template <typename LeftVector,
          typename RightVector,
          typename std::enable_if<is_vector_like_v<LeftVector> && is_vector_like_v<RightVector>,
                                  int>::type = 0>
LIBMESH_DEVICE_INLINE
auto contract(const LeftVector & left, const RightVector & right)
{
  return vector_dot(left, right);
}


// Forward declarations of the compound assignment operators, which the
// vector_ref member definitions below delegate to.

template <typename LeftVector, typename RightVector>
LIBMESH_DEVICE_INLINE
auto operator+=(LeftVector & left, const RightVector & right)
  -> std::enable_if_t<is_vector_like_v<LeftVector> && is_vector_like_v<RightVector> &&
                        (is_vector_ref_v<LeftVector> || is_vector_ref_v<RightVector>),
                      LeftVector &>;

template <typename LeftVector, typename RightVector>
LIBMESH_DEVICE_INLINE
auto operator-=(LeftVector & left, const RightVector & right)
  -> std::enable_if_t<is_vector_like_v<LeftVector> && is_vector_like_v<RightVector> &&
                        (is_vector_ref_v<LeftVector> || is_vector_ref_v<RightVector>),
                      LeftVector &>;

// Out-of-line vector_ref member definitions.  Several delegate to the
// generic operations above or the operators below, so the whole group
// lives here rather than in kokkos_linalg_base.h; see the class
// definition for their documentation.

template <typename ViewType>
template <typename RightVector>
LIBMESH_DEVICE_INLINE
void vector_ref<ViewType>::assign(const RightVector & right)
{
  detail::assign_vector_components(*this, right);
}

template <typename ViewType>
template <typename RightVector>
LIBMESH_DEVICE_INLINE
void vector_ref<ViewType>::add(const RightVector & right)
{
  libMesh::Kokkos::operator+=(*this, right);
}

template <typename ViewType>
template <typename RightVector>
LIBMESH_DEVICE_INLINE
void vector_ref<ViewType>::add_scaled(const RightVector & right, const value_type & factor)
{
  detail::update_vector_components(*this, right, factor);
}

template <typename ViewType>
template <typename RightVector>
LIBMESH_DEVICE_INLINE
void vector_ref<ViewType>::subtract(const RightVector & right)
{
  libMesh::Kokkos::operator-=(*this, right);
}

template <typename ViewType>
template <typename RightVector>
LIBMESH_DEVICE_INLINE
void vector_ref<ViewType>::subtract_scaled(const RightVector & right, const value_type & factor)
{
  detail::update_vector_components(*this, right, -factor);
}

template <typename ViewType>
LIBMESH_DEVICE_INLINE
void vector_ref<ViewType>::zero()
{
  detail::fill_vector_components(*this, value_type(0));
}

template <typename ViewType>
template <typename RightVector>
LIBMESH_DEVICE_INLINE
auto vector_ref<ViewType>::contract(const RightVector & right) const
{
  return libMesh::Kokkos::contract(*this, right);
}

template <typename ViewType>
LIBMESH_DEVICE_INLINE
auto vector_ref<ViewType>::norm() const
{
  using std::sqrt;
  return sqrt(this->norm_sq());
}

template <typename ViewType>
LIBMESH_DEVICE_INLINE
auto vector_ref<ViewType>::norm_sq() const
{
  // TensorTools::norm_sq handles complex-valued components correctly.
  using norm_type = detail::remove_cvref_t<decltype(libMesh::TensorTools::norm_sq((*this)(0)))>;

  norm_type sum = norm_type(0);
  for (unsigned int component = 0; component < LIBMESH_DIM; ++component)
    sum += libMesh::TensorTools::norm_sq((*this)(component));

  return sum;
}

template <typename ViewType>
LIBMESH_DEVICE_INLINE
auto vector_ref<ViewType>::l1_norm() const
{
  using std::abs;
  using norm_type = detail::remove_cvref_t<decltype(abs((*this)(0)))>;

  norm_type sum = norm_type(0);
  for (unsigned int component = 0; component < LIBMESH_DIM; ++component)
    sum += abs((*this)(component));

  return sum;
}

template <typename ViewType>
LIBMESH_DEVICE_INLINE
bool vector_ref<ViewType>::is_zero() const
{
  for (unsigned int component = 0; component < LIBMESH_DIM; ++component)
    if ((*this)(component) != value_type(0))
      return false;

  return true;
}

template <typename ViewType>
LIBMESH_DEVICE_INLINE
auto vector_ref<ViewType>::unit() const
{
  return vector_unit(*this);
}

template <typename ViewType>
template <typename RightVector>
LIBMESH_DEVICE_INLINE
auto vector_ref<ViewType>::cross(const RightVector & right) const
{
  return vector_cross(*this, right);
}

// Operator-compatible wrappers for storage-backed refs and mixed ref/owning math.
//
// Every overload requires at least one vector_ref operand (via
// is_vector_ref_v), so these never compete with the operators the
// owning types define for themselves.  Vector-valued results are
// returned as owning values (sums and differences use the left
// operand's semantic type); compound assignments update the left
// operand in place, writing through to storage when it is a ref.

/**
 * \returns The negative of the vector-like \p v, as an owning vector.
 */
template <typename VectorLike>
LIBMESH_DEVICE_INLINE
auto operator-(const VectorLike & v)
  -> std::enable_if_t<is_vector_like_v<VectorLike> && is_vector_ref_v<VectorLike>,
                      vector_semantic_type_t<VectorLike>>
{
  return detail::transformed_vector<vector_semantic_type_t<VectorLike>>(
    v,
    detail::negate_value<typename VectorLike::value_type>{});
}

/**
 * \returns The sum of \p left and \p right, as an owning vector.
 */
template <typename LeftVector, typename RightVector>
LIBMESH_DEVICE_INLINE
auto operator+(const LeftVector & left, const RightVector & right)
  -> std::enable_if_t<is_vector_like_v<LeftVector> && is_vector_like_v<RightVector> &&
                        (is_vector_ref_v<LeftVector> || is_vector_ref_v<RightVector>),
                      vector_semantic_type_t<LeftVector>>
{
  auto out = copy_vector<vector_semantic_type_t<LeftVector>>(left);
  out += right;
  return out;
}

/**
 * \returns The difference of \p left and \p right, as an owning vector.
 */
template <typename LeftVector, typename RightVector>
LIBMESH_DEVICE_INLINE
auto operator-(const LeftVector & left, const RightVector & right)
  -> std::enable_if_t<is_vector_like_v<LeftVector> && is_vector_like_v<RightVector> &&
                        (is_vector_ref_v<LeftVector> || is_vector_ref_v<RightVector>),
                      vector_semantic_type_t<LeftVector>>
{
  auto out = copy_vector<vector_semantic_type_t<LeftVector>>(left);
  out -= right;
  return out;
}

/**
 * \returns The dot product of \p left and \p right, matching the
 * TypeVector convention that vector * vector is a dot product.
 */
template <typename LeftVector,
          typename RightVector,
          typename std::enable_if<is_vector_like_v<LeftVector> && is_vector_like_v<RightVector> &&
                                    (is_vector_ref_v<LeftVector> || is_vector_ref_v<RightVector>),
                                  int>::type = 0>
LIBMESH_DEVICE_INLINE
auto operator*(const LeftVector & left, const RightVector & right)
{
  return vector_dot(left, right);
}

/**
 * \returns The vector-like \p v scaled by \p alpha, as an owning
 * vector.  The scalar may appear on either side of the product.
 */
template <typename Scalar,
          typename VectorLike,
          typename std::enable_if<!is_vector_like_v<Scalar> && !is_tensor_like_v<Scalar> &&
                                    is_vector_like_v<VectorLike> && is_vector_ref_v<VectorLike>,
                                  int>::type = 0>
LIBMESH_DEVICE_INLINE
auto operator*(const Scalar & alpha, const VectorLike & v)
{
  return v * alpha;
}

template <typename VectorLike,
          typename Scalar,
          typename std::enable_if<is_vector_like_v<VectorLike> && is_vector_ref_v<VectorLike> &&
                                    !is_vector_like_v<Scalar> && !is_tensor_like_v<Scalar>,
                                  int>::type = 0>
LIBMESH_DEVICE_INLINE
auto operator*(const VectorLike & v, const Scalar & alpha)
{
  return detail::transformed_vector<vector_semantic_type_t<VectorLike>>(
    v,
    detail::scale_value<Scalar>{alpha});
}

/**
 * \returns The vector-like \p v divided (componentwise) by the scalar
 * \p alpha, as an owning vector.
 */
template <typename VectorLike,
          typename Scalar,
          typename std::enable_if<is_vector_like_v<VectorLike> && is_vector_ref_v<VectorLike> &&
                                    !is_vector_like_v<Scalar> && !is_tensor_like_v<Scalar>,
                                  int>::type = 0>
LIBMESH_DEVICE_INLINE
auto operator/(const VectorLike & v, const Scalar & alpha)
{
  return detail::transformed_vector<vector_semantic_type_t<VectorLike>>(
    v,
    detail::divide_value<Scalar>{alpha});
}

/**
 * \returns \p true if \p left and \p right agree in every component,
 * \p false otherwise.  Each component is compared with exact
 * operator==, matching TypeVector::operator==; this is not a
 * tolerance-based floating-point comparison.
 */
template <typename LeftVector, typename RightVector>
LIBMESH_DEVICE_INLINE
auto operator==(const LeftVector & left, const RightVector & right)
  -> std::enable_if_t<is_vector_like_v<LeftVector> && is_vector_like_v<RightVector> &&
                        (is_vector_ref_v<LeftVector> || is_vector_ref_v<RightVector>),
                      bool>
{
  return detail::vector_equal_impl(left, right);
}

/**
 * \returns \p true if \p left and \p right differ in any component,
 * \p false otherwise.
 */
template <typename LeftVector, typename RightVector>
LIBMESH_DEVICE_INLINE
auto operator!=(const LeftVector & left, const RightVector & right)
  -> std::enable_if_t<is_vector_like_v<LeftVector> && is_vector_like_v<RightVector> &&
                        (is_vector_ref_v<LeftVector> || is_vector_ref_v<RightVector>),
                      bool>
{
  return !(left == right);
}

/**
 * Adds \p right to \p left in place; when \p left is a ref, the update
 * writes through to the underlying storage.
 */
template <typename LeftVector, typename RightVector>
LIBMESH_DEVICE_INLINE
auto operator+=(LeftVector & left, const RightVector & right)
  -> std::enable_if_t<is_vector_like_v<LeftVector> && is_vector_like_v<RightVector> &&
                        (is_vector_ref_v<LeftVector> || is_vector_ref_v<RightVector>),
                      LeftVector &>
{
  detail::update_vector_components(left, right, 1);
  return left;
}

/**
 * Subtracts \p right from \p left in place.
 */
template <typename LeftVector, typename RightVector>
LIBMESH_DEVICE_INLINE
auto operator-=(LeftVector & left, const RightVector & right)
  -> std::enable_if_t<is_vector_like_v<LeftVector> && is_vector_like_v<RightVector> &&
                        (is_vector_ref_v<LeftVector> || is_vector_ref_v<RightVector>),
                      LeftVector &>
{
  detail::update_vector_components(left, right, -1);
  return left;
}

/**
 * Scales \p left by \p alpha in place.
 */
template <typename LeftVector, typename Scalar>
LIBMESH_DEVICE_INLINE
auto operator*=(LeftVector & left, const Scalar & alpha)
  -> std::enable_if_t<is_vector_like_v<LeftVector> && is_vector_ref_v<LeftVector> &&
                        !is_vector_like_v<Scalar> && !is_tensor_like_v<Scalar>,
                      LeftVector &>
{
  detail::transform_vector_components(left, left, detail::scale_value<Scalar>{alpha});
  return left;
}

/**
 * Divides \p left by \p alpha in place.
 */
template <typename LeftVector, typename Scalar>
LIBMESH_DEVICE_INLINE
auto operator/=(LeftVector & left, const Scalar & alpha)
  -> std::enable_if_t<is_vector_like_v<LeftVector> && is_vector_ref_v<LeftVector> &&
                        !is_vector_like_v<Scalar> && !is_tensor_like_v<Scalar>,
                      LeftVector &>
{
  detail::transform_vector_components(left, left, detail::divide_value<Scalar>{alpha});
  return left;
}

} // namespace libMesh::Kokkos

#endif // LIBMESH_KOKKOS_VECTOR_OPS_H
