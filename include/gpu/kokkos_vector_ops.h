// libMesh Kokkos generic vector operations.
//
// These free functions build vector algebra on top of the primitive
// access/materialization layer in kokkos_linalg_base.h. They are written
// against vector-like inputs so both libMesh owning types and storage-backed
// refs can participate in the same math.

#ifndef LIBMESH_KOKKOS_VECTOR_OPS_H
#define LIBMESH_KOKKOS_VECTOR_OPS_H

#include "libmesh/kokkos_linalg_base.h"

#include "libmesh/tensor_tools.h"

#include <cmath>

namespace libMesh::Kokkos
{

// Construction and materialization

template <typename ResultVector>
LIBMESH_DEVICE_INLINE
ResultVector zero_vector_value()
{
  ResultVector out;
  out.zero();
  return out;
}

template <typename ResultVector, typename VectorLike>
LIBMESH_DEVICE_INLINE
ResultVector copy_vector(const VectorLike & v)
{
  return materialize_vector<ResultVector>(v);
}

template <typename Dummy = void,
          typename VectorLike,
          typename std::enable_if<std::is_void<Dummy>::value, int>::type = 0>
LIBMESH_DEVICE_INLINE
vector_semantic_type_t<VectorLike> copy_vector(const VectorLike & v)
{
  return copy_vector<vector_semantic_type_t<VectorLike>>(v);
}

namespace detail
{

template <typename LeftVector, typename RightVector>
LIBMESH_DEVICE_INLINE
void assign_vector_components(LeftVector & left, const RightVector & right)
{
  for (unsigned int component = 0; component < LIBMESH_DIM; ++component)
    vector_set_component(left, component, vector_get_component(right, component));
}

template <typename LeftVector, typename RightVector>
LIBMESH_DEVICE_INLINE
void add_vector_components(LeftVector & left, const RightVector & right)
{
  for (unsigned int component = 0; component < LIBMESH_DIM; ++component)
    vector_set_component(left,
                         component,
                         vector_get_component(left, component) + vector_get_component(right, component));
}

template <typename LeftVector, typename RightVector, typename Scalar>
LIBMESH_DEVICE_INLINE
void add_scaled_vector_components(LeftVector & left, const RightVector & right, const Scalar & factor)
{
  for (unsigned int component = 0; component < LIBMESH_DIM; ++component)
    vector_set_component(left,
                         component,
                         vector_get_component(left, component) +
                           factor * vector_get_component(right, component));
}

template <typename LeftVector, typename RightVector>
LIBMESH_DEVICE_INLINE
void subtract_vector_components(LeftVector & left, const RightVector & right)
{
  for (unsigned int component = 0; component < LIBMESH_DIM; ++component)
    vector_set_component(left,
                         component,
                         vector_get_component(left, component) - vector_get_component(right, component));
}

template <typename LeftVector, typename RightVector, typename Scalar>
LIBMESH_DEVICE_INLINE
void subtract_scaled_vector_components(LeftVector & left, const RightVector & right, const Scalar & factor)
{
  for (unsigned int component = 0; component < LIBMESH_DIM; ++component)
    vector_set_component(left,
                         component,
                         vector_get_component(left, component) -
                           factor * vector_get_component(right, component));
}

template <typename VectorLike>
LIBMESH_DEVICE_INLINE
void zero_vector_components(VectorLike & v)
{
  for (unsigned int component = 0; component < LIBMESH_DIM; ++component)
    vector_set_component(v, component, vector_value_type_t<VectorLike>(0));
}

template <typename VectorLike, typename Scalar>
LIBMESH_DEVICE_INLINE
void scale_vector_components(VectorLike & v, const Scalar & alpha)
{
  for (unsigned int component = 0; component < LIBMESH_DIM; ++component)
    vector_set_component(v, component, vector_get_component(v, component) * alpha);
}

template <typename VectorLike, typename Scalar>
LIBMESH_DEVICE_INLINE
void divide_vector_components(VectorLike & v, const Scalar & alpha)
{
  for (unsigned int component = 0; component < LIBMESH_DIM; ++component)
    vector_set_component(v, component, vector_get_component(v, component) / alpha);
}

} // namespace detail

// Reductions and predicates

template <typename LeftVector, typename RightVector>
LIBMESH_DEVICE_INLINE
auto vector_dot(const LeftVector & left, const RightVector & right)
{
  static_assert(is_vector_like_v<LeftVector>, "vector_dot() requires a vector-like left input");
  static_assert(is_vector_like_v<RightVector>, "vector_dot() requires a vector-like right input");

  using sum_type =
    detail::remove_cvref_t<decltype(vector_get_component(left, 0) * vector_get_component(right, 0))>;

  sum_type sum = sum_type(0);
  for (unsigned int component = 0; component < LIBMESH_DIM; ++component)
    sum += vector_get_component(left, component) * vector_get_component(right, component);

  return sum;
}

template <typename VectorLike>
LIBMESH_DEVICE_INLINE
auto vector_norm_sq(const VectorLike & v)
{
  static_assert(is_vector_like_v<VectorLike>, "vector_norm_sq() requires a vector-like input");

  using norm_type = detail::remove_cvref_t<decltype(libMesh::TensorTools::norm_sq(vector_get_component(v, 0)))>;

  norm_type sum = norm_type(0);
  for (unsigned int component = 0; component < LIBMESH_DIM; ++component)
    sum += libMesh::TensorTools::norm_sq(vector_get_component(v, component));

  return sum;
}

template <typename VectorLike>
LIBMESH_DEVICE_INLINE
auto vector_norm(const VectorLike & v)
{
  using std::sqrt;
  return sqrt(vector_norm_sq(v));
}

template <typename VectorLike>
LIBMESH_DEVICE_INLINE
auto vector_l1_norm(const VectorLike & v)
{
  static_assert(is_vector_like_v<VectorLike>, "vector_l1_norm() requires a vector-like input");

  using std::abs;
  using norm_type = detail::remove_cvref_t<decltype(abs(vector_get_component(v, 0)))>;

  norm_type sum = norm_type(0);
  for (unsigned int component = 0; component < LIBMESH_DIM; ++component)
    sum += abs(vector_get_component(v, component));

  return sum;
}

template <typename VectorLike>
LIBMESH_DEVICE_INLINE
bool vector_is_zero(const VectorLike & v)
{
  static_assert(is_vector_like_v<VectorLike>, "vector_is_zero() requires a vector-like input");

  for (unsigned int component = 0; component < LIBMESH_DIM; ++component)
    if (vector_get_component(v, component) != vector_value_type_t<VectorLike>(0))
      return false;

  return true;
}

// Arithmetic

template <typename ResultVector, typename LeftVector, typename RightVector>
LIBMESH_DEVICE_INLINE
ResultVector vector_add(const LeftVector & left, const RightVector & right)
{
  ResultVector out;
  out.zero();

  for (unsigned int component = 0; component < LIBMESH_DIM; ++component)
    vector_set_component(out,
                         component,
                         vector_get_component(left, component) + vector_get_component(right, component));

  return out;
}

template <typename Dummy = void,
          typename LeftVector,
          typename RightVector,
          typename std::enable_if<std::is_void<Dummy>::value, int>::type = 0>
LIBMESH_DEVICE_INLINE
vector_semantic_type_t<LeftVector> vector_add(const LeftVector & left, const RightVector & right)
{
  return vector_add<vector_semantic_type_t<LeftVector>>(left, right);
}

template <typename ResultVector, typename LeftVector, typename RightVector>
LIBMESH_DEVICE_INLINE
ResultVector vector_subtract(const LeftVector & left, const RightVector & right)
{
  ResultVector out;
  out.zero();

  for (unsigned int component = 0; component < LIBMESH_DIM; ++component)
    vector_set_component(out,
                         component,
                         vector_get_component(left, component) - vector_get_component(right, component));

  return out;
}

template <typename Dummy = void,
          typename LeftVector,
          typename RightVector,
          typename std::enable_if<std::is_void<Dummy>::value, int>::type = 0>
LIBMESH_DEVICE_INLINE
vector_semantic_type_t<LeftVector> vector_subtract(const LeftVector & left, const RightVector & right)
{
  return vector_subtract<vector_semantic_type_t<LeftVector>>(left, right);
}

template <typename ResultVector, typename Scalar, typename VectorLike>
LIBMESH_DEVICE_INLINE
ResultVector vector_scale(const Scalar & alpha, const VectorLike & v)
{
  ResultVector out;
  out.zero();

  for (unsigned int component = 0; component < LIBMESH_DIM; ++component)
    vector_set_component(out, component, alpha * vector_get_component(v, component));

  return out;
}

template <typename Dummy = void,
          typename Scalar,
          typename VectorLike,
          typename std::enable_if<std::is_void<Dummy>::value, int>::type = 0>
LIBMESH_DEVICE_INLINE
vector_semantic_type_t<VectorLike> vector_scale(const Scalar & alpha, const VectorLike & v)
{
  return vector_scale<vector_semantic_type_t<VectorLike>>(alpha, v);
}

template <typename ResultVector, typename VectorLike, typename Scalar>
LIBMESH_DEVICE_INLINE
ResultVector vector_divide(const VectorLike & v, const Scalar & alpha)
{
  ResultVector out;
  out.zero();

  for (unsigned int component = 0; component < LIBMESH_DIM; ++component)
    vector_set_component(out, component, vector_get_component(v, component) / alpha);

  return out;
}

template <typename Dummy = void,
          typename VectorLike,
          typename Scalar,
          typename std::enable_if<std::is_void<Dummy>::value, int>::type = 0>
LIBMESH_DEVICE_INLINE
vector_semantic_type_t<VectorLike> vector_divide(const VectorLike & v, const Scalar & alpha)
{
  return vector_divide<vector_semantic_type_t<VectorLike>>(v, alpha);
}

template <typename ResultVector, typename ScalarA, typename VectorA, typename ScalarB, typename VectorB>
LIBMESH_DEVICE_INLINE
ResultVector vector_linear_combination(const ScalarA & alpha,
                                       const VectorA & a,
                                       const ScalarB & beta,
                                       const VectorB & b)
{
  ResultVector out;
  out.zero();

  for (unsigned int component = 0; component < LIBMESH_DIM; ++component)
    vector_set_component(out,
                         component,
                         alpha * vector_get_component(a, component) +
                           beta * vector_get_component(b, component));

  return out;
}

template <typename Dummy = void,
          typename ScalarA,
          typename VectorA,
          typename ScalarB,
          typename VectorB,
          typename std::enable_if<std::is_void<Dummy>::value, int>::type = 0>
LIBMESH_DEVICE_INLINE
vector_semantic_type_t<VectorA> vector_linear_combination(const ScalarA & alpha,
                                                          const VectorA & a,
                                                          const ScalarB & beta,
                                                          const VectorB & b)
{
  return vector_linear_combination<vector_semantic_type_t<VectorA>>(alpha, a, beta, b);
}

template <typename ResultVector,
          typename ScalarA,
          typename VectorA,
          typename ScalarB,
          typename VectorB,
          typename ScalarC,
          typename VectorC>
LIBMESH_DEVICE_INLINE
ResultVector vector_linear_combination(const ScalarA & alpha,
                                       const VectorA & a,
                                       const ScalarB & beta,
                                       const VectorB & b,
                                       const ScalarC & gamma,
                                       const VectorC & c)
{
  ResultVector out;
  out.zero();

  for (unsigned int component = 0; component < LIBMESH_DIM; ++component)
    vector_set_component(out,
                         component,
                         alpha * vector_get_component(a, component) +
                           beta * vector_get_component(b, component) +
                           gamma * vector_get_component(c, component));

  return out;
}

template <typename Dummy = void,
          typename ScalarA,
          typename VectorA,
          typename ScalarB,
          typename VectorB,
          typename ScalarC,
          typename VectorC,
          typename std::enable_if<std::is_void<Dummy>::value, int>::type = 0>
LIBMESH_DEVICE_INLINE
vector_semantic_type_t<VectorA> vector_linear_combination(const ScalarA & alpha,
                                                          const VectorA & a,
                                                          const ScalarB & beta,
                                                          const VectorB & b,
                                                          const ScalarC & gamma,
                                                          const VectorC & c)
{
  return vector_linear_combination<vector_semantic_type_t<VectorA>>(alpha, a, beta, b, gamma, c);
}

template <typename ResultVector, typename VectorLike>
LIBMESH_DEVICE_INLINE
ResultVector vector_unit(const VectorLike & v)
{
  const auto length = vector_norm(v);
  libmesh_assert_not_equal_to(length, static_cast<Real>(0.));
  return vector_divide<ResultVector>(v, length);
}

template <typename Dummy = void,
          typename VectorLike,
          typename std::enable_if<std::is_void<Dummy>::value, int>::type = 0>
LIBMESH_DEVICE_INLINE
vector_semantic_type_t<VectorLike> vector_unit(const VectorLike & v)
{
  return vector_unit<vector_semantic_type_t<VectorLike>>(v);
}

// Geometry

template <typename ResultVector, typename LeftVector, typename RightVector>
LIBMESH_DEVICE_INLINE
ResultVector vector_cross(const LeftVector & left, const RightVector & right)
{
  ResultVector out;
  out.zero();

#if LIBMESH_DIM == 3
  vector_set_component(out,
                       0,
                       vector_get_component(left, 1) * vector_get_component(right, 2) -
                         vector_get_component(left, 2) * vector_get_component(right, 1));
  vector_set_component(out,
                       1,
                       -vector_get_component(left, 0) * vector_get_component(right, 2) +
                         vector_get_component(left, 2) * vector_get_component(right, 0));
  vector_set_component(out,
                       2,
                       vector_get_component(left, 0) * vector_get_component(right, 1) -
                         vector_get_component(left, 1) * vector_get_component(right, 0));
#else
  libmesh_ignore(left);
  libmesh_ignore(right);
#endif

  return out;
}

template <typename Dummy = void,
          typename LeftVector,
          typename RightVector,
          typename std::enable_if<std::is_void<Dummy>::value, int>::type = 0>
LIBMESH_DEVICE_INLINE
vector_semantic_type_t<LeftVector> vector_cross(const LeftVector & left, const RightVector & right)
{
  return vector_cross<vector_semantic_type_t<LeftVector>>(left, right);
}

template <typename LeftVector, typename MiddleVector, typename RightVector>
LIBMESH_DEVICE_INLINE
auto vector_triple_product(const LeftVector & left,
                           const MiddleVector & middle,
                           const RightVector & right)
{
#if LIBMESH_DIM == 3
  return vector_get_component(left, 0) *
           (vector_get_component(middle, 1) * vector_get_component(right, 2) -
            vector_get_component(middle, 2) * vector_get_component(right, 1)) -
         vector_get_component(left, 1) *
           (vector_get_component(middle, 0) * vector_get_component(right, 2) -
            vector_get_component(middle, 2) * vector_get_component(right, 0)) +
         vector_get_component(left, 2) *
           (vector_get_component(middle, 0) * vector_get_component(right, 1) -
            vector_get_component(middle, 1) * vector_get_component(right, 0));
#else
  libmesh_ignore(left, middle, right);
  using value_type =
    detail::remove_cvref_t<decltype(vector_get_component(left, 0) * vector_get_component(middle, 0))>;
  return value_type(0);
#endif
}

template <typename LeftVector, typename RightVector>
LIBMESH_DEVICE_INLINE
auto vector_cross_norm_sq(const LeftVector & left, const RightVector & right)
{
  const auto z = vector_get_component(left, 0) * vector_get_component(right, 1) -
                 vector_get_component(left, 1) * vector_get_component(right, 0);

#if LIBMESH_DIM == 3
  const auto x = vector_get_component(left, 1) * vector_get_component(right, 2) -
                 vector_get_component(left, 2) * vector_get_component(right, 1);
  const auto y = vector_get_component(left, 0) * vector_get_component(right, 2) -
                 vector_get_component(left, 2) * vector_get_component(right, 0);
  return x * x + y * y + z * z;
#else
  return z * z;
#endif
}

template <typename VectorA, typename VectorB, typename VectorC>
LIBMESH_DEVICE_INLINE
auto vector_solid_angle(const VectorA & v01, const VectorB & v02, const VectorC & v03)
{
  using std::atan;

  const auto norm01 = vector_norm(v01);
  const auto norm02 = vector_norm(v02);
  const auto norm03 = vector_norm(v03);
  const auto tan_half_angle =
    vector_triple_product(v01, v02, v03) /
    (vector_dot(v01, v02) * norm03 +
     vector_dot(v01, v03) * norm02 +
     vector_dot(v02, v03) * norm01 +
     norm01 * norm02 * norm03);

  return Real(2) * atan(tan_half_angle);
}

// libMesh-like convenience wrappers

template <typename LeftVector, typename RightVector>
LIBMESH_DEVICE_INLINE
auto contract(const LeftVector & left, const RightVector & right)
  -> std::enable_if_t<is_vector_like_v<LeftVector> && is_vector_like_v<RightVector>,
                      decltype(vector_dot(left, right))>
{
  return vector_dot(left, right);
}

template <typename VectorLike>
LIBMESH_DEVICE_INLINE
auto norm_sq(const VectorLike & v)
  -> std::enable_if_t<is_vector_like_v<VectorLike>, decltype(vector_norm_sq(v))>
{
  return vector_norm_sq(v);
}

template <typename VectorLike>
LIBMESH_DEVICE_INLINE
auto norm(const VectorLike & v)
  -> std::enable_if_t<is_vector_like_v<VectorLike>, decltype(vector_norm(v))>
{
  return vector_norm(v);
}

template <typename VectorLike>
LIBMESH_DEVICE_INLINE
auto is_zero(const VectorLike & v)
  -> std::enable_if_t<is_vector_like_v<VectorLike>, bool>
{
  return vector_is_zero(v);
}

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
  detail::add_vector_components(*this, right);
}

template <typename ViewType>
template <typename RightVector>
LIBMESH_DEVICE_INLINE
void vector_ref<ViewType>::add_scaled(const RightVector & right, const value_type & factor)
{
  detail::add_scaled_vector_components(*this, right, factor);
}

template <typename ViewType>
template <typename RightVector>
LIBMESH_DEVICE_INLINE
void vector_ref<ViewType>::subtract(const RightVector & right)
{
  detail::subtract_vector_components(*this, right);
}

template <typename ViewType>
template <typename RightVector>
LIBMESH_DEVICE_INLINE
void vector_ref<ViewType>::subtract_scaled(const RightVector & right, const value_type & factor)
{
  detail::subtract_scaled_vector_components(*this, right, factor);
}

template <typename ViewType>
LIBMESH_DEVICE_INLINE
void vector_ref<ViewType>::zero()
{
  detail::zero_vector_components(*this);
}

template <typename ViewType>
template <typename RightVector>
LIBMESH_DEVICE_INLINE
auto vector_ref<ViewType>::contract(const RightVector & right) const
{
  return vector_dot(*this, right);
}

template <typename ViewType>
LIBMESH_DEVICE_INLINE
auto vector_ref<ViewType>::norm() const
{
  return vector_norm(*this);
}

template <typename ViewType>
LIBMESH_DEVICE_INLINE
auto vector_ref<ViewType>::norm_sq() const
{
  return vector_norm_sq(*this);
}

template <typename ViewType>
LIBMESH_DEVICE_INLINE
auto vector_ref<ViewType>::l1_norm() const
{
  return vector_l1_norm(*this);
}

template <typename ViewType>
LIBMESH_DEVICE_INLINE
bool vector_ref<ViewType>::is_zero() const
{
  return vector_is_zero(*this);
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

template <typename VectorLike>
LIBMESH_DEVICE_INLINE
auto operator-(const VectorLike & v)
  -> std::enable_if_t<is_vector_like_v<VectorLike> && is_vector_ref_v<VectorLike>,
                      vector_semantic_type_t<VectorLike>>
{
  return vector_scale(vector_value_type_t<VectorLike>(-1), v);
}

template <typename LeftVector, typename RightVector>
LIBMESH_DEVICE_INLINE
auto operator+(const LeftVector & left, const RightVector & right)
  -> std::enable_if_t<is_vector_like_v<LeftVector> && is_vector_like_v<RightVector> &&
                        (is_vector_ref_v<LeftVector> || is_vector_ref_v<RightVector>),
                      vector_semantic_type_t<LeftVector>>
{
  return vector_add(left, right);
}

template <typename LeftVector, typename RightVector>
LIBMESH_DEVICE_INLINE
auto operator-(const LeftVector & left, const RightVector & right)
  -> std::enable_if_t<is_vector_like_v<LeftVector> && is_vector_like_v<RightVector> &&
                        (is_vector_ref_v<LeftVector> || is_vector_ref_v<RightVector>),
                      vector_semantic_type_t<LeftVector>>
{
  return vector_subtract(left, right);
}

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

template <typename Scalar,
          typename VectorLike,
          typename std::enable_if<!is_vector_like_v<Scalar> && !is_tensor_like_v<Scalar> &&
                                    is_vector_like_v<VectorLike> && is_vector_ref_v<VectorLike>,
                                  int>::type = 0>
LIBMESH_DEVICE_INLINE
auto operator*(const Scalar & alpha, const VectorLike & v)
{
  return vector_scale(alpha, v);
}

template <typename VectorLike,
          typename Scalar,
          typename std::enable_if<is_vector_like_v<VectorLike> && is_vector_ref_v<VectorLike> &&
                                    !is_vector_like_v<Scalar> && !is_tensor_like_v<Scalar>,
                                  int>::type = 0>
LIBMESH_DEVICE_INLINE
auto operator*(const VectorLike & v, const Scalar & alpha)
{
  return vector_scale(alpha, v);
}

template <typename VectorLike,
          typename Scalar,
          typename std::enable_if<is_vector_like_v<VectorLike> && is_vector_ref_v<VectorLike> &&
                                    !is_vector_like_v<Scalar> && !is_tensor_like_v<Scalar>,
                                  int>::type = 0>
LIBMESH_DEVICE_INLINE
auto operator/(const VectorLike & v, const Scalar & alpha)
{
  return vector_divide(v, alpha);
}

template <typename LeftVector, typename RightVector>
LIBMESH_DEVICE_INLINE
auto operator==(const LeftVector & left, const RightVector & right)
  -> std::enable_if_t<is_vector_like_v<LeftVector> && is_vector_like_v<RightVector> &&
                        (is_vector_ref_v<LeftVector> || is_vector_ref_v<RightVector>),
                      bool>
{
  for (unsigned int component = 0; component < LIBMESH_DIM; ++component)
    if (vector_get_component(left, component) != vector_get_component(right, component))
      return false;

  return true;
}

template <typename LeftVector, typename RightVector>
LIBMESH_DEVICE_INLINE
auto operator!=(const LeftVector & left, const RightVector & right)
  -> std::enable_if_t<is_vector_like_v<LeftVector> && is_vector_like_v<RightVector> &&
                        (is_vector_ref_v<LeftVector> || is_vector_ref_v<RightVector>),
                      bool>
{
  return !(left == right);
}

template <typename LeftVector, typename RightVector>
LIBMESH_DEVICE_INLINE
auto operator+=(LeftVector & left, const RightVector & right)
  -> std::enable_if_t<is_vector_like_v<LeftVector> && is_vector_like_v<RightVector> &&
                        (is_vector_ref_v<LeftVector> || is_vector_ref_v<RightVector>),
                      LeftVector &>
{
  detail::add_vector_components(left, right);
  return left;
}

template <typename LeftVector, typename RightVector>
LIBMESH_DEVICE_INLINE
auto operator-=(LeftVector & left, const RightVector & right)
  -> std::enable_if_t<is_vector_like_v<LeftVector> && is_vector_like_v<RightVector> &&
                        (is_vector_ref_v<LeftVector> || is_vector_ref_v<RightVector>),
                      LeftVector &>
{
  detail::subtract_vector_components(left, right);
  return left;
}

template <typename LeftVector, typename Scalar>
LIBMESH_DEVICE_INLINE
auto operator*=(LeftVector & left, const Scalar & alpha)
  -> std::enable_if_t<is_vector_like_v<LeftVector> && is_vector_ref_v<LeftVector> &&
                        !is_vector_like_v<Scalar> && !is_tensor_like_v<Scalar>,
                      LeftVector &>
{
  detail::scale_vector_components(left, alpha);
  return left;
}

template <typename LeftVector, typename Scalar>
LIBMESH_DEVICE_INLINE
auto operator/=(LeftVector & left, const Scalar & alpha)
  -> std::enable_if_t<is_vector_like_v<LeftVector> && is_vector_ref_v<LeftVector> &&
                        !is_vector_like_v<Scalar> && !is_tensor_like_v<Scalar>,
                      LeftVector &>
{
  detail::divide_vector_components(left, alpha);
  return left;
}

} // namespace libMesh::Kokkos

#endif // LIBMESH_KOKKOS_VECTOR_OPS_H
