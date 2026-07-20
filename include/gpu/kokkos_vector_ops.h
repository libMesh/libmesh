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

namespace detail
{

// These helpers are shared by the public functions and ref operators so
// Kokkos-backed refs use direct component access without extra materialization.

template <typename VectorLike>
constexpr bool is_vector_ref_unary_v =
  is_vector_like_v<VectorLike> && is_vector_ref_v<VectorLike>;

template <typename LeftVector, typename RightVector>
constexpr bool is_vector_ref_binary_v =
  is_vector_like_v<LeftVector> && is_vector_like_v<RightVector> &&
  (is_vector_ref_v<LeftVector> || is_vector_ref_v<RightVector>);

template <typename Scalar, typename VectorLike>
constexpr bool is_scalar_vector_ref_product_v =
  !is_vector_like_v<Scalar> && !is_tensor_like_v<Scalar> &&
  is_vector_like_v<VectorLike> && is_vector_ref_v<VectorLike>;

template <typename VectorLike, typename Scalar>
constexpr bool is_vector_ref_scalar_product_v =
  is_vector_like_v<VectorLike> && is_vector_ref_v<VectorLike> &&
  !is_vector_like_v<Scalar> && !is_tensor_like_v<Scalar>;

template <typename OutputVector, typename InputVector, typename TransformOp>
LIBMESH_DEVICE_INLINE
void transform_vector_components(OutputVector & out, const InputVector & in, const TransformOp & op)
{
  for (unsigned int component = 0; component < LIBMESH_DIM; ++component)
    out(component) = op(in(component));
}

template <typename ResultVector, typename VectorLike, typename TransformOp>
LIBMESH_DEVICE_INLINE
ResultVector transformed_vector(const VectorLike & v, const TransformOp & op)
{
  ResultVector out;
  out.zero();
  transform_vector_components(out, v, op);
  return out;
}

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

  template <typename ValueType>
  LIBMESH_DEVICE_INLINE
  auto operator()(const ValueType & value) const -> decltype(value / alpha)
  {
    return value / alpha;
  }
};

} // namespace detail

// Reductions and predicates

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

template <typename VectorLike>
LIBMESH_DEVICE_INLINE
auto vector_norm_sq(const VectorLike & v)
{
  static_assert(is_vector_like_v<VectorLike>, "vector_norm_sq() requires a vector-like input");

  using norm_type =
    detail::remove_cvref_t<decltype(libMesh::TensorTools::norm_sq(v(0)))>;

  norm_type sum = norm_type(0);
  for (unsigned int component = 0; component < LIBMESH_DIM; ++component)
    sum += libMesh::TensorTools::norm_sq(v(component));

  return sum;
}

template <typename VectorLike>
LIBMESH_DEVICE_INLINE
auto vector_norm(const VectorLike & v)
{
  using std::sqrt;
  return sqrt(vector_norm_sq(v));
}

template <typename ResultVector = void, typename VectorLike>
LIBMESH_DEVICE_INLINE
auto vector_unit(const VectorLike & v)
  -> std::conditional_t<std::is_void<ResultVector>::value,
                        vector_semantic_type_t<VectorLike>,
                        ResultVector>
{
  const auto length = vector_norm(v);
  libmesh_assert_not_equal_to(length, static_cast<Real>(0.));
  using output_type = std::conditional_t<std::is_void<ResultVector>::value,
                                         vector_semantic_type_t<VectorLike>,
                                         ResultVector>;
  return detail::transformed_vector<output_type>(v, detail::divide_value<decltype(length)>{length});
}

// Geometry

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

template <typename LeftVector,
          typename RightVector,
          typename std::enable_if<is_vector_like_v<LeftVector> && is_vector_like_v<RightVector>,
                                  int>::type = 0>
LIBMESH_DEVICE_INLINE
auto contract(const LeftVector & left, const RightVector & right)
{
  return vector_dot(left, right);
}


// Operator-compatible wrappers for storage-backed refs and mixed ref/owning math.

template <typename VectorLike>
LIBMESH_DEVICE_INLINE
auto operator-(const VectorLike & v)
  -> std::enable_if_t<detail::is_vector_ref_unary_v<VectorLike>,
                      vector_semantic_type_t<VectorLike>>
{
  return detail::transformed_vector<vector_semantic_type_t<VectorLike>>(
    v,
    detail::negate_value<typename VectorLike::value_type>{});
}

template <typename LeftVector, typename RightVector>
LIBMESH_DEVICE_INLINE
auto operator+(const LeftVector & left, const RightVector & right)
  -> std::enable_if_t<detail::is_vector_ref_binary_v<LeftVector, RightVector>,
                      vector_semantic_type_t<LeftVector>>
{
  auto out = materialize_vector<vector_semantic_type_t<LeftVector>>(left);
  out += right;
  return out;
}

template <typename LeftVector, typename RightVector>
LIBMESH_DEVICE_INLINE
auto operator-(const LeftVector & left, const RightVector & right)
  -> std::enable_if_t<detail::is_vector_ref_binary_v<LeftVector, RightVector>,
                      vector_semantic_type_t<LeftVector>>
{
  auto out = materialize_vector<vector_semantic_type_t<LeftVector>>(left);
  out -= right;
  return out;
}

template <typename LeftVector,
          typename RightVector,
          typename std::enable_if<detail::is_vector_ref_binary_v<LeftVector, RightVector>,
                                  int>::type = 0>
LIBMESH_DEVICE_INLINE
auto operator*(const LeftVector & left, const RightVector & right)
{
  return vector_dot(left, right);
}

template <typename Scalar,
          typename VectorLike,
          typename std::enable_if<detail::is_scalar_vector_ref_product_v<Scalar, VectorLike>,
                                  int>::type = 0>
LIBMESH_DEVICE_INLINE
auto operator*(const Scalar & alpha, const VectorLike & v)
{
  return v * alpha;
}

template <typename VectorLike,
          typename Scalar,
          typename std::enable_if<detail::is_vector_ref_scalar_product_v<VectorLike, Scalar>,
                                  int>::type = 0>
LIBMESH_DEVICE_INLINE
auto operator*(const VectorLike & v, const Scalar & alpha)
{
  return detail::transformed_vector<vector_semantic_type_t<VectorLike>>(
    v,
    detail::scale_value<Scalar>{alpha});
}

template <typename VectorLike,
          typename Scalar,
          typename std::enable_if<detail::is_vector_ref_scalar_product_v<VectorLike, Scalar>,
                                  int>::type = 0>
LIBMESH_DEVICE_INLINE
auto operator/(const VectorLike & v, const Scalar & alpha)
{
  return detail::transformed_vector<vector_semantic_type_t<VectorLike>>(
    v,
    detail::divide_value<Scalar>{alpha});
}

template <typename LeftVector, typename RightVector>
LIBMESH_DEVICE_INLINE
auto operator==(const LeftVector & left, const RightVector & right)
  -> std::enable_if_t<detail::is_vector_ref_binary_v<LeftVector, RightVector>,
                      bool>
{
  for (unsigned int component = 0; component < LIBMESH_DIM; ++component)
    if (left(component) != right(component))
      return false;

  return true;
}

template <typename LeftVector, typename RightVector>
LIBMESH_DEVICE_INLINE
auto operator!=(const LeftVector & left, const RightVector & right)
  -> std::enable_if_t<detail::is_vector_ref_binary_v<LeftVector, RightVector>,
                      bool>
{
  return !(left == right);
}

template <typename LeftVector, typename RightVector>
LIBMESH_DEVICE_INLINE
auto operator+=(LeftVector & left, const RightVector & right)
  -> std::enable_if_t<detail::is_vector_ref_binary_v<LeftVector, RightVector>,
                      LeftVector &>
{
  for (unsigned int component = 0; component < LIBMESH_DIM; ++component)
    left(component) = left(component) + right(component);

  return left;
}

template <typename LeftVector, typename RightVector>
LIBMESH_DEVICE_INLINE
auto operator-=(LeftVector & left, const RightVector & right)
  -> std::enable_if_t<detail::is_vector_ref_binary_v<LeftVector, RightVector>,
                      LeftVector &>
{
  for (unsigned int component = 0; component < LIBMESH_DIM; ++component)
    left(component) = left(component) - right(component);

  return left;
}

template <typename LeftVector, typename Scalar>
LIBMESH_DEVICE_INLINE
auto operator*=(LeftVector & left, const Scalar & alpha)
  -> std::enable_if_t<detail::is_vector_ref_scalar_product_v<LeftVector, Scalar>,
                      LeftVector &>
{
  detail::transform_vector_components(left, left, detail::scale_value<Scalar>{alpha});
  return left;
}

template <typename LeftVector, typename Scalar>
LIBMESH_DEVICE_INLINE
auto operator/=(LeftVector & left, const Scalar & alpha)
  -> std::enable_if_t<detail::is_vector_ref_scalar_product_v<LeftVector, Scalar>,
                      LeftVector &>
{
  detail::transform_vector_components(left, left, detail::divide_value<Scalar>{alpha});
  return left;
}

} // namespace libMesh::Kokkos

#endif // LIBMESH_KOKKOS_VECTOR_OPS_H
