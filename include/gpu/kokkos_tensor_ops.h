// libMesh Kokkos generic tensor operations.
//
// These free functions build tensor algebra on top of the primitive
// access/materialization layer in kokkos_linalg_base.h. They are written
// against tensor-like and vector-like inputs so both libMesh owning types and
// storage-backed refs can participate in the same math.

#ifndef LIBMESH_KOKKOS_TENSOR_OPS_H
#define LIBMESH_KOKKOS_TENSOR_OPS_H

#include "libmesh/kokkos_linalg_base.h"
#include "libmesh/kokkos_vector_ops.h"

#include "libmesh/tensor_tools.h"

#include <cmath>

namespace libMesh::Kokkos
{

// Construction and materialization

template <typename ResultTensor>
LIBMESH_DEVICE_INLINE
ResultTensor tensor_identity(const unsigned int dim = LIBMESH_DIM)
{
  ResultTensor out;
  out.zero();

  for (unsigned int i = 0; i < dim; ++i)
    out(i, i) = 1;

  return out;
}

namespace detail
{

// These helpers are shared by the public functions and ref operators so
// Kokkos-backed refs use direct component access without extra materialization.

template <typename TensorLike>
constexpr bool is_tensor_ref_unary_v =
  is_tensor_like_v<TensorLike> && is_tensor_ref_v<TensorLike>;

template <typename LeftTensor, typename RightTensor>
constexpr bool is_tensor_ref_binary_v =
  is_tensor_like_v<LeftTensor> && is_tensor_like_v<RightTensor> &&
  (is_tensor_ref_v<LeftTensor> || is_tensor_ref_v<RightTensor>);

template <typename Scalar, typename TensorLike>
constexpr bool is_scalar_tensor_ref_product_v =
  !is_vector_like_v<Scalar> && !is_tensor_like_v<Scalar> &&
  is_tensor_like_v<TensorLike> && is_tensor_ref_v<TensorLike>;

template <typename TensorLike, typename Scalar>
constexpr bool is_tensor_ref_scalar_product_v =
  is_tensor_like_v<TensorLike> && is_tensor_ref_v<TensorLike> &&
  !is_vector_like_v<Scalar> && !is_tensor_like_v<Scalar>;

template <typename TensorLike, typename VectorLike>
constexpr bool is_tensor_vector_ref_product_v =
  is_tensor_like_v<TensorLike> && is_vector_like_v<VectorLike> &&
  (is_tensor_ref_v<TensorLike> || is_vector_ref_v<VectorLike>);

template <typename VectorLike, typename TensorLike>
constexpr bool is_vector_tensor_ref_product_v =
  is_vector_like_v<VectorLike> && is_tensor_like_v<TensorLike> &&
  (is_vector_ref_v<VectorLike> || is_tensor_ref_v<TensorLike>);

template <typename TensorLike>
LIBMESH_DEVICE_INLINE
typename TensorLike::value_type
leading_determinant(const TensorLike & T_in, const unsigned int dim = LIBMESH_DIM)
{
  static_assert(is_tensor_like_v<TensorLike>,
                "detail::leading_determinant() requires a tensor-like input");

  if (dim == 0)
    return 1;

  if (dim == 1)
    return T_in(0, 0);

  if (dim == 2)
    return T_in(0, 0) * T_in(1, 1) -
           T_in(0, 1) * T_in(1, 0);

#if LIBMESH_DIM > 2
  const auto a00 = T_in(0, 0);
  const auto a01 = T_in(0, 1);
  const auto a02 = T_in(0, 2);
  const auto a10 = T_in(1, 0);
  const auto a11 = T_in(1, 1);
  const auto a12 = T_in(1, 2);
  const auto a20 = T_in(2, 0);
  const auto a21 = T_in(2, 1);
  const auto a22 = T_in(2, 2);

  return a00 * (a11 * a22 - a12 * a21) -
         a01 * (a10 * a22 - a12 * a20) +
         a02 * (a10 * a21 - a11 * a20);
#else
  libmesh_ignore(T_in);
  return 0;
#endif
}

template <typename ResultTensor, typename LeftVector, typename RightVector>
LIBMESH_DEVICE_INLINE
ResultTensor outer_product(const LeftVector & left, const RightVector & right)
{
  ResultTensor out;
  out.zero();

  for (unsigned int row = 0; row < LIBMESH_DIM; ++row)
    for (unsigned int col = 0; col < LIBMESH_DIM; ++col)
      out(row, col) = left(row) * libmesh_conj(right(col));

  return out;
}

template <typename ResultTensor, typename TensorLike>
LIBMESH_DEVICE_INLINE
ResultTensor inverse(const TensorLike & T_in, const unsigned int dim = LIBMESH_DIM)
{
  static_assert(is_tensor_like_v<TensorLike>, "detail::inverse() requires a tensor-like input");

  ResultTensor out;
  out.zero();

  if (dim == 0)
    return out;

  if (dim == 1)
  {
    using value_type = typename ResultTensor::value_type;
    out(0, 0) = value_type(1) / T_in(0, 0);
    return out;
  }

  const auto det = leading_determinant(T_in, dim);

  if (dim == 2)
  {
    out(0, 0) = T_in(1, 1) / det;
    out(0, 1) = -T_in(0, 1) / det;
    out(1, 0) = -T_in(1, 0) / det;
    out(1, 1) =  T_in(0, 0) / det;
    return out;
  }

#if LIBMESH_DIM > 2
  const auto a00 = T_in(0, 0);
  const auto a01 = T_in(0, 1);
  const auto a02 = T_in(0, 2);
  const auto a10 = T_in(1, 0);
  const auto a11 = T_in(1, 1);
  const auto a12 = T_in(1, 2);
  const auto a20 = T_in(2, 0);
  const auto a21 = T_in(2, 1);
  const auto a22 = T_in(2, 2);

  out(0, 0) = (a11 * a22 - a12 * a21) / det;
  out(0, 1) = (a02 * a21 - a01 * a22) / det;
  out(0, 2) = (a01 * a12 - a02 * a11) / det;
  out(1, 0) = (a12 * a20 - a10 * a22) / det;
  out(1, 1) = (a00 * a22 - a02 * a20) / det;
  out(1, 2) = (a02 * a10 - a00 * a12) / det;
  out(2, 0) = (a10 * a21 - a11 * a20) / det;
  out(2, 1) = (a01 * a20 - a00 * a21) / det;
  out(2, 2) = (a00 * a11 - a01 * a10) / det;
#else
  libmesh_ignore(T_in);
#endif

  return out;
}

template <typename ResultTensor, typename TensorLike>
LIBMESH_DEVICE_INLINE
ResultTensor transpose(const TensorLike & T_in)
{
  ResultTensor out;
  out.zero();

  for (unsigned int row = 0; row < LIBMESH_DIM; ++row)
    for (unsigned int col = 0; col < LIBMESH_DIM; ++col)
      out(row, col) = T_in(col, row);

  return out;
}

// Tensor/tensor product operators delegate here for the same direct-access path.

template <typename ResultTensor, typename LeftTensor, typename RightTensor>
LIBMESH_DEVICE_INLINE
ResultTensor multiply_tensors(const LeftTensor & left, const RightTensor & right)
{
  ResultTensor out;
  out.zero();

  for (unsigned int row = 0; row < LIBMESH_DIM; ++row)
    for (unsigned int col = 0; col < LIBMESH_DIM; ++col)
    {
      auto value = left(row, 0) * right(0, col);
      for (unsigned int k = 1; k < LIBMESH_DIM; ++k)
        value += left(row, k) * right(k, col);
      out(row, col) = value;
    }

  return out;
}

template <typename ResultVector, typename TensorLike>
LIBMESH_DEVICE_INLINE
ResultVector row(const TensorLike & T_in, const unsigned int row_index)
{
  ResultVector out;
  out.zero();

  for (unsigned int col = 0; col < LIBMESH_DIM; ++col)
    out(col) = T_in(row_index, col);

  return out;
}

template <typename ResultVector, typename TensorLike>
LIBMESH_DEVICE_INLINE
ResultVector column(const TensorLike & T_in, const unsigned int col_index)
{
  ResultVector out;
  out.zero();

  for (unsigned int row_index = 0; row_index < LIBMESH_DIM; ++row_index)
    out(row_index) = T_in(row_index, col_index);

  return out;
}

// Tensor/vector and vector/tensor product operators keep the direct-access path too.

template <typename ResultVector, typename TensorLike, typename VectorLike>
LIBMESH_DEVICE_INLINE
ResultVector multiply_tensor_vector(const TensorLike & T_in, const VectorLike & v)
{
  ResultVector out;
  out.zero();

  for (unsigned int row = 0; row < LIBMESH_DIM; ++row)
  {
    auto value = T_in(row, 0) * v(0);
    for (unsigned int col = 1; col < LIBMESH_DIM; ++col)
      value += T_in(row, col) * v(col);
    out(row) = value;
  }

  return out;
}

template <typename ResultVector, typename VectorLike, typename TensorLike>
LIBMESH_DEVICE_INLINE
ResultVector multiply_vector_tensor(const VectorLike & v, const TensorLike & T_in)
{
  ResultVector out;
  out.zero();

  for (unsigned int col = 0; col < LIBMESH_DIM; ++col)
  {
    auto value = v(0) * T_in(0, col);
    for (unsigned int row = 1; row < LIBMESH_DIM; ++row)
      value += v(row) * T_in(row, col);
    out(col) = value;
  }

  return out;
}

template <typename OutputTensor, typename InputTensor, typename TransformOp>
LIBMESH_DEVICE_INLINE
void transform_tensor_components(OutputTensor & out, const InputTensor & in, const TransformOp & op)
{
  for (unsigned int row = 0; row < LIBMESH_DIM; ++row)
    for (unsigned int col = 0; col < LIBMESH_DIM; ++col)
      out(row, col) = op(in(row, col));
}

template <typename ResultTensor, typename TensorLike, typename TransformOp>
LIBMESH_DEVICE_INLINE
ResultTensor transformed_tensor(const TensorLike & T_in, const TransformOp & op)
{
  ResultTensor out;
  out.zero();
  transform_tensor_components(out, T_in, op);
  return out;
}

// Tensor reductions and predicates

template <typename LeftTensor, typename RightTensor>
LIBMESH_DEVICE_INLINE
auto tensor_contract(const LeftTensor & left, const RightTensor & right)
{
  static_assert(is_tensor_like_v<LeftTensor>, "tensor_contract() requires a tensor-like left input");
  static_assert(is_tensor_like_v<RightTensor>, "tensor_contract() requires a tensor-like right input");

  using sum_type =
    detail::remove_cvref_t<decltype(left(0, 0) * right(0, 0))>;

  sum_type sum = sum_type(0);
  for (unsigned int row = 0; row < LIBMESH_DIM; ++row)
    for (unsigned int col = 0; col < LIBMESH_DIM; ++col)
      sum += left(row, col) * right(row, col);

  return sum;
}


template <typename TensorLike>
LIBMESH_DEVICE_INLINE
auto tensor_trace(const TensorLike & T_in)
{
  static_assert(is_tensor_like_v<TensorLike>, "tensor_trace() requires a tensor-like input");

  using trace_type = detail::remove_cvref_t<decltype(T_in(0, 0))>;
  trace_type sum = trace_type(0);
  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
    sum += T_in(i, i);

  return sum;
}

} // namespace detail

// libMesh-like convenience wrappers

template <typename LeftTensor,
          typename RightTensor,
          typename std::enable_if<is_tensor_like_v<LeftTensor> && is_tensor_like_v<RightTensor>,
                                  int>::type = 0>
LIBMESH_DEVICE_INLINE
auto contract(const LeftTensor & left, const RightTensor & right)
{
  return detail::tensor_contract(left, right);
}

template <typename TensorLike>
LIBMESH_DEVICE_INLINE
auto tensor_norm_sq(const TensorLike & T_in)
{
  static_assert(is_tensor_like_v<TensorLike>, "tensor_norm_sq() requires a tensor-like input");

  using norm_type =
    detail::remove_cvref_t<decltype(libMesh::TensorTools::norm_sq(T_in(0, 0)))>;

  norm_type sum = norm_type(0);
  for (unsigned int row_index = 0; row_index < LIBMESH_DIM; ++row_index)
    for (unsigned int col_index = 0; col_index < LIBMESH_DIM; ++col_index)
      sum += libMesh::TensorTools::norm_sq(T_in(row_index, col_index));

  return sum;
}

template <typename TensorLike>
LIBMESH_DEVICE_INLINE
auto tensor_norm(const TensorLike & T_in)
{
  using std::sqrt;
  return sqrt(tensor_norm_sq(T_in));
}

template <typename TensorLike>
LIBMESH_DEVICE_INLINE
bool tensor_is_zero(const TensorLike & T_in)
{
  static_assert(is_tensor_like_v<TensorLike>, "tensor_is_zero() requires a tensor-like input");

  using value_type = detail::remove_cvref_t<decltype(T_in(0, 0))>;
  for (unsigned int row_index = 0; row_index < LIBMESH_DIM; ++row_index)
    for (unsigned int col_index = 0; col_index < LIBMESH_DIM; ++col_index)
      if (T_in(row_index, col_index) != value_type(0))
        return false;

  return true;
}

template <typename TensorLike>
LIBMESH_DEVICE_INLINE
auto tr(const TensorLike & T_in)
  -> std::enable_if_t<is_tensor_like_v<TensorLike>,
                      decltype(detail::tensor_trace(T_in))>
{
  return detail::tensor_trace(T_in);
}

template <typename ResultTensor = void, typename LeftVector, typename RightVector>
LIBMESH_DEVICE_INLINE
auto outer_product(const LeftVector & left, const RightVector & right)
  -> std::enable_if_t<is_vector_like_v<LeftVector> && is_vector_like_v<RightVector>,
                      std::conditional_t<std::is_void<ResultTensor>::value,
                                         libMesh::TypeTensor<typename LeftVector::value_type>,
                                         ResultTensor>>
{
  using output_type = std::conditional_t<std::is_void<ResultTensor>::value,
                                         libMesh::TypeTensor<typename LeftVector::value_type>,
                                         ResultTensor>;
  return detail::outer_product<output_type>(left, right);
}

template <typename ResultTensor = void, typename TensorLike>
LIBMESH_DEVICE_INLINE
auto transpose(const TensorLike & T_in)
  -> std::enable_if_t<is_tensor_like_v<TensorLike>,
                      std::conditional_t<std::is_void<ResultTensor>::value,
                                         tensor_semantic_type_t<TensorLike>,
                                         ResultTensor>>
{
  using output_type = std::conditional_t<std::is_void<ResultTensor>::value,
                                         tensor_semantic_type_t<TensorLike>,
                                         ResultTensor>;
  return detail::transpose<output_type>(T_in);
}

template <typename TensorLike>
LIBMESH_DEVICE_INLINE
auto det(const TensorLike & T_in, const unsigned int dim = LIBMESH_DIM)
  -> std::enable_if_t<is_tensor_like_v<TensorLike>,
                      decltype(detail::leading_determinant(T_in, LIBMESH_DIM))>
{
  return detail::leading_determinant(T_in, dim);
}

template <typename ResultTensor = void, typename TensorLike>
LIBMESH_DEVICE_INLINE
auto inverse(const TensorLike & T_in, const unsigned int dim = LIBMESH_DIM)
  -> std::enable_if_t<is_tensor_like_v<TensorLike>,
                      std::conditional_t<std::is_void<ResultTensor>::value,
                                         tensor_semantic_type_t<TensorLike>,
                                         ResultTensor>>
{
  using output_type = std::conditional_t<std::is_void<ResultTensor>::value,
                                         tensor_semantic_type_t<TensorLike>,
                                         ResultTensor>;
  return detail::inverse<output_type>(T_in, dim);
}

template <typename ResultVector = void, typename TensorLike>
LIBMESH_DEVICE_INLINE
auto row(const TensorLike & T_in, const unsigned int i)
  -> std::enable_if_t<is_tensor_like_v<TensorLike>,
                      std::conditional_t<std::is_void<ResultVector>::value,
                                         libMesh::TypeVector<typename TensorLike::value_type>,
                                         ResultVector>>
{
  using output_type = std::conditional_t<std::is_void<ResultVector>::value,
                                         libMesh::TypeVector<typename TensorLike::value_type>,
                                         ResultVector>;
  return detail::row<output_type>(T_in, i);
}

template <typename ResultVector = void, typename TensorLike>
LIBMESH_DEVICE_INLINE
auto column(const TensorLike & T_in, const unsigned int i)
  -> std::enable_if_t<is_tensor_like_v<TensorLike>,
                      std::conditional_t<std::is_void<ResultVector>::value,
                                         libMesh::TypeVector<typename TensorLike::value_type>,
                                         ResultVector>>
{
  using output_type = std::conditional_t<std::is_void<ResultVector>::value,
                                         libMesh::TypeVector<typename TensorLike::value_type>,
                                         ResultVector>;
  return detail::column<output_type>(T_in, i);
}

// Operator-compatible wrappers for storage-backed refs and mixed ref/owning math.

template <typename TensorLike>
LIBMESH_DEVICE_INLINE
auto operator-(const TensorLike & T_in)
  -> std::enable_if_t<detail::is_tensor_ref_unary_v<TensorLike>,
                      tensor_semantic_type_t<TensorLike>>
{
  return detail::transformed_tensor<tensor_semantic_type_t<TensorLike>>(
    T_in,
    detail::negate_value<typename TensorLike::value_type>{});
}

template <typename LeftTensor, typename RightTensor>
LIBMESH_DEVICE_INLINE
auto operator+(const LeftTensor & left, const RightTensor & right)
  -> std::enable_if_t<detail::is_tensor_ref_binary_v<LeftTensor, RightTensor>,
                      tensor_semantic_type_t<LeftTensor>>
{
  auto out = materialize_tensor<tensor_semantic_type_t<LeftTensor>>(left);
  out += right;
  return out;
}

template <typename LeftTensor, typename RightTensor>
LIBMESH_DEVICE_INLINE
auto operator-(const LeftTensor & left, const RightTensor & right)
  -> std::enable_if_t<detail::is_tensor_ref_binary_v<LeftTensor, RightTensor>,
                      tensor_semantic_type_t<LeftTensor>>
{
  auto out = materialize_tensor<tensor_semantic_type_t<LeftTensor>>(left);
  out -= right;
  return out;
}

template <typename Scalar,
          typename TensorLike,
          typename std::enable_if<detail::is_scalar_tensor_ref_product_v<Scalar, TensorLike>,
                                  int>::type = 0>
LIBMESH_DEVICE_INLINE
auto operator*(const Scalar & alpha, const TensorLike & T_in)
{
  return T_in * alpha;
}

template <typename TensorLike,
          typename Scalar,
          typename std::enable_if<detail::is_tensor_ref_scalar_product_v<TensorLike, Scalar>,
                                  int>::type = 0>
LIBMESH_DEVICE_INLINE
auto operator*(const TensorLike & T_in, const Scalar & alpha)
{
  return detail::transformed_tensor<tensor_semantic_type_t<TensorLike>>(
    T_in,
    detail::scale_value<Scalar>{alpha});
}

template <typename TensorLike, typename Scalar>
LIBMESH_DEVICE_INLINE
auto operator/(const TensorLike & T_in, const Scalar & alpha)
  -> std::enable_if_t<detail::is_tensor_ref_scalar_product_v<TensorLike, Scalar>,
                      tensor_semantic_type_t<TensorLike>>
{
  return detail::transformed_tensor<tensor_semantic_type_t<TensorLike>>(
    T_in,
    detail::divide_value<Scalar>{alpha});
}

template <typename LeftTensor,
          typename RightTensor,
          typename std::enable_if<detail::is_tensor_ref_binary_v<LeftTensor, RightTensor>,
                                  int>::type = 0>
LIBMESH_DEVICE_INLINE
auto operator*(const LeftTensor & left, const RightTensor & right)
{
  return detail::multiply_tensors<tensor_semantic_type_t<LeftTensor>>(left, right);
}

template <typename TensorLike,
          typename VectorLike,
          typename std::enable_if<detail::is_tensor_vector_ref_product_v<TensorLike, VectorLike>,
                                  int>::type = 0>
LIBMESH_DEVICE_INLINE
auto operator*(const TensorLike & T_in, const VectorLike & v)
{
  return detail::multiply_tensor_vector<vector_semantic_type_t<VectorLike>>(T_in, v);
}

template <typename VectorLike,
          typename TensorLike,
          typename std::enable_if<detail::is_vector_tensor_ref_product_v<VectorLike, TensorLike>,
                                  int>::type = 0>
LIBMESH_DEVICE_INLINE
auto operator*(const VectorLike & v, const TensorLike & T_in)
{
  return detail::multiply_vector_tensor<vector_semantic_type_t<VectorLike>>(v, T_in);
}

template <typename LeftTensor, typename RightTensor>
LIBMESH_DEVICE_INLINE
auto operator==(const LeftTensor & left, const RightTensor & right)
  -> std::enable_if_t<detail::is_tensor_ref_binary_v<LeftTensor, RightTensor>,
                      bool>
{
  for (unsigned int row = 0; row < LIBMESH_DIM; ++row)
    for (unsigned int col = 0; col < LIBMESH_DIM; ++col)
      if (left(row, col) != right(row, col))
        return false;

  return true;
}

template <typename LeftTensor, typename RightTensor>
LIBMESH_DEVICE_INLINE
auto operator!=(const LeftTensor & left, const RightTensor & right)
  -> std::enable_if_t<detail::is_tensor_ref_binary_v<LeftTensor, RightTensor>,
                      bool>
{
  return !(left == right);
}

template <typename LeftTensor, typename RightTensor>
LIBMESH_DEVICE_INLINE
auto operator+=(LeftTensor & left, const RightTensor & right)
  -> std::enable_if_t<detail::is_tensor_ref_binary_v<LeftTensor, RightTensor>,
                      LeftTensor &>
{
  for (unsigned int row = 0; row < LIBMESH_DIM; ++row)
    for (unsigned int col = 0; col < LIBMESH_DIM; ++col)
      left(row, col) = left(row, col) + right(row, col);

  return left;
}

template <typename LeftTensor, typename RightTensor>
LIBMESH_DEVICE_INLINE
auto operator-=(LeftTensor & left, const RightTensor & right)
  -> std::enable_if_t<detail::is_tensor_ref_binary_v<LeftTensor, RightTensor>,
                      LeftTensor &>
{
  for (unsigned int row = 0; row < LIBMESH_DIM; ++row)
    for (unsigned int col = 0; col < LIBMESH_DIM; ++col)
      left(row, col) = left(row, col) - right(row, col);

  return left;
}

template <typename LeftTensor, typename Scalar>
LIBMESH_DEVICE_INLINE
auto operator*=(LeftTensor & left, const Scalar & alpha)
  -> std::enable_if_t<detail::is_tensor_ref_scalar_product_v<LeftTensor, Scalar>,
                      LeftTensor &>
{
  detail::transform_tensor_components(left, left, detail::scale_value<Scalar>{alpha});
  return left;
}

template <typename LeftTensor, typename Scalar>
LIBMESH_DEVICE_INLINE
auto operator/=(LeftTensor & left, const Scalar & alpha)
  -> std::enable_if_t<detail::is_tensor_ref_scalar_product_v<LeftTensor, Scalar>,
                      LeftTensor &>
{
  detail::transform_tensor_components(left, left, detail::divide_value<Scalar>{alpha});
  return left;
}

} // namespace libMesh::Kokkos

#endif // LIBMESH_KOKKOS_TENSOR_OPS_H
