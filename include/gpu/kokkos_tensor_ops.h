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
ResultTensor zero_tensor_value()
{
  ResultTensor out;
  out.zero();
  return out;
}

template <typename ResultTensor>
LIBMESH_DEVICE_INLINE
ResultTensor tensor_identity(const unsigned int dim = LIBMESH_DIM)
{
  ResultTensor out;
  out.zero();

  for (unsigned int i = 0; i < dim; ++i)
    tensor_set_component(out, i, i, tensor_value_type_t<ResultTensor>(1));

  return out;
}

template <typename ResultTensor, typename TensorLike>
LIBMESH_DEVICE_INLINE
ResultTensor copy_tensor(const TensorLike & T_in)
{
  return materialize_tensor<ResultTensor>(T_in);
}

template <typename Dummy = void,
          typename TensorLike,
          typename std::enable_if<std::is_void<Dummy>::value, int>::type = 0>
LIBMESH_DEVICE_INLINE
tensor_semantic_type_t<TensorLike> copy_tensor(const TensorLike & T_in)
{
  return copy_tensor<tensor_semantic_type_t<TensorLike>>(T_in);
}

namespace detail
{

template <typename TensorLike>
LIBMESH_DEVICE_INLINE
auto leading_determinant(const TensorLike & T_in, const unsigned int dim = LIBMESH_DIM)
{
  static_assert(is_tensor_like_v<TensorLike>,
                "detail::leading_determinant() requires a tensor-like input");

  if (dim == 0)
    return tensor_value_type_t<TensorLike>(1);

  if (dim == 1)
    return tensor_get_component(T_in, 0, 0);

  if (dim == 2)
    return tensor_get_component(T_in, 0, 0) * tensor_get_component(T_in, 1, 1) -
           tensor_get_component(T_in, 0, 1) * tensor_get_component(T_in, 1, 0);

#if LIBMESH_DIM > 2
  const auto a00 = tensor_get_component(T_in, 0, 0);
  const auto a01 = tensor_get_component(T_in, 0, 1);
  const auto a02 = tensor_get_component(T_in, 0, 2);
  const auto a10 = tensor_get_component(T_in, 1, 0);
  const auto a11 = tensor_get_component(T_in, 1, 1);
  const auto a12 = tensor_get_component(T_in, 1, 2);
  const auto a20 = tensor_get_component(T_in, 2, 0);
  const auto a21 = tensor_get_component(T_in, 2, 1);
  const auto a22 = tensor_get_component(T_in, 2, 2);

  return a00 * (a11 * a22 - a12 * a21) -
         a01 * (a10 * a22 - a12 * a20) +
         a02 * (a10 * a21 - a11 * a20);
#else
  libmesh_ignore(T_in);
  return tensor_value_type_t<TensorLike>(0);
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
      tensor_set_component(out,
                           row,
                           col,
                           vector_get_component(left, row) * libmesh_conj(vector_get_component(right, col)));

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
    tensor_set_component(out, 0, 0, tensor_value_type_t<ResultTensor>(1) / tensor_get_component(T_in, 0, 0));
    return out;
  }

  const auto det = leading_determinant(T_in, dim);

  if (dim == 2)
  {
    tensor_set_component(out, 0, 0,  tensor_get_component(T_in, 1, 1) / det);
    tensor_set_component(out, 0, 1, -tensor_get_component(T_in, 0, 1) / det);
    tensor_set_component(out, 1, 0, -tensor_get_component(T_in, 1, 0) / det);
    tensor_set_component(out, 1, 1,  tensor_get_component(T_in, 0, 0) / det);
    return out;
  }

#if LIBMESH_DIM > 2
  const auto a00 = tensor_get_component(T_in, 0, 0);
  const auto a01 = tensor_get_component(T_in, 0, 1);
  const auto a02 = tensor_get_component(T_in, 0, 2);
  const auto a10 = tensor_get_component(T_in, 1, 0);
  const auto a11 = tensor_get_component(T_in, 1, 1);
  const auto a12 = tensor_get_component(T_in, 1, 2);
  const auto a20 = tensor_get_component(T_in, 2, 0);
  const auto a21 = tensor_get_component(T_in, 2, 1);
  const auto a22 = tensor_get_component(T_in, 2, 2);

  tensor_set_component(out, 0, 0, (a11 * a22 - a12 * a21) / det);
  tensor_set_component(out, 0, 1, (a02 * a21 - a01 * a22) / det);
  tensor_set_component(out, 0, 2, (a01 * a12 - a02 * a11) / det);
  tensor_set_component(out, 1, 0, (a12 * a20 - a10 * a22) / det);
  tensor_set_component(out, 1, 1, (a00 * a22 - a02 * a20) / det);
  tensor_set_component(out, 1, 2, (a02 * a10 - a00 * a12) / det);
  tensor_set_component(out, 2, 0, (a10 * a21 - a11 * a20) / det);
  tensor_set_component(out, 2, 1, (a01 * a20 - a00 * a21) / det);
  tensor_set_component(out, 2, 2, (a00 * a11 - a01 * a10) / det);
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
      tensor_set_component(out, row, col, tensor_get_component(T_in, col, row));

  return out;
}

template <typename ResultTensor, typename LeftTensor, typename RightTensor>
LIBMESH_DEVICE_INLINE
ResultTensor multiply_tensors(const LeftTensor & left, const RightTensor & right)
{
  ResultTensor out;
  out.zero();

  for (unsigned int row = 0; row < LIBMESH_DIM; ++row)
    for (unsigned int col = 0; col < LIBMESH_DIM; ++col)
    {
      auto value = tensor_get_component(left, row, 0) * tensor_get_component(right, 0, col);
      for (unsigned int k = 1; k < LIBMESH_DIM; ++k)
        value += tensor_get_component(left, row, k) * tensor_get_component(right, k, col);
      tensor_set_component(out, row, col, value);
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
    vector_set_component(out, col, tensor_get_component(T_in, row_index, col));

  return out;
}

template <typename ResultVector, typename TensorLike>
LIBMESH_DEVICE_INLINE
ResultVector column(const TensorLike & T_in, const unsigned int col_index)
{
  ResultVector out;
  out.zero();

  for (unsigned int row_index = 0; row_index < LIBMESH_DIM; ++row_index)
    vector_set_component(out, row_index, tensor_get_component(T_in, row_index, col_index));

  return out;
}

template <typename ResultVector, typename TensorLike, typename VectorLike>
LIBMESH_DEVICE_INLINE
ResultVector multiply_tensor_vector(const TensorLike & T_in, const VectorLike & v)
{
  ResultVector out;
  out.zero();

  for (unsigned int row = 0; row < LIBMESH_DIM; ++row)
  {
    auto value = tensor_get_component(T_in, row, 0) * vector_get_component(v, 0);
    for (unsigned int col = 1; col < LIBMESH_DIM; ++col)
      value += tensor_get_component(T_in, row, col) * vector_get_component(v, col);
    vector_set_component(out, row, value);
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
    auto value = vector_get_component(v, 0) * tensor_get_component(T_in, 0, col);
    for (unsigned int row = 1; row < LIBMESH_DIM; ++row)
      value += vector_get_component(v, row) * tensor_get_component(T_in, row, col);
    vector_set_component(out, col, value);
  }

  return out;
}

template <typename LeftTensor, typename RightTensor>
LIBMESH_DEVICE_INLINE
void assign_tensor_components(LeftTensor & left, const RightTensor & right)
{
  for (unsigned int row = 0; row < LIBMESH_DIM; ++row)
    for (unsigned int col = 0; col < LIBMESH_DIM; ++col)
      tensor_set_component(left, row, col, tensor_get_component(right, row, col));
}

template <typename TensorLike, typename Scalar>
LIBMESH_DEVICE_INLINE
void fill_tensor_components(TensorLike & T_in, const Scalar & value)
{
  for (unsigned int row = 0; row < LIBMESH_DIM; ++row)
    for (unsigned int col = 0; col < LIBMESH_DIM; ++col)
      tensor_set_component(T_in, row, col, value);
}

template <typename LeftTensor, typename RightTensor, typename Scalar>
LIBMESH_DEVICE_INLINE
void update_tensor_components(LeftTensor & left, const RightTensor & right, const Scalar & factor)
{
  for (unsigned int row = 0; row < LIBMESH_DIM; ++row)
    for (unsigned int col = 0; col < LIBMESH_DIM; ++col)
      tensor_set_component(left,
                           row,
                           col,
                           tensor_get_component(left, row, col) +
                             factor * tensor_get_component(right, row, col));
}

template <typename ResultTensor, typename ScalarA, typename TensorA, typename ScalarB, typename TensorB>
LIBMESH_DEVICE_INLINE
ResultTensor linear_combination(const ScalarA & alpha,
                                const TensorA & A,
                                const ScalarB & beta,
                                const TensorB & B)
{
  ResultTensor out;
  out.zero();

  for (unsigned int row = 0; row < LIBMESH_DIM; ++row)
    for (unsigned int col = 0; col < LIBMESH_DIM; ++col)
      tensor_set_component(out,
                           row,
                           col,
                           alpha * tensor_get_component(A, row, col) +
                             beta * tensor_get_component(B, row, col));

  return out;
}

template <typename ResultTensor, typename Scalar, typename TensorLike>
LIBMESH_DEVICE_INLINE
ResultTensor scale_tensor(const Scalar & alpha, const TensorLike & T_in)
{
  ResultTensor out;
  out.zero();

  for (unsigned int row = 0; row < LIBMESH_DIM; ++row)
    for (unsigned int col = 0; col < LIBMESH_DIM; ++col)
      tensor_set_component(out, row, col, alpha * tensor_get_component(T_in, row, col));

  return out;
}

template <typename ResultTensor, typename TensorLike, typename Scalar>
LIBMESH_DEVICE_INLINE
ResultTensor divide_tensor(const TensorLike & T_in, const Scalar & alpha)
{
  ResultTensor out;
  out.zero();

  for (unsigned int row = 0; row < LIBMESH_DIM; ++row)
    for (unsigned int col = 0; col < LIBMESH_DIM; ++col)
      tensor_set_component(out, row, col, tensor_get_component(T_in, row, col) / alpha);

  return out;
}

template <typename TensorLike, typename Scalar>
LIBMESH_DEVICE_INLINE
void scale_tensor_components(TensorLike & T_in, const Scalar & alpha)
{
  for (unsigned int row = 0; row < LIBMESH_DIM; ++row)
    for (unsigned int col = 0; col < LIBMESH_DIM; ++col)
      tensor_set_component(T_in, row, col, tensor_get_component(T_in, row, col) * alpha);
}

template <typename TensorLike, typename Scalar>
LIBMESH_DEVICE_INLINE
void divide_tensor_components(TensorLike & T_in, const Scalar & alpha)
{
  for (unsigned int row = 0; row < LIBMESH_DIM; ++row)
    for (unsigned int col = 0; col < LIBMESH_DIM; ++col)
      tensor_set_component(T_in, row, col, tensor_get_component(T_in, row, col) / alpha);
}

// Tensor reductions and predicates

template <typename LeftTensor, typename RightTensor>
LIBMESH_DEVICE_INLINE
auto tensor_contract(const LeftTensor & left, const RightTensor & right)
{
  static_assert(is_tensor_like_v<LeftTensor>, "tensor_contract() requires a tensor-like left input");
  static_assert(is_tensor_like_v<RightTensor>, "tensor_contract() requires a tensor-like right input");

  using sum_type =
    detail::remove_cvref_t<decltype(tensor_get_component(left, 0, 0) * tensor_get_component(right, 0, 0))>;

  sum_type sum = sum_type(0);
  for (unsigned int row = 0; row < LIBMESH_DIM; ++row)
    for (unsigned int col = 0; col < LIBMESH_DIM; ++col)
      sum += tensor_get_component(left, row, col) * tensor_get_component(right, row, col);

  return sum;
}

template <typename TensorLike>
LIBMESH_DEVICE_INLINE
auto tensor_norm_sq(const TensorLike & T_in)
{
  static_assert(is_tensor_like_v<TensorLike>, "tensor_norm_sq() requires a tensor-like input");

  using norm_type = detail::remove_cvref_t<decltype(libMesh::TensorTools::norm_sq(tensor_get_component(T_in, 0, 0)))>;

  norm_type sum = norm_type(0);
  for (unsigned int row = 0; row < LIBMESH_DIM; ++row)
    for (unsigned int col = 0; col < LIBMESH_DIM; ++col)
      sum += libMesh::TensorTools::norm_sq(tensor_get_component(T_in, row, col));

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
auto tensor_trace(const TensorLike & T_in)
{
  static_assert(is_tensor_like_v<TensorLike>, "tensor_trace() requires a tensor-like input");

  using trace_type = detail::remove_cvref_t<decltype(tensor_get_component(T_in, 0, 0))>;
  trace_type sum = trace_type(0);
  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
    sum += tensor_get_component(T_in, i, i);

  return sum;
}

template <typename TensorLike>
LIBMESH_DEVICE_INLINE
bool tensor_is_zero(const TensorLike & T_in)
{
  static_assert(is_tensor_like_v<TensorLike>, "tensor_is_zero() requires a tensor-like input");

  for (unsigned int row = 0; row < LIBMESH_DIM; ++row)
    for (unsigned int col = 0; col < LIBMESH_DIM; ++col)
      if (tensor_get_component(T_in, row, col) != tensor_value_type_t<TensorLike>(0))
        return false;

  return true;
}

} // namespace detail

// libMesh-like convenience wrappers

template <typename LeftTensor, typename RightTensor>
LIBMESH_DEVICE_INLINE
auto contract(const LeftTensor & left, const RightTensor & right)
  -> std::enable_if_t<is_tensor_like_v<LeftTensor> && is_tensor_like_v<RightTensor>,
                      decltype(detail::tensor_contract(left, right))>
{
  return detail::tensor_contract(left, right);
}

template <typename TensorLike>
LIBMESH_DEVICE_INLINE
auto norm_sq(const TensorLike & T_in)
  -> std::enable_if_t<is_tensor_like_v<TensorLike>, decltype(detail::tensor_norm_sq(T_in))>
{
  return detail::tensor_norm_sq(T_in);
}

template <typename TensorLike>
LIBMESH_DEVICE_INLINE
auto norm(const TensorLike & T_in)
  -> std::enable_if_t<is_tensor_like_v<TensorLike>, decltype(detail::tensor_norm(T_in))>
{
  return detail::tensor_norm(T_in);
}

template <typename TensorLike>
LIBMESH_DEVICE_INLINE
auto is_zero(const TensorLike & T_in)
  -> std::enable_if_t<is_tensor_like_v<TensorLike>, bool>
{
  return detail::tensor_is_zero(T_in);
}

template <typename ResultTensor, typename LeftVector, typename RightVector>
LIBMESH_DEVICE_INLINE
auto outer_product(const LeftVector & left, const RightVector & right)
  -> std::enable_if_t<is_vector_like_v<LeftVector> && is_vector_like_v<RightVector>, ResultTensor>
{
  return detail::outer_product<ResultTensor>(left, right);
}

template <typename LeftVector, typename RightVector>
LIBMESH_DEVICE_INLINE
auto outer_product(const LeftVector & left, const RightVector & right)
  -> std::enable_if_t<is_vector_like_v<LeftVector> && is_vector_like_v<RightVector>,
                      libMesh::TypeTensor<vector_value_type_t<LeftVector>>>
{
  return outer_product<libMesh::TypeTensor<vector_value_type_t<LeftVector>>>(left, right);
}

template <typename ResultTensor, typename TensorLike>
LIBMESH_DEVICE_INLINE
auto transpose(const TensorLike & T_in)
  -> std::enable_if_t<is_tensor_like_v<TensorLike>, ResultTensor>
{
  return detail::transpose<ResultTensor>(T_in);
}

template <typename TensorLike>
LIBMESH_DEVICE_INLINE
auto transpose(const TensorLike & T_in)
  -> std::enable_if_t<is_tensor_like_v<TensorLike>, tensor_semantic_type_t<TensorLike>>
{
  return transpose<tensor_semantic_type_t<TensorLike>>(T_in);
}

template <typename TensorLike>
LIBMESH_DEVICE_INLINE
auto det(const TensorLike & T_in)
  -> std::enable_if_t<is_tensor_like_v<TensorLike>, decltype(T_in.det())>
{
  return T_in.det();
}

template <typename ResultTensor, typename TensorLike>
LIBMESH_DEVICE_INLINE
auto inverse(const TensorLike & T_in, const unsigned int dim = LIBMESH_DIM)
  -> std::enable_if_t<is_tensor_like_v<TensorLike>, ResultTensor>
{
  return detail::inverse<ResultTensor>(T_in, dim);
}

template <typename TensorLike>
LIBMESH_DEVICE_INLINE
auto inverse(const TensorLike & T_in, const unsigned int dim = LIBMESH_DIM)
  -> std::enable_if_t<is_tensor_like_v<TensorLike>, tensor_semantic_type_t<TensorLike>>
{
  return inverse<tensor_semantic_type_t<TensorLike>>(T_in, dim);
}

template <typename ResultVector, typename TensorLike>
LIBMESH_DEVICE_INLINE
auto row(const TensorLike & T_in, const unsigned int i)
  -> std::enable_if_t<is_tensor_like_v<TensorLike>, ResultVector>
{
  return detail::row<ResultVector>(T_in, i);
}

template <typename TensorLike>
LIBMESH_DEVICE_INLINE
auto row(const TensorLike & T_in, const unsigned int i)
  -> std::enable_if_t<is_tensor_like_v<TensorLike>, libMesh::TypeVector<tensor_value_type_t<TensorLike>>>
{
  return row<libMesh::TypeVector<tensor_value_type_t<TensorLike>>>(T_in, i);
}

template <typename ResultVector, typename TensorLike>
LIBMESH_DEVICE_INLINE
auto column(const TensorLike & T_in, const unsigned int i)
  -> std::enable_if_t<is_tensor_like_v<TensorLike>, ResultVector>
{
  return detail::column<ResultVector>(T_in, i);
}

template <typename TensorLike>
LIBMESH_DEVICE_INLINE
auto column(const TensorLike & T_in, const unsigned int i)
  -> std::enable_if_t<is_tensor_like_v<TensorLike>, libMesh::TypeVector<tensor_value_type_t<TensorLike>>>
{
  return column<libMesh::TypeVector<tensor_value_type_t<TensorLike>>>(T_in, i);
}

template <typename ViewType>
template <typename RightTensor>
LIBMESH_DEVICE_INLINE
void tensor_ref<ViewType>::assign(const RightTensor & right)
{
  detail::assign_tensor_components(*this, right);
}

template <typename ViewType>
template <typename RightTensor>
LIBMESH_DEVICE_INLINE
void tensor_ref<ViewType>::add(const RightTensor & right)
{
  detail::update_tensor_components(*this, right, value_type(1));
}

template <typename ViewType>
template <typename RightTensor>
LIBMESH_DEVICE_INLINE
void tensor_ref<ViewType>::add_scaled(const RightTensor & right, const value_type & factor)
{
  detail::update_tensor_components(*this, right, factor);
}

template <typename ViewType>
template <typename RightTensor>
LIBMESH_DEVICE_INLINE
void tensor_ref<ViewType>::subtract(const RightTensor & right)
{
  detail::update_tensor_components(*this, right, value_type(-1));
}

template <typename ViewType>
template <typename RightTensor>
LIBMESH_DEVICE_INLINE
void tensor_ref<ViewType>::subtract_scaled(const RightTensor & right, const value_type & factor)
{
  detail::update_tensor_components(*this, right, -factor);
}

template <typename ViewType>
LIBMESH_DEVICE_INLINE
void tensor_ref<ViewType>::zero()
{
  detail::fill_tensor_components(*this, value_type(0));
}

template <typename ViewType>
template <typename RightTensor>
LIBMESH_DEVICE_INLINE
auto tensor_ref<ViewType>::contract(const RightTensor & right) const
{
  return detail::tensor_contract(*this, right);
}

template <typename ViewType>
LIBMESH_DEVICE_INLINE
auto tensor_ref<ViewType>::norm() const
{
  return detail::tensor_norm(*this);
}

template <typename ViewType>
LIBMESH_DEVICE_INLINE
auto tensor_ref<ViewType>::norm_sq() const
{
  return detail::tensor_norm_sq(*this);
}

template <typename ViewType>
LIBMESH_DEVICE_INLINE
bool tensor_ref<ViewType>::is_zero() const
{
  return detail::tensor_is_zero(*this);
}

template <typename ViewType>
LIBMESH_DEVICE_INLINE
auto tensor_ref<ViewType>::transpose() const
{
  return libMesh::Kokkos::transpose(*this);
}

template <typename ViewType>
LIBMESH_DEVICE_INLINE
auto tensor_ref<ViewType>::det(const unsigned int dim) const
{
  return detail::leading_determinant(*this, dim);
}

template <typename ViewType>
LIBMESH_DEVICE_INLINE
auto tensor_ref<ViewType>::tr() const
{
  return detail::tensor_trace(*this);
}

template <typename ViewType>
LIBMESH_DEVICE_INLINE
auto tensor_ref<ViewType>::inverse(const unsigned int dim) const
{
  return libMesh::Kokkos::inverse(*this, dim);
}

template <typename ViewType>
template <typename VectorLike, typename ResultVector>
LIBMESH_DEVICE_INLINE
void tensor_ref<ViewType>::solve(const VectorLike & b, ResultVector & x) const
{
  const auto solution =
    detail::multiply_tensor_vector<vector_semantic_type_t<ResultVector>>(this->inverse(), b);
  for (unsigned int component = 0; component < LIBMESH_DIM; ++component)
    vector_set_component(x, component, vector_get_component(solution, component));
}

template <typename ViewType>
LIBMESH_DEVICE_INLINE
auto tensor_ref<ViewType>::row(const unsigned int i) const
{
  return libMesh::Kokkos::row(*this, i);
}

template <typename ViewType>
LIBMESH_DEVICE_INLINE
auto tensor_ref<ViewType>::column(const unsigned int i) const
{
  return libMesh::Kokkos::column(*this, i);
}

template <typename ViewType>
template <typename VectorLike>
LIBMESH_DEVICE_INLINE
auto tensor_ref<ViewType>::left_multiply(const VectorLike & v) const
{
  return v * *this;
}

// Operator-compatible wrappers for storage-backed refs and mixed ref/owning math.

template <typename TensorLike>
LIBMESH_DEVICE_INLINE
auto operator-(const TensorLike & T_in)
  -> std::enable_if_t<is_tensor_like_v<TensorLike> && is_tensor_ref_v<TensorLike>,
                      tensor_semantic_type_t<TensorLike>>
{
  return detail::scale_tensor<tensor_semantic_type_t<TensorLike>>(tensor_value_type_t<TensorLike>(-1), T_in);
}

template <typename LeftTensor, typename RightTensor>
LIBMESH_DEVICE_INLINE
auto operator+(const LeftTensor & left, const RightTensor & right)
  -> std::enable_if_t<is_tensor_like_v<LeftTensor> && is_tensor_like_v<RightTensor> &&
                        (is_tensor_ref_v<LeftTensor> || is_tensor_ref_v<RightTensor>),
                      tensor_semantic_type_t<LeftTensor>>
{
  return detail::linear_combination<tensor_semantic_type_t<LeftTensor>>(
    tensor_value_type_t<LeftTensor>(1), left, tensor_value_type_t<LeftTensor>(1), right);
}

template <typename LeftTensor, typename RightTensor>
LIBMESH_DEVICE_INLINE
auto operator-(const LeftTensor & left, const RightTensor & right)
  -> std::enable_if_t<is_tensor_like_v<LeftTensor> && is_tensor_like_v<RightTensor> &&
                        (is_tensor_ref_v<LeftTensor> || is_tensor_ref_v<RightTensor>),
                      tensor_semantic_type_t<LeftTensor>>
{
  return detail::linear_combination<tensor_semantic_type_t<LeftTensor>>(
    tensor_value_type_t<LeftTensor>(1), left, tensor_value_type_t<LeftTensor>(-1), right);
}

template <typename Scalar,
          typename TensorLike,
          typename std::enable_if<!is_vector_like_v<Scalar> && !is_tensor_like_v<Scalar> &&
                                    is_tensor_like_v<TensorLike> && is_tensor_ref_v<TensorLike>,
                                  int>::type = 0>
LIBMESH_DEVICE_INLINE
auto operator*(const Scalar & alpha, const TensorLike & T_in)
{
  return detail::scale_tensor<tensor_semantic_type_t<TensorLike>>(alpha, T_in);
}

template <typename TensorLike,
          typename Scalar,
          typename std::enable_if<is_tensor_like_v<TensorLike> && is_tensor_ref_v<TensorLike> &&
                                    !is_vector_like_v<Scalar> && !is_tensor_like_v<Scalar>,
                                  int>::type = 0>
LIBMESH_DEVICE_INLINE
auto operator*(const TensorLike & T_in, const Scalar & alpha)
{
  return detail::scale_tensor<tensor_semantic_type_t<TensorLike>>(alpha, T_in);
}

template <typename TensorLike, typename Scalar>
LIBMESH_DEVICE_INLINE
auto operator/(const TensorLike & T_in, const Scalar & alpha)
  -> std::enable_if_t<is_tensor_like_v<TensorLike> && is_tensor_ref_v<TensorLike> &&
                        !is_vector_like_v<Scalar> && !is_tensor_like_v<Scalar>,
                      tensor_semantic_type_t<TensorLike>>
{
  return detail::divide_tensor<tensor_semantic_type_t<TensorLike>>(T_in, alpha);
}

template <typename LeftTensor,
          typename RightTensor,
          typename std::enable_if<is_tensor_like_v<LeftTensor> && is_tensor_like_v<RightTensor> &&
                                    (is_tensor_ref_v<LeftTensor> || is_tensor_ref_v<RightTensor>),
                                  int>::type = 0>
LIBMESH_DEVICE_INLINE
auto operator*(const LeftTensor & left, const RightTensor & right)
{
  return detail::multiply_tensors<tensor_semantic_type_t<LeftTensor>>(left, right);
}

template <typename TensorLike,
          typename VectorLike,
          typename std::enable_if<is_tensor_like_v<TensorLike> && is_vector_like_v<VectorLike> &&
                                    (is_tensor_ref_v<TensorLike> || is_vector_ref_v<VectorLike>),
                                  int>::type = 0>
LIBMESH_DEVICE_INLINE
auto operator*(const TensorLike & T_in, const VectorLike & v)
{
  return detail::multiply_tensor_vector<vector_semantic_type_t<VectorLike>>(T_in, v);
}

template <typename VectorLike,
          typename TensorLike,
          typename std::enable_if<is_vector_like_v<VectorLike> && is_tensor_like_v<TensorLike> &&
                                    (is_vector_ref_v<VectorLike> || is_tensor_ref_v<TensorLike>),
                                  int>::type = 0>
LIBMESH_DEVICE_INLINE
auto operator*(const VectorLike & v, const TensorLike & T_in)
{
  return detail::multiply_vector_tensor<vector_semantic_type_t<VectorLike>>(v, T_in);
}

template <typename LeftTensor, typename RightTensor>
LIBMESH_DEVICE_INLINE
auto operator==(const LeftTensor & left, const RightTensor & right)
  -> std::enable_if_t<is_tensor_like_v<LeftTensor> && is_tensor_like_v<RightTensor> &&
                        (is_tensor_ref_v<LeftTensor> || is_tensor_ref_v<RightTensor>),
                      bool>
{
  for (unsigned int row = 0; row < LIBMESH_DIM; ++row)
    for (unsigned int col = 0; col < LIBMESH_DIM; ++col)
      if (tensor_get_component(left, row, col) != tensor_get_component(right, row, col))
        return false;

  return true;
}

template <typename LeftTensor, typename RightTensor>
LIBMESH_DEVICE_INLINE
auto operator!=(const LeftTensor & left, const RightTensor & right)
  -> std::enable_if_t<is_tensor_like_v<LeftTensor> && is_tensor_like_v<RightTensor> &&
                        (is_tensor_ref_v<LeftTensor> || is_tensor_ref_v<RightTensor>),
                      bool>
{
  return !(left == right);
}

template <typename LeftTensor, typename RightTensor>
LIBMESH_DEVICE_INLINE
auto operator+=(LeftTensor & left, const RightTensor & right)
  -> std::enable_if_t<is_tensor_like_v<LeftTensor> && is_tensor_like_v<RightTensor> &&
                        (is_tensor_ref_v<LeftTensor> || is_tensor_ref_v<RightTensor>),
                      LeftTensor &>
{
  detail::update_tensor_components(left, right, tensor_value_type_t<LeftTensor>(1));
  return left;
}

template <typename LeftTensor, typename RightTensor>
LIBMESH_DEVICE_INLINE
auto operator-=(LeftTensor & left, const RightTensor & right)
  -> std::enable_if_t<is_tensor_like_v<LeftTensor> && is_tensor_like_v<RightTensor> &&
                        (is_tensor_ref_v<LeftTensor> || is_tensor_ref_v<RightTensor>),
                      LeftTensor &>
{
  detail::update_tensor_components(left, right, tensor_value_type_t<LeftTensor>(-1));
  return left;
}

template <typename LeftTensor, typename Scalar>
LIBMESH_DEVICE_INLINE
auto operator*=(LeftTensor & left, const Scalar & alpha)
  -> std::enable_if_t<is_tensor_like_v<LeftTensor> && is_tensor_ref_v<LeftTensor> &&
                        !is_vector_like_v<Scalar> && !is_tensor_like_v<Scalar>,
                      LeftTensor &>
{
  detail::scale_tensor_components(left, alpha);
  return left;
}

template <typename LeftTensor, typename Scalar>
LIBMESH_DEVICE_INLINE
auto operator/=(LeftTensor & left, const Scalar & alpha)
  -> std::enable_if_t<is_tensor_like_v<LeftTensor> && is_tensor_ref_v<LeftTensor> &&
                        !is_vector_like_v<Scalar> && !is_tensor_like_v<Scalar>,
                      LeftTensor &>
{
  detail::divide_tensor_components(left, alpha);
  return left;
}

} // namespace libMesh::Kokkos

#endif // LIBMESH_KOKKOS_TENSOR_OPS_H
