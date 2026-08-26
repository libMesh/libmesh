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

#ifndef LIBMESH_KOKKOS_TENSOR_OPS_H
#define LIBMESH_KOKKOS_TENSOR_OPS_H

// libMesh Kokkos generic tensor operations.
//
// These free functions build tensor algebra on top of the primitive
// access/materialization layer in kokkos_linalg_base.h. They are written
// against tensor-like and vector-like inputs so both libMesh owning types and
// storage-backed refs can participate in the same math.
//
// The header is organized like kokkos_vector_ops.h:
//
// - construction/materialization helpers (tensor_identity, copy_tensor),
// - detail:: componentwise kernels shared by everything below,
// - libMesh-like wrappers (contract, transpose, det, inverse, ...), and
// - the out-of-line tensor_ref definitions and the operator overloads.
//
// As in the vector case, the operator overloads are constrained (via
// is_tensor_ref_v / is_vector_ref_v) to require at least one
// storage-backed ref operand, so expressions between owning types keep
// using the operators those types already define.
//
// Functions taking an optional ResultTensor/ResultVector template
// parameter pick a default owning result type - the input's semantic
// type for copy_tensor/transpose/inverse, a plain TypeTensor or
// TypeVector of the input's component type for outer_product/row/
// column - but a caller can request a different owning type
// explicitly, e.g. inverse<TensorValue<Real>>(T).

#include "libmesh/kokkos_linalg_base.h"
#include "libmesh/kokkos_vector_ops.h"

#include "libmesh/tensor_tools.h"

#include <cmath>

namespace libMesh::Kokkos
{

// Construction and materialization

/**
 * \returns An owning \p ResultTensor whose leading \p dim x \p dim
 * block is the identity and whose remaining components are zero.
 */
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

/**
 * \returns An owning copy of the tensor-like \p T_in: by default a
 * value of \p T_in 's semantic type, or of \p ResultTensor if that
 * parameter is given explicitly.  This is how a kernel snapshots a
 * storage-backed ref's current value for register-resident work.
 */
template <typename ResultTensor = void, typename TensorLike>
LIBMESH_DEVICE_INLINE
auto copy_tensor(const TensorLike & T_in)
  -> std::conditional_t<std::is_void<ResultTensor>::value,
                        tensor_semantic_type_t<TensorLike>,
                        ResultTensor>
{
  using output_type = std::conditional_t<std::is_void<ResultTensor>::value,
                                         tensor_semantic_type_t<TensorLike>,
                                         ResultTensor>;
  return materialize_tensor<output_type>(T_in);
}

namespace detail
{

// These helpers are shared by the public functions and ref operators so
// Kokkos-backed refs use direct component access without extra materialization.

/**
 * \returns The determinant of the leading \p dim x \p dim block of
 * \p T_in, expanded in closed form for dim <= 3.  A 0x0 determinant is
 * 1 by convention.
 */
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

/**
 * \returns The outer product left * conj(right)^T as an owning
 * \p ResultTensor, conjugating the right operand as the host-side
 * outer_product() free function in type_tensor.h does.
 */
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

/**
 * \returns The inverse of the leading \p dim x \p dim block of \p T_in
 * as an owning \p ResultTensor, computed from the closed-form adjugate
 * for dim <= 3.  Components outside the leading block are zero.
 */
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

/**
 * \returns The transpose of \p T_in as an owning \p ResultTensor.
 */
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

/**
 * \returns The matrix product left * right as an owning
 * \p ResultTensor.
 */
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

/**
 * \returns Row \p row_index of \p T_in as an owning \p ResultVector.
 */
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

/**
 * \returns Column \p col_index of \p T_in as an owning
 * \p ResultVector.
 */
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

/**
 * \returns The matrix-vector product T_in * v as an owning
 * \p ResultVector.
 */
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

/**
 * \returns The vector-matrix product v * T_in (i.e. T_in^T * v) as an
 * owning \p ResultVector.
 */
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

/**
 * Componentwise copy: left = right.
 */
template <typename LeftTensor, typename RightTensor>
LIBMESH_DEVICE_INLINE
void assign_tensor_components(LeftTensor & left, const RightTensor & right)
{
  for (unsigned int row = 0; row < LIBMESH_DIM; ++row)
    for (unsigned int col = 0; col < LIBMESH_DIM; ++col)
      left(row, col) = right(row, col);
}

/**
 * Sets every component of \p T_in to \p value.
 */
template <typename TensorLike, typename Scalar>
LIBMESH_DEVICE_INLINE
void fill_tensor_components(TensorLike & T_in, const Scalar & value)
{
  for (unsigned int row = 0; row < LIBMESH_DIM; ++row)
    for (unsigned int col = 0; col < LIBMESH_DIM; ++col)
      T_in(row, col) = value;
}

/**
 * Componentwise axpy-style update: left += factor * right.
 */
template <typename LeftTensor, typename RightTensor, typename Scalar>
LIBMESH_DEVICE_INLINE
void update_tensor_components(LeftTensor & left, const RightTensor & right, const Scalar & factor)
{
  for (unsigned int row = 0; row < LIBMESH_DIM; ++row)
    for (unsigned int col = 0; col < LIBMESH_DIM; ++col)
      left(row, col) = left(row, col) + factor * right(row, col);
}

/**
 * Componentwise transform: out(i,j) = op(in(i,j)).  \p out and \p in
 * may alias, as they do in the compound assignment operators below.
 */
template <typename OutputTensor, typename InputTensor, typename TransformOp>
LIBMESH_DEVICE_INLINE
void transform_tensor_components(OutputTensor & out, const InputTensor & in, const TransformOp & op)
{
  for (unsigned int row = 0; row < LIBMESH_DIM; ++row)
    for (unsigned int col = 0; col < LIBMESH_DIM; ++col)
      out(row, col) = op(in(row, col));
}

/**
 * \returns A new owning \p ResultTensor with components op(T_in(i,j)).
 */
template <typename ResultTensor, typename TensorLike, typename TransformOp>
LIBMESH_DEVICE_INLINE
ResultTensor transformed_tensor(const TensorLike & T_in, const TransformOp & op)
{
  ResultTensor out;
  out.zero();
  transform_tensor_components(out, T_in, op);
  return out;
}

/**
 * \returns \p true if \p left and \p right agree exactly in every
 * component, \p false otherwise.
 */
template <typename LeftTensor, typename RightTensor>
LIBMESH_DEVICE_INLINE
bool tensor_equal_impl(const LeftTensor & left, const RightTensor & right)
{
  for (unsigned int row = 0; row < LIBMESH_DIM; ++row)
    for (unsigned int col = 0; col < LIBMESH_DIM; ++col)
      if (left(row, col) != right(row, col))
        return false;

  return true;
}

// Tensor reductions

/**
 * \returns The full contraction (sum of componentwise products) of
 * \p left with \p right.
 */
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


/**
 * \returns The trace of \p T_in.
 */
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

/**
 * \returns The full contraction of \p left with \p right; the tensor
 * overload of the contract() name shared with kokkos_vector_ops.h.
 */
template <typename LeftTensor,
          typename RightTensor,
          typename std::enable_if<is_tensor_like_v<LeftTensor> && is_tensor_like_v<RightTensor>,
                                  int>::type = 0>
LIBMESH_DEVICE_INLINE
auto contract(const LeftTensor & left, const RightTensor & right)
{
  return detail::tensor_contract(left, right);
}


/**
 * \returns The outer product left * conj(right)^T of two vector-like
 * inputs, as an owning tensor: by default a TypeTensor of the left
 * operand's component type, or \p ResultTensor if given explicitly.
 */
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

/**
 * \returns The transpose of the tensor-like \p T_in, as an owning
 * tensor.
 */
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

/**
 * \returns The determinant of the tensor-like \p T_in, delegating to
 * its det() member so owning types use TypeTensor::det() while refs
 * use tensor_ref::det().
 */
template <typename TensorLike>
LIBMESH_DEVICE_INLINE
auto det(const TensorLike & T_in)
  -> std::enable_if_t<is_tensor_like_v<TensorLike>, decltype(T_in.det())>
{
  return T_in.det();
}

/**
 * \returns The inverse of the leading \p dim x \p dim block of the
 * tensor-like \p T_in, as an owning tensor.
 */
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

/**
 * \returns Row \p i of the tensor-like \p T_in, as an owning vector:
 * by default a TypeVector of \p T_in 's component type, or
 * \p ResultVector if given explicitly.
 */
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

/**
 * \returns Column \p i of the tensor-like \p T_in, as an owning
 * vector.
 */
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

// Forward declarations of the compound assignment operators, which the
// tensor_ref member definitions below delegate to.

template <typename LeftTensor, typename RightTensor>
LIBMESH_DEVICE_INLINE
auto operator+=(LeftTensor & left, const RightTensor & right)
  -> std::enable_if_t<is_tensor_like_v<LeftTensor> && is_tensor_like_v<RightTensor> &&
                        (is_tensor_ref_v<LeftTensor> || is_tensor_ref_v<RightTensor>),
                      LeftTensor &>;

template <typename LeftTensor, typename RightTensor>
LIBMESH_DEVICE_INLINE
auto operator-=(LeftTensor & left, const RightTensor & right)
  -> std::enable_if_t<is_tensor_like_v<LeftTensor> && is_tensor_like_v<RightTensor> &&
                        (is_tensor_ref_v<LeftTensor> || is_tensor_ref_v<RightTensor>),
                      LeftTensor &>;

// Out-of-line tensor_ref member definitions.  Several delegate to the
// generic operations above or the operators below, so the whole group
// lives here rather than in kokkos_linalg_base.h; see the class
// definition for their documentation.

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
  libMesh::Kokkos::operator+=(*this, right);
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
  libMesh::Kokkos::operator-=(*this, right);
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
  return libMesh::Kokkos::contract(*this, right);
}

template <typename ViewType>
LIBMESH_DEVICE_INLINE
auto tensor_ref<ViewType>::norm() const
{
  using std::sqrt;
  return sqrt(this->norm_sq());
}

template <typename ViewType>
LIBMESH_DEVICE_INLINE
auto tensor_ref<ViewType>::norm_sq() const
{
  // TensorTools::norm_sq handles complex-valued components correctly.
  using norm_type = detail::remove_cvref_t<decltype(libMesh::TensorTools::norm_sq((*this)(0, 0)))>;

  norm_type sum = norm_type(0);
  for (unsigned int row = 0; row < LIBMESH_DIM; ++row)
    for (unsigned int col = 0; col < LIBMESH_DIM; ++col)
      sum += libMesh::TensorTools::norm_sq((*this)(row, col));

  return sum;
}

template <typename ViewType>
LIBMESH_DEVICE_INLINE
bool tensor_ref<ViewType>::is_zero() const
{
  for (unsigned int row = 0; row < LIBMESH_DIM; ++row)
    for (unsigned int col = 0; col < LIBMESH_DIM; ++col)
      if ((*this)(row, col) != value_type(0))
        return false;

  return true;
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
    x(component) = solution(component);
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
//
// Every overload requires at least one tensor_ref or vector_ref
// operand, so these never compete with the operators the owning types
// define for themselves.  Tensor- and vector-valued results are
// returned as owning values; compound assignments update the left
// operand in place, writing through to storage when it is a ref.

/**
 * \returns The negative of the tensor-like \p T_in, as an owning
 * tensor.
 */
template <typename TensorLike>
LIBMESH_DEVICE_INLINE
auto operator-(const TensorLike & T_in)
  -> std::enable_if_t<is_tensor_like_v<TensorLike> && is_tensor_ref_v<TensorLike>,
                      tensor_semantic_type_t<TensorLike>>
{
  return detail::transformed_tensor<tensor_semantic_type_t<TensorLike>>(
    T_in,
    detail::negate_value<typename TensorLike::value_type>{});
}

/**
 * \returns The sum of \p left and \p right, as an owning tensor.
 */
template <typename LeftTensor, typename RightTensor>
LIBMESH_DEVICE_INLINE
auto operator+(const LeftTensor & left, const RightTensor & right)
  -> std::enable_if_t<is_tensor_like_v<LeftTensor> && is_tensor_like_v<RightTensor> &&
                        (is_tensor_ref_v<LeftTensor> || is_tensor_ref_v<RightTensor>),
                      tensor_semantic_type_t<LeftTensor>>
{
  auto out = copy_tensor<tensor_semantic_type_t<LeftTensor>>(left);
  out += right;
  return out;
}

/**
 * \returns The difference of \p left and \p right, as an owning
 * tensor.
 */
template <typename LeftTensor, typename RightTensor>
LIBMESH_DEVICE_INLINE
auto operator-(const LeftTensor & left, const RightTensor & right)
  -> std::enable_if_t<is_tensor_like_v<LeftTensor> && is_tensor_like_v<RightTensor> &&
                        (is_tensor_ref_v<LeftTensor> || is_tensor_ref_v<RightTensor>),
                      tensor_semantic_type_t<LeftTensor>>
{
  auto out = copy_tensor<tensor_semantic_type_t<LeftTensor>>(left);
  out -= right;
  return out;
}

/**
 * \returns The tensor-like \p T_in scaled by \p alpha, as an owning
 * tensor.  The scalar may appear on either side of the product.
 */
template <typename Scalar,
          typename TensorLike,
          typename std::enable_if<!is_vector_like_v<Scalar> && !is_tensor_like_v<Scalar> &&
                                    is_tensor_like_v<TensorLike> && is_tensor_ref_v<TensorLike>,
                                  int>::type = 0>
LIBMESH_DEVICE_INLINE
auto operator*(const Scalar & alpha, const TensorLike & T_in)
{
  return T_in * alpha;
}

template <typename TensorLike,
          typename Scalar,
          typename std::enable_if<is_tensor_like_v<TensorLike> && is_tensor_ref_v<TensorLike> &&
                                    !is_vector_like_v<Scalar> && !is_tensor_like_v<Scalar>,
                                  int>::type = 0>
LIBMESH_DEVICE_INLINE
auto operator*(const TensorLike & T_in, const Scalar & alpha)
{
  return detail::transformed_tensor<tensor_semantic_type_t<TensorLike>>(
    T_in,
    detail::scale_value<Scalar>{alpha});
}

/**
 * \returns The tensor-like \p T_in divided (componentwise) by the
 * scalar \p alpha, as an owning tensor.
 */
template <typename TensorLike, typename Scalar>
LIBMESH_DEVICE_INLINE
auto operator/(const TensorLike & T_in, const Scalar & alpha)
  -> std::enable_if_t<is_tensor_like_v<TensorLike> && is_tensor_ref_v<TensorLike> &&
                        !is_vector_like_v<Scalar> && !is_tensor_like_v<Scalar>,
                      tensor_semantic_type_t<TensorLike>>
{
  return detail::transformed_tensor<tensor_semantic_type_t<TensorLike>>(
    T_in,
    detail::divide_value<Scalar>{alpha});
}

/**
 * \returns The matrix product of \p left and \p right, as an owning
 * tensor.
 */
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

/**
 * \returns The matrix-vector product T_in * v, as an owning vector.
 */
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

/**
 * \returns The vector-matrix product v * T_in (i.e. T_in^T * v), as an
 * owning vector.
 */
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

/**
 * \returns \p true if \p left and \p right agree in every component,
 * \p false otherwise.  Each component is compared with exact
 * operator==.  \note This is stricter than TypeTensor::operator==,
 * which compares the summed componentwise differences against a
 * component-count multiple of TOLERANCE; two tensors that the host
 * operator calls equal may compare unequal here.
 */
template <typename LeftTensor, typename RightTensor>
LIBMESH_DEVICE_INLINE
auto operator==(const LeftTensor & left, const RightTensor & right)
  -> std::enable_if_t<is_tensor_like_v<LeftTensor> && is_tensor_like_v<RightTensor> &&
                        (is_tensor_ref_v<LeftTensor> || is_tensor_ref_v<RightTensor>),
                      bool>
{
  return detail::tensor_equal_impl(left, right);
}

/**
 * \returns \p true if \p left and \p right differ in any component,
 * \p false otherwise.
 */
template <typename LeftTensor, typename RightTensor>
LIBMESH_DEVICE_INLINE
auto operator!=(const LeftTensor & left, const RightTensor & right)
  -> std::enable_if_t<is_tensor_like_v<LeftTensor> && is_tensor_like_v<RightTensor> &&
                        (is_tensor_ref_v<LeftTensor> || is_tensor_ref_v<RightTensor>),
                      bool>
{
  return !(left == right);
}

/**
 * Adds \p right to \p left in place; when \p left is a ref, the update
 * writes through to the underlying storage.
 */
template <typename LeftTensor, typename RightTensor>
LIBMESH_DEVICE_INLINE
auto operator+=(LeftTensor & left, const RightTensor & right)
  -> std::enable_if_t<is_tensor_like_v<LeftTensor> && is_tensor_like_v<RightTensor> &&
                        (is_tensor_ref_v<LeftTensor> || is_tensor_ref_v<RightTensor>),
                      LeftTensor &>
{
  detail::update_tensor_components(left, right, 1);
  return left;
}

/**
 * Subtracts \p right from \p left in place.
 */
template <typename LeftTensor, typename RightTensor>
LIBMESH_DEVICE_INLINE
auto operator-=(LeftTensor & left, const RightTensor & right)
  -> std::enable_if_t<is_tensor_like_v<LeftTensor> && is_tensor_like_v<RightTensor> &&
                        (is_tensor_ref_v<LeftTensor> || is_tensor_ref_v<RightTensor>),
                      LeftTensor &>
{
  detail::update_tensor_components(left, right, -1);
  return left;
}

/**
 * Scales \p left by \p alpha in place.
 */
template <typename LeftTensor, typename Scalar>
LIBMESH_DEVICE_INLINE
auto operator*=(LeftTensor & left, const Scalar & alpha)
  -> std::enable_if_t<is_tensor_like_v<LeftTensor> && is_tensor_ref_v<LeftTensor> &&
                        !is_vector_like_v<Scalar> && !is_tensor_like_v<Scalar>,
                      LeftTensor &>
{
  detail::transform_tensor_components(left, left, detail::scale_value<Scalar>{alpha});
  return left;
}

/**
 * Divides \p left by \p alpha in place.
 */
template <typename LeftTensor, typename Scalar>
LIBMESH_DEVICE_INLINE
auto operator/=(LeftTensor & left, const Scalar & alpha)
  -> std::enable_if_t<is_tensor_like_v<LeftTensor> && is_tensor_ref_v<LeftTensor> &&
                        !is_vector_like_v<Scalar> && !is_tensor_like_v<Scalar>,
                      LeftTensor &>
{
  detail::transform_tensor_components(left, left, detail::divide_value<Scalar>{alpha});
  return left;
}

} // namespace libMesh::Kokkos

#endif // LIBMESH_KOKKOS_TENSOR_OPS_H
