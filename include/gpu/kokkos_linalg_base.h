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

#ifndef LIBMESH_KOKKOS_LINALG_BASE_H
#define LIBMESH_KOKKOS_LINALG_BASE_H

// libMesh Kokkos compile-time linalg foundation.
//
// This header defines the small access/materialization layer that sits
// underneath richer vector/tensor algebra. It is intentionally limited to
// component access, storage-backed references, and conversion between
// vector-like/tensor-like objects and libMesh semantic types.
//
// The key abstractions are:
//
// - "vector-like"/"tensor-like": the types registered with the
//   is_vector_like/is_tensor_like traits, whose components are read
//   through operator()(i) / operator()(i,j).  libMesh's owning types
//   (TypeVector, TypeTensor, Point, ...) and the non-owning refs below
//   are registered, so the algorithms in kokkos_vector_ops.h and
//   kokkos_tensor_ops.h can be written once against either kind.
//
// - vector_ref/tensor_ref: lightweight non-owning handles addressing
//   one entry of a Kokkos storage view (see kokkos_storage_policy.h).
//   Reads and writes through a ref go directly to the underlying view,
//   so kernels can update storage in place without materializing
//   temporaries.
//
// - materialize_vector()/materialize_tensor(): copy any vector-like or
//   tensor-like object into an owning libMesh value type, for use in
//   register-resident arithmetic.

#include "libmesh/libmesh_common.h"
#include "libmesh/libmesh_device.h"
#include "libmesh/point.h"
#include "libmesh/tensor_value.h"
#include "libmesh/type_tensor.h"
#include "libmesh/type_vector.h"

#include <type_traits>
#include <utility>

namespace libMesh::Kokkos
{

namespace detail
{

/**
 * Strips reference and const/volatile qualifiers from \p T, turning
 * e.g. \c const \c Real& into plain \c Real.  Equivalent of C++20's
 * \c std::remove_cvref_t, hand-rolled while we only require C++17
 * support.
 */
template <typename T>
using remove_cvref_t =
  typename std::remove_cv<typename std::remove_reference<T>::type>::type;

/**
 * The plain component type stored by a vector view.  Indexing a view
 * at an (entry, component) pair yields a reference into the view's
 * memory (e.g. \c Real& , or \c const \c Real& for a view of const
 * data), so this strips the reference and any cv-qualifiers to recover the
 * value type itself - the type needed to declare accumulators, make
 * zeros like \c value_type(0), or materialize owning results.
 */
template <typename ViewType>
using vector_view_value_t =
  remove_cvref_t<decltype(std::declval<ViewType &>()(0, 0))>;

/**
 * The plain component type stored by a tensor view, recovered from an
 * (entry, row, column) indexing expression the same way
 * vector_view_value_t is from its (entry, component) expression.
 */
template <typename ViewType>
using tensor_view_value_t =
  remove_cvref_t<decltype(std::declval<ViewType &>()(0, 0, 0))>;

} // namespace detail

/**
 * Traits mapping each vector-like type to its "semantic" owning libMesh
 * type - the type an algorithm should return when it needs to
 * materialize a result from that input.  Specialized below for the
 * libMesh owning types and for vector_ref.
 */
template <typename T>
struct vector_traits;

/**
 * Traits mapping each tensor-like type to its "semantic" owning libMesh
 * type.  Specialized below for the libMesh owning types and for
 * tensor_ref.
 */
template <typename T>
struct tensor_traits;

/**
 * Opt-in trait marking the types the generic vector operations in
 * kokkos_vector_ops.h (and their SFINAE constraints) accept, all of
 * which expose components as \c v(i): specialized below for the
 * libMesh owning vector types and vector_ref.  There is no structural
 * detection - a type with a suitable operator() still needs an
 * explicit specialization before the trait-constrained operations
 * accept it.
 */
template <typename T>
struct is_vector_like : std::false_type
{
};

/**
 * Opt-in trait marking the types the generic tensor operations in
 * kokkos_tensor_ops.h accept, all of which expose components as
 * \c T(i,j): specialized below for the libMesh owning tensor types
 * and tensor_ref.
 */
template <typename T>
struct is_tensor_like : std::false_type
{
};

/**
 * Type trait satisfied only by vector_ref instantiations.  The operator
 * overloads in kokkos_vector_ops.h require at least one ref operand, so
 * they never compete with the operators the owning types define for
 * themselves.
 */
template <typename T>
struct is_vector_ref : std::false_type
{
};

/**
 * Type trait satisfied only by tensor_ref instantiations; the
 * tensor-ref analogue of is_vector_ref.
 */
template <typename T>
struct is_tensor_ref : std::false_type
{
};

/**
 * Variable template shorthands for the traits above.  The traits are
 * specialized only for plain (unqualified) types, so each shorthand
 * strips reference and cv-qualifiers before consulting the trait;
 * callers can therefore apply them to qualified types (e.g.
 * is_vector_like_v<const TypeVector<Real> &> is still true) without
 * decaying first.
 */
template <typename T>
inline constexpr bool is_vector_like_v = is_vector_like<detail::remove_cvref_t<T>>::value;

template <typename T>
inline constexpr bool is_tensor_like_v = is_tensor_like<detail::remove_cvref_t<T>>::value;

template <typename T>
inline constexpr bool is_vector_ref_v = is_vector_ref<detail::remove_cvref_t<T>>::value;

template <typename T>
inline constexpr bool is_tensor_ref_v = is_tensor_ref<detail::remove_cvref_t<T>>::value;

/**
 * Non-owning handle to one vector entry of a Kokkos storage view.
 *
 * A vector_ref pairs a view with an entry index and exposes the entry's
 * components through \c operator(), so it can be used wherever a
 * vector-like object is expected.  Component writes go directly to the
 * underlying view; the arithmetic methods below (defined out-of-line in
 * kokkos_vector_ops.h) either update the view in place or return owning
 * results, mirroring the TypeVector interface.
 *
 * Refs are cheap to construct and copy (a view handle plus an index)
 * and are intended to be created on the fly inside kernels via
 * make_vector_ref().
 */
template <typename ViewType>
class vector_ref
{
public:
  // The underlying view type and the scalar type of its components.
  using view_type = ViewType;
  using value_type = detail::vector_view_value_t<ViewType>;

  /**
   * Constructs a reference to entry \p index of \p view.
   */
  LIBMESH_DEVICE_INLINE
  vector_ref(ViewType view, const unsigned int index) : _view(view), _index(index) {}

  /**
   * \returns A reference to component \p component of the referenced
   * entry.  Writes through it land directly in the underlying view;
   * since Kokkos view constness is shallow, even the const overload
   * yields a writable reference when the view's data type is
   * non-const.
   */
  LIBMESH_DEVICE_INLINE
  decltype(auto) operator()(const unsigned int component) const
  {
    return _view(_index, component);
  }

  LIBMESH_DEVICE_INLINE
  decltype(auto) operator()(const unsigned int component)
  {
    return _view(_index, component);
  }

  /**
   * Assigns \p right 's components to the referenced entry.
   */
  template <typename RightVector>
  LIBMESH_DEVICE_INLINE
  void assign(const RightVector & right);

  /**
   * Adds \p right to the referenced entry, in place.
   */
  template <typename RightVector>
  LIBMESH_DEVICE_INLINE
  void add(const RightVector & right);

  /**
   * Adds \p factor times \p right to the referenced entry, in place.
   */
  template <typename RightVector>
  LIBMESH_DEVICE_INLINE
  void add_scaled(const RightVector & right, const value_type & factor);

  /**
   * Subtracts \p right from the referenced entry, in place.
   */
  template <typename RightVector>
  LIBMESH_DEVICE_INLINE
  void subtract(const RightVector & right);

  /**
   * Subtracts \p factor times \p right from the referenced entry, in
   * place.
   */
  template <typename RightVector>
  LIBMESH_DEVICE_INLINE
  void subtract_scaled(const RightVector & right, const value_type & factor);

  /**
   * Sets every component of the referenced entry to zero.
   */
  LIBMESH_DEVICE_INLINE
  void zero();

  /**
   * \returns The dot product of the referenced entry with \p right.
   */
  template <typename RightVector>
  LIBMESH_DEVICE_INLINE
  auto contract(const RightVector & right) const;

  /**
   * \returns The magnitude of the referenced entry, i.e. the
   * square-root of the sum of the elements squared.
   */
  LIBMESH_DEVICE_INLINE
  auto norm() const;

  /**
   * \returns The magnitude squared of the referenced entry.
   */
  LIBMESH_DEVICE_INLINE
  auto norm_sq() const;

  /**
   * \returns The l1-norm of the referenced entry, i.e. the sum of the
   * absolute values of the components.
   */
  LIBMESH_DEVICE_INLINE
  auto l1_norm() const;

  /**
   * \returns \p true if every component of the referenced entry is
   * exactly zero, \p false otherwise.
   */
  LIBMESH_DEVICE_INLINE
  bool is_zero() const;

  /**
   * \returns An owning unit vector in the direction of the referenced
   * entry.
   */
  LIBMESH_DEVICE_INLINE
  auto unit() const;

  /**
   * \returns The cross product of the referenced entry with \p right,
   * as an owning vector.
   */
  template <typename RightVector>
  LIBMESH_DEVICE_INLINE
  auto cross(const RightVector & right) const;

  /**
   * \returns The index of the entry this ref addresses within its view.
   */
  LIBMESH_DEVICE_INLINE
  unsigned int index() const
  {
    return _index;
  }

private:
  ViewType _view;
  unsigned int _index;
};

/**
 * Non-owning handle to one tensor entry of a Kokkos storage view.
 *
 * The tensor analogue of vector_ref: it pairs a view with an entry
 * index, exposes the entry's components through \c operator()(row,col),
 * and provides the TypeTensor-like arithmetic methods below (defined
 * out-of-line in kokkos_tensor_ops.h).  Component writes go directly to
 * the underlying view.
 */
template <typename ViewType>
class tensor_ref
{
public:
  // The underlying view type and the scalar type of its components.
  using view_type = ViewType;
  using value_type = detail::tensor_view_value_t<ViewType>;

  /**
   * Constructs a reference to entry \p index of \p view.
   */
  LIBMESH_DEVICE_INLINE
  tensor_ref(ViewType view, const unsigned int index) : _view(view), _index(index) {}

  /**
   * \returns A reference to the (\p row, \p col) component of the
   * referenced entry, writable under the same shallow-constness rules
   * as vector_ref::operator().
   */
  LIBMESH_DEVICE_INLINE
  decltype(auto) operator()(const unsigned int row, const unsigned int col) const
  {
    return _view(_index, row, col);
  }

  LIBMESH_DEVICE_INLINE
  decltype(auto) operator()(const unsigned int row, const unsigned int col)
  {
    return _view(_index, row, col);
  }

  /**
   * Assigns \p right 's components to the referenced entry.
   */
  template <typename RightTensor>
  LIBMESH_DEVICE_INLINE
  void assign(const RightTensor & right);

  /**
   * Adds \p right to the referenced entry, in place.
   */
  template <typename RightTensor>
  LIBMESH_DEVICE_INLINE
  void add(const RightTensor & right);

  /**
   * Adds \p factor times \p right to the referenced entry, in place.
   */
  template <typename RightTensor>
  LIBMESH_DEVICE_INLINE
  void add_scaled(const RightTensor & right, const value_type & factor);

  /**
   * Subtracts \p right from the referenced entry, in place.
   */
  template <typename RightTensor>
  LIBMESH_DEVICE_INLINE
  void subtract(const RightTensor & right);

  /**
   * Subtracts \p factor times \p right from the referenced entry, in
   * place.
   */
  template <typename RightTensor>
  LIBMESH_DEVICE_INLINE
  void subtract_scaled(const RightTensor & right, const value_type & factor);

  /**
   * Sets every component of the referenced entry to zero.
   */
  LIBMESH_DEVICE_INLINE
  void zero();

  /**
   * \returns The full contraction (sum of componentwise products) of
   * the referenced entry with \p right.
   */
  template <typename RightTensor>
  LIBMESH_DEVICE_INLINE
  auto contract(const RightTensor & right) const;

  /**
   * \returns The Frobenius norm of the referenced entry, i.e. the
   * square-root of the sum of the elements squared.
   */
  LIBMESH_DEVICE_INLINE
  auto norm() const;

  /**
   * \returns The Frobenius norm squared of the referenced entry.
   */
  LIBMESH_DEVICE_INLINE
  auto norm_sq() const;

  /**
   * \returns \p true if every component of the referenced entry is
   * exactly zero, \p false otherwise.
   */
  LIBMESH_DEVICE_INLINE
  bool is_zero() const;

  /**
   * \returns The transpose of the referenced entry, as an owning
   * tensor.
   */
  LIBMESH_DEVICE_INLINE
  auto transpose() const;

  /**
   * \returns The determinant of the leading \p dim x \p dim block of
   * the referenced entry.
   */
  LIBMESH_DEVICE_INLINE
  auto det(const unsigned int dim = LIBMESH_DIM) const;

  /**
   * \returns The trace of the referenced entry.
   */
  LIBMESH_DEVICE_INLINE
  auto tr() const;

  /**
   * \returns The inverse of the leading \p dim x \p dim block of the
   * referenced entry, as an owning tensor.
   */
  LIBMESH_DEVICE_INLINE
  auto inverse(const unsigned int dim = LIBMESH_DIM) const;

  /**
   * Solves the linear system (*this)x = \p b for \p x, via the
   * explicit inverse.
   */
  template <typename VectorLike, typename ResultVector>
  LIBMESH_DEVICE_INLINE
  void solve(const VectorLike & b, ResultVector & x) const;

  /**
   * \returns Row \p i of the referenced entry, as an owning vector.
   */
  LIBMESH_DEVICE_INLINE
  auto row(const unsigned int i) const;

  /**
   * \returns Column \p i of the referenced entry, as an owning vector.
   */
  LIBMESH_DEVICE_INLINE
  auto column(const unsigned int i) const;

  /**
   * \returns The vector-tensor product \p v * (*this), as an owning
   * vector.
   */
  template <typename VectorLike>
  LIBMESH_DEVICE_INLINE
  auto left_multiply(const VectorLike & v) const;

  /**
   * \returns The index of the entry this ref addresses within its view.
   */
  LIBMESH_DEVICE_INLINE
  unsigned int index() const
  {
    return _index;
  }

private:
  ViewType _view;
  unsigned int _index;
};

// Owning vector types are their own semantic type; a vector_ref
// materializes to the TypeVector of its component type.

template <typename T>
struct vector_traits<libMesh::TypeVector<T>>
{
  using semantic_type = libMesh::TypeVector<T>;
};

template <typename T>
struct vector_traits<libMesh::VectorValue<T>>
{
  using semantic_type = libMesh::VectorValue<T>;
};

template <>
struct vector_traits<libMesh::Point>
{
  using semantic_type = libMesh::Point;
};

template <typename ViewType>
struct vector_traits<vector_ref<ViewType>>
{
  using value_type = typename vector_ref<ViewType>::value_type;
  using semantic_type = libMesh::TypeVector<value_type>;
};

// The vector-like family: libMesh owning vector types plus vector_ref.

template <typename T>
struct is_vector_like<libMesh::TypeVector<T>> : std::true_type
{
};

template <typename T>
struct is_vector_like<libMesh::VectorValue<T>> : std::true_type
{
};

template <>
struct is_vector_like<libMesh::Point> : std::true_type
{
};

template <typename ViewType>
struct is_vector_like<vector_ref<ViewType>> : std::true_type
{
};

template <typename ViewType>
struct is_vector_ref<vector_ref<ViewType>> : std::true_type
{
};

// Owning tensor types are their own semantic type; a tensor_ref
// materializes to the TypeTensor of its component type.

template <typename T>
struct tensor_traits<libMesh::TypeTensor<T>>
{
  using semantic_type = libMesh::TypeTensor<T>;
};

template <typename T>
struct tensor_traits<libMesh::TensorValue<T>>
{
  using semantic_type = libMesh::TensorValue<T>;
};

template <typename ViewType>
struct tensor_traits<tensor_ref<ViewType>>
{
  using value_type = typename tensor_ref<ViewType>::value_type;
  using semantic_type = libMesh::TypeTensor<value_type>;
};

// The tensor-like family: libMesh owning tensor types plus tensor_ref.

template <typename T>
struct is_tensor_like<libMesh::TypeTensor<T>> : std::true_type
{
};

template <typename T>
struct is_tensor_like<libMesh::TensorValue<T>> : std::true_type
{
};

template <typename ViewType>
struct is_tensor_like<tensor_ref<ViewType>> : std::true_type
{
};

template <typename ViewType>
struct is_tensor_ref<tensor_ref<ViewType>> : std::true_type
{
};

/**
 * The owning libMesh type a vector-like \p T materializes to, with any
 * cv/reference qualification on \p T stripped first.
 */
template <typename T>
using vector_semantic_type_t = typename vector_traits<detail::remove_cvref_t<T>>::semantic_type;

/**
 * The owning libMesh type a tensor-like \p T materializes to, with any
 * cv/reference qualification on \p T stripped first.
 */
template <typename T>
using tensor_semantic_type_t = typename tensor_traits<detail::remove_cvref_t<T>>::semantic_type;

/**
 * \returns A vector_ref addressing entry \p index of \p view, so
 * kernel code need not spell out the vector_ref's template argument.
 * The ref stores the view by value (Kokkos views are cheap
 * reference-counted handles), so the deduced reference qualification
 * on the forwarded \p view is stripped from the ref's type.
 */
template <typename ViewType>
LIBMESH_DEVICE_INLINE
vector_ref<typename std::remove_reference<ViewType>::type>
make_vector_ref(ViewType && view, const unsigned int index)
{
  return vector_ref<typename std::remove_reference<ViewType>::type>(std::forward<ViewType>(view),
                                                                    index);
}

/**
 * \returns A tensor_ref addressing entry \p index of \p view; the
 * tensor analogue of make_vector_ref(), with the same by-value view
 * semantics.
 */
template <typename ViewType>
LIBMESH_DEVICE_INLINE
tensor_ref<typename std::remove_reference<ViewType>::type>
make_tensor_ref(ViewType && view, const unsigned int index)
{
  return tensor_ref<typename std::remove_reference<ViewType>::type>(std::forward<ViewType>(view),
                                                                    index);
}

/**
 * \returns An owning vector of type \p OutputVector whose components
 * are copied from the vector-like \p v.
 *
 * This is the bridge from non-owning refs (or any other vector-like
 * input) to register-resident values a kernel can compute with.
 */
template <typename OutputVector, typename VectorLike>
LIBMESH_DEVICE_INLINE
OutputVector materialize_vector(const VectorLike & v)
{
  static_assert(is_vector_like_v<OutputVector>,
                "materialize_vector() requires a vector-like output type");
  static_assert(is_vector_like_v<VectorLike>,
                "materialize_vector() requires a vector-like input type");

  OutputVector out;
  out.zero();

  for (unsigned int component = 0; component < LIBMESH_DIM; ++component)
    out(component) = v(component);

  return out;
}

/**
 * \returns An owning tensor of type \p OutputTensor whose components
 * are copied from the tensor-like \p T_in.
 */
template <typename OutputTensor, typename TensorLike>
LIBMESH_DEVICE_INLINE
OutputTensor materialize_tensor(const TensorLike & T_in)
{
  static_assert(is_tensor_like_v<OutputTensor>,
                "materialize_tensor() requires a tensor-like output type");
  static_assert(is_tensor_like_v<TensorLike>,
                "materialize_tensor() requires a tensor-like input type");

  OutputTensor out;
  out.zero();

  for (unsigned int row = 0; row < LIBMESH_DIM; ++row)
    for (unsigned int col = 0; col < LIBMESH_DIM; ++col)
      out(row, col) = T_in(row, col);

  return out;
}

} // namespace libMesh::Kokkos

#endif // LIBMESH_KOKKOS_LINALG_BASE_H
