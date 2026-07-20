// libMesh Kokkos compile-time linalg foundation.
//
// This header defines the small access/materialization layer that sits
// underneath richer vector/tensor algebra. It is intentionally limited to
// component access, storage-backed references, and conversion between
// vector-like/tensor-like objects and libMesh semantic types.

#ifndef LIBMESH_KOKKOS_LINALG_BASE_H
#define LIBMESH_KOKKOS_LINALG_BASE_H

#include "libmesh/libmesh_common.h"
#include "libmesh/libmesh_device.h"
#include "libmesh/point.h"
#include "libmesh/tensor_value.h"
#include "libmesh/type_tensor.h"
#include "libmesh/type_vector.h"
#include "libmesh/vector_value.h"

#include <type_traits>
#include <utility>

namespace libMesh::Kokkos
{

namespace detail
{

template <typename T>
using remove_cvref_t =
  typename std::remove_cv<typename std::remove_reference<T>::type>::type;

template <typename T>
using remove_ref_t = typename std::remove_reference<T>::type;

template <typename ViewType>
using vector_view_value_t =
  remove_cvref_t<decltype(std::declval<ViewType &>()(0, 0))>;

template <typename ViewType>
using tensor_view_value_t =
  remove_cvref_t<decltype(std::declval<ViewType &>()(0, 0, 0))>;

} // namespace detail

template <typename T>
struct vector_traits;

template <typename T>
struct tensor_traits;

template <typename T>
struct is_vector_like : std::false_type
{
};

template <typename T>
struct is_tensor_like : std::false_type
{
};

template <typename T>
struct is_vector_ref : std::false_type
{
};

template <typename T>
struct is_tensor_ref : std::false_type
{
};

template <typename T>
inline constexpr bool is_vector_like_v = is_vector_like<detail::remove_cvref_t<T>>::value;

template <typename T>
inline constexpr bool is_tensor_like_v = is_tensor_like<detail::remove_cvref_t<T>>::value;

template <typename T>
inline constexpr bool is_vector_ref_v = is_vector_ref<detail::remove_cvref_t<T>>::value;

template <typename T>
inline constexpr bool is_tensor_ref_v = is_tensor_ref<detail::remove_cvref_t<T>>::value;

template <typename ViewType>
class vector_ref
{
public:
  using view_type = ViewType;
  using value_type = detail::vector_view_value_t<ViewType>;

  LIBMESH_DEVICE_INLINE
  vector_ref(ViewType view, const unsigned int index) : _view(view), _index(index) {}

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

  LIBMESH_DEVICE_INLINE
  unsigned int index() const
  {
    return _index;
  }

private:
  ViewType _view;
  unsigned int _index;
};

template <typename ViewType>
class tensor_ref
{
public:
  using view_type = ViewType;
  using value_type = detail::tensor_view_value_t<ViewType>;

  LIBMESH_DEVICE_INLINE
  tensor_ref(ViewType view, const unsigned int index) : _view(view), _index(index) {}

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

  LIBMESH_DEVICE_INLINE
  unsigned int index() const
  {
    return _index;
  }

private:
  ViewType _view;
  unsigned int _index;
};

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

template <typename T>
using vector_semantic_type_t = typename vector_traits<detail::remove_cvref_t<T>>::semantic_type;

template <typename T>
using tensor_semantic_type_t = typename tensor_traits<detail::remove_cvref_t<T>>::semantic_type;

template <typename ViewType>
LIBMESH_DEVICE_INLINE
vector_ref<typename detail::remove_ref_t<ViewType>>
make_vector_ref(ViewType && view, const unsigned int index)
{
  return vector_ref<typename detail::remove_ref_t<ViewType>>(std::forward<ViewType>(view), index);
}

template <typename ViewType>
LIBMESH_DEVICE_INLINE
tensor_ref<typename detail::remove_ref_t<ViewType>>
make_tensor_ref(ViewType && view, const unsigned int index)
{
  return tensor_ref<typename detail::remove_ref_t<ViewType>>(std::forward<ViewType>(view), index);
}

template <typename OutputVector, typename VectorLike>
LIBMESH_DEVICE_INLINE
OutputVector materialize_vector(const VectorLike & v)
{
  static_assert(is_vector_like<detail::remove_cvref_t<VectorLike>>::value,
                "materialize_vector() requires a vector-like input type");

  OutputVector out;
  out.zero();

  for (unsigned int component = 0; component < LIBMESH_DIM; ++component)
    out(component) = v(component);

  return out;
}

template <typename OutputTensor, typename TensorLike>
LIBMESH_DEVICE_INLINE
OutputTensor materialize_tensor(const TensorLike & T_in)
{
  static_assert(is_tensor_like<detail::remove_cvref_t<TensorLike>>::value,
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
