// libMesh Kokkos storage helpers for dimension-aware vector/tensor views.

#ifndef LIBMESH_KOKKOS_STORAGE_H
#define LIBMESH_KOKKOS_STORAGE_H

#include "libmesh/kokkos_linalg_base.h"

#include "libmesh/libmesh_common.h"
#include "libmesh/libmesh_device.h"
#include "libmesh/type_tensor.h"
#include "libmesh/type_vector.h"

namespace libMesh::Kokkos
{

template <typename VectorType, typename ViewType>
LIBMESH_DEVICE_INLINE
VectorType load_vector(const ViewType & view, const unsigned int i)
{
  return materialize_vector<VectorType>(make_vector_ref(view, i));
}

template <typename ViewType, typename VectorType>
LIBMESH_DEVICE_INLINE
void store_vector(const ViewType & view, const unsigned int i, const VectorType & v)
{
  auto out = make_vector_ref(view, i);

  for (unsigned int d = 0; d < LIBMESH_DIM; ++d)
    out(d) = v(d);
}

template <typename TensorType, typename ViewType>
LIBMESH_DEVICE_INLINE
TensorType load_tensor(const ViewType & view, const unsigned int i)
{
  return materialize_tensor<TensorType>(make_tensor_ref(view, i));
}

template <typename ViewType, typename TensorType>
LIBMESH_DEVICE_INLINE
void store_tensor(const ViewType & view, const unsigned int i, const TensorType & T)
{
  auto out = make_tensor_ref(view, i);

  for (unsigned int row = 0; row < LIBMESH_DIM; ++row)
    for (unsigned int col = 0; col < LIBMESH_DIM; ++col)
      out(row, col) = T(row, col);
}

} // namespace libMesh::Kokkos

#endif // LIBMESH_KOKKOS_STORAGE_H
