// libMesh Kokkos compile-time storage policies for fixed-dimension linalg data.
//
// These policies keep storage selection separate from the linalg algorithms:
// kernels operate on refs/materialized values, while the backend policy chooses
// the underlying Kokkos view layout.

#ifndef LIBMESH_KOKKOS_STORAGE_POLICY_H
#define LIBMESH_KOKKOS_STORAGE_POLICY_H

#include "libmesh/libmesh_common.h"

#define PETSC_SKIP_CXX_COMPLEX_FIX 1
#include "libmesh/ignore_warnings.h"
#include <Kokkos_Core.hpp>
#include "libmesh/restore_warnings.h"
#undef __CUDACC_VER__

#include <string>
#include <type_traits>
#include <vector>

namespace libMesh::Kokkos
{

template <typename Scalar, typename Layout>
struct static_dim_storage_policy
{
  using scalar_type = Scalar;
  using layout_type = Layout;
  using vector_view = ::Kokkos::View<scalar_type * [LIBMESH_DIM], layout_type>;
  using tensor_view = ::Kokkos::View<scalar_type * [LIBMESH_DIM][LIBMESH_DIM], layout_type>;

  static constexpr const char *
  name()
  {
    return std::is_same<layout_type, ::Kokkos::LayoutLeft>::value ? "layoutleft" :
           std::is_same<layout_type, ::Kokkos::LayoutRight>::value ? "layoutright" :
           "layoutcustom";
  }
};

using layout_left_storage_policy = static_dim_storage_policy<libMesh::Real, ::Kokkos::LayoutLeft>;
using layout_right_storage_policy = static_dim_storage_policy<libMesh::Real, ::Kokkos::LayoutRight>;
using default_storage_policy = layout_right_storage_policy;

template <typename StoragePolicy>
constexpr const char *
storage_policy_name()
{
  return StoragePolicy::name();
}

template <typename StoragePolicy>
inline typename StoragePolicy::vector_view
make_vector_storage(const char * label, const std::size_t n)
{
  return typename StoragePolicy::vector_view(std::string(label), n);
}

inline default_storage_policy::vector_view
make_vector_storage(const char * label, const std::size_t n)
{
  return make_vector_storage<default_storage_policy>(label, n);
}

template <typename StoragePolicy>
inline typename StoragePolicy::tensor_view
make_tensor_storage(const char * label, const std::size_t n)
{
  return typename StoragePolicy::tensor_view(std::string(label), n);
}

inline default_storage_policy::tensor_view
make_tensor_storage(const char * label, const std::size_t n)
{
  return make_tensor_storage<default_storage_policy>(label, n);
}

template <typename StoragePolicy, typename VectorType>
inline typename StoragePolicy::vector_view
upload_vector_storage(const std::vector<VectorType> & values, const char * label)
{
  auto d = make_vector_storage<StoragePolicy>(label, values.size());
  auto h = ::Kokkos::create_mirror_view(d);

  for (std::size_t i = 0; i < values.size(); ++i)
    for (unsigned int component = 0; component < LIBMESH_DIM; ++component)
      h(i, component) = values[i](component);

  ::Kokkos::deep_copy(d, h);
  return d;
}

template <typename VectorType>
inline default_storage_policy::vector_view
upload_vector_storage(const std::vector<VectorType> & values, const char * label)
{
  return upload_vector_storage<default_storage_policy>(values, label);
}

template <typename StoragePolicy, typename TensorType>
inline typename StoragePolicy::tensor_view
upload_tensor_storage(const std::vector<TensorType> & values, const char * label)
{
  auto d = make_tensor_storage<StoragePolicy>(label, values.size());
  auto h = ::Kokkos::create_mirror_view(d);

  for (std::size_t i = 0; i < values.size(); ++i)
    for (unsigned int row = 0; row < LIBMESH_DIM; ++row)
      for (unsigned int col = 0; col < LIBMESH_DIM; ++col)
        h(i, row, col) = values[i](row, col);

  ::Kokkos::deep_copy(d, h);
  return d;
}

template <typename TensorType>
inline default_storage_policy::tensor_view
upload_tensor_storage(const std::vector<TensorType> & values, const char * label)
{
  return upload_tensor_storage<default_storage_policy>(values, label);
}

} // namespace libMesh::Kokkos

#endif // LIBMESH_KOKKOS_STORAGE_POLICY_H
