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

#ifndef LIBMESH_KOKKOS_STORAGE_POLICY_H
#define LIBMESH_KOKKOS_STORAGE_POLICY_H

// libMesh Kokkos compile-time storage policies for fixed-dimension linalg data.
//
// These policies keep storage selection separate from the linalg algorithms:
// kernels operate on refs/materialized values, while the backend policy chooses
// the underlying Kokkos view layout.  Code written against a policy's
// vector_view/tensor_view typedefs can be instantiated with different memory
// layouts (or, in the future, different scalar types) without changing the
// math it performs.

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

/**
 * Storage policy describing Kokkos views that hold runtime-many
 * vectors or tensors whose per-entry dimensions are fixed at
 * LIBMESH_DIM.
 *
 * \tparam Scalar The component type stored in the views.
 * \tparam Layout The Kokkos memory layout (e.g. \c Kokkos::LayoutLeft
 * or \c Kokkos::LayoutRight) used to order components in memory.
 */
template <typename Scalar, typename Layout>
struct static_dim_storage_policy
{
  using scalar_type = Scalar;
  using layout_type = Layout;

  /**
   * View type holding \c n vectors: the runtime extent indexes the
   * entry, the static LIBMESH_DIM extent indexes the component.
   */
  using vector_view = ::Kokkos::View<scalar_type * [LIBMESH_DIM], layout_type>;

  /**
   * View type holding \c n tensors: the runtime extent indexes the
   * entry, the static LIBMESH_DIM x LIBMESH_DIM extents index the
   * components.
   */
  using tensor_view = ::Kokkos::View<scalar_type * [LIBMESH_DIM][LIBMESH_DIM], layout_type>;

  /**
   * \returns A short human-readable name for the policy's layout,
   * suitable for labeling tests or output.
   */
  static constexpr const char *
  name()
  {
    return std::is_same<layout_type, ::Kokkos::LayoutLeft>::value ? "layoutleft" :
           std::is_same<layout_type, ::Kokkos::LayoutRight>::value ? "layoutright" :
           "layoutcustom";
  }
};

/**
 * Convenience policies for libMesh's default \c Real scalar type with
 * the two standard Kokkos layouts, plus the default policy used by the
 * policy-free convenience overloads below.
 *
 * LayoutRight keeps each entry's components contiguous in memory, the
 * natural arrangement when one thread works on one whole entry;
 * LayoutLeft keeps a given component contiguous across entries, the
 * coalescing-friendly arrangement for GPU thread teams.  The Kokkos
 * unit tests run under both policies to keep the algebra
 * layout-agnostic.
 */
using layout_left_storage_policy = static_dim_storage_policy<libMesh::Real, ::Kokkos::LayoutLeft>;
using layout_right_storage_policy = static_dim_storage_policy<libMesh::Real, ::Kokkos::LayoutRight>;
using default_storage_policy = layout_right_storage_policy;

/**
 * \returns The name of \p StoragePolicy, as reported by
 * its static \c name() member; a free-function form convenient for
 * labeling per-policy results in generic code such as the Kokkos unit
 * tests.
 */
template <typename StoragePolicy>
constexpr const char *
storage_policy_name()
{
  return StoragePolicy::name();
}

/**
 * \returns A newly allocated (zero-initialized) device view holding
 * \p n vectors, laid out according to \p StoragePolicy and labeled
 * \p label for Kokkos profiling/debugging purposes.
 */
template <typename StoragePolicy>
inline typename StoragePolicy::vector_view
make_vector_storage(const char * label, const std::size_t n)
{
  return typename StoragePolicy::vector_view(std::string(label), n);
}

/**
 * \returns A newly allocated device view holding \p n vectors, using
 * the default storage policy.
 */
inline default_storage_policy::vector_view
make_vector_storage(const char * label, const std::size_t n)
{
  return make_vector_storage<default_storage_policy>(label, n);
}

/**
 * \returns A newly allocated (zero-initialized) device view holding
 * \p n tensors, laid out according to \p StoragePolicy and labeled
 * \p label for Kokkos profiling/debugging purposes.
 */
template <typename StoragePolicy>
inline typename StoragePolicy::tensor_view
make_tensor_storage(const char * label, const std::size_t n)
{
  return typename StoragePolicy::tensor_view(std::string(label), n);
}

/**
 * \returns A newly allocated device view holding \p n tensors, using
 * the default storage policy.
 */
inline default_storage_policy::tensor_view
make_tensor_storage(const char * label, const std::size_t n)
{
  return make_tensor_storage<default_storage_policy>(label, n);
}

/**
 * Allocates a device view sized for \p values, copies each vector's
 * components into a host mirror, and deep-copies the result to the
 * device.
 *
 * \returns The populated device view.
 *
 * \note Host-side setup helper, not something to call from within a
 * kernel.  With a device-resident default memory space the deep_copy
 * is a host-to-device transfer; on host-only builds the mirror
 * aliases the view and the copy is a no-op.
 */
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

/**
 * \returns A device view populated from \p values, using the default
 * storage policy.
 */
template <typename VectorType>
inline default_storage_policy::vector_view
upload_vector_storage(const std::vector<VectorType> & values, const char * label)
{
  return upload_vector_storage<default_storage_policy>(values, label);
}

/**
 * Allocates a device view sized for \p values, copies each tensor's
 * components into a host mirror, and deep-copies the result to the
 * device.
 *
 * \returns The populated device view.
 *
 * \note Host-side setup helper, not something to call from within a
 * kernel.  With a device-resident default memory space the deep_copy
 * is a host-to-device transfer; on host-only builds the mirror
 * aliases the view and the copy is a no-op.
 */
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

/**
 * \returns A device view populated from \p values, using the default
 * storage policy.
 */
template <typename TensorType>
inline default_storage_policy::tensor_view
upload_tensor_storage(const std::vector<TensorType> & values, const char * label)
{
  return upload_tensor_storage<default_storage_policy>(values, label);
}

} // namespace libMesh::Kokkos

#endif // LIBMESH_KOKKOS_STORAGE_POLICY_H
