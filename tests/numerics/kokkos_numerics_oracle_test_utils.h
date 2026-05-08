#ifndef KOKKOS_NUMERICS_ORACLE_TEST_UTILS_H
#define KOKKOS_NUMERICS_ORACLE_TEST_UTILS_H

#include "libmesh/libmesh.h"

// Avoid conflicting complex operators between CUDA and PETSc
#define PETSC_SKIP_CXX_COMPLEX_FIX 1
#include <Kokkos_Core.hpp>
#undef __CUDACC_VER__

#include <cmath>
#include <string>
#include <vector>

namespace libMeshTest
{
namespace KokkosOracle
{

using libMesh::Real;

template <typename T>
inline ::Kokkos::View<T *>
upload_objects(const std::vector<T> & values, const char * label)
{
  ::Kokkos::View<T *> d(std::string(label), values.size());
  auto h = ::Kokkos::create_mirror_view(d);
  for (std::size_t i = 0; i < values.size(); ++i)
    h(i) = values[i];
  ::Kokkos::deep_copy(d, h);
  return d;
}

inline int
compare_device_scalars(const ::Kokkos::View<double *> & d_values,
                       const std::vector<Real> & ref_values,
                       const double tol)
{
  auto h_values = ::Kokkos::create_mirror_view(d_values);
  ::Kokkos::deep_copy(h_values, d_values);

  int fail = 0;
  for (std::size_t i = 0; i < ref_values.size(); ++i)
    if (std::fabs(h_values(i) - ref_values[i]) > tol)
      ++fail;

  return fail;
}

template <typename ViewType, typename VectorType>
inline int
compare_device_vectors(const ViewType & d_values,
                       const std::vector<VectorType> & ref_values,
                       const double tol)
{
  auto h_values = ::Kokkos::create_mirror_view(d_values);
  ::Kokkos::deep_copy(h_values, d_values);

  int fail = 0;
  for (std::size_t i = 0; i < ref_values.size(); ++i)
    for (unsigned int d = 0; d < LIBMESH_DIM; ++d)
      if (std::fabs(h_values(i, d) - ref_values[i](d)) > tol)
        ++fail;

  return fail;
}

template <typename ViewType, typename TensorType>
inline int
compare_device_tensors(const ViewType & d_values,
                       const std::vector<TensorType> & ref_values,
                       const double tol)
{
  auto h_values = ::Kokkos::create_mirror_view(d_values);
  ::Kokkos::deep_copy(h_values, d_values);

  int fail = 0;
  for (std::size_t i = 0; i < ref_values.size(); ++i)
    for (unsigned int row = 0; row < LIBMESH_DIM; ++row)
      for (unsigned int col = 0; col < LIBMESH_DIM; ++col)
        if (std::fabs(h_values(i, row, col) - ref_values[i](row, col)) > tol)
          ++fail;

  return fail;
}

} // namespace KokkosOracle
} // namespace libMeshTest

#endif
