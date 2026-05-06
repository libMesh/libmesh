// libMesh Kokkos device-compatible scalar types.
//
// This header provides dimension-aware Kokkos aliases/helpers that mirror
// libMesh host numerics at LIBMESH_DIM=1/2/3.

#ifndef LIBMESH_KOKKOS_SCALAR_TYPES_H
#define LIBMESH_KOKKOS_SCALAR_TYPES_H

#include "libmesh/libmesh_common.h"
#include "libmesh/libmesh_device.h"
#include "libmesh/type_vector.h"
#include "libmesh/type_tensor.h"

namespace libMesh::Kokkos
{

using Real = libMesh::Real;
using RealVector = libMesh::TypeVector<Real>;
using RealTensor = libMesh::TypeTensor<Real>;

template <typename VectorType, typename ViewType>
LIBMESH_DEVICE_INLINE
VectorType load_vector(const ViewType & view, const unsigned int i)
{
  VectorType v;
  v.zero();

  for (unsigned int d = 0; d < LIBMESH_DIM; ++d)
    v(d) = view(i, d);

  return v;
}

template <typename ViewType, typename VectorType>
LIBMESH_DEVICE_INLINE
void store_vector(const ViewType & view, const unsigned int i, const VectorType & v)
{
  for (unsigned int d = 0; d < LIBMESH_DIM; ++d)
    view(i, d) = v(d);
}

template <typename ViewType>
LIBMESH_DEVICE_INLINE
Real vector_component(const ViewType & view, const unsigned int i, const unsigned int component)
{
  if (component < LIBMESH_DIM)
    return view(i, component);

  return Real(0);
}

template <typename TensorType, typename ViewType>
LIBMESH_DEVICE_INLINE
TensorType load_tensor(const ViewType & view, const unsigned int i)
{
  TensorType T;
  T.zero();

  for (unsigned int row = 0; row < LIBMESH_DIM; ++row)
    for (unsigned int col = 0; col < LIBMESH_DIM; ++col)
      T(row, col) = view(i, row, col);

  return T;
}

template <typename ViewType, typename TensorType>
LIBMESH_DEVICE_INLINE
void store_tensor(const ViewType & view, const unsigned int i, const TensorType & T)
{
  for (unsigned int row = 0; row < LIBMESH_DIM; ++row)
    for (unsigned int col = 0; col < LIBMESH_DIM; ++col)
      view(i, row, col) = T(row, col);
}

template <typename ViewType>
LIBMESH_DEVICE_INLINE
Real tensor_component(const ViewType & view,
                      const unsigned int i,
                      const unsigned int row,
                      const unsigned int col)
{
  if (row < LIBMESH_DIM && col < LIBMESH_DIM)
    return view(i, row, col);

  return Real(0);
}

LIBMESH_DEVICE_INLINE
RealVector zero_vector()
{
  RealVector v;
  v.zero();
  return v;
}

LIBMESH_DEVICE_INLINE
RealVector make_vector(const Real x, const Real y = 0, const Real z = 0)
{
  RealVector v = zero_vector();

  v(0) = x;

#if LIBMESH_DIM > 1
  v(1) = y;
#else
  libmesh_assert_equal_to(y, Real(0));
#endif

#if LIBMESH_DIM > 2
  v(2) = z;
#else
  libmesh_assert_equal_to(z, Real(0));
#endif

  return v;
}

LIBMESH_DEVICE_INLINE
RealTensor zero_tensor()
{
  RealTensor J;
  J.zero();
  return J;
}

LIBMESH_DEVICE_INLINE
RealTensor leading_identity(const unsigned int dim = LIBMESH_DIM)
{
  libmesh_assert_less_equal(dim, LIBMESH_DIM);

  RealTensor I = zero_tensor();
  for (unsigned int i = 0; i < dim; ++i)
    I(i, i) = Real(1);

  return I;
}

LIBMESH_DEVICE_INLINE
Real leading_determinant(const RealTensor & J, const unsigned int dim = LIBMESH_DIM)
{
  libmesh_assert_less_equal(dim, LIBMESH_DIM);

  if (dim == 0)
    return Real(1);

  if (dim == 1)
    return J(0, 0);

  if (dim == 2)
    return J(0, 0) * J(1, 1) - J(0, 1) * J(1, 0);

  return J.det();
}

LIBMESH_DEVICE_INLINE
RealTensor leading_inverse(const RealTensor & J, const unsigned int dim = LIBMESH_DIM)
{
  libmesh_assert_less_equal(dim, LIBMESH_DIM);

  if (dim == 0)
    return leading_identity(0);

  if (dim == 1)
  {
    RealTensor inv = zero_tensor();
    inv(0, 0) = Real(1) / J(0, 0);
    return inv;
  }

  if (dim == 2)
  {
    const Real inv_det = Real(1) / leading_determinant(J, dim);
    RealTensor inv = zero_tensor();
    inv(0, 0) =  J(1, 1) * inv_det;
    inv(0, 1) = -J(0, 1) * inv_det;
    inv(1, 0) = -J(1, 0) * inv_det;
    inv(1, 1) =  J(0, 0) * inv_det;
    return inv;
  }

  return J.inverse();
}

} // namespace libMesh::Kokkos

#endif // LIBMESH_KOKKOS_SCALAR_TYPES_H
