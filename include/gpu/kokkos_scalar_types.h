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
