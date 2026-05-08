#ifndef KOKKOS_TENSOR_OPS_ORACLE_FIXTURES_H
#define KOKKOS_TENSOR_OPS_ORACLE_FIXTURES_H

#include "libmesh/libmesh.h"
#include "libmesh/point.h"
#include "libmesh/tensor_value.h"
#include "libmesh/type_n_tensor.h"
#include "libmesh/vector_value.h"
#include "gpu/kokkos_tensor_ops.h"
#include "gpu/kokkos_storage.h"
#include "gpu/kokkos_storage_policy.h"

#include "kokkos_numerics_oracle_test_utils.h"

#include <type_traits>
#include <vector>

namespace libMeshTest
{
namespace KokkosTensorOracle
{

using libMesh::Real;

static constexpr double tol = 2.0e-13;

using oracle_vector = libMesh::TypeVector<Real>;
using oracle_tensor = libMesh::TypeTensor<Real>;

inline oracle_vector
make_host_vector(const Real x, const Real y = 0, const Real z = 0)
{
  oracle_vector v;
  v.zero();
  v(0) = x;
#if LIBMESH_DIM > 1
  v(1) = y;
#endif
#if LIBMESH_DIM > 2
  v(2) = z;
#endif
  return v;
}

inline oracle_tensor
make_host_tensor(const Real xx,
                 const Real xy = 0,
                 const Real xz = 0,
                 const Real yx = 0,
                 const Real yy = 0,
                 const Real yz = 0,
                 const Real zx = 0,
                 const Real zy = 0,
                 const Real zz = 0)
{
  oracle_tensor T;
  T.zero();
  T(0, 0) = xx;
#if LIBMESH_DIM > 1
  T(0, 1) = xy;
  T(1, 0) = yx;
  T(1, 1) = yy;
#endif
#if LIBMESH_DIM > 2
  T(0, 2) = xz;
  T(1, 2) = yz;
  T(2, 0) = zx;
  T(2, 1) = zy;
  T(2, 2) = zz;
#endif
  return T;
}

struct tensor_dim_case
{
  oracle_tensor J;
  unsigned int dim;
  const char * name;
};

static const tensor_dim_case dim_cases[] = {
  { make_host_tensor(1.7, -0.2, 0.5,
                     0.3,  1.1, -0.4,
                     -0.6, 0.8, 0.9),
    1,
    "leading_1d" },
#if LIBMESH_DIM > 1
  { make_host_tensor(2.5, -0.75, 0.4,
                     1.2,  1.8, -0.6,
                     -0.3, 0.9, 1.4),
    2,
    "leading_2d" },
#endif
#if LIBMESH_DIM > 2
  { make_host_tensor(9.08973348886179e-01, 3.36455579239923e-01, 5.16389236893863e-01,
                     9.44156071777472e-01, 1.35610910092516e-01, 1.49881119060538e-02,
                     1.15988384086146e-01, 6.79845197685518e-03, 3.77028969454745e-01),
    3,
    "leading_3d" }
#endif
};

inline oracle_tensor
build_identity_tensor(const unsigned int dim)
{
  oracle_tensor I;
  I.zero();
  for (unsigned int i = 0; i < dim; ++i)
    I(i, i) = Real(1);
  return I;
}

inline Real
host_leading_determinant(const oracle_tensor & J, const unsigned int dim)
{
  if (dim == 0)
    return Real(1);
  if (dim == 1)
    return J(0, 0);
  if (dim == 2)
    return J(0, 0) * J(1, 1) - J(0, 1) * J(1, 0);
#if LIBMESH_DIM > 2
  return J.det();
#else
  return Real(0);
#endif
}

inline oracle_tensor
host_leading_inverse(const oracle_tensor & J, const unsigned int dim)
{
  oracle_tensor inv;
  inv.zero();

  if (dim == 1)
  {
    inv(0, 0) = Real(1) / J(0, 0);
    return inv;
  }

  if (dim == 2)
  {
    const Real det = host_leading_determinant(J, dim);
    inv(0, 0) =  J(1, 1) / det;
    inv(0, 1) = -J(0, 1) / det;
    inv(1, 0) = -J(1, 0) / det;
    inv(1, 1) =  J(0, 0) / det;
    return inv;
  }

#if LIBMESH_DIM > 2
  return oracle_tensor(J.inverse());
#else
  return inv;
#endif
}

} // namespace KokkosTensorOracle
} // namespace libMeshTest

#endif
