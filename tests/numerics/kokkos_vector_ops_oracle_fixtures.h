#ifndef KOKKOS_VECTOR_OPS_ORACLE_FIXTURES_H
#define KOKKOS_VECTOR_OPS_ORACLE_FIXTURES_H

#include "libmesh/libmesh.h"
#include "libmesh/tensor_value.h"
#include "libmesh/type_vector.h"
#include "libmesh/vector_value.h"
#include "gpu/kokkos_vector_ops.h"
#include "gpu/kokkos_storage.h"
#include "gpu/kokkos_storage_policy.h"

#include "kokkos_numerics_oracle_test_utils.h"

#include <type_traits>
#include <vector>

namespace libMeshTest
{
namespace KokkosVectorOracle
{

using libMesh::Real;

static constexpr double tol = 2.0e-13;
static constexpr double unit_tol = 1.0e-14;
static constexpr Real golden_ratio = 1.6180339887498948482;
static constexpr unsigned int solid_angle_results =
  1 + ((LIBMESH_DIM > 1) ? 2u : 0u) + ((LIBMESH_DIM > 2) ? 4u : 0u);
static constexpr unsigned int vector_results =
  11 + ((LIBMESH_DIM > 2) ? 2u : 0u);
static constexpr unsigned int scalar_results = 12 + solid_angle_results;

template <typename Vec>
LIBMESH_DEVICE_INLINE
Vec
make_vector(const Real x, const Real y = 0, const Real z = 0)
{
  Vec v;
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

inline libMesh::TypeVector<Real>
as_type_vector(const libMesh::TypeVector<Real> & v)
{
  return v;
}

inline libMesh::TypeVector<Real>
as_type_vector(const libMesh::VectorValue<Real> & v)
{
  return make_vector<libMesh::TypeVector<Real>>(v(0)
#if LIBMESH_DIM > 1
                                                ,
                                                v(1)
#endif
#if LIBMESH_DIM > 2
                                                ,
                                                v(2)
#endif
  );
}

template <typename Vec>
struct host_oracle
{
  std::vector<Vec> vectors;
  std::vector<Real> scalars;
};

struct vector_case
{
  const char * name;
  Real ax, ay, az;
  Real bx, by, bz;
  Real cx, cy, cz;
};

static const vector_case cases[] = {
#if LIBMESH_DIM >= 1
  { "line_case_a", 2.0, 0.0, 0.0, -3.0, 0.0, 0.0, 0.5, 0.0, 0.0 },
  { "line_case_b", -1.25, 0.0, 0.0, 4.5, 0.0, 0.0, -2.0, 0.0, 0.0 },
#endif
#if LIBMESH_DIM >= 2
  { "plane_case_a", 2.0, 3.0, 0.0, 5.0, -6.0, 0.0, 1.25, -0.5, 0.0 },
  { "plane_case_b", -1.0, 4.0, 0.0, 0.5, 2.5, 0.0, -3.0, 1.5, 0.0 },
#endif
#if LIBMESH_DIM >= 3
  { "space_case_a", 2.0, 3.0, 4.0, 5.0, -6.0, 7.0, 1.25, -0.5, 2.0 },
  { "space_case_b", -1.0, 4.0, 0.75, 0.5, 2.5, -3.5, -3.0, 1.5, 2.25 },
#endif
};

template <typename Vec>
inline host_oracle<Vec>
build_host_oracle(const Vec & a, const Vec & b, const Vec & c)
{
  host_oracle<Vec> result;
  result.vectors.reserve(vector_results);
  result.scalars.reserve(scalar_results);

  const auto copied = a;

  Vec mix = a + b;
  mix -= c;

  Vec scaled = 1.25 * a;
  scaled += (-0.5) * b;
  scaled += (0.25) * c;

  Vec plus_assign = a;
  plus_assign += b;

  Vec minus_assign = a;
  minus_assign -= b;

  Vec accum;
  accum.zero();
  accum.add_scaled(a, 1.25);
  accum.add_scaled(b, -0.5);
  accum.subtract_scaled(c, -0.25);

  const auto divided = a / 5.0;
  const auto outer_right = libMesh::outer_product(a, 5.0);
  const auto outer_left = libMesh::outer_product(5.0, a);

  Vec mult_assign = a;
  mult_assign *= 5.0;

  Vec div_assign = a;
  div_assign /= 5.0;

  Vec assign_zero = a;
  assign_zero = 0.0;

  result.vectors.push_back(copied);
  result.vectors.push_back(mix);
  result.vectors.push_back(scaled);
  result.vectors.push_back(accum);
  result.vectors.push_back(plus_assign);
  result.vectors.push_back(minus_assign);
  result.vectors.push_back(divided);
  result.vectors.push_back(outer_right);
  result.vectors.push_back(outer_left);
  result.vectors.push_back(mult_assign);
  result.vectors.push_back(div_assign);

  result.scalars.push_back(a * b);
  result.scalars.push_back(a * b);
  result.scalars.push_back(a.contract(b));
  result.scalars.push_back(mix.norm());
  result.scalars.push_back(mix.norm_sq());
  result.scalars.push_back(make_vector<Vec>(0.0, 0.0, 0.0).is_zero() ? 1.0 : 0.0);
  result.scalars.push_back(mix.is_zero() ? 1.0 : 0.0);
  result.scalars.push_back((a == a) ? 1.0 : 0.0);
  result.scalars.push_back((a == b) ? 1.0 : 0.0);
  result.scalars.push_back((a != a) ? 1.0 : 0.0);
  result.scalars.push_back((a != b) ? 1.0 : 0.0);
  result.scalars.push_back(assign_zero.is_zero() ? 1.0 : 0.0);

  const auto xvec = make_vector<Vec>(1.3);
  result.scalars.push_back(libMesh::solid_angle(as_type_vector(xvec),
                                                as_type_vector(xvec),
                                                as_type_vector(xvec)));

#if LIBMESH_DIM > 1
  const auto yvec = make_vector<Vec>(0.0, 2.7);
  const auto xydiag = make_vector<Vec>(3.1, 3.1);
  result.scalars.push_back(libMesh::solid_angle(as_type_vector(xvec),
                                                as_type_vector(xvec),
                                                as_type_vector(yvec)));
  result.scalars.push_back(libMesh::solid_angle(as_type_vector(xvec),
                                                as_type_vector(yvec),
                                                as_type_vector(xydiag)));
#endif

#if LIBMESH_DIM > 2
  const auto xypdiag = make_vector<Vec>(0.8, -0.8);
  const auto zvec = make_vector<Vec>(0.0, 0.0, 1.1);
  const auto xzdiag = make_vector<Vec>(0.0, 0.7, 0.7);
  const auto icosa1 = make_vector<Vec>(1.0, golden_ratio, 0.0);
  const auto icosa2 = make_vector<Vec>(-1.0, golden_ratio, 0.0);
  const auto icosa3 = make_vector<Vec>(0.0, 1.0, golden_ratio);
  result.scalars.push_back(libMesh::solid_angle(as_type_vector(xydiag),
                                                as_type_vector(yvec),
                                                as_type_vector(zvec)));
  result.scalars.push_back(libMesh::solid_angle(as_type_vector(xvec),
                                                as_type_vector(yvec),
                                                as_type_vector(xzdiag)));
  result.scalars.push_back(libMesh::solid_angle(as_type_vector(xypdiag),
                                                as_type_vector(xydiag),
                                                as_type_vector(zvec)));
  result.scalars.push_back(libMesh::solid_angle(as_type_vector(icosa1),
                                                as_type_vector(icosa2),
                                                as_type_vector(icosa3)));
#endif

#if LIBMESH_DIM > 2
  const auto cross = a.cross(b);
  auto unit_cross = cross;
  if (cross.norm() > unit_tol)
    unit_cross = cross.unit();

  result.vectors.push_back(cross);
  result.vectors.push_back(unit_cross);
#endif

  libmesh_assert_equal_to(result.vectors.size(), vector_results);
  libmesh_assert_equal_to(result.scalars.size(), scalar_results);

  return result;
}

} // namespace KokkosVectorOracle
} // namespace libMeshTest

#endif
