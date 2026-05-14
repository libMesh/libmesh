// Kokkos kernel oracle tests for libMesh::TypeVector / VectorValue numerics.
//
// Compares device-side vector algebra
// against host libMesh TypeVector oracles built from direct hand-written
// fixtures.
//
// Returns 0 on success, non-zero on failure.

#include "libmesh/libmesh_config.h"

#include "libmesh/libmesh.h"
#include "libmesh/tensor_value.h"
#include "libmesh/type_vector.h"
#include "libmesh/vector_value.h"
#include "gpu/kokkos_storage.h"

// Avoid conflicting complex operators between CUDA and PETSc
#define PETSC_SKIP_CXX_COMPLEX_FIX 1
#include <Kokkos_Core.hpp>
#undef __CUDACC_VER__

#include <cmath>
#include <cstdio>
#include <string>
#include <type_traits>
#include <vector>

using libMesh::Real;

static constexpr double tol = 2.0e-13;
static constexpr double unit_tol = 1.0e-14;
static constexpr Real golden_ratio = 1.6180339887498948482;
static constexpr unsigned int solid_angle_results =
  1 + ((LIBMESH_DIM > 1) ? 2u : 0u) + ((LIBMESH_DIM > 2) ? 4u : 0u);
static constexpr unsigned int vector_results =
  11 + ((LIBMESH_DIM > 2) ? 2u : 0u);
static constexpr unsigned int scalar_results = 11 + solid_angle_results;

namespace
{

using vector_storage_view = Kokkos::View<Real * [LIBMESH_DIM]>;

inline int
compare_device_scalars(const Kokkos::View<double *> & d_values,
                       const std::vector<Real> & ref_values,
                       const double local_tol = tol)
{
  auto h_values = Kokkos::create_mirror_view(d_values);
  Kokkos::deep_copy(h_values, d_values);

  int fail = 0;
  for (std::size_t i = 0; i < ref_values.size(); ++i)
    if (std::fabs(h_values(i) - ref_values[i]) > local_tol)
      ++fail;

  return fail;
}

template <typename Vec>
inline vector_storage_view
upload_vectors(const std::vector<Vec> & values, const char * label)
{
  vector_storage_view d(std::string(label), values.size());
  auto h = Kokkos::create_mirror_view(d);
  for (std::size_t i = 0; i < values.size(); ++i)
    for (unsigned int component = 0; component < LIBMESH_DIM; ++component)
      h(i, component) = values[i](component);
  Kokkos::deep_copy(d, h);
  return d;
}

template <typename Vec>
inline int
compare_device_vectors(const vector_storage_view & d_values,
                       const std::vector<Vec> & ref_values,
                       const double local_tol = tol)
{
  auto h_values = Kokkos::create_mirror_view(d_values);
  Kokkos::deep_copy(h_values, d_values);

  int fail = 0;
  for (std::size_t i = 0; i < ref_values.size(); ++i)
    for (unsigned int d = 0; d < LIBMESH_DIM; ++d)
      if (std::fabs(h_values(i, d) - ref_values[i](d)) > local_tol)
        ++fail;

  return fail;
}

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

template <typename Vec>
struct vector_oracle
{
  std::vector<Vec> vectors;
  std::vector<Real> scalars;
};

struct vector_kernel_case
{
  const char * name;
  Real ax, ay, az;
  Real bx, by, bz;
  Real cx, cy, cz;
};

static const vector_kernel_case cases[] = {
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
inline vector_oracle<Vec>
make_manual_oracle(const Vec & a, const Vec & b, const Vec & c)
{
  vector_oracle<Vec> oracle;
  oracle.vectors.reserve(vector_results);
  oracle.scalars.reserve(scalar_results);

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

  oracle.vectors.push_back(copied);
  oracle.vectors.push_back(mix);
  oracle.vectors.push_back(scaled);
  oracle.vectors.push_back(accum);
  oracle.vectors.push_back(plus_assign);
  oracle.vectors.push_back(minus_assign);
  oracle.vectors.push_back(divided);
  oracle.vectors.push_back(outer_right);
  oracle.vectors.push_back(outer_left);
  oracle.vectors.push_back(mult_assign);
  oracle.vectors.push_back(div_assign);

  oracle.scalars.push_back(a * b);
  oracle.scalars.push_back(a.contract(b));
  oracle.scalars.push_back(mix.norm());
  oracle.scalars.push_back(mix.norm_sq());
  oracle.scalars.push_back(make_vector<Vec>(0.0, 0.0, 0.0).is_zero() ? 1.0 : 0.0);
  oracle.scalars.push_back(mix.is_zero() ? 1.0 : 0.0);
  oracle.scalars.push_back((a == a) ? 1.0 : 0.0);
  oracle.scalars.push_back((a == b) ? 1.0 : 0.0);
  oracle.scalars.push_back((a != a) ? 1.0 : 0.0);
  oracle.scalars.push_back((a != b) ? 1.0 : 0.0);
  oracle.scalars.push_back(assign_zero.is_zero() ? 1.0 : 0.0);

  const auto xvec = make_vector<Vec>(1.3);
  oracle.scalars.push_back(solid_angle(xvec, xvec, xvec));

#if LIBMESH_DIM > 1
  const auto yvec = make_vector<Vec>(0.0, 2.7);
  const auto xydiag = make_vector<Vec>(3.1, 3.1);
  oracle.scalars.push_back(solid_angle(xvec, xvec, yvec));
  oracle.scalars.push_back(solid_angle(xvec, yvec, xydiag));
#endif

#if LIBMESH_DIM > 2
  const auto xypdiag = make_vector<Vec>(0.8, -0.8);
  const auto zvec = make_vector<Vec>(0.0, 0.0, 1.1);
  const auto xzdiag = make_vector<Vec>(0.0, 0.7, 0.7);
  const auto icosa1 = make_vector<Vec>(1.0, golden_ratio, 0.0);
  const auto icosa2 = make_vector<Vec>(-1.0, golden_ratio, 0.0);
  const auto icosa3 = make_vector<Vec>(0.0, 1.0, golden_ratio);
  oracle.scalars.push_back(solid_angle(xydiag, yvec, zvec));
  oracle.scalars.push_back(solid_angle(xvec, yvec, xzdiag));
  oracle.scalars.push_back(solid_angle(xypdiag, xydiag, zvec));
  oracle.scalars.push_back(solid_angle(icosa1, icosa2, icosa3));
#endif

#if LIBMESH_DIM > 2
  const auto cross = a.cross(b);
  auto unit_cross = cross;
  if (cross.norm() > unit_tol)
    unit_cross = cross.unit();

  oracle.vectors.push_back(cross);
  oracle.vectors.push_back(unit_cross);
#endif

  libmesh_assert_equal_to(oracle.vectors.size(), vector_results);
  libmesh_assert_equal_to(oracle.scalars.size(), scalar_results);

  return oracle;
}

template <typename Vec>
static int
test_vector_ops_case(const vector_kernel_case & info)
{
  const auto a = make_vector<Vec>(info.ax, info.ay, info.az);
  const auto b = make_vector<Vec>(info.bx, info.by, info.bz);
  const auto c = make_vector<Vec>(info.cx, info.cy, info.cz);

  const auto oracle = make_manual_oracle(a, b, c);

  auto d_a = upload_vectors(std::vector<Vec>{a}, "vector_ops_a");
  auto d_b = upload_vectors(std::vector<Vec>{b}, "vector_ops_b");
  auto d_c = upload_vectors(std::vector<Vec>{c}, "vector_ops_c");
  vector_storage_view d_vectors("vector_ops_vectors", vector_results);
  Kokkos::View<double *> d_scalars("vector_ops_scalars", scalar_results);

  Kokkos::parallel_for(
    1,
    KOKKOS_LAMBDA(int) {
      const Vec a_d = libMesh::Kokkos::load_vector<Vec>(d_a, 0);
      const Vec b_d = libMesh::Kokkos::load_vector<Vec>(d_b, 0);
      const Vec c_d = libMesh::Kokkos::load_vector<Vec>(d_c, 0);

      const Vec copied = a_d;

      Vec mix = a_d + b_d;
      mix -= c_d;

      Vec scaled = 1.25 * a_d;
      scaled += (-0.5) * b_d;
      scaled += (0.25) * c_d;

      Vec plus_assign = a_d;
      plus_assign += b_d;

      Vec minus_assign = a_d;
      minus_assign -= b_d;

      Vec accum;
      accum.zero();
      accum.add_scaled(a_d, 1.25);
      accum.add_scaled(b_d, -0.5);
      accum.subtract_scaled(c_d, -0.25);

      const Vec divided = a_d / 5.0;
      const Vec outer_right = libMesh::outer_product(a_d, 5.0);
      const Vec outer_left = libMesh::outer_product(5.0, a_d);

      Vec mult_assign = a_d;
      mult_assign *= 5.0;

      Vec div_assign = a_d;
      div_assign /= 5.0;

      Vec assign_zero = a_d;
      assign_zero = 0.0;

      const Real dot = a_d * b_d;
      const Real contract = a_d.contract(b_d);
      const Real norm = mix.norm();
      const Real norm_sq = mix.norm_sq();
      const Vec zero = make_vector<Vec>(0.0, 0.0, 0.0);

      unsigned int vector_offset = 0;
      libMesh::Kokkos::store_vector(d_vectors, vector_offset++, copied);
      libMesh::Kokkos::store_vector(d_vectors, vector_offset++, mix);
      libMesh::Kokkos::store_vector(d_vectors, vector_offset++, scaled);
      libMesh::Kokkos::store_vector(d_vectors, vector_offset++, accum);
      libMesh::Kokkos::store_vector(d_vectors, vector_offset++, plus_assign);
      libMesh::Kokkos::store_vector(d_vectors, vector_offset++, minus_assign);
      libMesh::Kokkos::store_vector(d_vectors, vector_offset++, divided);
      libMesh::Kokkos::store_vector(d_vectors, vector_offset++, outer_right);
      libMesh::Kokkos::store_vector(d_vectors, vector_offset++, outer_left);
      libMesh::Kokkos::store_vector(d_vectors, vector_offset++, mult_assign);
      libMesh::Kokkos::store_vector(d_vectors, vector_offset++, div_assign);

      unsigned int scalar_offset = 0;
      d_scalars(scalar_offset++) = dot;
      d_scalars(scalar_offset++) = contract;
      d_scalars(scalar_offset++) = norm;
      d_scalars(scalar_offset++) = norm_sq;
      d_scalars(scalar_offset++) = zero.is_zero() ? 1.0 : 0.0;
      d_scalars(scalar_offset++) = mix.is_zero() ? 1.0 : 0.0;
      d_scalars(scalar_offset++) = (a_d == a_d) ? 1.0 : 0.0;
      d_scalars(scalar_offset++) = (a_d == b_d) ? 1.0 : 0.0;
      d_scalars(scalar_offset++) = (a_d != a_d) ? 1.0 : 0.0;
      d_scalars(scalar_offset++) = (a_d != b_d) ? 1.0 : 0.0;
      d_scalars(scalar_offset++) = assign_zero.is_zero() ? 1.0 : 0.0;

      const Vec xvec = make_vector<Vec>(1.3);
      d_scalars(scalar_offset++) = solid_angle(xvec, xvec, xvec);

#if LIBMESH_DIM > 1
      const Vec yvec = make_vector<Vec>(0.0, 2.7);
      const Vec xydiag = make_vector<Vec>(3.1, 3.1);
      d_scalars(scalar_offset++) = solid_angle(xvec, xvec, yvec);
      d_scalars(scalar_offset++) = solid_angle(xvec, yvec, xydiag);
#endif

#if LIBMESH_DIM > 2
      const Vec xypdiag = make_vector<Vec>(0.8, -0.8);
      const Vec zvec = make_vector<Vec>(0.0, 0.0, 1.1);
      const Vec xzdiag = make_vector<Vec>(0.0, 0.7, 0.7);
      const Vec icosa1 = make_vector<Vec>(1.0, golden_ratio, 0.0);
      const Vec icosa2 = make_vector<Vec>(-1.0, golden_ratio, 0.0);
      const Vec icosa3 = make_vector<Vec>(0.0, 1.0, golden_ratio);
      d_scalars(scalar_offset++) = solid_angle(xydiag, yvec, zvec);
      d_scalars(scalar_offset++) = solid_angle(xvec, yvec, xzdiag);
      d_scalars(scalar_offset++) = solid_angle(xypdiag, xydiag, zvec);
      d_scalars(scalar_offset++) = solid_angle(icosa1, icosa2, icosa3);
#endif

#if LIBMESH_DIM > 2
      const Vec cross = a_d.cross(b_d);
      Vec unit_cross = cross;
      if (cross.norm() > unit_tol)
        unit_cross = cross.unit();

      libMesh::Kokkos::store_vector(d_vectors, vector_offset++, cross);
      libMesh::Kokkos::store_vector(d_vectors, vector_offset++, unit_cross);
#endif

      libmesh_assert_equal_to(vector_offset, vector_results);
      libmesh_assert_equal_to(scalar_offset, scalar_results);
    });
  Kokkos::fence();

  return compare_device_vectors(d_vectors, oracle.vectors) +
         compare_device_scalars(d_scalars, oracle.scalars);
}

static int
test_vector_host_only_traits()
{
  int fail = 0;

#ifdef LIBMESH_HAVE_METAPHYSICL
  typedef typename MetaPhysicL::ReplaceAlgebraicType<
      std::vector<libMesh::TypeVector<double>>,
      typename libMesh::TensorTools::IncrementRank<
          typename MetaPhysicL::ValueType<std::vector<libMesh::TypeVector<double>>>::type>::type>::type
      ReplacedType;
  constexpr bool typevector_assertion =
    std::is_same<ReplacedType, std::vector<libMesh::TensorValue<double>>>::value;
  fail += typevector_assertion ? 0 : 1;

  typedef typename MetaPhysicL::ReplaceAlgebraicType<
      std::vector<libMesh::VectorValue<double>>,
      typename libMesh::TensorTools::IncrementRank<
          typename MetaPhysicL::ValueType<std::vector<libMesh::VectorValue<double>>>::type>::type>::type
      ReplacedValueType;
  constexpr bool vectorvalue_assertion =
    std::is_same<ReplacedValueType, std::vector<libMesh::TensorValue<double>>>::value;
  fail += vectorvalue_assertion ? 0 : 1;
#endif

  return fail;
}

} // anonymous namespace

int
main(int argc, char ** argv)
{
  Kokkos::initialize(argc, argv);
  libMesh::LibMeshInit init(argc, argv);

  int total_fail = 0;

  for (const auto & info : cases)
  {
    const int f = test_vector_ops_case<libMesh::TypeVector<Real>>(info);
    std::printf("[typevector_kernel_oracle] [%s] %s  (%d failures)\n",
                info.name, f ? "FAIL" : "PASS", f);
    total_fail += f;
  }

  for (const auto & info : cases)
  {
    const int f = test_vector_ops_case<libMesh::VectorValue<Real>>(info);
    std::printf("[vectorvalue_kernel_oracle] [%s] %s  (%d failures)\n",
                info.name, f ? "FAIL" : "PASS", f);
    total_fail += f;
  }

  const int host_fail = test_vector_host_only_traits();
  std::printf("[vector_host_traits_oracle] %s  (%d failures)\n",
              host_fail ? "FAIL" : "PASS",
              host_fail);
  total_fail += host_fail;

  Kokkos::finalize();

  if (total_fail == 0)
    std::printf("ALL TESTS PASSED\n");
  else
    std::printf("%d TEST(S) FAILED\n", total_fail);

  return total_fail ? 1 : 0;
}
