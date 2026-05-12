#ifndef KOKKOS_VECTOR_OPS_ORACLE_RUNNERS_H
#define KOKKOS_VECTOR_OPS_ORACLE_RUNNERS_H

#include "kokkos_vector_ops_oracle_fixtures.h"

#include <cstdio>

namespace libMeshTest
{
namespace KokkosVectorOracle
{

template <typename StoragePolicy, typename Vec>
static int
test_vector_ops_case(const vector_case & info)
{
  const auto a = make_vector<Vec>(info.ax, info.ay, info.az);
  const auto b = make_vector<Vec>(info.bx, info.by, info.bz);
  const auto c = make_vector<Vec>(info.cx, info.cy, info.cz);

  const auto expected = build_host_oracle(a, b, c);

  auto d_a = libMesh::Kokkos::upload_vector_storage<StoragePolicy>(std::vector<Vec>{a}, "vector_ops_a");
  auto d_b = libMesh::Kokkos::upload_vector_storage<StoragePolicy>(std::vector<Vec>{b}, "vector_ops_b");
  auto d_c = libMesh::Kokkos::upload_vector_storage<StoragePolicy>(std::vector<Vec>{c}, "vector_ops_c");
  auto d_vectors = libMesh::Kokkos::make_vector_storage<StoragePolicy>("vector_ops_vectors", vector_results);
  ::Kokkos::View<double *> d_scalars("vector_ops_scalars", scalar_results);

  ::Kokkos::parallel_for(
    1,
    KOKKOS_LAMBDA(int) {
      const auto a_ref = libMesh::Kokkos::make_vector_ref(d_a, 0);
      const auto b_ref = libMesh::Kokkos::make_vector_ref(d_b, 0);
      const auto c_ref = libMesh::Kokkos::make_vector_ref(d_c, 0);

      const Vec copied = libMesh::Kokkos::copy_vector<Vec>(a_ref);
      const Vec mix = a_ref + b_ref - c_ref;
      const Vec scaled = Real(1.25) * a_ref + Real(-0.5) * b_ref + Real(0.25) * c_ref;
      const Vec plus_assign = a_ref + b_ref;
      const Vec minus_assign = a_ref - b_ref;
      const Vec accum = Real(1.25) * a_ref + Real(-0.5) * b_ref + Real(0.25) * c_ref;
      const Vec divided = a_ref / Real(5.0);
      const Vec outer_right = Real(5.0) * a_ref;
      const Vec outer_left = a_ref * Real(5.0);
      const Vec mult_assign = a_ref * Real(5.0);
      const Vec div_assign = a_ref / Real(5.0);
      const Vec assign_zero = libMesh::Kokkos::zero_vector_value<Vec>();

      const Real dot = libMesh::Kokkos::vector_dot(a_ref, b_ref);
      const Real contract = a_ref.contract(b_ref);
      const Real norm = mix.norm();
      const Real norm_sq = mix.norm_sq();
      const Vec zero = libMesh::Kokkos::zero_vector_value<Vec>();

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
      d_scalars(scalar_offset++) = a_ref * b_ref;
      d_scalars(scalar_offset++) = contract;
      d_scalars(scalar_offset++) = norm;
      d_scalars(scalar_offset++) = norm_sq;
      d_scalars(scalar_offset++) = zero.is_zero() ? 1.0 : 0.0;
      d_scalars(scalar_offset++) = mix.is_zero() ? 1.0 : 0.0;
      d_scalars(scalar_offset++) = (a_ref == a_ref) ? 1.0 : 0.0;
      d_scalars(scalar_offset++) = (a_ref == b_ref) ? 1.0 : 0.0;
      d_scalars(scalar_offset++) = (a_ref != a_ref) ? 1.0 : 0.0;
      d_scalars(scalar_offset++) = (a_ref != b_ref) ? 1.0 : 0.0;
      d_scalars(scalar_offset++) = assign_zero.is_zero() ? 1.0 : 0.0;

      const Vec xvec = make_vector<Vec>(1.3);
      d_scalars(scalar_offset++) = libMesh::Kokkos::vector_solid_angle(xvec, xvec, xvec);

#if LIBMESH_DIM > 1
      const Vec yvec = make_vector<Vec>(0.0, 2.7);
      const Vec xydiag = make_vector<Vec>(3.1, 3.1);
      d_scalars(scalar_offset++) = libMesh::Kokkos::vector_solid_angle(xvec, xvec, yvec);
      d_scalars(scalar_offset++) = libMesh::Kokkos::vector_solid_angle(xvec, yvec, xydiag);
#endif

#if LIBMESH_DIM > 2
      const Vec xypdiag = make_vector<Vec>(0.8, -0.8);
      const Vec zvec = make_vector<Vec>(0.0, 0.0, 1.1);
      const Vec xzdiag = make_vector<Vec>(0.0, 0.7, 0.7);
      const Vec icosa1 = make_vector<Vec>(1.0, golden_ratio, 0.0);
      const Vec icosa2 = make_vector<Vec>(-1.0, golden_ratio, 0.0);
      const Vec icosa3 = make_vector<Vec>(0.0, 1.0, golden_ratio);
      d_scalars(scalar_offset++) = libMesh::Kokkos::vector_solid_angle(xydiag, yvec, zvec);
      d_scalars(scalar_offset++) = libMesh::Kokkos::vector_solid_angle(xvec, yvec, xzdiag);
      d_scalars(scalar_offset++) = libMesh::Kokkos::vector_solid_angle(xypdiag, xydiag, zvec);
      d_scalars(scalar_offset++) = libMesh::Kokkos::vector_solid_angle(icosa1, icosa2, icosa3);
#endif

#if LIBMESH_DIM > 2
      const Vec cross = a_ref.cross(b_ref);
      Vec unit_cross = cross;
      if (libMesh::Kokkos::vector_norm(cross) > unit_tol)
        unit_cross = libMesh::Kokkos::vector_unit<Vec>(cross);

      libMesh::Kokkos::store_vector(d_vectors, vector_offset++, cross);
      libMesh::Kokkos::store_vector(d_vectors, vector_offset++, unit_cross);
#endif

      libmesh_assert_equal_to(vector_offset, vector_results);
      libmesh_assert_equal_to(scalar_offset, scalar_results);
    });
  ::Kokkos::fence();

  return libMeshTest::KokkosOracle::compare_device_vectors(d_vectors, expected.vectors, tol) +
         libMeshTest::KokkosOracle::compare_device_scalars(d_scalars, expected.scalars, tol);
}

template <typename StoragePolicy, typename Vec>
int
run_vector_cases(const char * suite_name)
{
  int fail = 0;

  for (const auto & info : cases)
  {
    const int f = test_vector_ops_case<StoragePolicy, Vec>(info);
    std::printf("[%s] [%s] [%s] %s  (%d failures)\n",
                suite_name,
                libMesh::Kokkos::storage_policy_name<StoragePolicy>(),
                info.name,
                f ? "FAIL" : "PASS",
                f);
    fail += f;
  }

  return fail;
}

inline int
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

template <typename StoragePolicy, typename Vec>
static int
test_mixed_representation_ops()
{
  int fail = 0;

  const auto a = make_vector<Vec>(2.0, 3.0, 4.0);
  const auto b = make_vector<Vec>(5.0, -6.0, 7.0);
  const auto c = make_vector<Vec>(1.25, -0.5, 2.0);

  auto d_a = libMesh::Kokkos::upload_vector_storage<StoragePolicy>(std::vector<Vec>{a}, "mixed_vector_a");
  auto d_b = libMesh::Kokkos::upload_vector_storage<StoragePolicy>(std::vector<Vec>{b}, "mixed_vector_b");

  auto d_vectors =
    libMesh::Kokkos::make_vector_storage<StoragePolicy>("mixed_vector_vectors", (LIBMESH_DIM > 2) ? 5 : 3);
  ::Kokkos::View<double *> d_scalars("mixed_vector_scalars", (LIBMESH_DIM > 2) ? 7 : 5);

  const auto ref_sum = a + b;
  const auto ref_diff = a - b;
  const auto ref_scaled = 1.5 * a;
  const auto ref_dot = a * b;
  const auto ref_contract = a.contract(b);
  const auto ref_solid_angle =
    libMesh::solid_angle(as_type_vector(a), as_type_vector(b), as_type_vector(c));
  const auto ref_cross_norm_sq = libMesh::cross_norm_sq(as_type_vector(a), as_type_vector(b));

#if LIBMESH_DIM > 2
  const auto ref_cross = a.cross(b);
  auto ref_unit_cross = ref_cross;
  if (ref_cross.norm() > unit_tol)
    ref_unit_cross = ref_cross.unit();
#endif

  ::Kokkos::parallel_for(
    1,
    KOKKOS_LAMBDA(int) {
      const auto a_ref = libMesh::Kokkos::make_vector_ref(d_a, 0);
      const auto b_ref = libMesh::Kokkos::make_vector_ref(d_b, 0);

      const auto sum = a_ref + b;
      const auto diff = a - b_ref;
      const auto scaled = Real(1.5) * a_ref;

      libMesh::Kokkos::store_vector(d_vectors, 0, sum);
      libMesh::Kokkos::store_vector(d_vectors, 1, diff);
      libMesh::Kokkos::store_vector(d_vectors, 2, scaled);

      d_scalars(0) = a_ref * b;
      d_scalars(1) = b_ref.contract(a);
      d_scalars(2) = (a_ref == a) ? 1.0 : 0.0;
      d_scalars(3) = (a_ref != b) ? 1.0 : 0.0;
      d_scalars(4) = libMesh::Kokkos::vector_solid_angle(a_ref, b, c);

#if LIBMESH_DIM > 2
      const auto cross = a_ref.cross(b);
      Vec unit_cross = cross;
      if (libMesh::Kokkos::vector_norm(cross) > unit_tol)
        unit_cross = libMesh::Kokkos::vector_unit<Vec>(cross);

      libMesh::Kokkos::store_vector(d_vectors, 3, cross);
      libMesh::Kokkos::store_vector(d_vectors, 4, unit_cross);
      d_scalars(5) = libMesh::Kokkos::vector_cross_norm_sq(a_ref, b);
      d_scalars(6) = (cross == libMesh::Kokkos::vector_cross<Vec>(a, b_ref)) ? 1.0 : 0.0;
#endif
    });
  ::Kokkos::fence();

  fail += libMeshTest::KokkosOracle::compare_device_vectors(
    d_vectors,
    [&]() {
      std::vector<Vec> ref = {ref_sum, ref_diff, ref_scaled};
#if LIBMESH_DIM > 2
      ref.push_back(ref_cross);
      ref.push_back(ref_unit_cross);
#endif
      return ref;
    }(),
    tol);

  fail += libMeshTest::KokkosOracle::compare_device_scalars(
    d_scalars,
    [&]() {
      std::vector<Real> ref = {ref_dot, ref_contract, 1.0, 1.0, ref_solid_angle};
#if LIBMESH_DIM > 2
      ref.push_back(ref_cross_norm_sq);
      ref.push_back(1.0);
#endif
      return ref;
    }(),
    tol);

  return fail;
}

inline int
run_all_oracles()
{
  int total_fail = 0;

  total_fail += run_vector_cases<libMesh::Kokkos::layout_left_storage_policy, libMesh::TypeVector<Real>>(
    "typevector_kernel_oracle");
  total_fail += run_vector_cases<libMesh::Kokkos::layout_right_storage_policy, libMesh::TypeVector<Real>>(
    "typevector_kernel_oracle");
  total_fail += run_vector_cases<libMesh::Kokkos::layout_left_storage_policy, libMesh::VectorValue<Real>>(
    "vectorvalue_kernel_oracle");
  total_fail += run_vector_cases<libMesh::Kokkos::layout_right_storage_policy, libMesh::VectorValue<Real>>(
    "vectorvalue_kernel_oracle");

  const int mixed_typevector_left =
    test_mixed_representation_ops<libMesh::Kokkos::layout_left_storage_policy, libMesh::TypeVector<Real>>();
  std::printf("[vector_mixed_representation_oracle] [%s] [typevector] %s  (%d failures)\n",
              libMesh::Kokkos::storage_policy_name<libMesh::Kokkos::layout_left_storage_policy>(),
              mixed_typevector_left ? "FAIL" : "PASS",
              mixed_typevector_left);
  total_fail += mixed_typevector_left;

  const int mixed_typevector_right =
    test_mixed_representation_ops<libMesh::Kokkos::layout_right_storage_policy, libMesh::TypeVector<Real>>();
  std::printf("[vector_mixed_representation_oracle] [%s] [typevector] %s  (%d failures)\n",
              libMesh::Kokkos::storage_policy_name<libMesh::Kokkos::layout_right_storage_policy>(),
              mixed_typevector_right ? "FAIL" : "PASS",
              mixed_typevector_right);
  total_fail += mixed_typevector_right;

  const int mixed_vectorvalue_left =
    test_mixed_representation_ops<libMesh::Kokkos::layout_left_storage_policy, libMesh::VectorValue<Real>>();
  std::printf("[vector_mixed_representation_oracle] [%s] [vectorvalue] %s  (%d failures)\n",
              libMesh::Kokkos::storage_policy_name<libMesh::Kokkos::layout_left_storage_policy>(),
              mixed_vectorvalue_left ? "FAIL" : "PASS",
              mixed_vectorvalue_left);
  total_fail += mixed_vectorvalue_left;

  const int mixed_vectorvalue_right =
    test_mixed_representation_ops<libMesh::Kokkos::layout_right_storage_policy, libMesh::VectorValue<Real>>();
  std::printf("[vector_mixed_representation_oracle] [%s] [vectorvalue] %s  (%d failures)\n",
              libMesh::Kokkos::storage_policy_name<libMesh::Kokkos::layout_right_storage_policy>(),
              mixed_vectorvalue_right ? "FAIL" : "PASS",
              mixed_vectorvalue_right);
  total_fail += mixed_vectorvalue_right;

  const int host_fail = test_vector_host_only_traits();
  std::printf("[vector_host_traits_oracle] %s  (%d failures)\n",
              host_fail ? "FAIL" : "PASS",
              host_fail);
  total_fail += host_fail;

  return total_fail;
}

} // namespace KokkosVectorOracle
} // namespace libMeshTest

#endif
