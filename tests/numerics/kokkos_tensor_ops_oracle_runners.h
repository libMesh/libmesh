#ifndef KOKKOS_TENSOR_OPS_ORACLE_RUNNERS_H
#define KOKKOS_TENSOR_OPS_ORACLE_RUNNERS_H

#include "kokkos_tensor_ops_oracle_fixtures.h"

#include <cmath>
#include <cstdio>

namespace libMeshTest
{
namespace KokkosTensorOracle
{

template <typename StoragePolicy>
static int
test_dim_ops()
{
  const unsigned int ncases = sizeof(dim_cases) / sizeof(dim_cases[0]);

  std::vector<oracle_tensor> J_values(ncases);
  std::vector<unsigned int> dims(ncases);
  std::vector<Real> ref_det(ncases);
  std::vector<oracle_tensor> ref_inv(ncases);
  std::vector<oracle_tensor> ref_I(ncases);
  std::vector<oracle_tensor> ref_prod_left(ncases);
  std::vector<oracle_tensor> ref_prod_right(ncases);

  for (unsigned int c = 0; c < ncases; ++c)
  {
    const auto & info = dim_cases[c];
    J_values[c] = info.J;
    dims[c] = info.dim;

    ref_det[c] = host_leading_determinant(info.J, info.dim);
    ref_inv[c] = host_leading_inverse(info.J, info.dim);
    ref_I[c] = build_identity_tensor(info.dim);
    ref_prod_left[c] = info.J * ref_inv[c];
    ref_prod_right[c] = ref_inv[c] * info.J;
  }

  auto d_J = libMesh::Kokkos::upload_tensor_storage<StoragePolicy>(J_values, "tensor_dim_ops_J");
  auto d_dims = libMeshTest::KokkosOracle::upload_objects(dims, "tensor_dim_ops_dim");
  ::Kokkos::View<double *> d_det("tensor_dim_ops_det", ncases);
  auto d_inv = libMesh::Kokkos::make_tensor_storage<StoragePolicy>("tensor_dim_ops_inv", ncases);
  auto d_I = libMesh::Kokkos::make_tensor_storage<StoragePolicy>("tensor_dim_ops_I", ncases);
  auto d_prod_left = libMesh::Kokkos::make_tensor_storage<StoragePolicy>("tensor_dim_ops_prod_left", ncases);
  auto d_prod_right = libMesh::Kokkos::make_tensor_storage<StoragePolicy>("tensor_dim_ops_prod_right", ncases);

  ::Kokkos::parallel_for(
    static_cast<int>(ncases),
    KOKKOS_LAMBDA(int c) {
      const auto J_ref = libMesh::Kokkos::make_tensor_ref(d_J, c);
      const unsigned int dim = d_dims(c);
      const Real det = libMesh::Kokkos::tensor_determinant(J_ref, dim);
      const auto inv = J_ref.inverse(dim);
      const auto I = libMesh::Kokkos::tensor_identity<oracle_tensor>(dim);
      const auto prod_left = J_ref * inv;
      const auto prod_right = inv * J_ref;

      d_det(c) = det;
      libMesh::Kokkos::store_tensor(d_inv, c, inv);
      libMesh::Kokkos::store_tensor(d_I, c, I);
      libMesh::Kokkos::store_tensor(d_prod_left, c, prod_left);
      libMesh::Kokkos::store_tensor(d_prod_right, c, prod_right);
    });
  ::Kokkos::fence();

  return libMeshTest::KokkosOracle::compare_device_scalars(d_det, ref_det, tol) +
         libMeshTest::KokkosOracle::compare_device_tensors(d_inv, ref_inv, tol) +
         libMeshTest::KokkosOracle::compare_device_tensors(d_I, ref_I, tol) +
         libMeshTest::KokkosOracle::compare_device_tensors(d_prod_left, ref_prod_left, tol) +
         libMeshTest::KokkosOracle::compare_device_tensors(d_prod_right, ref_prod_right, tol);
}

template <typename StoragePolicy>
static int
test_tensor_ops()
{
  const auto A = make_host_tensor(1.1, -0.4, 0.7,
                                  0.3,  1.9, -1.2,
                                  -0.8, 0.5, 2.2);
  const auto a = make_host_vector(2.0, 3.0, 4.0);
  const auto b = make_host_vector(5.0, -6.0, 7.0);
  const auto c = make_host_vector(1.25, -0.5, 2.0);

  const auto outer = libMesh::outer_product(a, b);
  const auto transpose = A.transpose();
  const auto mix = 1.5 * A - 0.25 * outer;
  const auto right = A * c;
  const auto left = c * A;
  const Real contract = A.contract(outer);
  const Real norm = A.norm();
  const auto zero = libMesh::Kokkos::zero_tensor_value<oracle_tensor>();

  std::vector<oracle_tensor> ref_outer(1, outer);
  std::vector<oracle_tensor> ref_transpose(1, transpose);
  std::vector<oracle_tensor> ref_mix(1, mix);
  std::vector<oracle_vector> ref_rows(LIBMESH_DIM);
  std::vector<oracle_vector> ref_columns(LIBMESH_DIM);
  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
  {
    ref_rows[i] = A.row(i);
    ref_columns[i] = A.column(i);
  }
  std::vector<oracle_vector> ref_right(1, right);
  std::vector<oracle_vector> ref_left(1, left);
  std::vector<Real> ref_scalars = {contract, norm, zero.is_zero() ? 1.0 : 0.0, A.is_zero() ? 1.0 : 0.0};

  auto d_A = libMesh::Kokkos::upload_tensor_storage<StoragePolicy>(std::vector<oracle_tensor>{A}, "tensor_ops_A");
  auto d_a = libMesh::Kokkos::upload_vector_storage<StoragePolicy>(std::vector<oracle_vector>{a}, "tensor_ops_a");
  auto d_b = libMesh::Kokkos::upload_vector_storage<StoragePolicy>(std::vector<oracle_vector>{b}, "tensor_ops_b");
  auto d_c = libMesh::Kokkos::upload_vector_storage<StoragePolicy>(std::vector<oracle_vector>{c}, "tensor_ops_c");
  auto d_outer = libMesh::Kokkos::make_tensor_storage<StoragePolicy>("tensor_ops_outer", 1);
  auto d_transpose = libMesh::Kokkos::make_tensor_storage<StoragePolicy>("tensor_ops_transpose", 1);
  auto d_mix = libMesh::Kokkos::make_tensor_storage<StoragePolicy>("tensor_ops_mix", 1);
  auto d_rows = libMesh::Kokkos::make_vector_storage<StoragePolicy>("tensor_ops_rows", LIBMESH_DIM);
  auto d_columns = libMesh::Kokkos::make_vector_storage<StoragePolicy>("tensor_ops_columns", LIBMESH_DIM);
  auto d_right = libMesh::Kokkos::make_vector_storage<StoragePolicy>("tensor_ops_right", 1);
  auto d_left = libMesh::Kokkos::make_vector_storage<StoragePolicy>("tensor_ops_left", 1);
  ::Kokkos::View<double *> d_scalars("tensor_ops_scalars", 4);

  ::Kokkos::parallel_for(
    1,
    KOKKOS_LAMBDA(int) {
      const auto A_ref = libMesh::Kokkos::make_tensor_ref(d_A, 0);
      const auto a_ref = libMesh::Kokkos::make_vector_ref(d_a, 0);
      const auto b_ref = libMesh::Kokkos::make_vector_ref(d_b, 0);
      const auto c_ref = libMesh::Kokkos::make_vector_ref(d_c, 0);
      const auto outer_d = libMesh::Kokkos::tensor_outer_product<oracle_tensor>(a_ref, b_ref);
      const auto transpose_d = A_ref.transpose();
      const auto mix_d = Real(1.5) * A_ref - Real(0.25) * outer_d;
      const auto right_d = A_ref * c_ref;
      const auto left_d = c_ref * A_ref;
      const Real contract_d = A_ref.contract(outer_d);
      const Real norm_d = A_ref.norm();
      const bool zero_is_zero_d = libMesh::Kokkos::zero_tensor_value<oracle_tensor>().is_zero();
      const bool A_is_zero_d = A_ref.is_zero();

      for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
      {
        libMesh::Kokkos::store_vector(d_rows, i, A_ref.row(i));
        libMesh::Kokkos::store_vector(d_columns, i, A_ref.column(i));
      }

      libMesh::Kokkos::store_tensor(d_outer, 0, outer_d);
      libMesh::Kokkos::store_tensor(d_transpose, 0, transpose_d);
      libMesh::Kokkos::store_tensor(d_mix, 0, mix_d);
      libMesh::Kokkos::store_vector(d_right, 0, right_d);
      libMesh::Kokkos::store_vector(d_left, 0, left_d);
      d_scalars(0) = contract_d;
      d_scalars(1) = norm_d;
      d_scalars(2) = zero_is_zero_d ? 1.0 : 0.0;
      d_scalars(3) = A_is_zero_d ? 1.0 : 0.0;
    });
  ::Kokkos::fence();

  return libMeshTest::KokkosOracle::compare_device_tensors(d_outer, ref_outer, tol) +
         libMeshTest::KokkosOracle::compare_device_tensors(d_transpose, ref_transpose, tol) +
         libMeshTest::KokkosOracle::compare_device_tensors(d_mix, ref_mix, tol) +
         libMeshTest::KokkosOracle::compare_device_vectors(d_rows, ref_rows, tol) +
         libMeshTest::KokkosOracle::compare_device_vectors(d_columns, ref_columns, tol) +
         libMeshTest::KokkosOracle::compare_device_vectors(d_right, ref_right, tol) +
         libMeshTest::KokkosOracle::compare_device_vectors(d_left, ref_left, tol) +
         libMeshTest::KokkosOracle::compare_device_scalars(d_scalars, ref_scalars, tol);
}

inline int
test_tensor_host_only_ops()
{
  int fail = 0;

#if LIBMESH_DIM > 2
  {
    libMesh::TensorValue<double> tensor(2., 1., 0.,
                                        1., 2., 1.,
                                        0., 1., 2.);
    fail += tensor.is_hpd(/*rel_tol=*/0.) ? 0 : 1;
  }

  {
    libMesh::TensorValue<double> tensor(1., 0., 0.,
                                        0., 0., 1.,
                                        0., 1., 0.);
    fail += tensor.is_hpd() ? 1 : 0;
  }

  {
    const libMesh::Point x(1., 0., 0.);
    const auto R = libMesh::RealTensorValue::extrinsic_rotation_matrix(90., 0., 0.);
    const auto rotated = R * x;
    fail += (std::fabs(rotated(0)) <= tol) ? 0 : 1;
    fail += (std::fabs(rotated(1) - 1.) <= tol) ? 0 : 1;
    fail += (std::fabs(rotated(2)) <= tol) ? 0 : 1;

    const auto invR = libMesh::RealTensorValue::inverse_extrinsic_rotation_matrix(90., 0., 0.);
    const auto unrotated = invR * rotated;
    fail += (std::fabs(unrotated(0) - 1.) <= tol) ? 0 : 1;
    fail += (std::fabs(unrotated(1)) <= tol) ? 0 : 1;
    fail += (std::fabs(unrotated(2)) <= tol) ? 0 : 1;
  }

  {
    const libMesh::Point x(1., 1., 1.);
    const auto R = libMesh::RealTensorValue::extrinsic_rotation_matrix(90., 90., 90.);
    const auto rotated = R * x;
    fail += (std::fabs(rotated(0) - 1.) <= tol) ? 0 : 1;
    fail += (std::fabs(rotated(1) + 1.) <= tol) ? 0 : 1;
    fail += (std::fabs(rotated(2) - 1.) <= tol) ? 0 : 1;

    const auto invR = libMesh::RealTensorValue::inverse_extrinsic_rotation_matrix(90., 90., 90.);
    const auto unrotated = invR * rotated;
    fail += (std::fabs(unrotated(0) - 1.) <= tol) ? 0 : 1;
    fail += (std::fabs(unrotated(1) - 1.) <= tol) ? 0 : 1;
    fail += (std::fabs(unrotated(2) - 1.) <= tol) ? 0 : 1;
  }
#endif

#ifdef LIBMESH_HAVE_METAPHYSICL
  typedef typename MetaPhysicL::ReplaceAlgebraicType<
      std::vector<libMesh::TypeTensor<double>>,
      typename libMesh::TensorTools::IncrementRank<
          typename MetaPhysicL::ValueType<std::vector<libMesh::TypeTensor<double>>>::type>::type>::type
      ReplacedType;
  constexpr bool assertion =
    std::is_same<ReplacedType, std::vector<libMesh::TypeNTensor<3, double>>>::value;
  fail += assertion ? 0 : 1;
#endif

  return fail;
}

template <typename StoragePolicy>
static int
test_linalg_foundation_storage_roundtrip()
{
  int fail = 0;

  auto d_vector = libMesh::Kokkos::make_vector_storage<StoragePolicy>("foundation_vector", 1);
  auto d_tensor = libMesh::Kokkos::make_tensor_storage<StoragePolicy>("foundation_tensor", 1);

  {
    auto h_vector = ::Kokkos::create_mirror_view(d_vector);
    auto h_tensor = ::Kokkos::create_mirror_view(d_tensor);

    for (unsigned int d = 0; d < LIBMESH_DIM; ++d)
      h_vector(0, d) = Real(d + 1) * Real(0.5);

    for (unsigned int row = 0; row < LIBMESH_DIM; ++row)
      for (unsigned int col = 0; col < LIBMESH_DIM; ++col)
        h_tensor(0, row, col) = Real(10 * row + col + 1) * Real(0.25);

    ::Kokkos::deep_copy(d_vector, h_vector);
    ::Kokkos::deep_copy(d_tensor, h_tensor);
  }

  const auto vector_in = libMesh::Kokkos::make_vector_ref(d_vector, 0);
  const auto tensor_in = libMesh::Kokkos::make_tensor_ref(d_tensor, 0);

  const auto as_point = libMesh::Kokkos::materialize_vector<libMesh::Point>(vector_in);
  const auto as_vector_value =
    libMesh::Kokkos::materialize_vector<libMesh::VectorValue<Real>>(vector_in);
  const auto as_type_vector =
    libMesh::Kokkos::materialize_vector<libMesh::TypeVector<Real>>(vector_in);

  for (unsigned int d = 0; d < LIBMESH_DIM; ++d)
  {
    const Real expected = Real(d + 1) * Real(0.5);
    fail += (std::fabs(as_point(d) - expected) <= tol) ? 0 : 1;
    fail += (std::fabs(as_vector_value(d) - expected) <= tol) ? 0 : 1;
    fail += (std::fabs(as_type_vector(d) - expected) <= tol) ? 0 : 1;
  }

  const auto as_tensor_value =
    libMesh::Kokkos::materialize_tensor<libMesh::TensorValue<Real>>(tensor_in);
  const auto as_type_tensor =
    libMesh::Kokkos::materialize_tensor<libMesh::TypeTensor<Real>>(tensor_in);

  for (unsigned int row = 0; row < LIBMESH_DIM; ++row)
    for (unsigned int col = 0; col < LIBMESH_DIM; ++col)
    {
      const Real expected = Real(10 * row + col + 1) * Real(0.25);
      fail += (std::fabs(as_tensor_value(row, col) - expected) <= tol) ? 0 : 1;
      fail += (std::fabs(as_type_tensor(row, col) - expected) <= tol) ? 0 : 1;
    }

  auto d_vector_out = libMesh::Kokkos::make_vector_storage<StoragePolicy>("foundation_vector_out", 1);
  auto d_tensor_out = libMesh::Kokkos::make_tensor_storage<StoragePolicy>("foundation_tensor_out", 1);

  auto vector_out = libMesh::Kokkos::make_vector_ref(d_vector_out, 0);
  auto tensor_out = libMesh::Kokkos::make_tensor_ref(d_tensor_out, 0);

  vector_out.zero();
  vector_out.assign(as_vector_value);
  vector_out.add_scaled(as_type_vector, Real(0));
  vector_out.subtract_scaled(as_type_vector, Real(0));

  tensor_out.zero();
  tensor_out.assign(as_tensor_value);
  tensor_out.add_scaled(as_type_tensor, Real(0));
  tensor_out.subtract_scaled(as_type_tensor, Real(0));

  {
    auto h_vector_out = ::Kokkos::create_mirror_view(d_vector_out);
    auto h_tensor_out = ::Kokkos::create_mirror_view(d_tensor_out);
    ::Kokkos::deep_copy(h_vector_out, d_vector_out);
    ::Kokkos::deep_copy(h_tensor_out, d_tensor_out);

    for (unsigned int d = 0; d < LIBMESH_DIM; ++d)
      fail += (std::fabs(h_vector_out(0, d) - as_vector_value(d)) <= tol) ? 0 : 1;

    for (unsigned int row = 0; row < LIBMESH_DIM; ++row)
      for (unsigned int col = 0; col < LIBMESH_DIM; ++col)
        fail += (std::fabs(h_tensor_out(0, row, col) - as_tensor_value(row, col)) <= tol) ? 0 : 1;
  }

  return fail;
}

template <typename StoragePolicy>
static int
test_mixed_representation_ops()
{
  int fail = 0;

  const auto a = make_host_vector(2.0, 3.0, 4.0);
  const auto b = make_host_vector(5.0, -6.0, 7.0);
  const auto c = make_host_vector(1.25, -0.5, 2.0);
  const auto A = make_host_tensor(1.1, -0.4, 0.7,
                                  0.3,  1.9, -1.2,
                                  -0.8, 0.5, 2.2);

  auto d_a = libMesh::Kokkos::upload_vector_storage<StoragePolicy>(std::vector<oracle_vector>{a}, "mixed_ops_a");
  auto d_A = libMesh::Kokkos::upload_tensor_storage<StoragePolicy>(std::vector<oracle_tensor>{A}, "mixed_ops_A");

  ::Kokkos::View<double *> d_scalars("mixed_ops_scalars", 8);
  auto d_vectors = libMesh::Kokkos::make_vector_storage<StoragePolicy>("mixed_ops_vectors", 5);
  auto d_tensors = libMesh::Kokkos::make_tensor_storage<StoragePolicy>("mixed_ops_tensors", 4);

  const auto ref_dot = a * b;
  const auto ref_contract = A.contract(libMesh::outer_product(a, b));
  const auto ref_det = host_leading_determinant(A, LIBMESH_DIM);
  const auto ref_right = A * c;
  const auto ref_left = A.left_multiply(c);
  const auto ref_mix = a + b;
  const auto ref_row0 = A.row(0);
  const auto ref_col0 = A.column(0);
  const auto ref_transpose = A.transpose();
  const auto ref_inverse = host_leading_inverse(A, LIBMESH_DIM);
  const auto ref_add = A + ref_transpose;
  const auto ref_scaled = 0.5 * A;
  const auto ref_trace = A.tr();

  ::Kokkos::parallel_for(
    1,
    KOKKOS_LAMBDA(int) {
      const auto a_ref = libMesh::Kokkos::make_vector_ref(d_a, 0);
      const auto A_ref = libMesh::Kokkos::make_tensor_ref(d_A, 0);

      const auto mix = a_ref + b;
      const auto right = A_ref * c;
      const auto left = A_ref.left_multiply(c);
      const auto row0 = A_ref.row(0);
      const auto col0 = A_ref.column(0);
      const auto transpose = A_ref.transpose();
      const auto inverse = A_ref.inverse();
      const auto add = A_ref + ref_transpose;
      const auto scaled = Real(0.5) * A_ref;
      const auto outer = libMesh::Kokkos::tensor_outer_product<oracle_tensor>(a_ref, b);

      d_scalars(0) = a_ref * b;
      d_scalars(1) = A_ref.contract(outer);
      d_scalars(2) = A_ref.det();
      d_scalars(3) = (A_ref == A) ? 1.0 : 0.0;
      d_scalars(4) = (A_ref != inverse) ? 1.0 : 0.0;
      d_scalars(5) = libMesh::Kokkos::vector_equal(row0, ref_row0) ? 1.0 : 0.0;
      d_scalars(6) = libMesh::Kokkos::vector_equal(col0, ref_col0) ? 1.0 : 0.0;
      d_scalars(7) = A_ref.tr();

      libMesh::Kokkos::store_vector(d_vectors, 0, right);
      libMesh::Kokkos::store_vector(d_vectors, 1, left);
      libMesh::Kokkos::store_vector(d_vectors, 2, mix);
      libMesh::Kokkos::store_vector(d_vectors, 3, row0);
      libMesh::Kokkos::store_vector(d_vectors, 4, col0);
      libMesh::Kokkos::store_tensor(d_tensors, 0, transpose);
      libMesh::Kokkos::store_tensor(d_tensors, 1, inverse);
      libMesh::Kokkos::store_tensor(d_tensors, 2, add);
      libMesh::Kokkos::store_tensor(d_tensors, 3, scaled);
    });
  ::Kokkos::fence();

  fail += libMeshTest::KokkosOracle::compare_device_scalars(
    d_scalars,
    std::vector<Real>{ref_dot, ref_contract, ref_det, 1.0, 1.0, 1.0, 1.0, ref_trace},
    tol);
  fail += libMeshTest::KokkosOracle::compare_device_vectors(
    d_vectors,
    std::vector<oracle_vector>{ref_right, ref_left, ref_mix, ref_row0, ref_col0},
    tol);
  fail += libMeshTest::KokkosOracle::compare_device_tensors(
    d_tensors,
    std::vector<oracle_tensor>{ref_transpose, ref_inverse, ref_add, ref_scaled},
    tol);

  return fail;
}

inline int
run_all_oracles()
{
  int total_fail = 0;

  const int dim_fail_left = test_dim_ops<libMesh::Kokkos::layout_left_storage_policy>();
  std::printf("[tensor_dim_kernel_oracle] [%s] %s  (%d failures)\n",
              libMesh::Kokkos::storage_policy_name<libMesh::Kokkos::layout_left_storage_policy>(),
              dim_fail_left ? "FAIL" : "PASS",
              dim_fail_left);
  total_fail += dim_fail_left;

  const int dim_fail_right = test_dim_ops<libMesh::Kokkos::layout_right_storage_policy>();
  std::printf("[tensor_dim_kernel_oracle] [%s] %s  (%d failures)\n",
              libMesh::Kokkos::storage_policy_name<libMesh::Kokkos::layout_right_storage_policy>(),
              dim_fail_right ? "FAIL" : "PASS",
              dim_fail_right);
  total_fail += dim_fail_right;

  const int tensor_fail_left = test_tensor_ops<libMesh::Kokkos::layout_left_storage_policy>();
  std::printf("[tensor_ops_kernel_oracle] [%s] %s  (%d failures)\n",
              libMesh::Kokkos::storage_policy_name<libMesh::Kokkos::layout_left_storage_policy>(),
              tensor_fail_left ? "FAIL" : "PASS",
              tensor_fail_left);
  total_fail += tensor_fail_left;

  const int tensor_fail_right = test_tensor_ops<libMesh::Kokkos::layout_right_storage_policy>();
  std::printf("[tensor_ops_kernel_oracle] [%s] %s  (%d failures)\n",
              libMesh::Kokkos::storage_policy_name<libMesh::Kokkos::layout_right_storage_policy>(),
              tensor_fail_right ? "FAIL" : "PASS",
              tensor_fail_right);
  total_fail += tensor_fail_right;

  const int host_fail = test_tensor_host_only_ops();
  std::printf("[tensor_host_ops_oracle] %s  (%d failures)\n",
              host_fail ? "FAIL" : "PASS",
              host_fail);
  total_fail += host_fail;

  const int foundation_fail_left =
    test_linalg_foundation_storage_roundtrip<libMesh::Kokkos::layout_left_storage_policy>();
  std::printf("[kokkos_linalg_foundation_oracle] [%s] %s  (%d failures)\n",
              libMesh::Kokkos::storage_policy_name<libMesh::Kokkos::layout_left_storage_policy>(),
              foundation_fail_left ? "FAIL" : "PASS",
              foundation_fail_left);
  total_fail += foundation_fail_left;

  const int foundation_fail_right =
    test_linalg_foundation_storage_roundtrip<libMesh::Kokkos::layout_right_storage_policy>();
  std::printf("[kokkos_linalg_foundation_oracle] [%s] %s  (%d failures)\n",
              libMesh::Kokkos::storage_policy_name<libMesh::Kokkos::layout_right_storage_policy>(),
              foundation_fail_right ? "FAIL" : "PASS",
              foundation_fail_right);
  total_fail += foundation_fail_right;

  const int mixed_fail_left = test_mixed_representation_ops<libMesh::Kokkos::layout_left_storage_policy>();
  std::printf("[kokkos_linalg_mixed_representation_oracle] [%s] %s  (%d failures)\n",
              libMesh::Kokkos::storage_policy_name<libMesh::Kokkos::layout_left_storage_policy>(),
              mixed_fail_left ? "FAIL" : "PASS",
              mixed_fail_left);
  total_fail += mixed_fail_left;

  const int mixed_fail_right = test_mixed_representation_ops<libMesh::Kokkos::layout_right_storage_policy>();
  std::printf("[kokkos_linalg_mixed_representation_oracle] [%s] %s  (%d failures)\n",
              libMesh::Kokkos::storage_policy_name<libMesh::Kokkos::layout_right_storage_policy>(),
              mixed_fail_right ? "FAIL" : "PASS",
              mixed_fail_right);
  total_fail += mixed_fail_right;

  return total_fail;
}

} // namespace KokkosTensorOracle
} // namespace libMeshTest

#endif
