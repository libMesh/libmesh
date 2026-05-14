// Kokkos kernel oracle tests for libMesh tensor algebra used by Kokkos kernels.
//
// Compares device-side tensor/vector
// operations against host libMesh TypeTensor / TypeVector oracles built from
// direct hand-written fixtures.
//
// Returns 0 on success, non-zero on failure.

#include "libmesh/libmesh_config.h"

#include "libmesh/libmesh.h"
#include "libmesh/point.h"
#include "libmesh/tensor_value.h"
#include "libmesh/type_n_tensor.h"
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

using kernel_vector = libMesh::TypeVector<Real>;
using kernel_tensor = libMesh::TypeTensor<Real>;

namespace
{

using vector_storage_view = Kokkos::View<Real * [LIBMESH_DIM]>;
using tensor_storage_view = Kokkos::View<Real * [LIBMESH_DIM][LIBMESH_DIM]>;

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

template <typename T>
inline Kokkos::View<T *>
upload_objects(const std::vector<T> & values, const char * label)
{
  Kokkos::View<T *> d(std::string(label), values.size());
  auto h = Kokkos::create_mirror_view(d);
  for (std::size_t i = 0; i < values.size(); ++i)
    h(i) = values[i];
  Kokkos::deep_copy(d, h);
  return d;
}

inline int
compare_device_vectors(const vector_storage_view & d_values,
                       const std::vector<kernel_vector> & ref_values,
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

inline int
compare_device_tensors(const tensor_storage_view & d_values,
                       const std::vector<kernel_tensor> & ref_values,
                       const double local_tol = tol)
{
  auto h_values = Kokkos::create_mirror_view(d_values);
  Kokkos::deep_copy(h_values, d_values);

  int fail = 0;
  for (std::size_t i = 0; i < ref_values.size(); ++i)
    for (unsigned int row = 0; row < LIBMESH_DIM; ++row)
      for (unsigned int col = 0; col < LIBMESH_DIM; ++col)
        if (std::fabs(h_values(i, row, col) - ref_values[i](row, col)) > local_tol)
          ++fail;

  return fail;
}

inline vector_storage_view
upload_vectors(const std::vector<kernel_vector> & values, const char * label)
{
  vector_storage_view d(std::string(label), values.size());
  auto h = Kokkos::create_mirror_view(d);

  for (std::size_t i = 0; i < values.size(); ++i)
    for (unsigned int component = 0; component < LIBMESH_DIM; ++component)
      h(i, component) = values[i](component);

  Kokkos::deep_copy(d, h);
  return d;
}

inline tensor_storage_view
upload_tensors(const std::vector<kernel_tensor> & values, const char * label)
{
  tensor_storage_view d(std::string(label), values.size());
  auto h = Kokkos::create_mirror_view(d);

  for (std::size_t i = 0; i < values.size(); ++i)
    for (unsigned int row = 0; row < LIBMESH_DIM; ++row)
      for (unsigned int col = 0; col < LIBMESH_DIM; ++col)
        h(i, row, col) = values[i](row, col);

  Kokkos::deep_copy(d, h);
  return d;
}

inline libMesh::TypeVector<Real>
make_host_vector(const Real x, const Real y = 0, const Real z = 0)
{
  libMesh::TypeVector<Real> v;
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

inline libMesh::TypeTensor<Real>
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
  libMesh::TypeTensor<Real> T;
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
  libMesh::TypeTensor<Real> J;
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

inline libMesh::TypeTensor<Real>
host_identity(const unsigned int dim)
{
  libMesh::TypeTensor<Real> I;
  I.zero();
  for (unsigned int i = 0; i < dim; ++i)
    I(i, i) = Real(1);
  return I;
}

inline Real
host_leading_determinant(const libMesh::TypeTensor<Real> & J, const unsigned int dim)
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

inline libMesh::TypeTensor<Real>
host_leading_inverse(const libMesh::TypeTensor<Real> & J, const unsigned int dim)
{
  libMesh::TypeTensor<Real> inv;
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
  return libMesh::TypeTensor<Real>(J.inverse());
#else
  return inv;
#endif
}

LIBMESH_DEVICE_INLINE
kernel_vector
zero_kernel_vector()
{
  kernel_vector v;
  v.zero();
  return v;
}

LIBMESH_DEVICE_INLINE
kernel_tensor
zero_kernel_tensor()
{
  kernel_tensor T;
  T.zero();
  return T;
}

LIBMESH_DEVICE_INLINE
kernel_tensor
kernel_leading_identity(const unsigned int dim)
{
  kernel_tensor I = zero_kernel_tensor();
  for (unsigned int i = 0; i < dim; ++i)
    I(i, i) = Real(1);
  return I;
}

LIBMESH_DEVICE_INLINE
Real
kernel_leading_determinant(const kernel_tensor & J, const unsigned int dim)
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

LIBMESH_DEVICE_INLINE
kernel_tensor
kernel_leading_inverse(const kernel_tensor & J, const unsigned int dim)
{
  kernel_tensor inv = zero_kernel_tensor();

  if (dim == 1)
  {
    inv(0, 0) = Real(1) / J(0, 0);
    return inv;
  }

  if (dim == 2)
  {
    const Real det = kernel_leading_determinant(J, dim);
    inv(0, 0) =  J(1, 1) / det;
    inv(0, 1) = -J(0, 1) / det;
    inv(1, 0) = -J(1, 0) / det;
    inv(1, 1) =  J(0, 0) / det;
    return inv;
  }

#if LIBMESH_DIM > 2
  return J.inverse();
#else
  return inv;
#endif
}

static int
test_dim_ops()
{
  const unsigned int ncases = sizeof(dim_cases) / sizeof(dim_cases[0]);

  std::vector<kernel_tensor> J_values(ncases);
  std::vector<unsigned int> dims(ncases);
  std::vector<Real> ref_det(ncases);
  std::vector<kernel_tensor> ref_inv(ncases);
  std::vector<kernel_tensor> ref_I(ncases);
  std::vector<kernel_tensor> ref_prod_left(ncases);
  std::vector<kernel_tensor> ref_prod_right(ncases);

  for (unsigned int c = 0; c < ncases; ++c)
  {
    const auto & info = dim_cases[c];
    J_values[c] = info.J;
    dims[c] = info.dim;

    ref_det[c] = host_leading_determinant(info.J, info.dim);
    ref_inv[c] = host_leading_inverse(info.J, info.dim);
    ref_I[c] = host_identity(info.dim);
    ref_prod_left[c] = info.J * ref_inv[c];
    ref_prod_right[c] = ref_inv[c] * info.J;
  }

  auto d_J = upload_tensors(J_values, "tensor_dim_ops_J");
  auto d_dims = upload_objects(dims, "tensor_dim_ops_dim");
  Kokkos::View<double *> d_det("tensor_dim_ops_det", ncases);
  tensor_storage_view d_inv("tensor_dim_ops_inv", ncases);
  tensor_storage_view d_I("tensor_dim_ops_I", ncases);
  tensor_storage_view d_prod_left("tensor_dim_ops_prod_left", ncases);
  tensor_storage_view d_prod_right("tensor_dim_ops_prod_right", ncases);

  Kokkos::parallel_for(
    static_cast<int>(ncases),
    KOKKOS_LAMBDA(int c) {
      const auto J = libMesh::Kokkos::load_tensor<kernel_tensor>(d_J, c);
      const unsigned int dim = d_dims(c);
      const Real det = kernel_leading_determinant(J, dim);
      const auto inv = kernel_leading_inverse(J, dim);
      const auto I = kernel_leading_identity(dim);
      const auto prod_left = J * inv;
      const auto prod_right = inv * J;

      d_det(c) = det;
      libMesh::Kokkos::store_tensor(d_inv, c, inv);
      libMesh::Kokkos::store_tensor(d_I, c, I);
      libMesh::Kokkos::store_tensor(d_prod_left, c, prod_left);
      libMesh::Kokkos::store_tensor(d_prod_right, c, prod_right);
    });
  Kokkos::fence();

  return compare_device_scalars(d_det, ref_det) +
         compare_device_tensors(d_inv, ref_inv) +
         compare_device_tensors(d_I, ref_I) +
         compare_device_tensors(d_prod_left, ref_prod_left) +
         compare_device_tensors(d_prod_right, ref_prod_right);
}

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
  const auto zero = zero_kernel_tensor();

  std::vector<kernel_tensor> ref_outer(1, outer);
  std::vector<kernel_tensor> ref_transpose(1, transpose);
  std::vector<kernel_tensor> ref_mix(1, mix);
  std::vector<kernel_vector> ref_rows(LIBMESH_DIM);
  std::vector<kernel_vector> ref_columns(LIBMESH_DIM);
  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
  {
    ref_rows[i] = A.row(i);
    ref_columns[i] = A.column(i);
  }
  std::vector<kernel_vector> ref_right(1, right);
  std::vector<kernel_vector> ref_left(1, left);
  std::vector<Real> ref_scalars = {contract, norm, zero.is_zero() ? 1.0 : 0.0, A.is_zero() ? 1.0 : 0.0};

  auto d_A = upload_tensors(std::vector<kernel_tensor>{A}, "tensor_ops_A");
  auto d_a = upload_vectors(std::vector<kernel_vector>{a}, "tensor_ops_a");
  auto d_b = upload_vectors(std::vector<kernel_vector>{b}, "tensor_ops_b");
  auto d_c = upload_vectors(std::vector<kernel_vector>{c}, "tensor_ops_c");
  tensor_storage_view d_outer("tensor_ops_outer", 1);
  tensor_storage_view d_transpose("tensor_ops_transpose", 1);
  tensor_storage_view d_mix("tensor_ops_mix", 1);
  vector_storage_view d_rows("tensor_ops_rows", LIBMESH_DIM);
  vector_storage_view d_columns("tensor_ops_columns", LIBMESH_DIM);
  vector_storage_view d_right("tensor_ops_right", 1);
  vector_storage_view d_left("tensor_ops_left", 1);
  Kokkos::View<double *> d_scalars("tensor_ops_scalars", 4);

  Kokkos::parallel_for(
    1,
    KOKKOS_LAMBDA(int) {
      const auto A_d = libMesh::Kokkos::load_tensor<kernel_tensor>(d_A, 0);
      const auto a_d = libMesh::Kokkos::load_vector<kernel_vector>(d_a, 0);
      const auto b_d = libMesh::Kokkos::load_vector<kernel_vector>(d_b, 0);
      const auto c_d = libMesh::Kokkos::load_vector<kernel_vector>(d_c, 0);
      const auto outer_d = libMesh::outer_product(a_d, b_d);
      const auto transpose_d = A_d.transpose();
      const auto mix_d = 1.5 * A_d - 0.25 * outer_d;
      const auto right_d = A_d * c_d;
      const auto left_d = c_d * A_d;
      const Real contract_d = A_d.contract(outer_d);
      const Real norm_d = A_d.norm();
      const bool zero_is_zero_d = zero_kernel_tensor().is_zero();
      const bool A_is_zero_d = A_d.is_zero();

      for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
      {
        libMesh::Kokkos::store_vector(d_rows, i, A_d.row(i));
        libMesh::Kokkos::store_vector(d_columns, i, A_d.column(i));
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
  Kokkos::fence();

  return compare_device_tensors(d_outer, ref_outer) +
         compare_device_tensors(d_transpose, ref_transpose) +
         compare_device_tensors(d_mix, ref_mix) +
         compare_device_vectors(d_rows, ref_rows) +
         compare_device_vectors(d_columns, ref_columns) +
         compare_device_vectors(d_right, ref_right) +
         compare_device_vectors(d_left, ref_left) +
         compare_device_scalars(d_scalars, ref_scalars);
}

static int
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

} // anonymous namespace

int
main(int argc, char ** argv)
{
  Kokkos::initialize(argc, argv);
  libMesh::LibMeshInit init(argc, argv);

  int total_fail = 0;

  const int dim_fail = test_dim_ops();
  std::printf("[tensor_dim_kernel_oracle] %s  (%d failures)\n",
              dim_fail ? "FAIL" : "PASS",
              dim_fail);
  total_fail += dim_fail;

  const int tensor_fail = test_tensor_ops();
  std::printf("[tensor_ops_kernel_oracle] %s  (%d failures)\n",
              tensor_fail ? "FAIL" : "PASS",
              tensor_fail);
  total_fail += tensor_fail;

  const int host_fail = test_tensor_host_only_ops();
  std::printf("[tensor_host_ops_oracle] %s  (%d failures)\n",
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
