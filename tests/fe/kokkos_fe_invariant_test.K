// Kokkos kernel regression tests for libMesh::Kokkos FE invariants and quadrature exactness.
// The test suite covers:
//   A. Partition of unity for all implemented LAGRANGE map topologies.
//   B. Zero-sum gradients for the same map topologies.
//   C. Nodal Kronecker-delta behavior at master nodes for the same map topologies.
//   D. Quadrature exactness sweeps for the canonical tensor-product and simplex
//      reference topologies using analytic monomial integrals on the reference
//      element.
//
// Returns 0 on success, non-zero on failure.

#include "gpu/kokkos_fe_evaluator.h"
#include "gpu/kokkos_quadrature.h"
#include "../quadrature/quadrature_exactness.h"

#include "libmesh/elem.h"
#include "libmesh/libmesh.h"
#include "libmesh/quadrature_gauss.h"

// Avoid conflicting complex operators between CUDA and PETSc
#define PETSC_SKIP_CXX_COMPLEX_FIX 1
#include <Kokkos_Core.hpp>
#undef __CUDACC_VER__

#include "kokkos_fe_oracle_test_utils.h"

#include <cmath>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

using libMesh::Kokkos::GaussQuadrature;
using libMesh::Kokkos::Real;
using libMesh::Kokkos::RealVector;
using kokkos_test_utils::build_qps;
using kokkos_test_utils::build_reference_elem;
using kokkos_test_utils::compare_device_values;
using kokkos_test_utils::dispatch_supported_lagrange_map_topology;
using kokkos_test_utils::upload_real;

static constexpr double invariant_tol = 1.0e-13;
static constexpr double exactness_tol = 2.0e-12;
static constexpr unsigned int quadrature_order = 4;

namespace
{

struct map_elem_info
{
  libMesh::ElemType topo;
  unsigned int      dim;
  unsigned int      n_dofs;
  const char *      name;
};

struct quadrature_exactness_case
{
  libMesh::ElemType topo;
  unsigned int      dim;
  unsigned int      max_order;
  const char *      name;
};

static const map_elem_info map_elems[] = {
  { libMesh::EDGE2,  1,  2, "EDGE2"  },
  { libMesh::EDGE3,  1,  3, "EDGE3"  },
  { libMesh::TRI3,   2,  3, "TRI3"   },
  { libMesh::TRI6,   2,  6, "TRI6"   },
  { libMesh::QUAD4,  2,  4, "QUAD4"  },
  { libMesh::QUAD8,  2,  8, "QUAD8"  },
  { libMesh::QUAD9,  2,  9, "QUAD9"  },
  { libMesh::TET4,   3,  4, "TET4"   },
  { libMesh::TET10,  3, 10, "TET10"  },
  { libMesh::HEX8,   3,  8, "HEX8"   },
  { libMesh::HEX20,  3, 20, "HEX20"  },
  { libMesh::HEX27,  3, 27, "HEX27"  }
};

// These cases sweep the full exactness range provided by the current Kokkos
// QGauss helper:
//   - tensor-product rules: through order 13 (7-point 1D Gauss-Legendre)
//   - simplex rules: through order 6 (highest explicit triangle/tet tables)
static const quadrature_exactness_case quadrature_cases[] = {
  { libMesh::EDGE2,  1, 13, "EDGE2"  },
  { libMesh::TRI3,   2,  6, "TRI3"   },
  { libMesh::QUAD4,  2, 13, "QUAD4"  },
  { libMesh::TET4,   3,  6, "TET4"   },
  { libMesh::HEX8,   3, 13, "HEX8"   }
};

} // anonymous namespace

LIBMESH_DEVICE_INLINE Real
int_pow(Real x, unsigned int p)
{
  Real result = 1.0;
  for (unsigned int i = 0; i < p; ++i)
    result *= x;
  return result;
}

template <libMesh::ElemType Topo>
static int
test_partition_of_unity_impl(const map_elem_info & info)
{
  std::vector<Real> xi_h, eta_h, zeta_h;
  const unsigned int nqp = build_qps(info.topo, info.dim, quadrature_order, xi_h, eta_h, zeta_h);

  auto d_xi = upload_real(xi_h, "unity_xi");
  auto d_eta = upload_real(eta_h, "unity_eta");
  auto d_zeta = upload_real(zeta_h, "unity_zeta");

  Kokkos::View<double *> d_sum(std::string("unity_sum"), nqp);

  const unsigned int n_dofs = info.n_dofs;
  Kokkos::parallel_for(
    nqp,
    KOKKOS_LAMBDA(int q) {
      Real sum = 0.0;
      for (unsigned int i = 0; i < n_dofs; ++i)
        sum += libMesh::Kokkos::map_shape<libMesh::LAGRANGE, Topo>(i, d_xi(q), d_eta(q), d_zeta(q));
      d_sum(q) = sum;
    });
  Kokkos::fence();

  std::vector<Real> ref_values(nqp, 1.0);
  return compare_device_values(d_sum, ref_values, invariant_tol);
}

struct partition_of_unity_dispatch
{
  explicit partition_of_unity_dispatch(const map_elem_info & in_info) : info(in_info) {}

  template <libMesh::ElemType Topo>
  int operator()() const
  {
    return test_partition_of_unity_impl<Topo>(info);
  }

  int unsupported_topology(libMesh::ElemType topo) const
  {
    std::printf("  unexpected unsupported partition topology: %s type=%d\n",
                info.name,
                static_cast<int>(topo));
    return 1;
  }

  const map_elem_info & info;
};

static int
test_partition_of_unity(const map_elem_info & info)
{
  const partition_of_unity_dispatch dispatch(info);
  return dispatch_supported_lagrange_map_topology(info.topo, dispatch);
}

template <libMesh::ElemType Topo>
static int
test_zero_sum_gradients_impl(const map_elem_info & info)
{
  std::vector<Real> xi_h, eta_h, zeta_h;
  const unsigned int nqp = build_qps(info.topo, info.dim, quadrature_order, xi_h, eta_h, zeta_h);

  auto d_xi = upload_real(xi_h, "gradsum_xi");
  auto d_eta = upload_real(eta_h, "gradsum_eta");
  auto d_zeta = upload_real(zeta_h, "gradsum_zeta");

  Kokkos::View<double *> d_sum(std::string("gradsum"), info.dim * nqp);

  const unsigned int n_dofs = info.n_dofs;
  const unsigned int dim = info.dim;
  Kokkos::parallel_for(
    nqp,
    KOKKOS_LAMBDA(int q) {
      RealVector sum = libMesh::Kokkos::zero_vector();
      for (unsigned int i = 0; i < n_dofs; ++i)
        sum += libMesh::Kokkos::grad_map_shape<libMesh::LAGRANGE, Topo>(
          i, d_xi(q), d_eta(q), d_zeta(q));

      for (unsigned int d = 0; d < dim; ++d)
        d_sum(dim * q + d) = sum(d);
    });
  Kokkos::fence();

  std::vector<Real> ref_values(info.dim * nqp, 0.0);
  return compare_device_values(d_sum, ref_values, invariant_tol);
}

struct zero_sum_gradients_dispatch
{
  explicit zero_sum_gradients_dispatch(const map_elem_info & in_info) : info(in_info) {}

  template <libMesh::ElemType Topo>
  int operator()() const
  {
    return test_zero_sum_gradients_impl<Topo>(info);
  }

  int unsupported_topology(libMesh::ElemType topo) const
  {
    std::printf("  unexpected unsupported zero-sum gradient topology: %s type=%d\n",
                info.name,
                static_cast<int>(topo));
    return 1;
  }

  const map_elem_info & info;
};

static int
test_zero_sum_gradients(const map_elem_info & info)
{
  const zero_sum_gradients_dispatch dispatch(info);
  return dispatch_supported_lagrange_map_topology(info.topo, dispatch);
}

template <libMesh::ElemType Topo>
static int
test_kronecker_delta_impl(const map_elem_info & info)
{
  std::vector<Real> xi_h(info.n_dofs), eta_h(info.n_dofs), zeta_h(info.n_dofs);
  for (unsigned int j = 0; j < info.n_dofs; ++j)
  {
    libMesh::Point p;
    libmesh_error_msg_if(!libMesh::try_reference_node(Topo, j, p),
                         "test_kronecker_delta_impl(): unsupported reference-node lookup");
    xi_h[j] = p(0);
    eta_h[j] = p(1);
    zeta_h[j] = p(2);
  }

  auto d_xi = upload_real(xi_h, "delta_xi");
  auto d_eta = upload_real(eta_h, "delta_eta");
  auto d_zeta = upload_real(zeta_h, "delta_zeta");

  Kokkos::View<double *> d_values(std::string("delta_values"), info.n_dofs * info.n_dofs);

  const unsigned int n_dofs = info.n_dofs;
  Kokkos::parallel_for(
    n_dofs * n_dofs,
    KOKKOS_LAMBDA(int idx) {
      const unsigned int i = static_cast<unsigned int>(idx) / n_dofs;
      const unsigned int j = static_cast<unsigned int>(idx) % n_dofs;
      d_values(idx) = libMesh::Kokkos::map_shape<libMesh::LAGRANGE, Topo>(
        i, d_xi(j), d_eta(j), d_zeta(j));
    });
  Kokkos::fence();

  std::vector<Real> ref_values(info.n_dofs * info.n_dofs, 0.0);
  for (unsigned int i = 0; i < info.n_dofs; ++i)
    ref_values[i * info.n_dofs + i] = 1.0;

  return compare_device_values(d_values, ref_values, invariant_tol);
}

struct kronecker_delta_dispatch
{
  explicit kronecker_delta_dispatch(const map_elem_info & in_info) : info(in_info) {}

  template <libMesh::ElemType Topo>
  int operator()() const
  {
    return test_kronecker_delta_impl<Topo>(info);
  }

  int unsupported_topology(libMesh::ElemType topo) const
  {
    std::printf("  unexpected unsupported kronecker topology: %s type=%d\n",
                info.name,
                static_cast<int>(topo));
    return 1;
  }

  const map_elem_info & info;
};

static int
test_kronecker_delta(const map_elem_info & info)
{
  const kronecker_delta_dispatch dispatch(info);
  return dispatch_supported_lagrange_map_topology(info.topo, dispatch);
}

static double
integrate_monomial_on_device(const quadrature_exactness_case & info,
                             unsigned int order,
                             unsigned int a,
                             unsigned int b,
                             unsigned int c)
{
  const unsigned int nqp = GaussQuadrature::n_points(info.topo, order);
  const libMesh::ElemType topo = info.topo;

  double integral = 0.0;
  Kokkos::parallel_reduce(
    nqp,
    KOKKOS_LAMBDA(int qp, double & local_sum) {
      const RealVector pt = GaussQuadrature::point(topo, order, qp);
      const Real weight = GaussQuadrature::weight(topo, order, qp);
      Real monomial = int_pow(pt(0), a);

#if LIBMESH_DIM > 1
      monomial *= int_pow(pt(1), b);
#else
      libmesh_assert_equal_to(b, 0);
#endif

#if LIBMESH_DIM > 2
      monomial *= int_pow(pt(2), c);
#else
      libmesh_assert_equal_to(c, 0);
#endif

      local_sum += static_cast<double>(weight) * static_cast<double>(monomial);
    },
    integral);

  return integral;
}

static int
test_quadrature_exactness(const quadrature_exactness_case & info)
{
  int fail = 0;

  for (unsigned int order = 0; order <= info.max_order; ++order)
    switch (info.dim)
    {
      case 1:
        for (unsigned int a = 0; a <= order; ++a)
        {
          const double actual = integrate_monomial_on_device(info, order, a, 0, 0);
          const double expected = quadrature_exactness::monomial_integral(info.topo, a, 0, 0);
          if (std::fabs(actual - expected) > exactness_tol)
            ++fail;
        }
        break;

      case 2:
        for (unsigned int a = 0; a <= order; ++a)
          for (unsigned int b = 0; a + b <= order; ++b)
          {
            const double actual = integrate_monomial_on_device(info, order, a, b, 0);
            const double expected = quadrature_exactness::monomial_integral(info.topo, a, b, 0);
            if (std::fabs(actual - expected) > exactness_tol)
              ++fail;
          }
        break;

      case 3:
        for (unsigned int a = 0; a <= order; ++a)
          for (unsigned int b = 0; a + b <= order; ++b)
            for (unsigned int c = 0; a + b + c <= order; ++c)
            {
              const double actual = integrate_monomial_on_device(info, order, a, b, c);
              const double expected = quadrature_exactness::monomial_integral(info.topo, a, b, c);
              if (std::fabs(actual - expected) > exactness_tol)
                ++fail;
            }
        break;

      default:
        ++fail;
        break;
    }

  return fail;
}

int
main(int argc, char ** argv)
{
  Kokkos::initialize(argc, argv);
  libMesh::LibMeshInit init(argc, argv);

  int total_fail = 0;

  for (const auto & info : map_elems)
  {
    {
      const int f = test_partition_of_unity(info);
      std::printf("[partition_of_unity]     [%s] %s  (%d failures)\n",
                  info.name, f ? "FAIL" : "PASS", f);
      total_fail += f;
    }
    {
      const int f = test_zero_sum_gradients(info);
      std::printf("[zero_sum_gradients]     [%s] %s  (%d failures)\n",
                  info.name, f ? "FAIL" : "PASS", f);
      total_fail += f;
    }
    {
      const int f = test_kronecker_delta(info);
      std::printf("[kronecker_delta]        [%s] %s  (%d failures)\n",
                  info.name, f ? "FAIL" : "PASS", f);
      total_fail += f;
    }
  }

  for (const auto & info : quadrature_cases)
  {
    const int f = test_quadrature_exactness(info);
    std::printf("[quadrature_exactness]   [%s] %s  (%d failures)\n",
                info.name, f ? "FAIL" : "PASS", f);
    total_fail += f;
  }

  Kokkos::finalize();

  if (total_fail == 0)
    std::printf("ALL TESTS PASSED\n");
  else
    std::printf("%d TEST(S) FAILED\n", total_fail);

  return total_fail ? 1 : 0;
}
