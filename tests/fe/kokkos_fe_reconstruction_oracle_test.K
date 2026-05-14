// GPU kernel tests for libMesh::Kokkos FE reconstruction on physical elements.
//
// Standalone executable (no CppUnit). Uses libMesh::LibMeshInit so that
// FEBase provides the host physical-space oracle values and gradients.
//
// The test suite covers:
//   A. Reconstruction of solution values on distorted physical elements for
//      every exact LAGRANGE key currently supported by the Kokkos evaluator.
//   B. Reconstruction of physical-space gradients on the same elements by
//      pulling reference-space gradients through the element Jacobian.
//   C. The same value/gradient reconstruction parity for representative exact
//      MONOMIAL keys across all supported dimensions and orders.
//
// Returns 0 on success, non-zero on failure.

#include "libmesh/libmesh_config.h"

#include "gpu/kokkos_fe_evaluator.h"
#include "gpu/kokkos_fe_map.h"
#include "gpu/kokkos_fe_types.h"
#include "gpu/kokkos_scalar_types.h"

#include "libmesh/elem.h"
#include "libmesh/fe_base.h"
#include "libmesh/libmesh.h"
#include "libmesh/quadrature_gauss.h"

// Avoid conflicting complex operators between CUDA and PETSc
#define PETSC_SKIP_CXX_COMPLEX_FIX 1
#include <Kokkos_Core.hpp>
#undef __CUDACC_VER__

#include "kokkos_fe_oracle_test_utils.h"

#include <cstdio>
#include <string>
#include <vector>

using libMesh::Kokkos::FEShapeKey;
using libMesh::Kokkos::Real;
using libMesh::Kokkos::RealTensor;
using libMesh::Kokkos::RealVector;
using kokkos_test_utils::build_flat_reference_fixture;
using kokkos_test_utils::compare_device_values;
using kokkos_test_utils::dispatch_supported_shape_key;
using kokkos_test_utils::dispatch_supported_shape_key_with_lagrange_map;
using kokkos_test_utils::element_fixture;
using kokkos_test_utils::grad_shape_for_key;
using kokkos_test_utils::shape_for_key;
using kokkos_test_utils::upload_point_coordinates;
using kokkos_test_utils::upload_real;
using kokkos_test_utils::vector_component;

static constexpr double value_tol = 5.0e-13;
static constexpr double grad_tol = 5.0e-12;
static constexpr unsigned int quad_order = 4;

namespace
{

struct reconstruction_case
{
  FEShapeKey   key;
  unsigned int dim;
  unsigned int n_dofs;
  const char * name;
};

constexpr unsigned int
monomial_n_dofs_for_dim(unsigned int dim, libMesh::Order order)
{
  const unsigned int p = static_cast<unsigned int>(order);

  switch (dim)
  {
    case 1:
      return p + 1;
    case 2:
      return (p + 1) * (p + 2) / 2;
    case 3:
      return (p + 1) * (p + 2) * (p + 3) / 6;
    default:
      return 0;
  }
}

static const reconstruction_case lagrange_cases[] = {
  { { libMesh::LAGRANGE, libMesh::EDGE2,  libMesh::FIRST  }, 1,  2, "LAGRANGE/EDGE2/FIRST"   },
  { { libMesh::LAGRANGE, libMesh::EDGE3,  libMesh::FIRST  }, 1,  2, "LAGRANGE/EDGE3/FIRST"   },
  { { libMesh::LAGRANGE, libMesh::EDGE3,  libMesh::SECOND }, 1,  3, "LAGRANGE/EDGE3/SECOND"  },

  { { libMesh::LAGRANGE, libMesh::TRI3,   libMesh::FIRST  }, 2,  3, "LAGRANGE/TRI3/FIRST"    },
  { { libMesh::LAGRANGE, libMesh::TRI6,   libMesh::FIRST  }, 2,  3, "LAGRANGE/TRI6/FIRST"    },
  { { libMesh::LAGRANGE, libMesh::TRI6,   libMesh::SECOND }, 2,  6, "LAGRANGE/TRI6/SECOND"   },

  { { libMesh::LAGRANGE, libMesh::QUAD4,  libMesh::FIRST  }, 2,  4, "LAGRANGE/QUAD4/FIRST"   },
  { { libMesh::LAGRANGE, libMesh::QUAD8,  libMesh::FIRST  }, 2,  4, "LAGRANGE/QUAD8/FIRST"   },
  { { libMesh::LAGRANGE, libMesh::QUAD8,  libMesh::SECOND }, 2,  8, "LAGRANGE/QUAD8/SECOND"  },
  { { libMesh::LAGRANGE, libMesh::QUAD9,  libMesh::FIRST  }, 2,  4, "LAGRANGE/QUAD9/FIRST"   },
  { { libMesh::LAGRANGE, libMesh::QUAD9,  libMesh::SECOND }, 2,  9, "LAGRANGE/QUAD9/SECOND"  },

  { { libMesh::LAGRANGE, libMesh::TET4,   libMesh::FIRST  }, 3,  4, "LAGRANGE/TET4/FIRST"    },
  { { libMesh::LAGRANGE, libMesh::TET10,  libMesh::FIRST  }, 3,  4, "LAGRANGE/TET10/FIRST"   },
  { { libMesh::LAGRANGE, libMesh::TET10,  libMesh::SECOND }, 3, 10, "LAGRANGE/TET10/SECOND"  },

  { { libMesh::LAGRANGE, libMesh::HEX8,   libMesh::FIRST  }, 3,  8, "LAGRANGE/HEX8/FIRST"    },
  { { libMesh::LAGRANGE, libMesh::HEX20,  libMesh::FIRST  }, 3,  8, "LAGRANGE/HEX20/FIRST"   },
  { { libMesh::LAGRANGE, libMesh::HEX20,  libMesh::SECOND }, 3, 20, "LAGRANGE/HEX20/SECOND"  },
  { { libMesh::LAGRANGE, libMesh::HEX27,  libMesh::FIRST  }, 3,  8, "LAGRANGE/HEX27/FIRST"   },
  { { libMesh::LAGRANGE, libMesh::HEX27,  libMesh::SECOND }, 3, 27, "LAGRANGE/HEX27/SECOND"  }
};

static const reconstruction_case monomial_cases[] = {
  { { libMesh::MONOMIAL, libMesh::EDGE2,     libMesh::CONSTANT }, 1, monomial_n_dofs_for_dim(1, libMesh::CONSTANT), "MONOMIAL/EDGE2/CONSTANT"  },
  { { libMesh::MONOMIAL, libMesh::EDGE2,     libMesh::FIRST    }, 1, monomial_n_dofs_for_dim(1, libMesh::FIRST),    "MONOMIAL/EDGE2/FIRST"     },
  { { libMesh::MONOMIAL, libMesh::EDGE3,     libMesh::SECOND   }, 1, monomial_n_dofs_for_dim(1, libMesh::SECOND),   "MONOMIAL/EDGE3/SECOND"    },
  { { libMesh::MONOMIAL, libMesh::EDGE3,     libMesh::FOURTH   }, 1, monomial_n_dofs_for_dim(1, libMesh::FOURTH),   "MONOMIAL/EDGE3/FOURTH"    },
  { { libMesh::MONOMIAL, libMesh::EDGE3,     libMesh::FIFTH    }, 1, monomial_n_dofs_for_dim(1, libMesh::FIFTH),    "MONOMIAL/EDGE3/FIFTH"     },

  { { libMesh::MONOMIAL, libMesh::TRI3,      libMesh::CONSTANT }, 2, monomial_n_dofs_for_dim(2, libMesh::CONSTANT), "MONOMIAL/TRI3/CONSTANT"   },
  { { libMesh::MONOMIAL, libMesh::TRI3,      libMesh::FIRST    }, 2, monomial_n_dofs_for_dim(2, libMesh::FIRST),    "MONOMIAL/TRI3/FIRST"      },
  { { libMesh::MONOMIAL, libMesh::TRI6,      libMesh::SECOND   }, 2, monomial_n_dofs_for_dim(2, libMesh::SECOND),   "MONOMIAL/TRI6/SECOND"     },
  { { libMesh::MONOMIAL, libMesh::QUAD4,     libMesh::FIRST    }, 2, monomial_n_dofs_for_dim(2, libMesh::FIRST),    "MONOMIAL/QUAD4/FIRST"     },
  { { libMesh::MONOMIAL, libMesh::QUAD9,     libMesh::SECOND   }, 2, monomial_n_dofs_for_dim(2, libMesh::SECOND),   "MONOMIAL/QUAD9/SECOND"    },
  { { libMesh::MONOMIAL, libMesh::QUAD9,     libMesh::FIFTH    }, 2, monomial_n_dofs_for_dim(2, libMesh::FIFTH),    "MONOMIAL/QUAD9/FIFTH"     },

  { { libMesh::MONOMIAL, libMesh::TET4,      libMesh::CONSTANT }, 3, monomial_n_dofs_for_dim(3, libMesh::CONSTANT), "MONOMIAL/TET4/CONSTANT"   },
  { { libMesh::MONOMIAL, libMesh::TET4,      libMesh::FIRST    }, 3, monomial_n_dofs_for_dim(3, libMesh::FIRST),    "MONOMIAL/TET4/FIRST"      },
  { { libMesh::MONOMIAL, libMesh::TET10,     libMesh::SECOND   }, 3, monomial_n_dofs_for_dim(3, libMesh::SECOND),   "MONOMIAL/TET10/SECOND"    },
  { { libMesh::MONOMIAL, libMesh::HEX8,      libMesh::FIRST    }, 3, monomial_n_dofs_for_dim(3, libMesh::FIRST),    "MONOMIAL/HEX8/FIRST"      },
  { { libMesh::MONOMIAL, libMesh::HEX27,     libMesh::SECOND   }, 3, monomial_n_dofs_for_dim(3, libMesh::SECOND),   "MONOMIAL/HEX27/SECOND"    },
  { { libMesh::MONOMIAL, libMesh::HEX27,     libMesh::FIFTH    }, 3, monomial_n_dofs_for_dim(3, libMesh::FIFTH),    "MONOMIAL/HEX27/FIFTH"     }
};

} // anonymous namespace

static std::vector<Real>
build_coefficients(const reconstruction_case & info)
{
  std::vector<Real> coeffs(info.n_dofs);
  const Real family_bias = (info.key.family == libMesh::MONOMIAL) ? Real(0.19) : Real(0.07);
  const Real order_bias = Real(static_cast<unsigned int>(info.key.order) + 1u) * Real(0.013);

  for (unsigned int i = 0; i < info.n_dofs; ++i)
  {
    const Real sign = (i % 2) ? Real(-1.0) : Real(1.0);
    coeffs[i] = sign * (Real(0.17) + Real(0.041) * Real(i + 1u) + family_bias + order_bias);
  }

  return coeffs;
}

template <libMesh::FEFamily Family, libMesh::ElemType ExactTopo, libMesh::Order ExactOrder>
static int
test_reconstruction_case_impl(const reconstruction_case & info)
{
  constexpr unsigned int max_geom_nodes = 27;

  element_fixture fixture = build_flat_reference_fixture(ExactTopo);
  const unsigned int n_geom_nodes = fixture.elem->n_nodes();

  const libMesh::FEType fe_type(info.key.order, info.key.family);
  auto fe = libMesh::FEBase::build(info.dim, fe_type);

  libMesh::QGauss qr(info.dim, static_cast<libMesh::Order>(quad_order));
  qr.allow_rules_with_negative_weights = true;
  qr.init(info.key.elem_type);

  fe->attach_quadrature_rule(&qr);
  fe->get_phi();
  fe->get_dphi();
  fe->reinit(fixture.elem.get());

  const auto & phi = fe->get_phi();
  const auto & dphi = fe->get_dphi();
  const unsigned int nqp = qr.n_points();

  const std::vector<Real> coeffs = build_coefficients(info);
  std::vector<Real> ref_u(nqp, 0.0);
  std::vector<Real> ref_gx(nqp, 0.0);
  std::vector<Real> ref_gy(nqp, 0.0);
  std::vector<Real> ref_gz(nqp, 0.0);
  std::vector<Real> xi_h(nqp), eta_h(nqp), zeta_h(nqp);

  for (unsigned int q = 0; q < nqp; ++q)
  {
    for (unsigned int i = 0; i < info.n_dofs; ++i)
    {
      ref_u[q] += phi[i][q] * coeffs[i];
      ref_gx[q] += dphi[i][q](0) * coeffs[i];
      ref_gy[q] += dphi[i][q](1) * coeffs[i];
      ref_gz[q] += dphi[i][q](2) * coeffs[i];
    }

    xi_h[q] = qr.qp(q)(0);
    eta_h[q] = (info.dim >= 2) ? qr.qp(q)(1) : Real(0);
    zeta_h[q] = (info.dim >= 3) ? qr.qp(q)(2) : Real(0);
  }

  auto d_coords = upload_point_coordinates(*fixture.elem, "recon_coords");
  auto d_coeffs = upload_real(coeffs, "recon_coeffs");
  auto d_xi = upload_real(xi_h, "recon_xi");
  auto d_eta = upload_real(eta_h, "recon_eta");
  auto d_zeta = upload_real(zeta_h, "recon_zeta");

  Kokkos::View<double *> d_u(std::string("recon_u"), nqp);
  Kokkos::View<double *> d_gx(std::string("recon_gx"), nqp);
  Kokkos::View<double *> d_gy(std::string("recon_gy"), nqp);
  Kokkos::View<double *> d_gz(std::string("recon_gz"), nqp);

  const unsigned int dim = info.dim;
  const unsigned int n_dofs = info.n_dofs;
  const unsigned int n_geom_nodes_ = n_geom_nodes;

  Kokkos::parallel_for(
    nqp,
    KOKKOS_LAMBDA(int q) {
      RealVector geom_nodes[max_geom_nodes];
      for (unsigned int i = 0; i < n_geom_nodes_; ++i)
        geom_nodes[i] = libMesh::Kokkos::make_vector(
          d_coords(3 * i + 0), d_coords(3 * i + 1), d_coords(3 * i + 2));

      const RealTensor J = libMesh::Kokkos::jacobian<libMesh::LAGRANGE, ExactTopo>(
        geom_nodes, n_geom_nodes_, d_xi(q), d_eta(q), d_zeta(q));

      Real u = 0.0;
      RealVector grad_ref_sum = libMesh::Kokkos::zero_vector();
      for (unsigned int i = 0; i < n_dofs; ++i)
      {
        const Real coeff = d_coeffs(i);
        u += coeff * shape_for_key<Family, ExactTopo, ExactOrder>(i, d_xi(q), d_eta(q), d_zeta(q));
        grad_ref_sum +=
          coeff * grad_shape_for_key<Family, ExactTopo, ExactOrder>(i, d_xi(q), d_eta(q), d_zeta(q));
      }

      const RealTensor invJ = libMesh::Kokkos::leading_inverse(J, dim);
      const RealVector grad_phys = invJ * grad_ref_sum;

      d_u(q) = u;
      d_gx(q) = vector_component(grad_phys, 0);
      d_gy(q) = vector_component(grad_phys, 1);
      d_gz(q) = vector_component(grad_phys, 2);
    });
  Kokkos::fence();

  int fail = 0;
  fail += compare_device_values(d_u, ref_u, value_tol);
  fail += compare_device_values(d_gx, ref_gx, grad_tol);
  fail += compare_device_values(d_gy, ref_gy, grad_tol);
  fail += compare_device_values(d_gz, ref_gz, grad_tol);

  if (fail)
    std::printf("  reconstruction mismatch: %s (%d failures)\n", info.name, fail);

  return fail;
}

struct reconstruction_dispatch
{
  explicit reconstruction_dispatch(const reconstruction_case & in_info) : info(in_info) {}

  template <libMesh::FEFamily Family, libMesh::ElemType ExactTopo, libMesh::Order ExactOrder>
  int operator()() const
  {
    return test_reconstruction_case_impl<Family, ExactTopo, ExactOrder>(info);
  }

  int unsupported_key(FEShapeKey key) const
  {
    std::printf("  unexpected unsupported reconstruction key: %s family=%d elem_type=%d order=%d\n",
                info.name,
                static_cast<int>(key.family),
                static_cast<int>(key.elem_type),
                static_cast<int>(key.order));
    return 1;
  }

  const reconstruction_case & info;
};

static int
test_reconstruction_case(const reconstruction_case & info)
{
  const reconstruction_dispatch dispatch(info);
  return dispatch_supported_shape_key_with_lagrange_map(info.key, dispatch);
}

int
main(int argc, char ** argv)
{
  Kokkos::initialize(argc, argv);
  libMesh::LibMeshInit init(argc, argv);

  int total_fail = 0;

  for (const auto & info : lagrange_cases)
  {
    const int f = test_reconstruction_case(info);
    std::printf("[reconstruction_lagrange] [%s] %s  (%d failures)\n",
                info.name, f ? "FAIL" : "PASS", f);
    total_fail += f;
  }

  for (const auto & info : monomial_cases)
  {
    const int f = test_reconstruction_case(info);
    std::printf("[reconstruction_monomial] [%s] %s  (%d failures)\n",
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
