// GPU kernel tests for libMesh::Kokkos side-trace FE parity.
//
// Standalone executable (no CppUnit). Uses libMesh::LibMeshInit so that
// FEBase side reinit provides host side-trace oracle values.
//
// The test suite covers:
//   A. Side-restricted shape values for supported exact LAGRANGE keys on
//      distorted physical elements.
//   B. Tangential physical gradients on those same sides, using host FEBase
//      side traces as the oracle.
//   C. The same side value and tangential-gradient parity for representative
//      exact MONOMIAL keys whose parent geometry topology is supported by the
//      current Kokkos map implementation.
//
// Returns 0 on success, non-zero on failure.

#include "libmesh/libmesh_config.h"

#include "gpu/kokkos_fe_base.h"
#include "gpu/kokkos_fe_evaluator.h"
#include "gpu/kokkos_fe_face_map.h"
#include "gpu/kokkos_fe_map.h"
#include "gpu/kokkos_fe_types.h"

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
using kokkos_test_utils::grad_shape_for_key;
using kokkos_test_utils::shape_for_key;
using kokkos_test_utils::upload_point_coordinates;
using kokkos_test_utils::upload_real;
using kokkos_test_utils::vector_component;

static constexpr double value_tol = 5.0e-13;
static constexpr double grad_tol = 5.0e-12;

namespace
{

struct side_trace_case
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

// Restrict to parent geometries currently handled by the Kokkos map path.
static const side_trace_case lagrange_cases[] = {
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

static const side_trace_case monomial_cases[] = {
  { { libMesh::MONOMIAL, libMesh::TRI3,   libMesh::CONSTANT }, 2, monomial_n_dofs_for_dim(2, libMesh::CONSTANT), "MONOMIAL/TRI3/CONSTANT" },
  { { libMesh::MONOMIAL, libMesh::TRI3,   libMesh::FIRST    }, 2, monomial_n_dofs_for_dim(2, libMesh::FIRST),    "MONOMIAL/TRI3/FIRST"    },
  { { libMesh::MONOMIAL, libMesh::TRI6,   libMesh::SECOND   }, 2, monomial_n_dofs_for_dim(2, libMesh::SECOND),   "MONOMIAL/TRI6/SECOND"   },
  { { libMesh::MONOMIAL, libMesh::QUAD4,  libMesh::FIRST    }, 2, monomial_n_dofs_for_dim(2, libMesh::FIRST),    "MONOMIAL/QUAD4/FIRST"   },
  { { libMesh::MONOMIAL, libMesh::QUAD9,  libMesh::SECOND   }, 2, monomial_n_dofs_for_dim(2, libMesh::SECOND),   "MONOMIAL/QUAD9/SECOND"  },
  { { libMesh::MONOMIAL, libMesh::QUAD9,  libMesh::FIFTH    }, 2, monomial_n_dofs_for_dim(2, libMesh::FIFTH),    "MONOMIAL/QUAD9/FIFTH"   },

  { { libMesh::MONOMIAL, libMesh::TET4,   libMesh::CONSTANT }, 3, monomial_n_dofs_for_dim(3, libMesh::CONSTANT), "MONOMIAL/TET4/CONSTANT" },
  { { libMesh::MONOMIAL, libMesh::TET4,   libMesh::FIRST    }, 3, monomial_n_dofs_for_dim(3, libMesh::FIRST),    "MONOMIAL/TET4/FIRST"    },
  { { libMesh::MONOMIAL, libMesh::TET10,  libMesh::SECOND   }, 3, monomial_n_dofs_for_dim(3, libMesh::SECOND),   "MONOMIAL/TET10/SECOND"  },
  { { libMesh::MONOMIAL, libMesh::HEX8,   libMesh::FIRST    }, 3, monomial_n_dofs_for_dim(3, libMesh::FIRST),    "MONOMIAL/HEX8/FIRST"    },
  { { libMesh::MONOMIAL, libMesh::HEX27,  libMesh::SECOND   }, 3, monomial_n_dofs_for_dim(3, libMesh::SECOND),   "MONOMIAL/HEX27/SECOND"  },
  { { libMesh::MONOMIAL, libMesh::HEX27,  libMesh::FIFTH    }, 3, monomial_n_dofs_for_dim(3, libMesh::FIFTH),    "MONOMIAL/HEX27/FIFTH"   }
};

} // anonymous namespace

LIBMESH_DEVICE_INLINE RealVector
tangential_component(const RealVector & v, const RealVector & normal)
{
  return v - ((v * normal) * normal);
}

template <libMesh::FEFamily Family, libMesh::ElemType ExactTopo, libMesh::Order ExactOrder>
static int
test_side_trace_case_impl(const side_trace_case & info)
{
  constexpr unsigned int max_geom_nodes = 27;

  auto fixture = build_flat_reference_fixture(ExactTopo);
  const unsigned int n_geom_nodes = fixture.elem->n_nodes();
  const unsigned int parent_dim = fixture.elem->dim();

  const libMesh::FEType fe_type(info.key.order, info.key.family);
  auto side_fe = libMesh::FEBase::build(parent_dim, fe_type);

  int fail = 0;

  for (unsigned int side_id = 0; side_id < fixture.elem->n_sides(); ++side_id)
  {
    auto side = fixture.elem->build_side_ptr(side_id);

    libMesh::QGauss qr(parent_dim - 1, libMesh::FOURTH);
    qr.allow_rules_with_negative_weights = true;
    qr.init(side->type());

    side_fe->attach_quadrature_rule(&qr);
    side_fe->get_phi();
    side_fe->get_dphi();
    side_fe->get_normals();
    side_fe->reinit(fixture.elem.get(), side_id);

    const auto & phi = side_fe->get_phi();
    const auto & dphi = side_fe->get_dphi();
    const auto & normals = side_fe->get_normals();

    if (phi.size() != info.n_dofs || dphi.size() != info.n_dofs)
    {
      std::printf("  unexpected host side FE size: %s side_id=%u phi=%llu dphi=%llu expected=%u\n",
                  info.name,
                  side_id,
                  libMesh::cast_int<unsigned long long>(phi.size()),
                  libMesh::cast_int<unsigned long long>(dphi.size()),
                  info.n_dofs);
      fail += 1;
      continue;
    }

    const unsigned int nqp = qr.n_points();
    std::vector<Real> ref_phi(info.n_dofs * nqp);
    std::vector<Real> ref_tgx(info.n_dofs * nqp);
    std::vector<Real> ref_tgy(info.n_dofs * nqp);
    std::vector<Real> ref_tgz(info.n_dofs * nqp);
    std::vector<Real> parent_xi_h(nqp), parent_eta_h(nqp), parent_zeta_h(nqp);
    std::vector<Real> normal_x_h(nqp), normal_y_h(nqp), normal_z_h(nqp);

    for (unsigned int q = 0; q < nqp; ++q)
    {
      const RealVector face_qpt = libMesh::Kokkos::make_vector(
        qr.qp(q)(0),
        (side->dim() >= 2) ? qr.qp(q)(1) : Real(0),
        (side->dim() >= 3) ? qr.qp(q)(2) : Real(0));
      const RealVector parent_qpt =
        libMesh::Kokkos::map_face_qp_to_parent(*side, libMesh::LAGRANGE_MAP, side->type(), face_qpt);

      parent_xi_h[q] = vector_component(parent_qpt, 0);
      parent_eta_h[q] = vector_component(parent_qpt, 1);
      parent_zeta_h[q] = vector_component(parent_qpt, 2);

      const auto & n = normals[q];
      normal_x_h[q] = n(0);
      normal_y_h[q] = n(1);
      normal_z_h[q] = n(2);

      for (unsigned int i = 0; i < info.n_dofs; ++i)
      {
        const unsigned int idx = q * info.n_dofs + i;
        const RealVector host_tg = tangential_component(
          libMesh::Kokkos::make_vector(dphi[i][q](0), dphi[i][q](1), dphi[i][q](2)),
          libMesh::Kokkos::make_vector(n(0), n(1), n(2)));

        ref_phi[idx] = phi[i][q];
        ref_tgx[idx] = vector_component(host_tg, 0);
        ref_tgy[idx] = vector_component(host_tg, 1);
        ref_tgz[idx] = vector_component(host_tg, 2);
      }
    }

    auto d_coords = upload_point_coordinates(*fixture.elem, "side_trace_coords");
    auto d_parent_xi = upload_real(parent_xi_h, "side_trace_parent_xi");
    auto d_parent_eta = upload_real(parent_eta_h, "side_trace_parent_eta");
    auto d_parent_zeta = upload_real(parent_zeta_h, "side_trace_parent_zeta");
    auto d_normal_x = upload_real(normal_x_h, "side_trace_normal_x");
    auto d_normal_y = upload_real(normal_y_h, "side_trace_normal_y");
    auto d_normal_z = upload_real(normal_z_h, "side_trace_normal_z");

    Kokkos::View<double *> d_phi(std::string("side_trace_phi"), ref_phi.size());
    Kokkos::View<double *> d_tgx(std::string("side_trace_tgx"), ref_tgx.size());
    Kokkos::View<double *> d_tgy(std::string("side_trace_tgy"), ref_tgy.size());
    Kokkos::View<double *> d_tgz(std::string("side_trace_tgz"), ref_tgz.size());

    const unsigned int n_dofs = info.n_dofs;
    const unsigned int n_geom_nodes_ = n_geom_nodes;

    Kokkos::parallel_for(
      static_cast<int>(ref_phi.size()),
      KOKKOS_LAMBDA(int idx) {
        const unsigned int q = static_cast<unsigned int>(idx) / n_dofs;
        const unsigned int i = static_cast<unsigned int>(idx) % n_dofs;

        RealVector geom_nodes[max_geom_nodes];
        for (unsigned int k = 0; k < n_geom_nodes_; ++k)
          geom_nodes[k] = libMesh::Kokkos::make_vector(
            d_coords(3 * k + 0), d_coords(3 * k + 1), d_coords(3 * k + 2));

        const Real xi = d_parent_xi(q);
        const Real eta = d_parent_eta(q);
        const Real zeta = d_parent_zeta(q);
        const RealTensor J =
          libMesh::Kokkos::jacobian<libMesh::LAGRANGE, ExactTopo>(geom_nodes, n_geom_nodes_, xi, eta, zeta);
        const RealVector grad_ref = grad_shape_for_key<Family, ExactTopo, ExactOrder>(i, xi, eta, zeta);
        const RealVector grad_phys = libMesh::Kokkos::inverse<RealTensor>(J, parent_dim) * grad_ref;
        const RealVector normal = libMesh::Kokkos::make_vector(d_normal_x(q), d_normal_y(q), d_normal_z(q));
        const RealVector tangential_grad = tangential_component(grad_phys, normal);

        d_phi(idx) = shape_for_key<Family, ExactTopo, ExactOrder>(i, xi, eta, zeta);
        d_tgx(idx) = vector_component(tangential_grad, 0);
        d_tgy(idx) = vector_component(tangential_grad, 1);
        d_tgz(idx) = vector_component(tangential_grad, 2);
      });
    Kokkos::fence();

    const int side_fail =
      compare_device_values(d_phi, ref_phi, value_tol) +
      compare_device_values(d_tgx, ref_tgx, grad_tol) +
      compare_device_values(d_tgy, ref_tgy, grad_tol) +
      compare_device_values(d_tgz, ref_tgz, grad_tol);

    if (side_fail)
      std::printf("  side-trace mismatch: %s side_id=%u side_type=%d (%d failures)\n",
                  info.name,
                  side_id,
                  static_cast<int>(side->type()),
                  side_fail);

    fail += side_fail;
  }

  return fail;
}

struct side_trace_dispatch
{
  explicit side_trace_dispatch(const side_trace_case & in_info) : info(in_info) {}

  template <libMesh::FEFamily Family, libMesh::ElemType ExactTopo, libMesh::Order ExactOrder>
  int operator()() const
  {
    return test_side_trace_case_impl<Family, ExactTopo, ExactOrder>(info);
  }

  int unsupported_key(FEShapeKey key) const
  {
    std::printf("  unexpected unsupported side-trace key: %s family=%d elem_type=%d order=%d\n",
                info.name,
                static_cast<int>(key.family),
                static_cast<int>(key.elem_type),
                static_cast<int>(key.order));
    return 1;
  }

  const side_trace_case & info;
};

static int
test_side_trace_case(const side_trace_case & info)
{
  const side_trace_dispatch dispatch(info);
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
    const int f = test_side_trace_case(info);
    std::printf("[side_trace_lagrange]    [%s] %s  (%d failures)\n",
                info.name, f ? "FAIL" : "PASS", f);
    total_fail += f;
  }

  for (const auto & info : monomial_cases)
  {
    const int f = test_side_trace_case(info);
    std::printf("[side_trace_monomial]    [%s] %s  (%d failures)\n",
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
