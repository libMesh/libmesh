// GPU kernel oracle tests for libMesh::Kokkos FE shape functions.
// The test suite covers:
//   A. Geometry-map shape parity against the libMesh FE map oracle for the
//      12 implemented LAGRANGE map topologies.
//   B. Geometry-map gradient parity against the libMesh FE map oracle for the
//      same topologies.
//   C. Physics FE parity for exact libMesh LAGRANGE keys that the Kokkos
//      evaluator currently supports.
//   D. Physics FE gradient parity for the same exact LAGRANGE keys.
//   E. Physics FE parity for Kokkos-supported exact MONOMIAL keys.
//   F. Physics FE gradient parity for the same exact MONOMIAL keys.
//
// Unsupported exact keys are expected to hard-abort in the Kokkos path and
// are therefore intentionally not invoked here.

#include "libmesh/elem.h"
#include "gpu/kokkos_fe_base.h"
#include "gpu/kokkos_fe_evaluator.h"
#include "gpu/kokkos_fe_types.h"

#include "libmesh/fe.h"
#include "libmesh/fe_interface.h"
#include "libmesh/fe_map.h"
#include "libmesh/libmesh.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/enum_elem_type.h"
#include "libmesh/enum_order.h"

// Avoid conflicting complex operators between CUDA and PETSc
#define PETSC_SKIP_CXX_COMPLEX_FIX 1
#include <Kokkos_Core.hpp>
#undef __CUDACC_VER__

#include "kokkos_fe_oracle_test_utils.h"
#include <cstdio>
#include <cmath>
#include <string>
#include <vector>

using libMesh::Kokkos::Real;
using libMesh::Kokkos::RealVector;
using kokkos_test_utils::build_qps;
using kokkos_test_utils::build_reference_elem;
using kokkos_test_utils::compare_device_values;
using kokkos_test_utils::dispatch_supported_lagrange_map_topology;
using kokkos_test_utils::dispatch_supported_shape_key;
using kokkos_test_utils::grad_shape_for_key;
using kokkos_test_utils::shape_for_key;
using kokkos_test_utils::upload_real;
using kokkos_test_utils::vector_component;

static constexpr double tol = 1.0e-13;
static constexpr unsigned int quad_order = 4;

namespace
{

struct map_elem_info
{
  libMesh::ElemType topo;
  unsigned int      dim;
  unsigned int      n_dofs;
  const char *      name;
};

struct physics_shape_info
{
  libMesh::Kokkos::FEShapeKey key;
  unsigned int                dim;
  unsigned int                n_dofs;
  const char *                name;
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
  { libMesh::HEX27,  3, 27, "HEX27"  },
};
static constexpr unsigned int n_map_elems = sizeof(map_elems) / sizeof(map_elems[0]);

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

// Only exact libMesh LAGRANGE keys whose evaluator topology is implemented in
// the current Kokkos path are included here.
static const physics_shape_info lagrange_physics_cases[] = {
  { { libMesh::LAGRANGE, libMesh::EDGE2,  libMesh::FIRST  }, 1,  2, "EDGE2/FIRST"   },
  { { libMesh::LAGRANGE, libMesh::EDGE3,  libMesh::FIRST  }, 1,  2, "EDGE3/FIRST"   },
  { { libMesh::LAGRANGE, libMesh::EDGE3,  libMesh::SECOND }, 1,  3, "EDGE3/SECOND"  },
  { { libMesh::LAGRANGE, libMesh::EDGE4,  libMesh::FIRST  }, 1,  2, "EDGE4/FIRST"   },
  { { libMesh::LAGRANGE, libMesh::EDGE4,  libMesh::THIRD  }, 1,  4, "EDGE4/THIRD"   },

  { { libMesh::LAGRANGE, libMesh::TRI3,   libMesh::FIRST  }, 2,  3, "TRI3/FIRST"    },
  { { libMesh::LAGRANGE, libMesh::TRI6,   libMesh::FIRST  }, 2,  3, "TRI6/FIRST"    },
  { { libMesh::LAGRANGE, libMesh::TRI6,   libMesh::SECOND }, 2,  6, "TRI6/SECOND"   },
  { { libMesh::LAGRANGE, libMesh::TRI7,   libMesh::FIRST  }, 2,  3, "TRI7/FIRST"    },
  { { libMesh::LAGRANGE, libMesh::TRI7,   libMesh::SECOND }, 2,  6, "TRI7/SECOND"   },

  { { libMesh::LAGRANGE, libMesh::QUAD4,  libMesh::FIRST  }, 2,  4, "QUAD4/FIRST"   },
  { { libMesh::LAGRANGE, libMesh::QUAD8,  libMesh::FIRST  }, 2,  4, "QUAD8/FIRST"   },
  { { libMesh::LAGRANGE, libMesh::QUAD8,  libMesh::SECOND }, 2,  8, "QUAD8/SECOND"  },
  { { libMesh::LAGRANGE, libMesh::QUAD9,  libMesh::FIRST  }, 2,  4, "QUAD9/FIRST"   },
  { { libMesh::LAGRANGE, libMesh::QUAD9,  libMesh::SECOND }, 2,  9, "QUAD9/SECOND"  },

  { { libMesh::LAGRANGE, libMesh::TET4,   libMesh::FIRST  }, 3,  4, "TET4/FIRST"    },
  { { libMesh::LAGRANGE, libMesh::TET10,  libMesh::FIRST  }, 3,  4, "TET10/FIRST"   },
  { { libMesh::LAGRANGE, libMesh::TET10,  libMesh::SECOND }, 3, 10, "TET10/SECOND"  },
  { { libMesh::LAGRANGE, libMesh::TET14,  libMesh::FIRST  }, 3,  4, "TET14/FIRST"   },
  { { libMesh::LAGRANGE, libMesh::TET14,  libMesh::SECOND }, 3, 10, "TET14/SECOND"  },

  { { libMesh::LAGRANGE, libMesh::HEX8,   libMesh::FIRST  }, 3,  8, "HEX8/FIRST"    },
  { { libMesh::LAGRANGE, libMesh::HEX20,  libMesh::FIRST  }, 3,  8, "HEX20/FIRST"   },
  { { libMesh::LAGRANGE, libMesh::HEX20,  libMesh::SECOND }, 3, 20, "HEX20/SECOND"  },
  { { libMesh::LAGRANGE, libMesh::HEX27,  libMesh::FIRST  }, 3,  8, "HEX27/FIRST"   },
  { { libMesh::LAGRANGE, libMesh::HEX27,  libMesh::SECOND }, 3, 27, "HEX27/SECOND"  },
};
static constexpr unsigned int n_lagrange_physics_cases =
  sizeof(lagrange_physics_cases) / sizeof(lagrange_physics_cases[0]);

// These MONOMIAL cases cover all implemented MonomialImpl<dim, order> paths
// for orders 0..5, plus representative non-simplex exact keys.
static const physics_shape_info monomial_physics_cases[] = {
  { { libMesh::MONOMIAL, libMesh::EDGE2,     libMesh::CONSTANT }, 1, monomial_n_dofs_for_dim(1, libMesh::CONSTANT), "MONOMIAL/EDGE2/CONSTANT"  },
  { { libMesh::MONOMIAL, libMesh::EDGE2,     libMesh::FIRST    }, 1, monomial_n_dofs_for_dim(1, libMesh::FIRST),    "MONOMIAL/EDGE2/FIRST"     },
  { { libMesh::MONOMIAL, libMesh::EDGE3,     libMesh::SECOND   }, 1, monomial_n_dofs_for_dim(1, libMesh::SECOND),   "MONOMIAL/EDGE3/SECOND"    },
  { { libMesh::MONOMIAL, libMesh::EDGE4,     libMesh::THIRD    }, 1, monomial_n_dofs_for_dim(1, libMesh::THIRD),    "MONOMIAL/EDGE4/THIRD"     },
  { { libMesh::MONOMIAL, libMesh::EDGE3,     libMesh::FOURTH   }, 1, monomial_n_dofs_for_dim(1, libMesh::FOURTH),   "MONOMIAL/EDGE3/FOURTH"    },
  { { libMesh::MONOMIAL, libMesh::EDGE3,     libMesh::FIFTH    }, 1, monomial_n_dofs_for_dim(1, libMesh::FIFTH),    "MONOMIAL/EDGE3/FIFTH"     },

  { { libMesh::MONOMIAL, libMesh::TRI3,      libMesh::CONSTANT }, 2, monomial_n_dofs_for_dim(2, libMesh::CONSTANT), "MONOMIAL/TRI3/CONSTANT"   },
  { { libMesh::MONOMIAL, libMesh::TRI3,      libMesh::FIRST    }, 2, monomial_n_dofs_for_dim(2, libMesh::FIRST),    "MONOMIAL/TRI3/FIRST"      },
  { { libMesh::MONOMIAL, libMesh::TRI6,      libMesh::SECOND   }, 2, monomial_n_dofs_for_dim(2, libMesh::SECOND),   "MONOMIAL/TRI6/SECOND"     },
  { { libMesh::MONOMIAL, libMesh::TRI7,      libMesh::THIRD    }, 2, monomial_n_dofs_for_dim(2, libMesh::THIRD),    "MONOMIAL/TRI7/THIRD"      },
  { { libMesh::MONOMIAL, libMesh::TRI7,      libMesh::FOURTH   }, 2, monomial_n_dofs_for_dim(2, libMesh::FOURTH),   "MONOMIAL/TRI7/FOURTH"     },
  { { libMesh::MONOMIAL, libMesh::TRI7,      libMesh::FIFTH    }, 2, monomial_n_dofs_for_dim(2, libMesh::FIFTH),    "MONOMIAL/TRI7/FIFTH"      },
  { { libMesh::MONOMIAL, libMesh::QUAD4,     libMesh::FIRST    }, 2, monomial_n_dofs_for_dim(2, libMesh::FIRST),    "MONOMIAL/QUAD4/FIRST"     },
  { { libMesh::MONOMIAL, libMesh::QUAD9,     libMesh::SECOND   }, 2, monomial_n_dofs_for_dim(2, libMesh::SECOND),   "MONOMIAL/QUAD9/SECOND"    },
  { { libMesh::MONOMIAL, libMesh::QUAD9,     libMesh::FIFTH    }, 2, monomial_n_dofs_for_dim(2, libMesh::FIFTH),    "MONOMIAL/QUAD9/FIFTH"     },

  { { libMesh::MONOMIAL, libMesh::TET4,      libMesh::CONSTANT }, 3, monomial_n_dofs_for_dim(3, libMesh::CONSTANT), "MONOMIAL/TET4/CONSTANT"   },
  { { libMesh::MONOMIAL, libMesh::TET4,      libMesh::FIRST    }, 3, monomial_n_dofs_for_dim(3, libMesh::FIRST),    "MONOMIAL/TET4/FIRST"      },
  { { libMesh::MONOMIAL, libMesh::TET10,     libMesh::SECOND   }, 3, monomial_n_dofs_for_dim(3, libMesh::SECOND),   "MONOMIAL/TET10/SECOND"    },
  { { libMesh::MONOMIAL, libMesh::TET14,     libMesh::THIRD    }, 3, monomial_n_dofs_for_dim(3, libMesh::THIRD),    "MONOMIAL/TET14/THIRD"     },
  { { libMesh::MONOMIAL, libMesh::TET14,     libMesh::FOURTH   }, 3, monomial_n_dofs_for_dim(3, libMesh::FOURTH),   "MONOMIAL/TET14/FOURTH"    },
  { { libMesh::MONOMIAL, libMesh::TET14,     libMesh::FIFTH    }, 3, monomial_n_dofs_for_dim(3, libMesh::FIFTH),    "MONOMIAL/TET14/FIFTH"     },
  { { libMesh::MONOMIAL, libMesh::HEX8,      libMesh::FIRST    }, 3, monomial_n_dofs_for_dim(3, libMesh::FIRST),    "MONOMIAL/HEX8/FIRST"      },
  { { libMesh::MONOMIAL, libMesh::HEX27,     libMesh::SECOND   }, 3, monomial_n_dofs_for_dim(3, libMesh::SECOND),   "MONOMIAL/HEX27/SECOND"    },
  { { libMesh::MONOMIAL, libMesh::HEX27,     libMesh::FIFTH    }, 3, monomial_n_dofs_for_dim(3, libMesh::FIFTH),    "MONOMIAL/HEX27/FIFTH"     },
  { { libMesh::MONOMIAL, libMesh::PRISM6,    libMesh::FIRST    }, 3, monomial_n_dofs_for_dim(3, libMesh::FIRST),    "MONOMIAL/PRISM6/FIRST"    },
  { { libMesh::MONOMIAL, libMesh::PRISM21,   libMesh::FIFTH    }, 3, monomial_n_dofs_for_dim(3, libMesh::FIFTH),    "MONOMIAL/PRISM21/FIFTH"   },
  { { libMesh::MONOMIAL, libMesh::PYRAMID5,  libMesh::FIRST    }, 3, monomial_n_dofs_for_dim(3, libMesh::FIRST),    "MONOMIAL/PYRAMID5/FIRST"  },
  { { libMesh::MONOMIAL, libMesh::PYRAMID14, libMesh::FIFTH    }, 3, monomial_n_dofs_for_dim(3, libMesh::FIFTH),    "MONOMIAL/PYRAMID14/FIFTH" },
};
static constexpr unsigned int n_monomial_physics_cases =
  sizeof(monomial_physics_cases) / sizeof(monomial_physics_cases[0]);

} // anonymous namespace

static Real
host_map_shape(const libMesh::Elem & elem,
               const libMesh::FEType & fe_type,
               unsigned int i,
               Real xi,
               Real eta,
               Real zeta)
{
  return libMesh::FEInterface::shape(fe_type, 0, &elem, i, libMesh::Point(xi, eta, zeta));
}

static RealVector
host_map_grad(const libMesh::Elem & elem,
              const libMesh::FEType & fe_type,
              unsigned int dim,
              unsigned int i,
              Real xi,
              Real eta,
              Real zeta)
{
  const libMesh::Point p(xi, eta, zeta);
  const Real gx = libMesh::FEInterface::shape_deriv(fe_type, 0, &elem, i, 0, p);
  const Real gy = (dim >= 2) ? libMesh::FEInterface::shape_deriv(fe_type, 0, &elem, i, 1, p)
                             : Real(0);
  const Real gz = (dim >= 3) ? libMesh::FEInterface::shape_deriv(fe_type, 0, &elem, i, 2, p)
                             : Real(0);
  return libMesh::Kokkos::make_vector(gx, gy, gz);
}

static Real
host_physics_shape(const libMesh::Elem & elem,
                   const libMesh::FEType & fe_type,
                   unsigned int i,
                   Real xi,
                   Real eta,
                   Real zeta)
{
  return libMesh::FEInterface::shape(fe_type, 0, &elem, i, libMesh::Point(xi, eta, zeta));
}

static RealVector
host_physics_grad(const libMesh::Elem & elem,
                  const libMesh::FEType & fe_type,
                  unsigned int dim,
                  unsigned int i,
                  Real xi,
                  Real eta,
                  Real zeta)
{
  const libMesh::Point p(xi, eta, zeta);
  const Real gx = libMesh::FEInterface::shape_deriv(fe_type, 0, &elem, i, 0, p);
  const Real gy = (dim >= 2) ? libMesh::FEInterface::shape_deriv(fe_type, 0, &elem, i, 1, p)
                             : Real(0);
  const Real gz = (dim >= 3) ? libMesh::FEInterface::shape_deriv(fe_type, 0, &elem, i, 2, p)
                             : Real(0);
  return libMesh::Kokkos::make_vector(gx, gy, gz);
}

// ---------------------------------------------------------------------------
// Test A: Geometry-map shape parity against libMesh FE map dispatch.
// ---------------------------------------------------------------------------
template <libMesh::ElemType Topo>
static int
test_map_shape_parity_impl(const map_elem_info & e)
{
  auto elem = build_reference_elem(Topo);
  const libMesh::FEType fe_type(elem->default_order(), libMesh::FEMap::map_fe_type(*elem));

  std::vector<Real> xi_h, eta_h, zeta_h;
  const unsigned int nqp = build_qps(e.topo, e.dim, quad_order, xi_h, eta_h, zeta_h);
  const unsigned int nd  = e.n_dofs;

  std::vector<Real> ref_phi(nd * nqp);
  for (unsigned int i = 0; i < nd; ++i)
    for (unsigned int q = 0; q < nqp; ++q)
      ref_phi[i * nqp + q] =
        host_map_shape(*elem, fe_type, i, xi_h[q], eta_h[q], zeta_h[q]);

  auto d_xi   = upload_real(xi_h,   "xi");
  auto d_eta  = upload_real(eta_h,  "eta");
  auto d_zeta = upload_real(zeta_h, "zeta");

  Kokkos::View<double *> d_phi(std::string("dev_phi"), nd * nqp);

  const unsigned int      nd_  = nd;
  const unsigned int      nqp_ = nqp;

  Kokkos::parallel_for(
    nd_ * nqp_,
    KOKKOS_LAMBDA(int idx) {
      const int i = idx / static_cast<int>(nqp_);
      const int q = idx % static_cast<int>(nqp_);
      d_phi(idx) =
        libMesh::Kokkos::map_shape<libMesh::LAGRANGE, Topo>(static_cast<unsigned int>(i),
                                                            d_xi(q),
                                                            d_eta(q),
                                                            d_zeta(q));
    });
  Kokkos::fence();

  return compare_device_values(d_phi, ref_phi, tol);
}

struct map_shape_dispatch
{
  explicit map_shape_dispatch(const map_elem_info & in_info) : info(in_info) {}

  template <libMesh::ElemType Topo>
  int operator()() const
  {
    return test_map_shape_parity_impl<Topo>(info);
  }

  int unsupported_topology(libMesh::ElemType topo) const
  {
    std::printf("  unexpected unsupported map-shape topology: %s type=%d\n",
                info.name,
                static_cast<int>(topo));
    return 1;
  }

  const map_elem_info & info;
};

static int
test_map_shape_parity(const map_elem_info & e)
{
  const map_shape_dispatch dispatch(e);
  return dispatch_supported_lagrange_map_topology(e.topo, dispatch);
}

// ---------------------------------------------------------------------------
// Test B: Geometry-map gradient parity against libMesh FE map dispatch.
// ---------------------------------------------------------------------------
template <libMesh::ElemType Topo>
static int
test_map_grad_parity_impl(const map_elem_info & e)
{
  auto elem = build_reference_elem(Topo);
  const libMesh::FEType fe_type(elem->default_order(), libMesh::FEMap::map_fe_type(*elem));

  std::vector<Real> xi_h, eta_h, zeta_h;
  const unsigned int nqp = build_qps(e.topo, e.dim, quad_order, xi_h, eta_h, zeta_h);
  const unsigned int nd  = e.n_dofs;
  const unsigned int dim = e.dim;

  std::vector<Real> ref_gx(nd * nqp);
  std::vector<Real> ref_gy(nd * nqp);
  std::vector<Real> ref_gz(nd * nqp);
  for (unsigned int i = 0; i < nd; ++i)
    for (unsigned int q = 0; q < nqp; ++q)
    {
      const RealVector g = host_map_grad(*elem, fe_type, dim, i, xi_h[q], eta_h[q], zeta_h[q]);
      ref_gx[i * nqp + q] = vector_component(g, 0);
      ref_gy[i * nqp + q] = vector_component(g, 1);
      ref_gz[i * nqp + q] = vector_component(g, 2);
    }

  auto d_xi   = upload_real(xi_h,   "xi");
  auto d_eta  = upload_real(eta_h,  "eta");
  auto d_zeta = upload_real(zeta_h, "zeta");

  Kokkos::View<double *> d_gx(std::string("map_gx"), nd * nqp);
  Kokkos::View<double *> d_gy(std::string("map_gy"), nd * nqp);
  Kokkos::View<double *> d_gz(std::string("map_gz"), nd * nqp);

  const unsigned int      nd_  = nd;
  const unsigned int      nqp_ = nqp;

  Kokkos::parallel_for(
    nd_ * nqp_,
    KOKKOS_LAMBDA(int idx) {
      const int i = idx / static_cast<int>(nqp_);
      const int q = idx % static_cast<int>(nqp_);
      const RealVector g =
        libMesh::Kokkos::grad_map_shape<libMesh::LAGRANGE, Topo>(static_cast<unsigned int>(i),
                                                                 d_xi(q),
                                                                 d_eta(q),
                                                                 d_zeta(q));
      d_gx(idx) = vector_component(g, 0);
      d_gy(idx) = vector_component(g, 1);
      d_gz(idx) = vector_component(g, 2);
    });
  Kokkos::fence();

  return compare_device_values(d_gx, ref_gx, tol) +
         compare_device_values(d_gy, ref_gy, tol) +
         compare_device_values(d_gz, ref_gz, tol);
}

struct map_grad_dispatch
{
  explicit map_grad_dispatch(const map_elem_info & in_info) : info(in_info) {}

  template <libMesh::ElemType Topo>
  int operator()() const
  {
    return test_map_grad_parity_impl<Topo>(info);
  }

  int unsupported_topology(libMesh::ElemType topo) const
  {
    std::printf("  unexpected unsupported map-grad topology: %s type=%d\n",
                info.name,
                static_cast<int>(topo));
    return 1;
  }

  const map_elem_info & info;
};

static int
test_map_grad_parity(const map_elem_info & e)
{
  const map_grad_dispatch dispatch(e);
  return dispatch_supported_lagrange_map_topology(e.topo, dispatch);
}

// ---------------------------------------------------------------------------
// Test C: Exact-key physics shape parity against libMesh CPU FE dispatch.
// ---------------------------------------------------------------------------
template <libMesh::FEFamily Family, libMesh::ElemType ExactTopo, libMesh::Order ExactOrder>
static int
test_shape_parity_impl(const physics_shape_info & info)
{
  auto elem = build_reference_elem(info.key.elem_type);
  const libMesh::FEType fe_type(info.key.order, info.key.family);

  std::vector<Real> xi_h, eta_h, zeta_h;
  const unsigned int nqp = build_qps(info.key.elem_type, info.dim, quad_order, xi_h, eta_h, zeta_h);
  const unsigned int nd  = info.n_dofs;

  std::vector<Real> ref_phi(nd * nqp);
  for (unsigned int i = 0; i < nd; ++i)
    for (unsigned int q = 0; q < nqp; ++q)
      ref_phi[i * nqp + q] =
        host_physics_shape(*elem, fe_type, i, xi_h[q], eta_h[q], zeta_h[q]);

  auto d_xi   = upload_real(xi_h,   "phys_xi");
  auto d_eta  = upload_real(eta_h,  "phys_eta");
  auto d_zeta = upload_real(zeta_h, "phys_zeta");

  Kokkos::View<double *> d_phi(std::string("phys_phi"), nd * nqp);

  const unsigned int                nqp_ = nqp;
  const unsigned int                nd_  = nd;

  Kokkos::parallel_for(
    nd_ * nqp_,
    KOKKOS_LAMBDA(int idx) {
      const int i = idx / static_cast<int>(nqp_);
      const int q = idx % static_cast<int>(nqp_);
      d_phi(idx) = shape_for_key<Family, ExactTopo, ExactOrder>(
        static_cast<unsigned int>(i), d_xi(q), d_eta(q), d_zeta(q));
    });
  Kokkos::fence();

  return compare_device_values(d_phi, ref_phi, tol);
}

// ---------------------------------------------------------------------------
// Test D: Exact-key physics gradient parity against libMesh CPU FE dispatch.
// ---------------------------------------------------------------------------
template <libMesh::FEFamily Family, libMesh::ElemType ExactTopo, libMesh::Order ExactOrder>
static int
test_grad_shape_parity_impl(const physics_shape_info & info)
{
  auto elem = build_reference_elem(info.key.elem_type);
  const libMesh::FEType fe_type(info.key.order, info.key.family);

  std::vector<Real> xi_h, eta_h, zeta_h;
  const unsigned int nqp = build_qps(info.key.elem_type, info.dim, quad_order, xi_h, eta_h, zeta_h);
  const unsigned int nd  = info.n_dofs;

  std::vector<Real> ref_gx(nd * nqp);
  std::vector<Real> ref_gy(nd * nqp);
  std::vector<Real> ref_gz(nd * nqp);
  for (unsigned int i = 0; i < nd; ++i)
    for (unsigned int q = 0; q < nqp; ++q)
    {
      const RealVector g =
        host_physics_grad(*elem, fe_type, info.dim, i, xi_h[q], eta_h[q], zeta_h[q]);
      ref_gx[i * nqp + q] = vector_component(g, 0);
      ref_gy[i * nqp + q] = vector_component(g, 1);
      ref_gz[i * nqp + q] = vector_component(g, 2);
    }

  auto d_xi   = upload_real(xi_h,   "grad_xi");
  auto d_eta  = upload_real(eta_h,  "grad_eta");
  auto d_zeta = upload_real(zeta_h, "grad_zeta");

  Kokkos::View<double *> d_gx(std::string("phys_gx"), nd * nqp);
  Kokkos::View<double *> d_gy(std::string("phys_gy"), nd * nqp);
  Kokkos::View<double *> d_gz(std::string("phys_gz"), nd * nqp);

  const unsigned int                nqp_ = nqp;
  const unsigned int                nd_  = nd;

  Kokkos::parallel_for(
    nd_ * nqp_,
    KOKKOS_LAMBDA(int idx) {
      const int i = idx / static_cast<int>(nqp_);
      const int q = idx % static_cast<int>(nqp_);
      const RealVector g = grad_shape_for_key<Family, ExactTopo, ExactOrder>(
        static_cast<unsigned int>(i), d_xi(q), d_eta(q), d_zeta(q));
      d_gx(idx) = vector_component(g, 0);
      d_gy(idx) = vector_component(g, 1);
      d_gz(idx) = vector_component(g, 2);
    });
  Kokkos::fence();

  return compare_device_values(d_gx, ref_gx, tol) +
         compare_device_values(d_gy, ref_gy, tol) +
         compare_device_values(d_gz, ref_gz, tol);
}

struct shape_dispatch
{
  explicit shape_dispatch(const physics_shape_info & in_info) : info(in_info) {}

  template <libMesh::FEFamily Family, libMesh::ElemType ExactTopo, libMesh::Order ExactOrder>
  int operator()() const
  {
    return test_shape_parity_impl<Family, ExactTopo, ExactOrder>(info);
  }

  int unsupported_key(libMesh::Kokkos::FEShapeKey key) const
  {
    std::printf("  unexpected unsupported shape key: %s family=%d elem_type=%d order=%d\n",
                info.name,
                static_cast<int>(key.family),
                static_cast<int>(key.elem_type),
                static_cast<int>(key.order));
    return 1;
  }

  const physics_shape_info & info;
};

struct grad_shape_dispatch
{
  explicit grad_shape_dispatch(const physics_shape_info & in_info) : info(in_info) {}

  template <libMesh::FEFamily Family, libMesh::ElemType ExactTopo, libMesh::Order ExactOrder>
  int operator()() const
  {
    return test_grad_shape_parity_impl<Family, ExactTopo, ExactOrder>(info);
  }

  int unsupported_key(libMesh::Kokkos::FEShapeKey key) const
  {
    std::printf("  unexpected unsupported grad-shape key: %s family=%d elem_type=%d order=%d\n",
                info.name,
                static_cast<int>(key.family),
                static_cast<int>(key.elem_type),
                static_cast<int>(key.order));
    return 1;
  }

  const physics_shape_info & info;
};

static int
test_shape_parity(const physics_shape_info & info)
{
  const shape_dispatch dispatch(info);
  return dispatch_supported_shape_key(info.key, dispatch);
}

static int
test_grad_shape_parity(const physics_shape_info & info)
{
  const grad_shape_dispatch dispatch(info);
  return dispatch_supported_shape_key(info.key, dispatch);
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------
int main(int argc, char ** argv)
{
  Kokkos::initialize(argc, argv);
  libMesh::LibMeshInit init(argc, argv);

  int total_fail = 0;

  for (unsigned int e = 0; e < n_map_elems; ++e)
  {
    const map_elem_info & info = map_elems[e];

    {
      const int f = test_map_shape_parity(info);
      std::printf("[map_shape_parity     ][%s] %s (%d)\n",
                  info.name, f ? "FAIL" : "PASS", f);
      total_fail += f;
    }
    {
      const int f = test_map_grad_parity(info);
      std::printf("[map_grad_parity      ][%s] %s (%d)\n",
                  info.name, f ? "FAIL" : "PASS", f);
      total_fail += f;
    }
  }

  for (unsigned int c = 0; c < n_lagrange_physics_cases; ++c)
  {
    const physics_shape_info & info = lagrange_physics_cases[c];

    {
      const int f = test_shape_parity(info);
      std::printf("[shape_parity         ][%s] %s (%d)\n",
                  info.name, f ? "FAIL" : "PASS", f);
      total_fail += f;
    }
    {
      const int f = test_grad_shape_parity(info);
      std::printf("[grad_shape_parity    ][%s] %s (%d)\n",
                  info.name, f ? "FAIL" : "PASS", f);
      total_fail += f;
    }
  }

  for (unsigned int c = 0; c < n_monomial_physics_cases; ++c)
  {
    const physics_shape_info & info = monomial_physics_cases[c];

    {
      const int f = test_shape_parity(info);
      std::printf("[shape_parity         ][%s] %s (%d)\n",
                  info.name, f ? "FAIL" : "PASS", f);
      total_fail += f;
    }
    {
      const int f = test_grad_shape_parity(info);
      std::printf("[grad_shape_parity    ][%s] %s (%d)\n",
                  info.name, f ? "FAIL" : "PASS", f);
      total_fail += f;
    }
  }

  if (total_fail == 0)
    std::printf("ALL TESTS PASSED\n");
  else
    std::printf("%d TEST(S) FAILED\n", total_fail);

  Kokkos::finalize();
  return total_fail ? 1 : 0;
}
