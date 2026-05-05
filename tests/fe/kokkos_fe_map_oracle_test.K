// GPU kernel tests for libMesh::Kokkos map helpers across broad topology coverage.
//
// Standalone executable (no CppUnit). Uses libMesh::LibMeshInit so that
// FEMap, FEBase, and FEBase::side_map are available for oracle values.
//
// The test suite covers:
//   A. physical_point_and_jacobian() and volume_jxw() against libMesh FEBase
//      for all implemented LAGRANGE map topologies.
//   B. face_jacobian(), face_jxw(), face_normal(), and
//      edge_normal_on_parent_surface() against libMesh FE oracles for all
//      sides of the implemented 2D and 3D parent topologies.
//   C. map_face_qp_to_parent() against libMesh FEBase::side_map() for all
//      sides and multiple side quadrature points, including supported
//      mixed-face prism and pyramid element types.
//
// Returns 0 on success, non-zero on failure.

#include "libmesh/libmesh_config.h"

#include "gpu/kokkos_fe_face_map.h"
#include "gpu/kokkos_fe_map.h"

#include "libmesh/elem.h"
#include "libmesh/fe_base.h"
#include "libmesh/fe_map.h"
#include "libmesh/libmesh.h"
#include "libmesh/node.h"
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

using libMesh::Kokkos::Real;
using libMesh::Kokkos::RealVector;
using kokkos_test_utils::build_face_helper_context;
using kokkos_test_utils::build_map_helper_context;
using kokkos_test_utils::build_reference_fixture;
using kokkos_test_utils::dispatch_supported_lagrange_face_map_topology;
using kokkos_test_utils::dispatch_supported_lagrange_map_topology;
using kokkos_test_utils::element_fixture;
using kokkos_test_utils::evaluate_face_helper_context_2d;
using kokkos_test_utils::evaluate_face_helper_context_3d;
using kokkos_test_utils::evaluate_map_helper_context;
using kokkos_test_utils::face_helper_context;
using kokkos_test_utils::is_supported_lagrange_face_map_topology;
using kokkos_test_utils::vector_component;

static constexpr double tol = 1.0e-13;

namespace
{

struct map_helper_case
{
  libMesh::ElemType topo;
  const char *      name;
};

struct face_parent_case
{
  libMesh::ElemType topo;
  const char *      name;
};

struct face_qp_parent_case
{
  libMesh::ElemType topo;
  const char *      name;
};

} // anonymous namespace

template <libMesh::ElemType Topo>
static int
test_map_helpers_case_impl(const map_helper_case & info)
{
  auto fixture = build_reference_fixture(Topo);
  const auto context = build_map_helper_context(fixture, info.topo, "map_oracle");
  const int fail = evaluate_map_helper_context<Topo>(context, "map_oracle_results", tol);
  if (fail)
    std::printf("  map helper mismatch: topo=%s (%d failures)\n", info.name, fail);

  return fail;
}

struct map_helper_dispatch
{
  explicit map_helper_dispatch(const map_helper_case & in_info) : info(in_info) {}

  template <libMesh::ElemType Topo>
  int operator()() const
  {
    return test_map_helpers_case_impl<Topo>(info);
  }

  int unsupported_topology(libMesh::ElemType topo) const
  {
    std::printf("  unexpected unsupported map-helper topology: topo=%s type=%d\n",
                info.name,
                static_cast<int>(topo));
    return 1;
  }

  const map_helper_case & info;
};

static int
test_map_helpers_case(const map_helper_case & info)
{
  const map_helper_dispatch dispatch(info);
  return dispatch_supported_lagrange_map_topology(info.topo, dispatch);
}

template <libMesh::ElemType SideTopo>
static int
test_face_helper_side_case_3d_impl(const face_helper_context & context,
                                   unsigned int side_id,
                                   const char * parent_name,
                                   libMesh::ElemType side_topo)
{
  const int fail = evaluate_face_helper_context_3d<SideTopo>(context, "face_oracle_results", tol);
  if (fail)
    std::printf("  face helper mismatch: parent=%s side_id=%u side_type=%d (%d failures)\n",
                parent_name,
                side_id,
                static_cast<int>(side_topo),
                fail);

  return fail;
}

template <libMesh::ElemType ParentTopo, libMesh::ElemType SideTopo>
static int
test_face_helper_side_case_2d_impl(const face_helper_context & context,
                                   unsigned int side_id,
                                   const char * parent_name,
                                   libMesh::ElemType side_topo)
{
  const int fail =
    evaluate_face_helper_context_2d<ParentTopo, SideTopo>(context, "face_oracle_results", tol);
  if (fail)
    std::printf("  face helper mismatch: parent=%s side_id=%u side_type=%d (%d failures)\n",
                parent_name,
                side_id,
                static_cast<int>(side_topo),
                fail);

  return fail;
}

struct face_helper_side_dispatch_3d
{
  face_helper_side_dispatch_3d(const face_helper_context & in_context,
                               unsigned int in_side_id,
                               const char * in_parent_name,
                               libMesh::ElemType in_side_topo)
    : context(in_context),
      side_id(in_side_id),
      parent_name(in_parent_name),
      side_topo(in_side_topo)
  {
  }

  template <libMesh::ElemType SideTopo>
  int operator()() const
  {
    return test_face_helper_side_case_3d_impl<SideTopo>(context, side_id, parent_name, side_topo);
  }

  int unsupported_topology(libMesh::ElemType topo) const
  {
    std::printf("  unexpected unsupported face-helper side: parent=%s side_id=%u side_type=%d\n",
                parent_name,
                side_id,
                static_cast<int>(topo));
    return 1;
  }

  const face_helper_context & context;
  unsigned int                side_id;
  const char *                parent_name;
  libMesh::ElemType           side_topo;
};

template <libMesh::ElemType ParentTopo>
struct face_helper_side_dispatch_2d
{
  face_helper_side_dispatch_2d(const face_helper_context & in_context,
                               unsigned int in_side_id,
                               const char * in_parent_name,
                               libMesh::ElemType in_side_topo)
    : context(in_context),
      side_id(in_side_id),
      parent_name(in_parent_name),
      side_topo(in_side_topo)
  {
  }

  template <libMesh::ElemType SideTopo>
  int operator()() const
  {
    return test_face_helper_side_case_2d_impl<ParentTopo, SideTopo>(
      context, side_id, parent_name, side_topo);
  }

  int unsupported_topology(libMesh::ElemType topo) const
  {
    std::printf("  unexpected unsupported face-helper side: parent=%s side_id=%u side_type=%d\n",
                parent_name,
                side_id,
                static_cast<int>(topo));
    return 1;
  }

  const face_helper_context & context;
  unsigned int                side_id;
  const char *                parent_name;
  libMesh::ElemType           side_topo;
};

struct face_helper_parent_dispatch_2d
{
  face_helper_parent_dispatch_2d(const face_helper_context & in_context,
                                 unsigned int in_side_id,
                                 const char * in_parent_name,
                                 libMesh::ElemType in_side_topo)
    : context(in_context),
      side_id(in_side_id),
      parent_name(in_parent_name),
      side_topo(in_side_topo)
  {
  }

  template <libMesh::ElemType ParentTopo>
  int operator()() const
  {
    const face_helper_side_dispatch_2d<ParentTopo> dispatch(
      context, side_id, parent_name, side_topo);
    return dispatch_supported_lagrange_face_map_topology(side_topo, dispatch);
  }

  int unsupported_topology(libMesh::ElemType topo) const
  {
    std::printf("  unexpected unsupported face-helper parent: parent=%s parent_type=%d side_id=%u\n",
                parent_name,
                static_cast<int>(topo),
                side_id);
    return 1;
  }

  const face_helper_context & context;
  unsigned int                side_id;
  const char *                parent_name;
  libMesh::ElemType           side_topo;
};

static int
test_face_helper_side_case(const element_fixture & fixture,
                           unsigned int side_id,
                           const char * parent_name)
{
  auto side = fixture.elem->build_side_ptr(side_id);
  const face_helper_context context =
    build_face_helper_context(fixture, *side, side_id, "face_oracle");

  if (context.parent_dim == 3)
  {
    const face_helper_side_dispatch_3d dispatch(context, side_id, parent_name, side->type());
    return dispatch_supported_lagrange_face_map_topology(side->type(), dispatch);
  }

  if (context.parent_dim == 2)
  {
    const face_helper_parent_dispatch_2d dispatch(context, side_id, parent_name, side->type());
    return dispatch_supported_lagrange_map_topology(fixture.elem->type(), dispatch);
  }

  std::printf("  unexpected unsupported face-helper parent dimension: parent=%s dim=%u side_id=%u\n",
              parent_name,
              context.parent_dim,
              side_id);
  return 1;
}

static int
test_face_helpers_for_parent(const face_parent_case & info)
{
  auto fixture = build_reference_fixture(info.topo);

  int fail = 0;
  for (unsigned int side_id = 0; side_id < fixture.elem->n_sides(); ++side_id)
    fail += test_face_helper_side_case(fixture, side_id, info.name);

  return fail;
}

static RealVector
host_face_qp_to_parent_oracle(const libMesh::Elem & parent,
                              const libMesh::Elem & side,
                              unsigned int side_id,
                              RealVector face_qpt)
{
  const libMesh::FEType fe_type(parent.default_order(), libMesh::FEMap::map_fe_type(parent));
  auto fe = libMesh::FEBase::build(parent.dim(), fe_type);

  // FE::side_map() relies on FEMap::psi_map, which is only populated after
  // some mapping quantity (e.g. xyz) has been requested on the FE object.
  fe->get_xyz();

  std::vector<libMesh::Point> ref_side_points(1);
  ref_side_points[0] = libMesh::Point(
    vector_component(face_qpt, 0), vector_component(face_qpt, 1), vector_component(face_qpt, 2));

  std::vector<libMesh::Point> ref_points;
  fe->side_map(&parent, &side, side_id, ref_side_points, ref_points);

  return libMesh::Kokkos::make_vector(ref_points[0](0), ref_points[0](1), ref_points[0](2));
}

static int
check_face_qp_to_parent_case(const char * parent_name,
                             const libMesh::Elem & parent,
                             const libMesh::Elem & side,
                             unsigned int side_id,
                             RealVector face_qpt)
{
  using libMesh::Kokkos::map_face_qp_to_parent;

  const RealVector host = host_face_qp_to_parent_oracle(parent, side, side_id, face_qpt);
  const RealVector kokkos =
    map_face_qp_to_parent(side, libMesh::LAGRANGE_MAP, side.type(), face_qpt);

  int fail = 0;
  for (unsigned int d = 0; d < 3; ++d)
    if (std::fabs(vector_component(kokkos, d) - vector_component(host, d)) > tol)
      ++fail;

  if (fail)
  {
    std::vector<libMesh::Point> refspace_nodes;
    libMesh::FEBase::get_refspace_nodes(parent.type(), refspace_nodes);

    std::printf("  face_qp mismatch: parent=%s side_id=%u parent_type=%d side_type=%d\n",
                parent_name,
                side_id,
                static_cast<int>(parent.type()),
                static_cast<int>(side.type()));
    std::printf("    face_qpt=(%.17g, %.17g, %.17g)\n",
                vector_component(face_qpt, 0), vector_component(face_qpt, 1), vector_component(face_qpt, 2));
    std::printf("    host   =(%.17g, %.17g, %.17g)\n",
                vector_component(host, 0), vector_component(host, 1), vector_component(host, 2));
    std::printf("    kokkos =(%.17g, %.17g, %.17g)\n",
                vector_component(kokkos, 0), vector_component(kokkos, 1), vector_component(kokkos, 2));
    std::printf("    side nodes / parent refspace nodes:\n");

    for (unsigned int k = 0; k < side.n_nodes(); ++k)
    {
      const unsigned int parent_node = parent.local_side_node(side_id, k);
      const libMesh::Point parent_refspace = refspace_nodes[parent_node];
      std::printf("      k=%u side_node_id=%llu parent_node=%u parent_refspace=(%.17g, %.17g, %.17g)\n",
                  k,
                  libMesh::cast_int<unsigned long long>(side.node_id(k)),
                  parent_node,
                  parent_refspace(0),
                  parent_refspace(1),
                  parent_refspace(2));
    }
  }

  return fail;
}

static int
test_face_qp_to_parent_for_parent(const face_qp_parent_case & info)
{
  auto fixture = build_reference_fixture(info.topo);

  int fail = 0;
  for (unsigned int side_id = 0; side_id < fixture.elem->n_sides(); ++side_id)
  {
    auto side = fixture.elem->build_side_ptr(side_id);

    if (side->n_nodes() == 1)
    {
      fail +=
        check_face_qp_to_parent_case(info.name, *fixture.elem, *side, side_id, libMesh::Kokkos::zero_vector());
      continue;
    }

    if (!is_supported_lagrange_face_map_topology(side->type()))
    {
      std::printf("  unexpected unsupported face_qp side: parent=%s side_id=%u side_type=%d\n",
                  info.name,
                  side_id,
                  static_cast<int>(side->type()));
      ++fail;
      continue;
    }

    libMesh::QGauss qr(side->dim(), libMesh::FOURTH);
    qr.allow_rules_with_negative_weights = true;
    qr.init(side->type());

    for (unsigned int q = 0; q < qr.n_points(); ++q)
    {
      const RealVector face_qpt = libMesh::Kokkos::make_vector(
        qr.qp(q)(0),
        (side->dim() >= 2) ? qr.qp(q)(1) : Real(0),
        (side->dim() >= 3) ? qr.qp(q)(2) : Real(0));
      fail += check_face_qp_to_parent_case(info.name, *fixture.elem, *side, side_id, face_qpt);
    }
  }

  return fail;
}

int
main(int argc, char ** argv)
{
  Kokkos::initialize(argc, argv);
  libMesh::LibMeshInit init(argc, argv);

  int total_fail = 0;

  {
    const map_helper_case cases[] = {
      { libMesh::EDGE2,  "EDGE2"  },
      { libMesh::EDGE3,  "EDGE3"  },
      { libMesh::TRI3,   "TRI3"   },
      { libMesh::TRI6,   "TRI6"   },
      { libMesh::QUAD4,  "QUAD4"  },
      { libMesh::QUAD8,  "QUAD8"  },
      { libMesh::QUAD9,  "QUAD9"  },
      { libMesh::TET4,   "TET4"   },
      { libMesh::TET10,  "TET10"  },
      { libMesh::HEX8,   "HEX8"   },
      { libMesh::HEX20,  "HEX20"  },
      { libMesh::HEX27,  "HEX27"  }
    };

    for (const auto & info : cases)
    {
      const int f = test_map_helpers_case(info);
      std::printf("[map_helper_breadth]     [%s] %s  (%d failures)\n",
                  info.name, f ? "FAIL" : "PASS", f);
      total_fail += f;
    }
  }

  {
    const face_parent_case cases[] = {
      { libMesh::TRI3,  "TRI3"  },
      { libMesh::TRI6,  "TRI6"  },
      { libMesh::QUAD4, "QUAD4" },
      { libMesh::QUAD8, "QUAD8" },
      { libMesh::QUAD9, "QUAD9" },
      { libMesh::TET4,  "TET4"  },
      { libMesh::TET10, "TET10" },
      { libMesh::HEX8,  "HEX8"  },
      { libMesh::HEX20, "HEX20" },
      { libMesh::HEX27, "HEX27" },
      { libMesh::PRISM6, "PRISM6" },
      { libMesh::PRISM15, "PRISM15" },
      { libMesh::PRISM18, "PRISM18" },
      { libMesh::PYRAMID5, "PYRAMID5" },
      { libMesh::PYRAMID13, "PYRAMID13" },
      { libMesh::PYRAMID14, "PYRAMID14" }
    };

    for (const auto & info : cases)
    {
      const int f = test_face_helpers_for_parent(info);
      std::printf("[face_helper_breadth]    [%s] %s  (%d failures)\n",
                  info.name, f ? "FAIL" : "PASS", f);
      total_fail += f;
    }
  }

  {
    const face_qp_parent_case cases[] = {
      { libMesh::EDGE2,     "EDGE2"     },
      { libMesh::EDGE3,     "EDGE3"     },
      { libMesh::EDGE4,     "EDGE4"     },
      { libMesh::TRI3,      "TRI3"      },
      { libMesh::TRI6,      "TRI6"      },
      { libMesh::QUAD4,     "QUAD4"     },
      { libMesh::QUAD8,     "QUAD8"     },
      { libMesh::QUAD9,     "QUAD9"     },
      { libMesh::TET4,      "TET4"      },
      { libMesh::TET10,     "TET10"     },
      { libMesh::HEX8,      "HEX8"      },
      { libMesh::HEX20,     "HEX20"     },
      { libMesh::HEX27,     "HEX27"     },
      { libMesh::PRISM6,    "PRISM6"    },
      { libMesh::PRISM15,   "PRISM15"   },
      { libMesh::PRISM18,   "PRISM18"   },
      { libMesh::PYRAMID5,  "PYRAMID5"  },
      { libMesh::PYRAMID13, "PYRAMID13" },
      { libMesh::PYRAMID14, "PYRAMID14" }
    };

    for (const auto & info : cases)
    {
      const int f = test_face_qp_to_parent_for_parent(info);
      std::printf("[face_qp_parent_breadth] [%s] %s  (%d failures)\n",
                  info.name, f ? "FAIL" : "PASS", f);
      total_fail += f;
    }
  }

  Kokkos::finalize();

  if (total_fail == 0)
    std::printf("ALL TESTS PASSED\n");
  else
    std::printf("%d TEST(S) FAILED\n", total_fail);

  return total_fail ? 1 : 0;
}
