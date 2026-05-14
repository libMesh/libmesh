// GPU kernel tests for permuted libMesh::Kokkos map helpers.
//
// Standalone executable (no CppUnit). Uses libMesh::LibMeshInit so that
// FEMap, FEBase, and FEBase::side_map are available for oracle values.
//
// The test suite covers:
//   A. physical_point_and_jacobian() and volume_jxw() against libMesh FEBase
//      for every non-identity element permutation of the implemented
//      LAGRANGE map topologies.
//   B. face_jacobian(), face_jxw(), face_normal(), and
//      edge_normal_on_parent_surface() against libMesh FE oracles for every
//      non-identity parent permutation of the supported 2D and 3D parent
//      topologies, including mixed-face prism and pyramid parents.
//   C. map_face_qp_to_parent() against libMesh FEBase::side_map() for every
//      non-identity permutation of those same parent topologies.
//
// Returns 0 on success, non-zero on failure.

#include "libmesh/libmesh_config.h"

#include "gpu/kokkos_fe_face_map.h"
#include "gpu/kokkos_fe_map.h"

#include "libmesh/elem.h"
#include "libmesh/fe_base.h"
#include "libmesh/fe_map.h"
#include "libmesh/libmesh.h"
#include "libmesh/quadrature_gauss.h"

// Avoid conflicting complex operators between CUDA and PETSc
#define PETSC_SKIP_CXX_COMPLEX_FIX 1
#include <Kokkos_Core.hpp>
#undef __CUDACC_VER__

#include "kokkos_fe_oracle_test_utils.h"

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

using libMesh::Kokkos::Real;
using libMesh::Kokkos::RealVector;
using kokkos_test_utils::build_face_helper_context;
using kokkos_test_utils::build_map_helper_context;
using kokkos_test_utils::build_permuted_reference_fixture;
using kokkos_test_utils::build_reference_elem;
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

static const map_helper_case map_cases[] = {
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

static const face_parent_case face_cases[] = {
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

} // anonymous namespace

template <libMesh::ElemType Topo>
static int
test_permuted_map_case_impl(const map_helper_case & info, unsigned int perm_num)
{
  auto fixture = build_permuted_reference_fixture(Topo, perm_num);
  const auto context = build_map_helper_context(fixture, info.topo, "perm_map");
  const int fail = evaluate_map_helper_context<Topo>(context, "perm_map_results", tol);
  if (fail)
    std::printf("  permuted map mismatch: topo=%s perm=%u (%d failures)\n",
                info.name, perm_num, fail);

  return fail;
}

struct permuted_map_dispatch
{
  permuted_map_dispatch(const map_helper_case & in_info, unsigned int in_perm_num)
    : info(in_info), perm_num(in_perm_num)
  {
  }

  template <libMesh::ElemType Topo>
  int operator()() const
  {
    return test_permuted_map_case_impl<Topo>(info, perm_num);
  }

  int unsupported_topology(libMesh::ElemType topo) const
  {
    std::printf("  unexpected unsupported permuted map topology: topo=%s perm=%u type=%d\n",
                info.name,
                perm_num,
                static_cast<int>(topo));
    return 1;
  }

  const map_helper_case & info;
  unsigned int            perm_num;
};

static int
test_permuted_map_case(const map_helper_case & info, unsigned int perm_num)
{
  const permuted_map_dispatch dispatch(info, perm_num);
  return dispatch_supported_lagrange_map_topology(info.topo, dispatch);
}

template <libMesh::ElemType SideTopo>
static int
test_permuted_face_helper_side_case_3d_impl(const face_helper_context & context,
                                            unsigned int side_id,
                                            const char * parent_name,
                                            unsigned int perm_num,
                                            libMesh::ElemType side_topo)
{
  const int fail = evaluate_face_helper_context_3d<SideTopo>(context, "perm_face_results", tol);
  if (fail)
    std::printf("  permuted face mismatch: parent=%s perm=%u side_id=%u side_type=%d (%d failures)\n",
                parent_name,
                perm_num,
                side_id,
                static_cast<int>(side_topo),
                fail);

  return fail;
}

template <libMesh::ElemType ParentTopo, libMesh::ElemType SideTopo>
static int
test_permuted_face_helper_side_case_2d_impl(const face_helper_context & context,
                                            unsigned int side_id,
                                            const char * parent_name,
                                            unsigned int perm_num,
                                            libMesh::ElemType side_topo)
{
  const int fail =
    evaluate_face_helper_context_2d<ParentTopo, SideTopo>(context, "perm_face_results", tol);
  if (fail)
    std::printf("  permuted face mismatch: parent=%s perm=%u side_id=%u side_type=%d (%d failures)\n",
                parent_name,
                perm_num,
                side_id,
                static_cast<int>(side_topo),
                fail);

  return fail;
}

struct permuted_face_side_dispatch_3d
{
  permuted_face_side_dispatch_3d(const face_helper_context & in_context,
                                 unsigned int in_side_id,
                                 const char * in_parent_name,
                                 unsigned int in_perm_num,
                                 libMesh::ElemType in_side_topo)
    : context(in_context),
      side_id(in_side_id),
      parent_name(in_parent_name),
      perm_num(in_perm_num),
      side_topo(in_side_topo)
  {
  }

  template <libMesh::ElemType SideTopo>
  int operator()() const
  {
    return test_permuted_face_helper_side_case_3d_impl<SideTopo>(
      context, side_id, parent_name, perm_num, side_topo);
  }

  int unsupported_topology(libMesh::ElemType topo) const
  {
    std::printf("  unexpected unsupported permuted face-helper side: parent=%s perm=%u side_id=%u side_type=%d\n",
                parent_name,
                perm_num,
                side_id,
                static_cast<int>(topo));
    return 1;
  }

  const face_helper_context & context;
  unsigned int                side_id;
  const char *                parent_name;
  unsigned int                perm_num;
  libMesh::ElemType           side_topo;
};

template <libMesh::ElemType ParentTopo>
struct permuted_face_side_dispatch_2d
{
  permuted_face_side_dispatch_2d(const face_helper_context & in_context,
                                 unsigned int in_side_id,
                                 const char * in_parent_name,
                                 unsigned int in_perm_num,
                                 libMesh::ElemType in_side_topo)
    : context(in_context),
      side_id(in_side_id),
      parent_name(in_parent_name),
      perm_num(in_perm_num),
      side_topo(in_side_topo)
  {
  }

  template <libMesh::ElemType SideTopo>
  int operator()() const
  {
    return test_permuted_face_helper_side_case_2d_impl<ParentTopo, SideTopo>(
      context, side_id, parent_name, perm_num, side_topo);
  }

  int unsupported_topology(libMesh::ElemType topo) const
  {
    std::printf("  unexpected unsupported permuted face-helper side: parent=%s perm=%u side_id=%u side_type=%d\n",
                parent_name,
                perm_num,
                side_id,
                static_cast<int>(topo));
    return 1;
  }

  const face_helper_context & context;
  unsigned int                side_id;
  const char *                parent_name;
  unsigned int                perm_num;
  libMesh::ElemType           side_topo;
};

struct permuted_face_parent_dispatch_2d
{
  permuted_face_parent_dispatch_2d(const face_helper_context & in_context,
                                   unsigned int in_side_id,
                                   const char * in_parent_name,
                                   unsigned int in_perm_num,
                                   libMesh::ElemType in_side_topo)
    : context(in_context),
      side_id(in_side_id),
      parent_name(in_parent_name),
      perm_num(in_perm_num),
      side_topo(in_side_topo)
  {
  }

  template <libMesh::ElemType ParentTopo>
  int operator()() const
  {
    const permuted_face_side_dispatch_2d<ParentTopo> dispatch(
      context, side_id, parent_name, perm_num, side_topo);
    return dispatch_supported_lagrange_face_map_topology(side_topo, dispatch);
  }

  int unsupported_topology(libMesh::ElemType topo) const
  {
    std::printf("  unexpected unsupported permuted face-helper parent: parent=%s perm=%u side_id=%u parent_type=%d\n",
                parent_name,
                perm_num,
                side_id,
                static_cast<int>(topo));
    return 1;
  }

  const face_helper_context & context;
  unsigned int                side_id;
  const char *                parent_name;
  unsigned int                perm_num;
  libMesh::ElemType           side_topo;
};

static int
test_permuted_face_helper_side_case(const element_fixture & fixture,
                                    unsigned int side_id,
                                    const char * parent_name,
                                    unsigned int perm_num)
{
  auto side = fixture.elem->build_side_ptr(side_id);
  const face_helper_context context =
    build_face_helper_context(fixture, *side, side_id, "perm_face");

  if (context.parent_dim == 3)
  {
    const permuted_face_side_dispatch_3d dispatch(
      context, side_id, parent_name, perm_num, side->type());
    return dispatch_supported_lagrange_face_map_topology(side->type(), dispatch);
  }

  if (context.parent_dim == 2)
  {
    const permuted_face_parent_dispatch_2d dispatch(
      context, side_id, parent_name, perm_num, side->type());
    return dispatch_supported_lagrange_map_topology(fixture.elem->type(), dispatch);
  }

  std::printf("  unexpected unsupported permuted face-helper parent dimension: parent=%s perm=%u side_id=%u dim=%u\n",
              parent_name,
              perm_num,
              side_id,
              context.parent_dim);
  return 1;
}

static RealVector
host_face_qp_to_parent_oracle(const libMesh::Elem & parent,
                              const libMesh::Elem & side,
                              unsigned int side_id,
                              RealVector face_qpt)
{
  const libMesh::FEType fe_type(parent.default_order(), libMesh::FEMap::map_fe_type(parent));
  auto fe = libMesh::FEBase::build(parent.dim(), fe_type);
  fe->get_xyz();

  std::vector<libMesh::Point> ref_side_points(1);
  ref_side_points[0] = libMesh::Point(
    vector_component(face_qpt, 0), vector_component(face_qpt, 1), vector_component(face_qpt, 2));

  std::vector<libMesh::Point> ref_points;
  fe->side_map(&parent, &side, side_id, ref_side_points, ref_points);

  return libMesh::Kokkos::make_vector(ref_points[0](0), ref_points[0](1), ref_points[0](2));
}

static int
check_permuted_face_qp_case(const char * parent_name,
                            const libMesh::Elem & parent,
                            const libMesh::Elem & side,
                            unsigned int side_id,
                            unsigned int perm_num,
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
    std::printf("  permuted face_qp mismatch: parent=%s perm=%u side_id=%u side_type=%d\n",
                parent_name,
                perm_num,
                side_id,
                static_cast<int>(side.type()));

  return fail;
}

static int
test_map_helpers_for_all_permutations()
{
  int fail = 0;
  for (const auto & info : map_cases)
  {
    const auto elem = build_reference_elem(info.topo);
    for (unsigned int perm = 1; perm < elem->n_permutations(); ++perm)
      fail += test_permuted_map_case(info, perm);
  }
  return fail;
}

static int
test_face_helpers_for_all_permutations()
{
  int fail = 0;
  for (const auto & info : face_cases)
  {
    const auto elem = build_reference_elem(info.topo);
    for (unsigned int perm = 1; perm < elem->n_permutations(); ++perm)
    {
      auto fixture = build_permuted_reference_fixture(info.topo, perm);
      for (unsigned int side_id = 0; side_id < fixture.elem->n_sides(); ++side_id)
      {
        auto side = fixture.elem->build_side_ptr(side_id);
        if (!is_supported_lagrange_face_map_topology(side->type()))
        {
          std::printf("  unexpected unsupported permuted face-helper side: parent=%s perm=%u side_id=%u side_type=%d\n",
                      info.name,
                      perm,
                      side_id,
                      static_cast<int>(side->type()));
          ++fail;
          continue;
        }

        fail += test_permuted_face_helper_side_case(fixture, side_id, info.name, perm);
      }
    }
  }
  return fail;
}

static int
test_face_qp_to_parent_for_all_permutations()
{
  int fail = 0;
  for (const auto & info : face_cases)
  {
    const auto elem = build_reference_elem(info.topo);
    for (unsigned int perm = 1; perm < elem->n_permutations(); ++perm)
    {
      auto fixture = build_permuted_reference_fixture(info.topo, perm);
      for (unsigned int side_id = 0; side_id < fixture.elem->n_sides(); ++side_id)
      {
        auto side = fixture.elem->build_side_ptr(side_id);
        if (side->n_nodes() == 1)
        {
          fail += check_permuted_face_qp_case(
            info.name, *fixture.elem, *side, side_id, perm, libMesh::Kokkos::zero_vector());
          continue;
        }

        if (!is_supported_lagrange_face_map_topology(side->type()))
        {
          std::printf("  unexpected unsupported permuted face_qp side: parent=%s perm=%u side_id=%u side_type=%d\n",
                      info.name,
                      perm,
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
          fail += check_permuted_face_qp_case(info.name, *fixture.elem, *side, side_id, perm, face_qpt);
        }
      }
    }
  }
  return fail;
}

int
main(int argc, char ** argv)
{
  Kokkos::initialize(argc, argv);
  libMesh::LibMeshInit init(argc, argv);

  int fail = 0;

  const int map_fail = test_map_helpers_for_all_permutations();
  fail += map_fail;
  std::printf("[permuted_map_helpers]    %s  (%d failures)\n", map_fail ? "FAIL" : "OK", map_fail);

  const int face_fail = test_face_helpers_for_all_permutations();
  fail += face_fail;
  std::printf("[permuted_face_helpers]   %s  (%d failures)\n", face_fail ? "FAIL" : "OK", face_fail);

  const int face_qp_fail = test_face_qp_to_parent_for_all_permutations();
  fail += face_qp_fail;
  std::printf("[permuted_face_qp]        %s  (%d failures)\n", face_qp_fail ? "FAIL" : "OK", face_qp_fail);

  Kokkos::finalize();
  return fail ? 1 : 0;
}
