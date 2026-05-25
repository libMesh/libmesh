// GPU kernel oracle tests for libMesh::Kokkos quadrature and map helpers.
//
// Standalone executable (no CppUnit). Uses libMesh::LibMeshInit so that
// QGauss, FEMap, and FEBase::side_map are available for oracle values.
//
// The test suite covers:
//   A. GaussQuadrature point and weight tables against libMesh QGauss.
//   B. physical_point() and jacobian() against libMesh FEMap::map() and
//      FEMap::map_deriv().
//   C. physical_point_and_jacobian() and volume_jxw() against libMesh FEBase.
//   D. face_jacobian(), face_jxw(), face_normal(), and
//      edge_normal_on_parent_surface() against libMesh FE oracles.
//   E. map_face_qp_to_parent() against libMesh FEBase::side_map().
//
// Returns 0 on success, non-zero on failure.

#include "libmesh/libmesh_config.h"

#include "gpu/kokkos_fe_face_map.h"
#include "gpu/kokkos_fe_map.h"
#include "gpu/kokkos_quadrature.h"

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

namespace libMeshTest
{
namespace KokkosFEQuadratureOracle
{

using libMesh::Kokkos::GaussQuadrature;
using libMesh::Kokkos::Real;
using libMesh::Kokkos::RealTensor;
using libMesh::Kokkos::RealVector;
using kokkos_test_utils::build_face_helper_context;
using kokkos_test_utils::build_map_helper_context;
using kokkos_test_utils::build_host_qgauss;
using kokkos_test_utils::build_reference_fixture;
using kokkos_test_utils::compare_device_values;
using kokkos_test_utils::dispatch_supported_lagrange_face_map_topology;
using kokkos_test_utils::dispatch_supported_lagrange_map_topology;
using kokkos_test_utils::element_fixture;
using kokkos_test_utils::evaluate_face_helper_context_2d;
using kokkos_test_utils::evaluate_face_helper_context_3d;
using kokkos_test_utils::evaluate_map_helper_context;
using kokkos_test_utils::face_helper_context;
using kokkos_test_utils::tensor_component;
using kokkos_test_utils::upload_point_coordinates;
using kokkos_test_utils::vector_component;

static constexpr double tol = 1.0e-13;

namespace
{

struct quadrature_case
{
  libMesh::ElemType topo;
  unsigned int      dim;
  unsigned int      order;
};

struct map_helper_case
{
  libMesh::ElemType topo;
  const char *      name;
};

struct face_helper_case
{
  libMesh::ElemType parent_topo;
  unsigned int      side_id;
  const char *      name;
};

} // anonymous namespace

static int
test_quadrature_case(const quadrature_case & info)
{
  std::vector<Real> x_ref, y_ref, z_ref, w_ref;
  const unsigned int host_nqp =
    build_host_qgauss(info.topo, info.dim, info.order, x_ref, y_ref, z_ref, w_ref);

  Kokkos::View<unsigned int> d_nqp(std::string("nqp"));
  Kokkos::parallel_for(
    1,
    KOKKOS_LAMBDA(int) {
      d_nqp() = GaussQuadrature::n_points(info.topo, info.order);
    });
  Kokkos::fence();

  unsigned int device_nqp = 0;
  Kokkos::deep_copy(device_nqp, d_nqp);

  int fail = 0;
  if (device_nqp != host_nqp)
    ++fail;

  Kokkos::View<double *> d_x(std::string("qx"), host_nqp);
  Kokkos::View<double *> d_y(std::string("qy"), host_nqp);
  Kokkos::View<double *> d_z(std::string("qz"), host_nqp);
  Kokkos::View<double *> d_w(std::string("qw"), host_nqp);

  Kokkos::parallel_for(
    host_nqp,
    KOKKOS_LAMBDA(int qp) {
      const RealVector pt = GaussQuadrature::point(info.topo, info.order, qp);
      d_x(qp) = vector_component(pt, 0);
      d_y(qp) = vector_component(pt, 1);
      d_z(qp) = vector_component(pt, 2);
      d_w(qp) = GaussQuadrature::weight(info.topo, info.order, qp);
    });
  Kokkos::fence();

  fail += compare_device_values(d_x, x_ref);
  fail += compare_device_values(d_y, y_ref);
  fail += compare_device_values(d_z, z_ref);
  fail += compare_device_values(d_w, w_ref);

  if (fail)
    std::printf("  quadrature mismatch: topo=%d dim=%u order=%u (%d failures)\n",
                static_cast<int>(info.topo), info.dim, info.order, fail);

  return fail;
}

static int
test_quadrature_against_qgauss()
{
  int fail = 0;

  const libMesh::ElemType edge_topos[] = { libMesh::EDGE2, libMesh::EDGE3 };
  const libMesh::ElemType quad_topos[] = { libMesh::QUAD4, libMesh::QUAD8, libMesh::QUAD9 };
  const libMesh::ElemType hex_topos[]  = { libMesh::HEX8,  libMesh::HEX20, libMesh::HEX27 };
  const libMesh::ElemType tri_topos[]  = { libMesh::TRI3,  libMesh::TRI6 };
  const libMesh::ElemType tet_topos[]  = { libMesh::TET4,  libMesh::TET10 };

  for (const auto topo : edge_topos)
    for (unsigned int order = 0; order <= 12; ++order)
      fail += test_quadrature_case({ topo, 1, order });

  for (const auto topo : quad_topos)
    for (unsigned int order = 0; order <= 12; ++order)
      fail += test_quadrature_case({ topo, 2, order });

  for (const auto topo : hex_topos)
    for (unsigned int order = 0; order <= 12; ++order)
      fail += test_quadrature_case({ topo, 3, order });

  for (const auto topo : tri_topos)
    for (unsigned int order = 0; order <= 6; ++order)
      fail += test_quadrature_case({ topo, 2, order });

  for (const auto topo : tet_topos)
    for (unsigned int order = 0; order <= 6; ++order)
      fail += test_quadrature_case({ topo, 3, order });

  return fail;
}


static element_fixture
build_hex8_fixture()
{
  element_fixture fixture;
  fixture.elem = libMesh::Elem::build(libMesh::HEX8);
  fixture.elem->set_mapping_type(libMesh::LAGRANGE_MAP);
  fixture.nodes.reserve(8);

  static const double coords[8][3] = {
    {0.0, 0.0, 0.0},
    {1.0, 0.0, 0.0},
    {1.0, 1.0, 0.0},
    {0.0, 1.0, 0.0},
    {0.0, 0.0, 1.0},
    {1.0, 0.0, 1.0},
    {1.0, 1.0, 1.0},
    {0.0, 1.0, 1.0}
  };

  for (unsigned int i = 0; i < 8; ++i)
  {
    fixture.nodes.push_back(libMesh::Node::build(coords[i][0], coords[i][1], coords[i][2], i));
    fixture.elem->set_node(i, fixture.nodes.back().get());
  }

  return fixture;
}

static element_fixture
build_tri3_fixture()
{
  element_fixture fixture;
  fixture.elem = libMesh::Elem::build(libMesh::TRI3);
  fixture.elem->set_mapping_type(libMesh::LAGRANGE_MAP);
  fixture.nodes.reserve(3);

  static const double coords[3][3] = {
    {0.0, 0.0, 0.0},
    {1.0, 0.0, 0.0},
    {0.0, 1.0, 0.0}
  };

  for (unsigned int i = 0; i < 3; ++i)
  {
    fixture.nodes.push_back(libMesh::Node::build(coords[i][0], coords[i][1], coords[i][2], i));
    fixture.elem->set_node(i, fixture.nodes.back().get());
  }

  return fixture;
}

static int
test_physical_map_hex8()
{
  auto fixture = build_hex8_fixture();

  const libMesh::Point ref_center(0.0, 0.0, 0.0);
  const libMesh::Point ref_corner(-1.0, -1.0, -1.0);

  const libMesh::Point host_center = libMesh::FEMap::map(3, fixture.elem.get(), ref_center);
  const libMesh::Point host_corner = libMesh::FEMap::map(3, fixture.elem.get(), ref_corner);
  const libMesh::Point host_dxi = libMesh::FEMap::map_deriv(3, fixture.elem.get(), 0, ref_center);
  const libMesh::Point host_deta = libMesh::FEMap::map_deriv(3, fixture.elem.get(), 1, ref_center);
  const libMesh::Point host_dzeta = libMesh::FEMap::map_deriv(3, fixture.elem.get(), 2, ref_center);

  std::vector<Real> ref_values = {
    host_center(0), host_center(1), host_center(2),
    host_dxi(0), host_dxi(1), host_dxi(2),
    host_deta(0), host_deta(1), host_deta(2),
    host_dzeta(0), host_dzeta(1), host_dzeta(2),
    host_corner(0), host_corner(1), host_corner(2)
  };

  auto d_coords = upload_point_coordinates(*fixture.elem, "hex_coords");

  Kokkos::View<double *> d_results(std::string("hex_results"), ref_values.size());
  Kokkos::parallel_for(
    1,
    KOKKOS_LAMBDA(int) {
      const RealVector xyz_center =
        libMesh::Kokkos::physical_point<libMesh::LAGRANGE, libMesh::HEX8>(d_coords, 8, 0.0, 0.0, 0.0);
      const RealTensor J_center =
        libMesh::Kokkos::jacobian<libMesh::LAGRANGE, libMesh::HEX8>(d_coords, 8, 0.0, 0.0, 0.0);
      const RealVector xyz_corner = libMesh::Kokkos::physical_point<libMesh::LAGRANGE, libMesh::HEX8>(
        d_coords, 8, -1.0, -1.0, -1.0);

      d_results(0) = vector_component(xyz_center, 0);
      d_results(1) = vector_component(xyz_center, 1);
      d_results(2) = vector_component(xyz_center, 2);
      d_results(3) = tensor_component(J_center, 0, 0);
      d_results(4) = tensor_component(J_center, 0, 1);
      d_results(5) = tensor_component(J_center, 0, 2);
      d_results(6) = tensor_component(J_center, 1, 0);
      d_results(7) = tensor_component(J_center, 1, 1);
      d_results(8) = tensor_component(J_center, 1, 2);
      d_results(9) = tensor_component(J_center, 2, 0);
      d_results(10) = tensor_component(J_center, 2, 1);
      d_results(11) = tensor_component(J_center, 2, 2);
      d_results(12) = vector_component(xyz_corner, 0);
      d_results(13) = vector_component(xyz_corner, 1);
      d_results(14) = vector_component(xyz_corner, 2);
    });
  Kokkos::fence();

  return compare_device_values(d_results, ref_values);
}

static int
test_physical_map_tri3()
{
  auto fixture = build_tri3_fixture();

  const libMesh::Point ref_pt(1.0 / 3.0, 1.0 / 3.0, 0.0);

  const libMesh::Point host_xyz = libMesh::FEMap::map(2, fixture.elem.get(), ref_pt);
  const libMesh::Point host_dxi = libMesh::FEMap::map_deriv(2, fixture.elem.get(), 0, ref_pt);
  const libMesh::Point host_deta = libMesh::FEMap::map_deriv(2, fixture.elem.get(), 1, ref_pt);

  std::vector<Real> ref_values = {
    host_xyz(0), host_xyz(1), host_xyz(2),
    host_dxi(0), host_dxi(1), host_dxi(2),
    host_deta(0), host_deta(1), host_deta(2)
  };

  auto d_coords = upload_point_coordinates(*fixture.elem, "tri_coords");

  Kokkos::View<double *> d_results(std::string("tri_results"), ref_values.size());
  Kokkos::parallel_for(
    1,
    KOKKOS_LAMBDA(int) {
      const RealVector xyz = libMesh::Kokkos::physical_point<libMesh::LAGRANGE, libMesh::TRI3>(
        d_coords, 3, 1.0 / 3.0, 1.0 / 3.0, 0.0);
      const RealTensor J = libMesh::Kokkos::jacobian<libMesh::LAGRANGE, libMesh::TRI3>(
        d_coords, 3, 1.0 / 3.0, 1.0 / 3.0, 0.0);

      d_results(0) = vector_component(xyz, 0);
      d_results(1) = vector_component(xyz, 1);
      d_results(2) = vector_component(xyz, 2);
      d_results(3) = tensor_component(J, 0, 0);
      d_results(4) = tensor_component(J, 0, 1);
      d_results(5) = tensor_component(J, 0, 2);
      d_results(6) = tensor_component(J, 1, 0);
      d_results(7) = tensor_component(J, 1, 1);
      d_results(8) = tensor_component(J, 1, 2);
    });
  Kokkos::fence();

  return compare_device_values(d_results, ref_values);
}

template <libMesh::ElemType Topo>
static int
test_map_helpers_case_impl(const map_helper_case & info)
{
  auto fixture = build_reference_fixture(Topo);
  const auto context = build_map_helper_context(fixture, info.topo, "map_helper");
  const int fail = evaluate_map_helper_context<Topo>(context, "map_helper_results", tol);
  if (fail)
    std::printf("  quadrature map-helper mismatch: %s (%d failures)\n",
                info.name,
                fail);
  return fail;
}

struct quadrature_map_helper_dispatch
{
  explicit quadrature_map_helper_dispatch(const map_helper_case & in_info) : info(in_info) {}

  template <libMesh::ElemType Topo>
  int operator()() const
  {
    return test_map_helpers_case_impl<Topo>(info);
  }

  int unsupported_topology(libMesh::ElemType topo) const
  {
    std::printf("  unexpected unsupported quadrature map-helper topology: %s type=%d\n",
                info.name,
                static_cast<int>(topo));
    return 1;
  }

  const map_helper_case & info;
};

static int
test_map_helpers_case(const map_helper_case & info)
{
  const quadrature_map_helper_dispatch dispatch(info);
  return dispatch_supported_lagrange_map_topology(info.topo, dispatch);
}

template <libMesh::ElemType SideTopo>
static int
test_face_map_helpers_case_3d_impl(const face_helper_context & context,
                                   const face_helper_case & info,
                                   libMesh::ElemType side_topo)
{
  const int fail = evaluate_face_helper_context_3d<SideTopo>(context, "face_helper_results", tol);
  if (fail)
    std::printf("  quadrature face-helper mismatch: %s parent_type=%d side_id=%u side_type=%d (%d failures)\n",
                info.name,
                static_cast<int>(info.parent_topo),
                info.side_id,
                static_cast<int>(side_topo),
                fail);

  return fail;
}

template <libMesh::ElemType ParentTopo, libMesh::ElemType SideTopo>
static int
test_face_map_helpers_case_2d_impl(const face_helper_context & context,
                                   const face_helper_case & info,
                                   libMesh::ElemType side_topo)
{
  const int fail =
    evaluate_face_helper_context_2d<ParentTopo, SideTopo>(context, "face_helper_results", tol);
  if (fail)
    std::printf("  quadrature face-helper mismatch: %s parent_type=%d side_id=%u side_type=%d (%d failures)\n",
                info.name,
                static_cast<int>(info.parent_topo),
                info.side_id,
                static_cast<int>(side_topo),
                fail);

  return fail;
}

struct quadrature_face_side_dispatch_3d
{
  quadrature_face_side_dispatch_3d(const face_helper_context & in_context,
                                   const face_helper_case & in_info,
                                   libMesh::ElemType in_side_topo)
    : context(in_context), info(in_info), side_topo(in_side_topo)
  {
  }

  template <libMesh::ElemType SideTopo>
  int operator()() const
  {
    return test_face_map_helpers_case_3d_impl<SideTopo>(context, info, side_topo);
  }

  int unsupported_topology(libMesh::ElemType topo) const
  {
    std::printf("  unexpected unsupported quadrature face-helper side: %s parent_type=%d side_id=%u side_type=%d\n",
                info.name,
                static_cast<int>(info.parent_topo),
                info.side_id,
                static_cast<int>(topo));
    return 1;
  }

  const face_helper_context & context;
  const face_helper_case &    info;
  libMesh::ElemType           side_topo;
};

template <libMesh::ElemType ParentTopo>
struct quadrature_face_side_dispatch_2d
{
  quadrature_face_side_dispatch_2d(const face_helper_context & in_context,
                                   const face_helper_case & in_info,
                                   libMesh::ElemType in_side_topo)
    : context(in_context), info(in_info), side_topo(in_side_topo)
  {
  }

  template <libMesh::ElemType SideTopo>
  int operator()() const
  {
    return test_face_map_helpers_case_2d_impl<ParentTopo, SideTopo>(context, info, side_topo);
  }

  int unsupported_topology(libMesh::ElemType topo) const
  {
    std::printf("  unexpected unsupported quadrature face-helper side: %s parent_type=%d side_id=%u side_type=%d\n",
                info.name,
                static_cast<int>(info.parent_topo),
                info.side_id,
                static_cast<int>(topo));
    return 1;
  }

  const face_helper_context & context;
  const face_helper_case &    info;
  libMesh::ElemType           side_topo;
};

struct quadrature_face_parent_dispatch_2d
{
  quadrature_face_parent_dispatch_2d(const face_helper_context & in_context,
                                     const face_helper_case & in_info,
                                     libMesh::ElemType in_side_topo)
    : context(in_context), info(in_info), side_topo(in_side_topo)
  {
  }

  template <libMesh::ElemType ParentTopo>
  int operator()() const
  {
    const quadrature_face_side_dispatch_2d<ParentTopo> dispatch(context, info, side_topo);
    return dispatch_supported_lagrange_face_map_topology(side_topo, dispatch);
  }

  int unsupported_topology(libMesh::ElemType topo) const
  {
    std::printf("  unexpected unsupported quadrature face-helper parent: %s parent_type=%d side_id=%u\n",
                info.name,
                static_cast<int>(topo),
                info.side_id);
    return 1;
  }

  const face_helper_context & context;
  const face_helper_case &    info;
  libMesh::ElemType           side_topo;
};

static int
test_face_map_helpers_case(const face_helper_case & info)
{
  auto fixture = build_reference_fixture(info.parent_topo);
  auto side = fixture.elem->build_side_ptr(info.side_id);
  const face_helper_context context =
    build_face_helper_context(fixture, *side, info.side_id, "face_helper");

  if (context.parent_dim == 3)
  {
    const quadrature_face_side_dispatch_3d dispatch(context, info, side->type());
    return dispatch_supported_lagrange_face_map_topology(side->type(), dispatch);
  }

  if (context.parent_dim == 2)
  {
    const quadrature_face_parent_dispatch_2d dispatch(context, info, side->type());
    return dispatch_supported_lagrange_map_topology(fixture.elem->type(), dispatch);
  }

  std::printf("  unexpected unsupported quadrature face-helper parent dimension: %s parent_type=%d side_id=%u dim=%u\n",
              info.name,
              static_cast<int>(info.parent_topo),
              info.side_id,
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
check_face_qp_to_parent_case(const char * case_name,
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
    std::printf("  face_qp mismatch: case=%s parent_type=%d side_type=%d side_id=%u\n",
                case_name,
                static_cast<int>(parent.type()),
                static_cast<int>(side.type()),
                side_id);
    std::printf("    face_qpt=(%.17g, %.17g, %.17g)\n",
                vector_component(face_qpt, 0), vector_component(face_qpt, 1), vector_component(face_qpt, 2));
    std::printf("    host   =(%.17g, %.17g, %.17g)\n",
                vector_component(host, 0), vector_component(host, 1), vector_component(host, 2));
    std::printf("    kokkos =(%.17g, %.17g, %.17g)\n",
                vector_component(kokkos, 0), vector_component(kokkos, 1), vector_component(kokkos, 2));
    std::printf("    diff   =(%.17g, %.17g, %.17g)\n",
                vector_component(kokkos, 0) - vector_component(host, 0),
                vector_component(kokkos, 1) - vector_component(host, 1),
                vector_component(kokkos, 2) - vector_component(host, 2));
    std::printf("    side nodes / parent refspace nodes:\n");

    for (unsigned int k = 0; k < side.n_nodes(); ++k)
    {
      libMesh::Point parent_refspace;
      libmesh_error_msg_if(!libMesh::try_reference_side_node(parent.type(), side_id, k, parent_refspace),
                           "check_face_qp_to_parent_case(): unsupported parent side-node lookup");
      std::printf("      k=%u side_node_id=%llu parent_refspace=(%.17g, %.17g, %.17g)\n",
                  k,
                  libMesh::cast_int<unsigned long long>(side.node_id(k)),
                  parent_refspace(0),
                  parent_refspace(1),
                  parent_refspace(2));
    }
  }

  return fail;
}

static int
test_face_qp_to_parent_ref_coords()
{
  using libMesh::Elem;
  using libMesh::Node;

  int fail = 0;

  {
    auto edge = Elem::build(libMesh::EDGE2);
    edge->set_mapping_type(libMesh::LAGRANGE_MAP);
    auto n0 = Node::build(3.25, -2.0, 5.0, 0);
    auto n1 = Node::build(9.50, 4.0, -1.0, 1);
    edge->set_node(0, n0.get());
    edge->set_node(1, n1.get());

    auto side0 = edge->build_side_ptr(0);
    auto side1 = edge->build_side_ptr(1);

    fail += check_face_qp_to_parent_case("edge2_side0", *edge, *side0, 0, libMesh::Kokkos::zero_vector());
    fail += check_face_qp_to_parent_case("edge2_side1", *edge, *side1, 1, libMesh::Kokkos::zero_vector());
  }

  {
    auto tri3 = Elem::build(libMesh::TRI3);
    tri3->set_mapping_type(libMesh::LAGRANGE_MAP);
    auto n0 = Node::build(10.0, 20.0, 0.0, 0);
    auto n1 = Node::build(14.0, 20.0, 0.0, 1);
    auto n2 = Node::build(10.0, 23.0, 0.0, 2);
    tri3->set_node(0, n0.get());
    tri3->set_node(1, n1.get());
    tri3->set_node(2, n2.get());

    auto side0 = tri3->build_side_ptr(0);

    fail += check_face_qp_to_parent_case("tri3_side0", *tri3, *side0, 0, libMesh::Kokkos::zero_vector());
  }

  {
    auto tri6 = Elem::build(libMesh::TRI6);
    tri6->set_mapping_type(libMesh::LAGRANGE_MAP);
    auto n0 = Node::build(4.0, 1.0, 0.0, 0);
    auto n1 = Node::build(9.0, 2.0, 0.0, 1);
    auto n2 = Node::build(3.0, 8.0, 0.0, 2);
    auto n3 = Node::build(42.0, -17.0, 5.0, 3);
    auto n4 = Node::build(11.0, 11.0, 1.0, 4);
    auto n5 = Node::build(-7.0, 4.0, 2.0, 5);
    tri6->set_node(0, n0.get());
    tri6->set_node(1, n1.get());
    tri6->set_node(2, n2.get());
    tri6->set_node(3, n3.get());
    tri6->set_node(4, n4.get());
    tri6->set_node(5, n5.get());

    auto side0 = tri6->build_side_ptr(0);

    fail += check_face_qp_to_parent_case("tri6_side0", *tri6, *side0, 0, libMesh::Kokkos::zero_vector());
  }

  return fail;
}

inline int
run_all_oracles()
{
  int total_fail = 0;

  {
    const int f = test_quadrature_against_qgauss();
    std::printf("[quadrature_qgauss_oracle] %s  (%d failures)\n", f ? "FAIL" : "PASS", f);
    total_fail += f;
  }
  {
    const int f = test_physical_map_hex8();
    std::printf("[physical_map_hex8]       %s  (%d failures)\n", f ? "FAIL" : "PASS", f);
    total_fail += f;
  }
  {
    const int f = test_physical_map_tri3();
    std::printf("[physical_map_tri3]       %s  (%d failures)\n", f ? "FAIL" : "PASS", f);
    total_fail += f;
  }
  {
    const map_helper_case cases[] = {
      { libMesh::EDGE3, "EDGE3" },
      { libMesh::TRI6,  "TRI6"  },
      { libMesh::QUAD9, "QUAD9" },
      { libMesh::TET10, "TET10" },
      { libMesh::HEX20, "HEX20" }
    };

    for (const auto & info : cases)
    {
      const int f = test_map_helpers_case(info);
      std::printf("[map_helper_oracle]      [%s] %s  (%d failures)\n",
                  info.name, f ? "FAIL" : "PASS", f);
      total_fail += f;
    }
  }
  {
    const face_helper_case cases[] = {
      { libMesh::TRI6,  0, "TRI6/side0"  },
      { libMesh::TET10, 0, "TET10/side0" },
      { libMesh::HEX20, 0, "HEX20/side0" }
    };

    for (const auto & info : cases)
    {
      const int f = test_face_map_helpers_case(info);
      std::printf("[face_helper_oracle]     [%s] %s  (%d failures)\n",
                  info.name, f ? "FAIL" : "PASS", f);
      total_fail += f;
    }
  }
  {
    const int f = test_face_qp_to_parent_ref_coords();
    std::printf("[face_qp_parent_oracle]   %s  (%d failures)\n", f ? "FAIL" : "PASS", f);
    total_fail += f;
  }

  return total_fail ? 1 : 0;
}

} // namespace KokkosFEQuadratureOracle
} // namespace libMeshTest
