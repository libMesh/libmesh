// Host-side contract tests for libMesh::Kokkos hard-fail paths.
//
// This executable self-spawns child processes that intentionally invoke
// unsupported Kokkos FE entry points. A child succeeds only if the call
// returns normally; the parent test expects those calls to terminate with a
// non-zero exit status instead.
//
// Returns 0 on success, non-zero on failure.

#include "gpu/kokkos_fe_evaluator.h"
#include "gpu/kokkos_fe_face_map.h"
#include "gpu/kokkos_fe_map.h"
#include "gpu/kokkos_fe_types.h"

#include "libmesh/elem.h"
#include "libmesh/fe_base.h"
#include "libmesh/node.h"

#include <cstdio>
#include <cstdlib>
#include <memory>
#include <string>
#include <sys/wait.h>
#include <vector>

using libMesh::Kokkos::RealTensor;
using libMesh::Kokkos::RealVector;

namespace
{

struct contract_case
{
  const char * name;
  bool expect_abort;
};

struct element_fixture
{
  std::unique_ptr<libMesh::Elem>              elem;
  std::vector<std::unique_ptr<libMesh::Node>> nodes;
};

} // anonymous namespace

static element_fixture
build_master_fixture(libMesh::ElemType elem_type)
{
  element_fixture fixture;
  fixture.elem = libMesh::Elem::build(elem_type);
  fixture.elem->set_mapping_type(libMesh::LAGRANGE_MAP);
  fixture.nodes.reserve(fixture.elem->n_nodes());

  std::vector<libMesh::Point> refspace_nodes;
  libMesh::FEBase::get_refspace_nodes(elem_type, refspace_nodes);

  for (unsigned int i = 0; i < fixture.elem->n_nodes(); ++i)
  {
    const libMesh::Point refspace = refspace_nodes[i];
    fixture.nodes.push_back(libMesh::Node::build(refspace(0), refspace(1), refspace(2), i));
    fixture.elem->set_node(i, fixture.nodes.back().get());
  }

  return fixture;
}

static std::unique_ptr<const libMesh::Elem>
find_side_of_type(const libMesh::Elem & parent,
                  libMesh::ElemType desired_type,
                  unsigned int & side_id)
{
  for (unsigned int s = 0; s < parent.n_sides(); ++s)
  {
    auto side = parent.build_side_ptr(s);
    if (side->type() == desired_type)
    {
      side_id = s;
      return side;
    }
  }

  side_id = libMesh::invalid_uint;
  return nullptr;
}

static void
invoke_face_jacobian_on_side(const libMesh::Elem & side)
{
  constexpr unsigned int max_face_nodes = 9;

  RealVector face_nodes[max_face_nodes];
  for (unsigned int i = 0; i < side.n_nodes(); ++i)
    face_nodes[i] = libMesh::Kokkos::point_to_real_vector(side.point(i));

  (void)libMesh::Kokkos::face_jacobian(
    libMesh::LAGRANGE_MAP, side.type(), face_nodes, side.n_nodes(), 0.0, 0.0, 0.0);
}

static int
run_child_case(const std::string & case_name)
{
  using libMesh::Kokkos::FEShapeKey;

  if (case_name == "noop")
    return 0;

  if (case_name == "get_side_topology_prism6")
  {
    (void)libMesh::Kokkos::get_side_topology(libMesh::PRISM6);
    return 0;
  }

  if (case_name == "get_side_topology_pyramid5")
  {
    (void)libMesh::Kokkos::get_side_topology(libMesh::PYRAMID5);
    return 0;
  }

  if (case_name == "shape_lagrange_edge4_third")
  {
    (void)libMesh::Kokkos::shape(
      FEShapeKey{libMesh::LAGRANGE, libMesh::EDGE4, libMesh::THIRD}, 0, 0.0, 0.0, 0.0);
    return 0;
  }

  if (case_name == "grad_shape_lagrange_prism6_first")
  {
    (void)libMesh::Kokkos::grad_shape(
      FEShapeKey{libMesh::LAGRANGE, libMesh::PRISM6, libMesh::FIRST}, 0, 0.0, 0.0, 0.0);
    return 0;
  }

  if (case_name == "shape_monomial_hex27_sixth")
  {
    (void)libMesh::Kokkos::shape(
      FEShapeKey{libMesh::MONOMIAL, libMesh::HEX27, libMesh::SIXTH}, 0, 0.0, 0.0, 0.0);
    return 0;
  }

  if (case_name == "grad_shape_monomial_tri7_sixth")
  {
    (void)libMesh::Kokkos::grad_shape(
      FEShapeKey{libMesh::MONOMIAL, libMesh::TRI7, libMesh::SIXTH}, 0, 0.0, 0.0, 0.0);
    return 0;
  }

  if (case_name == "ndofs_lagrange_prism6_first")
  {
    (void)libMesh::Kokkos::n_dofs(
      FEShapeKey{libMesh::LAGRANGE, libMesh::PRISM6, libMesh::FIRST});
    return 0;
  }

  if (case_name == "map_shape_rational")
  {
    (void)libMesh::Kokkos::map_shape(
      libMesh::RATIONAL_BERNSTEIN_MAP, libMesh::QUAD4, 0, 0.0, 0.0, 0.0);
    return 0;
  }

  if (case_name == "grad_map_shape_rational")
  {
    (void)libMesh::Kokkos::grad_map_shape(
      libMesh::RATIONAL_BERNSTEIN_MAP, libMesh::QUAD4, 0, 0.0, 0.0, 0.0);
    return 0;
  }

  if (case_name == "face_normal_parent_dim2")
  {
    RealTensor J = libMesh::Kokkos::zero_tensor();
    J(0, 0) = 1.0;
    J(0, 1) = 2.0;
    (void)libMesh::Kokkos::face_normal(J, 2);
    return 0;
  }

  if (case_name == "map_face_qp_to_parent_prism20_tri7")
  {
    auto fixture = build_master_fixture(libMesh::PRISM20);
    unsigned int side_id = libMesh::invalid_uint;
    auto side = find_side_of_type(*fixture.elem, libMesh::TRI7, side_id);
    if (!side)
      return 2;
    (void)libMesh::Kokkos::map_face_qp_to_parent(
      *side, libMesh::LAGRANGE_MAP, side->type(), libMesh::Kokkos::zero_vector());
    return 0;
  }

  if (case_name == "map_face_qp_to_parent_prism21_tri7")
  {
    auto fixture = build_master_fixture(libMesh::PRISM21);
    unsigned int side_id = libMesh::invalid_uint;
    auto side = find_side_of_type(*fixture.elem, libMesh::TRI7, side_id);
    if (!side)
      return 2;
    (void)libMesh::Kokkos::map_face_qp_to_parent(
      *side, libMesh::LAGRANGE_MAP, side->type(), libMesh::Kokkos::zero_vector());
    return 0;
  }

  if (case_name == "map_face_qp_to_parent_pyramid18_tri7")
  {
    auto fixture = build_master_fixture(libMesh::PYRAMID18);
    unsigned int side_id = libMesh::invalid_uint;
    auto side = find_side_of_type(*fixture.elem, libMesh::TRI7, side_id);
    if (!side)
      return 2;
    (void)libMesh::Kokkos::map_face_qp_to_parent(
      *side, libMesh::LAGRANGE_MAP, side->type(), libMesh::Kokkos::zero_vector());
    return 0;
  }

  if (case_name == "face_jacobian_prism20_tri7")
  {
    auto fixture = build_master_fixture(libMesh::PRISM20);
    unsigned int side_id = libMesh::invalid_uint;
    auto side = find_side_of_type(*fixture.elem, libMesh::TRI7, side_id);
    if (!side)
      return 2;
    invoke_face_jacobian_on_side(*side);
    return 0;
  }

  if (case_name == "face_jacobian_prism21_tri7")
  {
    auto fixture = build_master_fixture(libMesh::PRISM21);
    unsigned int side_id = libMesh::invalid_uint;
    auto side = find_side_of_type(*fixture.elem, libMesh::TRI7, side_id);
    if (!side)
      return 2;
    invoke_face_jacobian_on_side(*side);
    return 0;
  }

  if (case_name == "face_jacobian_pyramid18_tri7")
  {
    auto fixture = build_master_fixture(libMesh::PYRAMID18);
    unsigned int side_id = libMesh::invalid_uint;
    auto side = find_side_of_type(*fixture.elem, libMesh::TRI7, side_id);
    if (!side)
      return 2;
    invoke_face_jacobian_on_side(*side);
    return 0;
  }

  std::fprintf(stderr, "Unknown child case: %s\n", case_name.c_str());
  return 3;
}

static int
run_command(const std::string & command)
{
  std::fflush(nullptr);
  return std::system(command.c_str());
}

static bool
expect_child_success(const char * argv0, const char * case_name)
{
  const std::string command =
    std::string(argv0) + " --child " + case_name + " >/dev/null 2>&1";
  const int status = run_command(command);
  return status == 0;
}

static bool
expect_child_abort(const char * argv0, const char * case_name)
{
  const std::string command =
    std::string(argv0) + " --child " + case_name + " >/dev/null 2>&1";
  const int status = run_command(command);
  if (status == -1)
    return false;

  if (WIFSIGNALED(status))
    return true;

  if (!WIFEXITED(status))
    return false;

  const int exit_code = WEXITSTATUS(status);
  return exit_code != 0 && exit_code != 2 && exit_code != 3;
}

int
main(int argc, char ** argv)
{
  if (argc == 3 && std::string(argv[1]) == "--child")
    return run_child_case(argv[2]);

  if (!expect_child_success(argv[0], "noop"))
  {
    std::printf("[contract_spawn]          FAIL  (could not successfully respawn test executable)\n");
    return 1;
  }

  const contract_case cases[] = {
    { "get_side_topology_prism6", true },
    { "get_side_topology_pyramid5", true },
    { "shape_lagrange_edge4_third", false },
    { "grad_shape_lagrange_prism6_first", true },
    { "shape_monomial_hex27_sixth", true },
    { "grad_shape_monomial_tri7_sixth", true },
    { "ndofs_lagrange_prism6_first", true },
    { "map_shape_rational", true },
    { "grad_map_shape_rational", true },
    { "face_normal_parent_dim2", true },
    { "face_jacobian_prism20_tri7", true },
    { "face_jacobian_prism21_tri7", true },
    { "face_jacobian_pyramid18_tri7", true },
    { "map_face_qp_to_parent_prism20_tri7", true },
    { "map_face_qp_to_parent_prism21_tri7", true },
    { "map_face_qp_to_parent_pyramid18_tri7", true }
  };

  int total_fail = 0;
  for (const auto & info : cases)
  {
    const bool passed = info.expect_abort ? expect_child_abort(argv[0], info.name)
                                          : expect_child_success(argv[0], info.name);
    const int fail = passed ? 0 : 1;
    std::printf("[%s]          [%s] %s  (%d failures)\n",
                info.expect_abort ? "contract_abort" : "contract_success",
                info.name,
                passed ? "PASS" : "FAIL",
                fail);
    total_fail += fail;
  }

  if (total_fail == 0)
    std::printf("ALL TESTS PASSED\n");
  else
    std::printf("%d TEST(S) FAILED\n", total_fail);

  return total_fail ? 1 : 0;
}
