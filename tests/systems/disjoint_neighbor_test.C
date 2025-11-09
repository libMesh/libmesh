#include <libmesh/dirichlet_boundaries.h>
#include <libmesh/dof_map.h>
#include <libmesh/elem.h>
#include <libmesh/equation_systems.h>
#include <libmesh/exodusII_io.h>
#include <libmesh/function_base.h>
#include <libmesh/linear_implicit_system.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/numeric_vector.h>
#include <libmesh/quadrature_gauss.h>
#include <libmesh/sparse_matrix.h>
#include <libmesh/wrapped_function.h>

#include "test_comm.h"
#include "libmesh_cppunit.h"

#include "libmesh/mesh_refinement.h"
#include "libmesh/periodic_boundaries.h"
#include "libmesh/periodic_boundary_base.h"

using namespace libMesh;

static const Real b = 1.0;
static const Real Cgap = 1.0;

// Boundary IDs (counterclockwise): bottom=0, right=1, top=2, left=3
static const boundary_id_type left_id = 3;
static const boundary_id_type right_id = 1;
static const boundary_id_type interface_left_id = 5;
static const boundary_id_type interface_right_id = 6;

static const boundary_id_type interface1_left_id  = 5;
static const boundary_id_type interface1_right_id = 6;
static const boundary_id_type interface2_left_id  = 7;
static const boundary_id_type interface2_right_id = 8;
static const boundary_id_type interface3_left_id  = 9;
static const boundary_id_type interface3_right_id = 10;

static const int Nx = 2, Ny = 2;

Number
heat_exact (const Point & p,
            const Parameters &,
            const std::string &,
            const std::string &)
{
  return (p(0) < 0.5) ? p(0) : p(0) + b;
}

Number left_solution_fn(const Point &,
                        const Parameters &,
                        const std::string &,
                        const std::string &)
{
  return 0.0;
}

Number right_solution_fn(const Point &,
                         const Parameters & params,
                         const std::string &,
                         const std::string &)
{
  const Real b = params.get<Real>("b");
  return 1.0 + b;
}

// Assemble system with temperature jump interface conditions
void assemble_temperature_jump(EquationSystems &es,
                               const std::string & /*system_name*/)
  {
    const MeshBase &mesh = es.get_mesh();
    LinearImplicitSystem &system = es.get_system<LinearImplicitSystem>("TempJump");
    const DofMap &dof_map = system.get_dof_map();
    FEType fe_type = dof_map.variable_type(0);

    // FE objects for volume and face integration
    auto fe = FEBase::build(mesh.mesh_dimension(), fe_type);
    auto fe_face_L = FEBase::build(mesh.mesh_dimension(), fe_type);
    auto fe_face_R = FEBase::build(mesh.mesh_dimension(), fe_type);

    // Quadrature rules
    QGauss qrule(mesh.mesh_dimension(), fe_type.default_quadrature_order());
    QGauss qface(mesh.mesh_dimension() - 1, fe_type.default_quadrature_order());
    fe->attach_quadrature_rule(&qrule);
    fe_face_L->attach_quadrature_rule(&qface);
    fe_face_R->attach_quadrature_rule(&qface);

    const auto &JxW_vol = fe->get_JxW();
    const auto &dphi_vol = fe->get_dphi();
    const auto &phi_L = fe_face_L->get_phi();
    const auto &phi_R = fe_face_R->get_phi();
    const auto &JxW_face = fe_face_L->get_JxW();

    DenseMatrix<Number> Ke;
    DenseVector<Number> Fe;
    std::vector<dof_id_type> dof_indices;

    const BoundaryInfo &boundary = mesh.get_boundary_info();
    const double conductance = Cgap * b; // interface conductance

    for (const Elem *elem : mesh.active_local_element_ptr_range())
    {
      dof_map.dof_indices(elem, dof_indices);
      const unsigned int n_dofs = dof_indices.size();

      Ke.resize(n_dofs, n_dofs);
      Fe.resize(n_dofs);
      Ke.zero();
      Fe.zero();

      // --- Volume contribution ---
      fe->reinit(elem);
      for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
        for (unsigned int i = 0; i < n_dofs; i++)
          for (unsigned int j = 0; j < n_dofs; j++)
            Ke(i, j) += conductance * JxW_vol[qp] * (dphi_vol[i][qp] * dphi_vol[j][qp]);

      // --- Left-side interface ---
      unsigned int side = boundary.side_with_boundary_id(elem, interface_left_id);
      if (side != libMesh::invalid_uint)
        {
          fe_face_L->reinit(elem, side);
          const Elem *neighbor = elem->neighbor_ptr(side);
          libmesh_assert(neighbor);
          unsigned int n_side = neighbor->which_neighbor_am_i(elem);
          fe_face_R->reinit(neighbor, n_side);

          std::vector<dof_id_type> dofs_L, dofs_R;
          dof_map.dof_indices(elem, dofs_L);
          dof_map.dof_indices(neighbor, dofs_R);

          DenseMatrix<Number> Ke_ll(n_dofs, n_dofs);
          DenseMatrix<Number> Ke_lr(n_dofs, dofs_R.size());

          for (unsigned int qp = 0; qp < qface.n_points(); qp++)
            {
              for (unsigned int i = 0; i < n_dofs; i++)
                for (unsigned int j = 0; j < n_dofs; j++)
                  Ke_ll(i, j) += Cgap * phi_L[i][qp] * phi_L[j][qp] * JxW_face[qp];

              for (unsigned int i = 0; i < n_dofs; i++)
                for (unsigned int j = 0; j < dofs_R.size(); j++)
                  Ke_lr(i, j) += -Cgap * phi_L[i][qp] * phi_R[j][qp] * JxW_face[qp];
            }

          Ke += Ke_ll;
          system.matrix->add_matrix(Ke_lr, dofs_L, dofs_R);
        }

      // --- Right-side interface ---
      side = boundary.side_with_boundary_id(elem, interface_right_id);
      if (side != libMesh::invalid_uint)
        {
          fe_face_R->reinit(elem, side);
          const Elem *neighbor = elem->neighbor_ptr(side);
          libmesh_assert(neighbor);
          unsigned int n_side = neighbor->which_neighbor_am_i(elem);
          fe_face_L->reinit(neighbor, n_side);

          std::vector<dof_id_type> dofs_R, dofs_L;
          dof_map.dof_indices(elem, dofs_R);
          dof_map.dof_indices(neighbor, dofs_L);

          DenseMatrix<Number> Ke_rr(n_dofs, n_dofs);
          DenseMatrix<Number> Ke_rl(n_dofs, dofs_L.size());

          for (unsigned int qp = 0; qp < qface.n_points(); qp++)
            {
              for (unsigned int i = 0; i < n_dofs; i++)
                for (unsigned int j = 0; j < n_dofs; j++)
                  Ke_rr(i, j) += Cgap * phi_R[i][qp] * phi_R[j][qp] * JxW_face[qp];

              for (unsigned int i = 0; i < n_dofs; i++)
                for (unsigned int j = 0; j < dofs_L.size(); j++)
                  Ke_rl(i, j) += -Cgap * phi_R[i][qp] * phi_L[j][qp] * JxW_face[qp];
            }

          Ke += Ke_rr;
          system.matrix->add_matrix(Ke_rl, dofs_R, dofs_L);
        }

      // Apply constraints and add local contributions
      dof_map.heterogenously_constrain_element_matrix_and_vector(Ke, Fe, dof_indices);
      system.matrix->add_matrix(Ke, dof_indices);
      system.rhs->add_vector(Fe, dof_indices);
    }

    system.matrix->close();
    system.rhs->close();
  }

class DisjointNeighborTest : public CppUnit::TestCase {
public:
  LIBMESH_CPPUNIT_TEST_SUITE( DisjointNeighborTest );
#ifdef LIBMESH_ENABLE_PERIODIC
#if defined(LIBMESH_HAVE_SOLVER)
  CPPUNIT_TEST( testTempJump );
  CPPUNIT_TEST( testTempJumpRefine );
#endif
#ifdef LIBMESH_ENABLE_AMR
  // // This test intentionally triggers find_neighbors() consistency check
  // // failure after AMR refinement across disjoint interfaces.
  // // Expected: libmesh_assert_valid_neighbors() fails.
  CPPUNIT_TEST(testTempJumpLocalRefineFail);
#endif
  CPPUNIT_TEST( testCloneEquality );
  CPPUNIT_TEST( testStitchingDiscontinuousBoundaries );
  CPPUNIT_TEST( testPreserveDisjointNeighborPairsAfterStitch );
  CPPUNIT_TEST( testStitchCrossMesh );
  CPPUNIT_TEST( testDisjointNeighborConflictError );
#endif

  CPPUNIT_TEST_SUITE_END();

private:

  void testTempJump()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld, 2);
    build_two_disjoint_elems(mesh);

    // Check elem 0 (the left element)
    if (const Elem * elem_left = mesh.query_elem_ptr(0))
      {
        // This processor owns elem 0, so we can safely check its neighbor
        const Elem * elem_right_neighbor = elem_left->neighbor_ptr(1);

        // The neighbor relationship should have been set up by prepare_for_use()
        libmesh_assert(elem_right_neighbor);

        // Verify the neighbor is indeed elem 1
        LIBMESH_ASSERT_NUMBERS_EQUAL(elem_right_neighbor->id(), 1, 1e-15);
      }

    // Check elem 1 (the right element)
    if (const Elem * elem_right = mesh.query_elem_ptr(1))
      {
        // This processor owns elem 1, so we can safely check its neighbor
        const Elem * elem_left_neighbor = elem_right->neighbor_ptr(3);

        // The neighbor relationship should be valid
        libmesh_assert(elem_left_neighbor);

        // Verify the neighbor is indeed elem 0
        LIBMESH_ASSERT_NUMBERS_EQUAL(elem_left_neighbor->id(), 0, 1e-15);
      }

    EquationSystems es(mesh);
    LinearImplicitSystem & sys =
      es.add_system<LinearImplicitSystem>("TempJump");

    const unsigned int u_var = sys.add_variable("u", FEType(FIRST, LAGRANGE));

    std::set<boundary_id_type> left_bdy  { left_id };
    std::set<boundary_id_type> right_bdy { right_id };
    std::vector<unsigned int> vars (1, u_var);

    Parameters params;
    params.set<Real>("b") = b;
    WrappedFunction<Number> left_val(sys, &left_solution_fn);
    WrappedFunction<Number> right_val(sys, &right_solution_fn, &params);
    DirichletBoundary bc_left (left_bdy,  vars, left_val);
    DirichletBoundary bc_right(right_bdy, vars, right_val);

    sys.get_dof_map().add_dirichlet_boundary(bc_left);
    sys.get_dof_map().add_dirichlet_boundary(bc_right);

    sys.attach_assemble_function(assemble_temperature_jump);

    // Ensure the DofMap creates sparsity entries between neighboring elements.
    // Without this, PETSc may report "New nonzero at (a,b) caused a malloc."
    sys.get_dof_map().set_implicit_neighbor_dofs(true);

    es.init();
    sys.solve();

    // ExodusII_IO(mesh).write_equation_systems("temperature_jump.e", es);

    for (Real x=0.; x<=1.; x+=0.05)
      {
        if (std::abs(x - 0.5) < 1e-12) continue; // skip interface
        for (Real y=0.; y<=1.; y+=0.05)
          {
            Point p(x,y);
            const Number exact = heat_exact(p,params,"","");
            const Number approx = sys.point_value(0,p);
            LIBMESH_ASSERT_NUMBERS_EQUAL(exact, approx, 1e-2);
          }
      }
  }


  void testCloneEquality()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld, 2);
    build_two_disjoint_elems(mesh);

    std::unique_ptr<MeshBase> mesh_clone = mesh.clone();
    CPPUNIT_ASSERT(*mesh_clone == mesh);
  }

  void testStitchingDiscontinuousBoundaries()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld, 2);
    build_two_disjoint_elems(mesh);

    CPPUNIT_ASSERT(mesh.parallel_n_nodes() == 8);

    mesh.stitch_surfaces(interface_left_id, interface_right_id, 1e-8 /*tolerance*/,
                        true/*clear_stitched_boundary_ids*/, false /*verbose*/,false/*use_binary_search*/,
                        false /*enforce_all_nodes_match_on_boundaries*/, false /*merge_boundary_nodes_all_or_nothing*/,
                        false /*prepare_after_stitching*/);

    // ExodusII_IO(mesh).write("stitched_disjoint_neighbor_boundary_pairs.e");
    CPPUNIT_ASSERT(mesh.parallel_n_nodes() == 6);

    Mesh mesh2(*TestCommWorld, 2);
    build_two_disjoint_elems(mesh2);

    CPPUNIT_ASSERT(mesh2.parallel_n_nodes() == 8);

    mesh2.stitch_surfaces(interface_left_id, interface_right_id, 1e-8 /*tolerance*/,
                        true/*clear_stitched_boundary_ids*/, false /*verbose*/,false/*use_binary_search*/,
                        false /*enforce_all_nodes_match_on_boundaries*/, false /*merge_boundary_nodes_all_or_nothing*/,
                        true /*prepare_after_stitching*/);

    CPPUNIT_ASSERT(mesh2.parallel_n_nodes() == 6);

    // ExodusII_IO(mesh2).write("stitched_disjoint_neighbor_boundary_pairs_2.e");
  }

  void testTempJumpRefine()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld, 2);
    build_split_mesh_with_interface(mesh, Nx, Ny, 1.0, 1.0);

    EquationSystems es(mesh);
    LinearImplicitSystem & sys =
      es.add_system<LinearImplicitSystem>("TempJump");

    const unsigned int u_var = sys.add_variable("u", FEType(FIRST, LAGRANGE));

    std::set<boundary_id_type> left_bdy  { left_id };
    std::set<boundary_id_type> right_bdy { right_id };
    std::vector<unsigned int> vars (1, u_var);

    Parameters params;
    params.set<Real>("b") = b;
    WrappedFunction<Number> left_val(sys, &left_solution_fn);
    WrappedFunction<Number> right_val(sys, &right_solution_fn, &params);
    DirichletBoundary bc_left (left_bdy,  vars, left_val);
    DirichletBoundary bc_right(right_bdy, vars, right_val);

    sys.get_dof_map().add_dirichlet_boundary(bc_left);
    sys.get_dof_map().add_dirichlet_boundary(bc_right);

    sys.attach_assemble_function(assemble_temperature_jump);

    // Ensure the DofMap creates sparsity entries between neighboring elements.
    // Without this, PETSc may report "New nonzero at (a,b) caused a malloc."
    sys.get_dof_map().set_implicit_neighbor_dofs(true);

    es.init();
    sys.solve();

    // ExodusII_IO(mesh).write_equation_systems("temperature_jump_refine.e", es);

    for (Real x=0.; x<=1.; x+=0.05)
      {
        if (std::abs(x - 0.5) < 1e-12) continue; // skip interface
        for (Real y=0.; y<=1.; y+=0.05)
          {
            Point p(x,y);
            const Number exact = heat_exact(p,params,"","");
            const Number approx = sys.point_value(0,p);
            LIBMESH_ASSERT_NUMBERS_EQUAL(exact, approx, 1e-2);
          }
      }
  }


  void testTempJumpLocalRefineFail()
  {
    LOG_UNIT_TEST;

    try
    {
      Mesh mesh(*TestCommWorld, 2);
      build_split_mesh_with_interface(mesh, Nx, Ny, 1.0, 1.0, true);
    }
    catch (const std::exception &e)
    {
      CPPUNIT_ASSERT(std::string(e.what()).find("Mesh refinement is not yet implemented") != std::string::npos);
    }
  }

  void testPreserveDisjointNeighborPairsAfterStitch()
  {
    LOG_UNIT_TEST;
    Mesh mesh(*TestCommWorld, 2);
    build_four_disjoint_elems(mesh);
    // ExodusII_IO(mesh).write("four_elem_one_row.e");
    CPPUNIT_ASSERT(mesh.parallel_n_nodes() == 16);

    mesh.stitch_surfaces(interface2_left_id,
                        interface2_right_id,
                        1e-8 /* tolerance */,
                        true  /* clear_stitched_boundary_ids */,
                        false /* verbose */,
                        false /* use_binary_search */,
                        false /* enforce_all_nodes_match_on_boundaries */,
                        false /* merge_boundary_nodes_all_or_nothing */,
                        true  /* prepare_after_stitching */);

    // ExodusII_IO(mesh).write("four_elem_one_row_after_stitching.e");
    CPPUNIT_ASSERT(mesh.parallel_n_nodes() == 14);
    const auto * disjoint_pairs = mesh.get_disjoint_neighbor_boundary_pairs();

    // make sure that both the interface1 and the interface3 pairing still exist afterward.
    CPPUNIT_ASSERT(disjoint_pairs->size() == 4);
    for (const auto & [bid, pb_ptr] : *disjoint_pairs)
      {
          const auto & pb = *pb_ptr;
          const boundary_id_type a = pb.myboundary;
          const boundary_id_type b = pb.pairedboundary;
          if ((a == interface1_left_id && b == interface1_right_id) ||
              (a == interface3_left_id && b == interface3_right_id) ||
              (a == interface1_right_id && b == interface1_left_id) ||
              (a == interface3_right_id && b == interface3_left_id))
            continue;
          else
            CPPUNIT_FAIL("Disjoint neighbor boundary pair missing after stitching!");
      }

  }

  void testStitchCrossMesh()
  {
    LOG_UNIT_TEST;

    // Create a mesh with interface boundaries but without pairing
    Mesh mesh_Y(*TestCommWorld, 2);
    build_two_disjoint_elems(mesh_Y);

    // Create another mesh that defines the same pair properly
    auto left_id_x = left_id + 10; // to avoid boundary ID conflict with mesh_Y
    auto right_id_x = right_id + 10;
    auto interface_left_id_x = interface_left_id + 10;
    auto interface_right_id_x = interface_right_id + 10;

    Mesh mesh_X(*TestCommWorld, 2);
    const auto translate_vec = Point(1.0, 0.0, 0.0);
    // Left subdomain nodes
    mesh_X.add_point(Point(0.0, 0.0) + translate_vec, 0);   // bottom-left corner
    mesh_X.add_point(Point(0.5, 0.0) + translate_vec, 1);   // bottom-right corner of left element (interface node)
    mesh_X.add_point(Point(0.0, 1.0) + translate_vec, 2);   // top-left corner
    mesh_X.add_point(Point(0.5, 1.0) + translate_vec, 3);   // top-right corner of left element (interface node)

    // Right subdomain nodes (duplicated interface points)
    mesh_X.add_point(Point(0.5, 0.0) + translate_vec, 4);   // bottom-left corner of right element (same coords as node 1) (interface node)
    mesh_X.add_point(Point(1.0, 0.0) + translate_vec, 5);   // bottom-right corner
    mesh_X.add_point(Point(0.5, 1.0) + translate_vec, 6);   // top-left corner of right element (same coords as node 3) (interface node)
    mesh_X.add_point(Point(1.0, 1.0) + translate_vec, 7);   // top-right corner


    // ---- Define elements ----
    BoundaryInfo & boundary = mesh_X.get_boundary_info();
    // Left element (element ID = 0)
    {
      Elem * elem = mesh_X.add_elem(Elem::build_with_id(QUAD4, 0));
      elem->set_node(0, mesh_X.node_ptr(0)); // bottom-left  (0,0)
      elem->set_node(1, mesh_X.node_ptr(1)); // bottom-right (0.5,0)
      elem->set_node(2, mesh_X.node_ptr(3)); // top-right    (0.5,1)
      elem->set_node(3, mesh_X.node_ptr(2)); // top-left     (0,1)
      boundary.add_side(elem, 3, left_id_x); // left boundary
      boundary.add_side(elem, 1, interface_left_id_x);
      boundary.sideset_name(left_id_x) = "left_boundary_x";
      boundary.sideset_name(interface_left_id_x) = "interface_left_x";
    }

    // Right element (element ID = 1)
    {
      Elem * elem = mesh_X.add_elem(Elem::build_with_id(QUAD4, 1));
      elem->set_node(0, mesh_X.node_ptr(4)); // bottom-left  (0.5,0)_R
      elem->set_node(1, mesh_X.node_ptr(5)); // bottom-right (1,0)
      elem->set_node(2, mesh_X.node_ptr(7)); // top-right    (1,1)
      elem->set_node(3, mesh_X.node_ptr(6)); // top-left     (0.5,1)_R
      boundary.add_side(elem, 1, right_id_x); // right boundary
      boundary.add_side(elem, 3, interface_right_id_x);
      boundary.sideset_name(right_id_x) = "right_boundary_x";
      boundary.sideset_name(interface_right_id_x) = "interface_right_x";
    }

    mesh_X.add_disjoint_neighbor_boundary_pairs(interface_left_id_x, interface_right_id_x, RealVectorValue(0.0, 0.0, 0.0));

    // libMesh shouldn't renumber, or our based-on-initial-id
    // assertions later may fail.
    mesh_X.allow_renumbering(false);

    mesh_X.prepare_for_use();

    // Stitching should trigger an error since Y already has the same
    // boundary IDs but not registered as a pair
    mesh_Y.stitch_meshes(mesh_X,
                        right_id, left_id_x,
                        1e-8,
                        true,  // clear_stitched_boundary_ids
                        false, // verbose
                        false, // use_binary_search
                        false, // enforce_all_nodes_match_on_boundaries
                        false, // merge_boundary_nodes_all_or_nothing
                        false); // prepare_after_stitching

    CPPUNIT_ASSERT(mesh_Y.parallel_n_nodes() == 14);
    const auto * disjoint_pairs = mesh_Y.get_disjoint_neighbor_boundary_pairs();

    // ExodusII_IO(mesh_Y).write("mesh_Y_after_stitching.e");
    // make sure that both the interface and the interface_x pairing still exist afterward.
    CPPUNIT_ASSERT(disjoint_pairs->size() == 4);
    for (const auto & [bid, pb_ptr] : *disjoint_pairs)
      {
          const auto & pb = *pb_ptr;
          const boundary_id_type a = pb.myboundary;
          const boundary_id_type b = pb.pairedboundary;
          if ((a == interface_left_id_x && b == interface_right_id_x) ||
              (a == interface_right_id_x && b == interface_left_id_x) ||
              (a == interface_left_id && b == interface_right_id) ||
              (a == interface_right_id && b == interface_left_id))
            continue;
          else
            CPPUNIT_FAIL("Disjoint neighbor boundary pair missing after stitching!");
      }


  }


  void testDisjointNeighborConflictError()
  {
    LOG_UNIT_TEST;

    try
    {
      Mesh mesh_Y(*TestCommWorld, 2);
      build_two_disjoint_elems(mesh_Y, false);
      mesh_Y.add_disjoint_neighbor_boundary_pairs(right_id, left_id, RealVectorValue(-1.0, 0.0, 0.0));

      // Create another mesh that defines the same pair properly
      auto right_id_x = right_id + 10;
      auto interface_left_id_x = interface_left_id + 10;
      auto interface_right_id_x = interface_right_id + 10;

      Mesh mesh_X(*TestCommWorld, 2);
      const auto translate_vec = Point(1.0, 0.0, 0.0);
      // Left subdomain nodes
      mesh_X.add_point(Point(0.0, 0.0) + translate_vec, 0);   // bottom-left corner
      mesh_X.add_point(Point(0.5, 0.0) + translate_vec, 1);   // bottom-right corner of left element (interface node)
      mesh_X.add_point(Point(0.0, 1.0) + translate_vec, 2);   // top-left corner
      mesh_X.add_point(Point(0.5, 1.0) + translate_vec, 3);   // top-right corner of left element (interface node)

      // Right subdomain nodes (duplicated interface points)
      mesh_X.add_point(Point(0.5, 0.0) + translate_vec, 4);   // bottom-left corner of right element (same coords as node 1) (interface node)
      mesh_X.add_point(Point(1.0, 0.0) + translate_vec, 5);   // bottom-right corner
      mesh_X.add_point(Point(0.5, 1.0) + translate_vec, 6);   // top-left corner of right element (same coords as node 3) (interface node)
      mesh_X.add_point(Point(1.0, 1.0) + translate_vec, 7);   // top-right corner


      // ---- Define elements ----
      BoundaryInfo & boundary = mesh_X.get_boundary_info();
      // Left element (element ID = 0)
      {
        Elem * elem = mesh_X.add_elem(Elem::build_with_id(QUAD4, 0));
        elem->set_node(0, mesh_X.node_ptr(0)); // bottom-left  (0,0)
        elem->set_node(1, mesh_X.node_ptr(1)); // bottom-right (0.5,0)
        elem->set_node(2, mesh_X.node_ptr(3)); // top-right    (0.5,1)
        elem->set_node(3, mesh_X.node_ptr(2)); // top-left     (0,1)
        boundary.add_side(elem, 3, right_id); // left boundary -> use the same ID as mesh_Y to trigger conflict
        boundary.add_side(elem, 1, interface_left_id_x);
        boundary.sideset_name(right_id) = "left_boundary_x";
        boundary.sideset_name(interface_left_id_x) = "interface_left_x";
      }

      // Right element (element ID = 1)
      {
        Elem * elem = mesh_X.add_elem(Elem::build_with_id(QUAD4, 1));
        elem->set_node(0, mesh_X.node_ptr(4)); // bottom-left  (0.5,0)_R
        elem->set_node(1, mesh_X.node_ptr(5)); // bottom-right (1,0)
        elem->set_node(2, mesh_X.node_ptr(7)); // top-right    (1,1)
        elem->set_node(3, mesh_X.node_ptr(6)); // top-left     (0.5,1)_R
        boundary.add_side(elem, 1, right_id_x); // right boundary
        boundary.add_side(elem, 3, interface_right_id_x);
        boundary.sideset_name(right_id_x) = "right_boundary_x";
        boundary.sideset_name(interface_right_id_x) = "interface_right_x";
      }

      mesh_X.add_disjoint_neighbor_boundary_pairs(right_id, right_id_x, RealVectorValue(1.0, 0.0, 0.0));

      // libMesh shouldn't renumber, or our based-on-initial-id
      // assertions later may fail.
      mesh_X.allow_renumbering(false);

      mesh_X.prepare_for_use();

      // Stitching should trigger an error since Y already has the same
      // boundary IDs but not registered as a pair
      mesh_Y.stitch_meshes(mesh_X,
                          right_id, right_id_x,
                          1e-8,
                          true,  // clear_stitched_boundary_ids
                          false, // verbose
                          false, // use_binary_search
                          false, // enforce_all_nodes_match_on_boundaries
                          false, // merge_boundary_nodes_all_or_nothing
                          false); // prepare_after_stitching
    }
    catch (const std::exception &e)
    {
      CPPUNIT_ASSERT(std::string(e.what()).find("Disjoint boundary pairing mismatch") != std::string::npos);
    }
  }



  void build_two_disjoint_elems(Mesh &mesh, bool add_pair = true, const Point &translate_vec = Point(0.0, 0.0, 0.0))
  {
    // Domain: x in (0, 1), y in (0, 1)
    // Split into two subdomains:
    //   Left subdomain:  0 <= x <= 0.5
    //   Right subdomain: 0.5 <= x <= 1
    //
    // Note: Points at x = 0.5 are duplicated (same coordinates but different node IDs)
    //       to represent an interface or discontinuity between the two subdomains.
    //
    // Coordinates layout:
    //
    //  (0,1)   (0.5,1)_L   (0.5,1)_R   (1,1)
    //    x--------x           x--------x
    //    |        |           |        |
    //    |  Left  | Interface |  Right |
    //    |        |           |        |
    //    x--------x           x--------x
    //  (0,0)   (0.5,0)_L   (0.5,0)_R   (1,0)


    // ---- Define points ----

    // Left subdomain nodes
    mesh.add_point(Point(0.0, 0.0) + translate_vec, 0);   // bottom-left corner
    mesh.add_point(Point(0.5, 0.0) + translate_vec, 1);   // bottom-right corner of left element (interface node)
    mesh.add_point(Point(0.0, 1.0) + translate_vec, 2);   // top-left corner
    mesh.add_point(Point(0.5, 1.0) + translate_vec, 3);   // top-right corner of left element (interface node)

    // Right subdomain nodes (duplicated interface points)
    mesh.add_point(Point(0.5, 0.0) + translate_vec, 4);   // bottom-left corner of right element (same coords as node 1) (interface node)
    mesh.add_point(Point(1.0, 0.0) + translate_vec, 5);   // bottom-right corner
    mesh.add_point(Point(0.5, 1.0) + translate_vec, 6);   // top-left corner of right element (same coords as node 3) (interface node)
    mesh.add_point(Point(1.0, 1.0) + translate_vec, 7);   // top-right corner


    // ---- Define elements ----
    BoundaryInfo & boundary = mesh.get_boundary_info();
    // Left element (element ID = 0)
    {
      Elem * elem = mesh.add_elem(Elem::build_with_id(QUAD4, 0));
      elem->set_node(0, mesh.node_ptr(0)); // bottom-left  (0,0)
      elem->set_node(1, mesh.node_ptr(1)); // bottom-right (0.5,0)
      elem->set_node(2, mesh.node_ptr(3)); // top-right    (0.5,1)
      elem->set_node(3, mesh.node_ptr(2)); // top-left     (0,1)
      boundary.add_side(elem, 3, left_id); // left boundary
      boundary.add_side(elem, 1, interface_left_id);
      boundary.sideset_name(left_id) = "left_boundary";
      boundary.sideset_name(interface_left_id) = "interface_left";
    }

    // Right element (element ID = 1)
    {
      Elem * elem = mesh.add_elem(Elem::build_with_id(QUAD4, 1));
      elem->set_node(0, mesh.node_ptr(4)); // bottom-left  (0.5,0)_R
      elem->set_node(1, mesh.node_ptr(5)); // bottom-right (1,0)
      elem->set_node(2, mesh.node_ptr(7)); // top-right    (1,1)
      elem->set_node(3, mesh.node_ptr(6)); // top-left     (0.5,1)_R
      boundary.add_side(elem, 1, right_id); // right boundary
      boundary.add_side(elem, 3, interface_right_id);
      boundary.sideset_name(right_id) = "right_boundary";
      boundary.sideset_name(interface_right_id) = "interface_right";
    }

    // This is the key testing step: inform libMesh about the disjoint neighbor boundary pairs
    // And, in `prepare_for_use()`, libMesh will set up the disjoint neighbor relationships.
    if (add_pair)
      mesh.add_disjoint_neighbor_boundary_pairs(interface_left_id, interface_right_id, RealVectorValue(0.0, 0.0, 0.0));

    // libMesh shouldn't renumber, or our based-on-initial-id
    // assertions later may fail.
    mesh.allow_renumbering(false);

    mesh.prepare_for_use();
  }

  // The interface is located at x = 0.5; nx must be even (split evenly into left and right subdomains)
  void build_split_mesh_with_interface(Mesh &mesh, unsigned nx, unsigned ny, double Lx, double Ly, bool test_local_refinement = false)
  {
    // Ensure nx is even so the interface aligns with element boundaries
    libmesh_error_msg_if(nx % 2 != 0, "nx must be even!");

    mesh.clear();
    mesh.set_mesh_dimension(2);

    const unsigned nxL = nx / 2;       // Number of elements in x-direction (left half)
    const unsigned nxR = nx / 2;       // Number of elements in x-direction (right half)
    const double dx = 1.0 / static_cast<double>(nx);
    const double dy = 1.0 / static_cast<double>(ny);

    // --- Generate points for the left subdomain [0, 0.5] ---
    // Each row of the left subdomain has (nxL + 1) nodes.
    // Total nodes = (nxL + 1) * (ny + 1).
    auto nidL = [&](unsigned i, unsigned j) {
      return static_cast<dof_id_type>(j * (nxL + 1) + i);
    };

    for (unsigned j = 0; j <= ny; ++j)
    {
      const double y = j * dy;
      for (unsigned i = 0; i <= nxL; ++i)
      {
        const double x = i * dx; // At i = nxL, x = 0.5 (interface)
        mesh.add_point(Point(x, y), nidL(i, j));
      }
    }

    // --- Generate points for the right subdomain [0.5, 1.0] ---
    // Interface nodes at x = 0.5 are duplicated with new node IDs.
    const dof_id_type baseR = static_cast<dof_id_type>((nxL + 1) * (ny + 1));

    auto nidR = [&](unsigned i, unsigned j) {
      return static_cast<dof_id_type>(baseR + j * (nxR + 1) + i);
    };

    for (unsigned j = 0; j <= ny; ++j)
    {
      const double y = j * dy;
      for (unsigned i = 0; i <= nxR; ++i)
      {
        const double x = 0.5 + i * dx; // At i = 0, x = 0.5 (same coords as left interface, different ID)
        mesh.add_point(Point(x, y), nidR(i, j));
      }
    }

    BoundaryInfo &boundary = mesh.get_boundary_info();

    // --- Create left subdomain elements (QUAD4) ---
    // Node order: 0 = BL, 1 = BR, 2 = TR, 3 = TL
    // Side order: 0 = bottom(0-1), 1 = right(1-2), 2 = top(2-3), 3 = left(3-0)
    dof_id_type eid = 0;
    for (unsigned j = 0; j < ny; ++j)
    {
      for (unsigned i = 0; i < nxL; ++i)
      {
        Elem *e = mesh.add_elem(Elem::build_with_id(QUAD4, eid++));
        e->set_node(0, mesh.node_ptr(nidL(i,   j)));
        e->set_node(1, mesh.node_ptr(nidL(i+1, j)));
        e->set_node(2, mesh.node_ptr(nidL(i+1, j+1)));
        e->set_node(3, mesh.node_ptr(nidL(i,   j+1)));

        // Left outer boundary (i == 0 -> side 3)
        if (i == 0)
          boundary.add_side(e, 3, left_id);

        // Interface (left side) boundary (i == nxL - 1 -> side 1)
        if (i == nxL - 1)
          boundary.add_side(e, 1, interface_left_id);
      }
    }

    // --- Create right subdomain elements (QUAD4) ---
    for (unsigned j = 0; j < ny; ++j)
    {
      for (unsigned i = 0; i < nxR; ++i)
      {
        Elem *e = mesh.add_elem(Elem::build_with_id(QUAD4, eid++));
        e->set_node(0, mesh.node_ptr(nidR(i,   j)));
        e->set_node(1, mesh.node_ptr(nidR(i+1, j)));
        e->set_node(2, mesh.node_ptr(nidR(i+1, j+1)));
        e->set_node(3, mesh.node_ptr(nidR(i,   j+1)));

        // Interface (right side) boundary (i == 0 -> side 3)
        if (i == 0)
          boundary.add_side(e, 3, interface_right_id);

        // Right outer boundary (i == nxR - 1 -> side 1)
        if (i == nxR - 1)
          boundary.add_side(e, 1, right_id);
      }
    }

    // --- Assign human-readable boundary names ---
    boundary.sideset_name(left_id)            = "left_boundary";
    boundary.sideset_name(right_id)           = "right_boundary";
    boundary.sideset_name(interface_left_id)  = "interface_left";
    boundary.sideset_name(interface_right_id) = "interface_right";


    // This is the key testing step: inform libMesh about the disjoint neighbor boundary pairs
    // And, in `prepare_for_use()`, libMesh will set up the disjoint neighbor relationships.
    mesh.add_disjoint_neighbor_boundary_pairs(interface_left_id, interface_right_id, RealVectorValue(0.0, 0.0, 0.0));

    mesh.prepare_for_use();

#if LIBMESH_ENABLE_AMR
    if (test_local_refinement)
      {
      // set refinement flags
      for (auto & elem : mesh.element_ptr_range())
        {
            const Real xmid = elem->vertex_average()(0);
            if ((xmid < 0.5 && xmid > 0.5 - dx) ||
                (xmid > 0.5 && xmid < 0.5 + dx))
                elem->set_refinement_flag(Elem::REFINE);
        }

      // refine
      MeshRefinement(mesh).refine_elements();
      }
#endif
  }

void build_four_disjoint_elems(Mesh &mesh)
{
  // Clear any existing data and initialize the mesh
  mesh.clear();
  mesh.set_mesh_dimension(2);

  BoundaryInfo &boundary = mesh.get_boundary_info();

  // --------------------------------------------------------------------
  // Domain layout:
  //   4 elements aligned in a single row, each of width = 0.25
  //
  //   y=1:  x0---x1_L---x1_R---x2_L---x2_R---x3_L---x3_R---x4
  //
  //   y=0:  x0---x1_L---x1_R---x2_L---x2_R---x3_L---x3_R---x4
  //
  //   The interfaces between adjacent elements are duplicated:
  //   for each internal interface, left and right sides have distinct nodes
  //   (geometrically coincident but topologically disconnected).
  // --------------------------------------------------------------------
  const Real dx = 0.25;

  // ---- Add nodal coordinates ----
  // Each interface has duplicated nodes (left and right) to represent discontinuity.
  unsigned id = 0;
  for (unsigned i = 0; i < 4; ++i)
    {
      Real xL = i * dx;
      Real xR = (i + 1) * dx;

      // Left nodes of the element
      mesh.add_point(Point(xL, 0.0), id++); // bottom-left
      mesh.add_point(Point(xL, 1.0), id++); // top-left

      // Right nodes (duplicated interface nodes)
      mesh.add_point(Point(xR, 0.0), id++); // bottom-right
      mesh.add_point(Point(xR, 1.0), id++); // top-right
    }

  // Add one final pair of rightmost boundary nodes
  mesh.add_point(Point(1.0, 0.0), id++);
  mesh.add_point(Point(1.0, 1.0), id++);

  // ---- Create four QUAD4 elements ----
  // Node ordering: 0 = BL, 1 = BR, 2 = TR, 3 = TL
  unsigned elem_id = 0;
  for (unsigned i = 0; i < 4; ++i)
    {
      Elem *elem = mesh.add_elem(Elem::build_with_id(QUAD4, elem_id++));

      unsigned base = i * 4;
      elem->set_node(0, mesh.node_ptr(base));     // bottom-left
      elem->set_node(1, mesh.node_ptr(base + 2)); // bottom-right
      elem->set_node(2, mesh.node_ptr(base + 3)); // top-right
      elem->set_node(3, mesh.node_ptr(base + 1)); // top-left

      // Add outer boundaries
      if (i == 0)
        boundary.add_side(elem, 3, left_id);  // left outer boundary
      if (i == 3)
        boundary.add_side(elem, 1, right_id); // right outer boundary
    }

  // ---- Define internal interface sides ----
  // Element side 1 (right) connects to the next element’s side 3 (left)
  // Each interface is treated as a disjoint (non-conforming) boundary pair
  boundary.add_side(mesh.elem_ptr(0), 1, interface1_left_id);
  boundary.add_side(mesh.elem_ptr(1), 3, interface1_right_id);

  boundary.add_side(mesh.elem_ptr(1), 1, interface2_left_id);
  boundary.add_side(mesh.elem_ptr(2), 3, interface2_right_id);

  boundary.add_side(mesh.elem_ptr(2), 1, interface3_left_id);
  boundary.add_side(mesh.elem_ptr(3), 3, interface3_right_id);

  // ---- Register disjoint neighbor boundary pairs ----
  // These pairs inform libMesh that the specified sides are geometrically
  // adjacent but topologically disconnected.
  mesh.add_disjoint_neighbor_boundary_pairs(interface1_left_id, interface1_right_id, RealVectorValue(0.0, 0.0, 0.0));
  mesh.add_disjoint_neighbor_boundary_pairs(interface2_left_id, interface2_right_id, RealVectorValue(0.0, 0.0, 0.0));
  mesh.add_disjoint_neighbor_boundary_pairs(interface3_left_id, interface3_right_id, RealVectorValue(0.0, 0.0, 0.0));

  // ---- Assign descriptive names to boundary IDs ----
  boundary.sideset_name(left_id)  = "left_boundary";
  boundary.sideset_name(right_id) = "right_boundary";

  boundary.sideset_name(interface1_left_id)  = "interface1_left";
  boundary.sideset_name(interface1_right_id) = "interface1_right";
  boundary.sideset_name(interface2_left_id)  = "interface2_left";
  boundary.sideset_name(interface2_right_id) = "interface2_right";
  boundary.sideset_name(interface3_left_id)  = "interface3_left";
  boundary.sideset_name(interface3_right_id) = "interface3_right";

  // ---- Finalize mesh setup ----
  mesh.allow_renumbering(false);
  mesh.prepare_for_use();
}

};

CPPUNIT_TEST_SUITE_REGISTRATION( DisjointNeighborTest );
