#include <libmesh/distributed_mesh.h>
#include <libmesh/dof_map.h>
#include <libmesh/equation_systems.h>
#include <libmesh/linear_implicit_system.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_communication.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/numeric_vector.h>
#include <libmesh/replicated_mesh.h>
#include <libmesh/enum_norm_type.h>

#include <libmesh/dyna_io.h>
#include <libmesh/exodusII_io.h>
#include <libmesh/nemesis_io.h>

#include "test_comm.h"
#include "libmesh_cppunit.h"


using namespace libMesh;


Number six_x_plus_sixty_y (const Point& p,
                           const Parameters&,
                           const std::string&,
                           const std::string&)
{
  const Real & x = p(0);
  const Real & y = p(1);

  return 6*x + 60*y;
}


Number sin_x_plus_cos_y (const Point& p,
                         const Parameters&,
                         const std::string&,
                         const std::string&)
{
  const Real & x = p(0);
  const Real & y = p(1);

  return sin(x) + cos(y);
}


class MeshInputTest : public CppUnit::TestCase {
public:
  CPPUNIT_TEST_SUITE( MeshInputTest );

#if LIBMESH_DIM > 1
#ifdef LIBMESH_HAVE_EXODUS_API
  CPPUNIT_TEST( testExodusCopyNodalSolutionDistributed );
  CPPUNIT_TEST( testExodusCopyElementSolutionDistributed );
  CPPUNIT_TEST( testExodusCopyNodalSolutionReplicated );
  CPPUNIT_TEST( testExodusCopyElementSolutionReplicated );
  CPPUNIT_TEST( testExodusReadHeader );
#ifndef LIBMESH_USE_COMPLEX_NUMBERS
  // Eventually this will support complex numbers.
  CPPUNIT_TEST( testExodusWriteElementDataFromDiscontinuousNodalData );
#endif // LIBMESH_USE_COMPLEX_NUMBERS
#endif // LIBMESH_HAVE_EXODUS_API

#if defined(LIBMESH_HAVE_EXODUS_API) && defined(LIBMESH_HAVE_NEMESIS_API)
  CPPUNIT_TEST( testNemesisReadReplicated );
  CPPUNIT_TEST( testNemesisReadDistributed );

  CPPUNIT_TEST( testNemesisCopyNodalSolutionDistributed );
  CPPUNIT_TEST( testNemesisCopyNodalSolutionReplicated );
  CPPUNIT_TEST( testNemesisCopyElementSolutionDistributed );
  CPPUNIT_TEST( testNemesisCopyElementSolutionReplicated );
#endif

#ifdef LIBMESH_HAVE_GZSTREAM
  CPPUNIT_TEST( testDynaReadElem );
  CPPUNIT_TEST( testDynaReadPatch );
  CPPUNIT_TEST( testDynaFileMappingsFEMEx5);
  CPPUNIT_TEST( testDynaFileMappingsBlockWithHole);
  CPPUNIT_TEST( testDynaFileMappingsPlateWithHole);
  CPPUNIT_TEST( testDynaFileMappingsCyl3d);
#endif // LIBMESH_HAVE_GZSTREAM
#endif // LIBMESH_DIM > 1

  CPPUNIT_TEST_SUITE_END();

private:

public:
  void setUp()
  {}

  void tearDown()
  {}

#ifdef LIBMESH_HAVE_EXODUS_API
  void testExodusReadHeader ()
  {
    // first scope: write file
    {
      ReplicatedMesh mesh(*TestCommWorld);
      MeshTools::Generation::build_square (mesh, 3, 3, 0., 1., 0., 1.);
      ExodusII_IO exii(mesh);
      mesh.write("read_header_test.e");
    }

    // Make sure that the writing is done before the reading starts.
    TestCommWorld->barrier();

    // second scope: read header
    // Note: The header information is read from file on processor 0
    // and then broadcast to the other procs, so with this test we are
    // checking both that the header information is read correctly and
    // that it is correctly communicated to other procs.
    {
      ReplicatedMesh mesh(*TestCommWorld);
      ExodusII_IO exii(mesh);
      ExodusHeaderInfo header_info = exii.read_header("read_header_test.e");

      // Make sure the header information is as expected.
      CPPUNIT_ASSERT_EQUAL(std::string(header_info.title.data()), std::string("read_header_test.e"));
      CPPUNIT_ASSERT_EQUAL(header_info.num_dim, 2);
      CPPUNIT_ASSERT_EQUAL(header_info.num_elem, 9);
      CPPUNIT_ASSERT_EQUAL(header_info.num_elem_blk, 1);
      CPPUNIT_ASSERT_EQUAL(header_info.num_node_sets, 0);
      CPPUNIT_ASSERT_EQUAL(header_info.num_side_sets, 4);
      CPPUNIT_ASSERT_EQUAL(header_info.num_edge_blk, 0);
      CPPUNIT_ASSERT_EQUAL(header_info.num_edge, 0);
    }
  }


  template <typename MeshType, typename IOType>
  void testCopyNodalSolutionImpl (const std::string & filename)
  {
    {
      MeshType mesh(*TestCommWorld);

      EquationSystems es(mesh);
      System &sys = es.add_system<System> ("SimpleSystem");
      sys.add_variable("n", FIRST, LAGRANGE);

      MeshTools::Generation::build_square (mesh,
                                           3, 3,
                                           0., 1., 0., 1.);

      es.init();
      sys.project_solution(six_x_plus_sixty_y, nullptr, es.parameters);

      IOType meshoutput(mesh);

      meshoutput.write_equation_systems(filename, es);
    }

    {
      MeshType mesh(*TestCommWorld);
      IOType meshinput(mesh);

      // Avoid getting Nemesis solution values mixed up
      if (meshinput.is_parallel_format())
        {
          mesh.allow_renumbering(false);
          mesh.skip_noncritical_partitioning(true);
        }

      EquationSystems es(mesh);
      System &sys = es.add_system<System> ("SimpleSystem");
      sys.add_variable("testn", FIRST, LAGRANGE);

      if (mesh.processor_id() == 0 || meshinput.is_parallel_format())
        meshinput.read(filename);
      if (!meshinput.is_parallel_format())
        MeshCommunication().broadcast(mesh);
      mesh.prepare_for_use();

      es.init();

      // Read the solution e into variable teste.
      //
      // With complex numbers, we'll only bother reading the real
      // part.
#ifdef LIBMESH_USE_COMPLEX_NUMBERS
      meshinput.copy_nodal_solution(sys, "testn", "r_n");
#else
      meshinput.copy_nodal_solution(sys, "testn", "n");
#endif

      // Exodus only handles double precision
      Real exotol = std::max(TOLERANCE*TOLERANCE, Real(1e-12));

      for (Real x = 0; x < 1 + TOLERANCE; x += Real(1.L/3.L))
        for (Real y = 0; y < 1 + TOLERANCE; y += Real(1.L/3.L))
          {
            Point p(x,y);
            LIBMESH_ASSERT_FP_EQUAL(libmesh_real(sys.point_value(0,p)),
                                    libmesh_real(6*x+60*y),
                                    exotol);
          }
    }
  }


  void testExodusCopyNodalSolutionReplicated ()
  { testCopyNodalSolutionImpl<ReplicatedMesh,ExodusII_IO>("repl_with_nodal_soln.e"); }

  void testExodusCopyNodalSolutionDistributed ()
  { testCopyNodalSolutionImpl<DistributedMesh,ExodusII_IO>("dist_with_nodal_soln.e"); }

#if defined(LIBMESH_HAVE_NEMESIS_API)
  void testNemesisCopyNodalSolutionReplicated ()
  { testCopyNodalSolutionImpl<ReplicatedMesh,Nemesis_IO>("repl_with_nodal_soln.nem"); }

  void testNemesisCopyNodalSolutionDistributed ()
  { testCopyNodalSolutionImpl<DistributedMesh,Nemesis_IO>("dist_with_nodal_soln.nem"); }
#endif


  template <typename MeshType, typename IOType>
  void testCopyElementSolutionImpl (const std::string & filename)
  {
    {
      MeshType mesh(*TestCommWorld);

      EquationSystems es(mesh);
      System &sys = es.add_system<System> ("SimpleSystem");
      sys.add_variable("e", CONSTANT, MONOMIAL);

      MeshTools::Generation::build_square (mesh,
                                           3, 3,
                                           0., 1., 0., 1.);

      es.init();
      sys.project_solution(six_x_plus_sixty_y, nullptr, es.parameters);

      IOType meshinput(mesh);

      // Don't try to write element data as nodal data
      std::set<std::string> sys_list;
      meshinput.write_equation_systems(filename, es, &sys_list);

      // Just write it as element data
      meshinput.write_element_data(es);
    }

    {
      MeshType mesh(*TestCommWorld);
      IOType meshinput(mesh);

      // Avoid getting Nemesis solution values mixed up
      if (meshinput.is_parallel_format())
        {
          mesh.allow_renumbering(false);
          mesh.skip_noncritical_partitioning(true);
        }

      EquationSystems es(mesh);
      System &sys = es.add_system<System> ("SimpleSystem");
      sys.add_variable("teste", CONSTANT, MONOMIAL);

      if (mesh.processor_id() == 0 || meshinput.is_parallel_format())
        meshinput.read(filename);
      if (!meshinput.is_parallel_format())
        MeshCommunication().broadcast(mesh);
      mesh.prepare_for_use();

      es.init();

      // Read the solution e into variable teste.
      //
      // With complex numbers, we'll only bother reading the real
      // part.
#ifdef LIBMESH_USE_COMPLEX_NUMBERS
      meshinput.copy_elemental_solution(sys, "teste", "r_e");
#else
      meshinput.copy_elemental_solution(sys, "teste", "e");
#endif

      // Exodus only handles double precision
      Real exotol = std::max(TOLERANCE*TOLERANCE, Real(1e-12));

      for (Real x = Real(1.L/6.L); x < 1; x += Real(1.L/3.L))
        for (Real y = Real(1.L/6.L); y < 1; y += Real(1.L/3.L))
          {
            Point p(x,y);
            LIBMESH_ASSERT_FP_EQUAL(libmesh_real(sys.point_value(0,p)),
                                    libmesh_real(6*x+60*y),
                                    exotol);
          }
    }
  }


  void testExodusCopyElementSolutionReplicated ()
  { testCopyElementSolutionImpl<ReplicatedMesh,ExodusII_IO>("repl_with_elem_soln.e"); }

  void testExodusCopyElementSolutionDistributed ()
  { testCopyElementSolutionImpl<DistributedMesh,ExodusII_IO>("dist_with_elem_soln.e"); }

#if defined(LIBMESH_HAVE_NEMESIS_API)
  void testNemesisCopyElementSolutionReplicated ()
  { testCopyElementSolutionImpl<ReplicatedMesh,Nemesis_IO>("repl_with_elem_soln.nem"); }

  void testNemesisCopyElementSolutionDistributed ()
  { testCopyElementSolutionImpl<DistributedMesh,Nemesis_IO>("dist_with_elem_soln.nem"); }
#endif


#ifndef LIBMESH_USE_COMPLEX_NUMBERS
  void testExodusWriteElementDataFromDiscontinuousNodalData()
  {
    // first scope: write file
    {
      Mesh mesh(*TestCommWorld);

      EquationSystems es(mesh);
      System & sys = es.add_system<System> ("SimpleSystem");
      sys.add_variable("u", FIRST, L2_LAGRANGE);

      MeshTools::Generation::build_cube
        (mesh, 2, 2, 2, 0., 1., 0., 1., 0., 1., HEX8);

      es.init();

      // Set solution u^e_i = i, for the ith vertex of a given element e.
      const DofMap & dof_map = sys.get_dof_map();
      std::vector<dof_id_type> dof_indices;
      for (const auto & elem : mesh.element_ptr_range())
        {
          dof_map.dof_indices(elem, dof_indices, /*var_id=*/0);
          for (unsigned int i=0; i<dof_indices.size(); ++i)
            sys.solution->set(dof_indices[i], i);
        }
      sys.solution->close();

      // Now write to file.
      ExodusII_IO exii(mesh);

      // Don't try to write element data as averaged nodal data.
      std::set<std::string> sys_list;
      exii.write_equation_systems("elemental_from_nodal.e", es, &sys_list);

      // Write one elemental data field per vertex value.
      sys_list = {"SimpleSystem"};

      exii.write_element_data_from_discontinuous_nodal_data
        (es, &sys_list, /*var_suffix=*/"_elem_corner_");
    } // end first scope

    // second scope: read values back in, verify they are correct.
    {
      std::vector<std::string> file_var_names =
        {"u_elem_corner_0",
         "u_elem_corner_1",
         "u_elem_corner_2",
         "u_elem_corner_3"};
      std::vector<Real> expected_values = {0., 1., 2., 3.};

      // copy_elemental_solution currently requires ReplicatedMesh
      ReplicatedMesh mesh(*TestCommWorld);

      EquationSystems es(mesh);
      System & sys = es.add_system<System> ("SimpleSystem");
      for (auto i : index_range(file_var_names))
        sys.add_variable(file_var_names[i], CONSTANT, MONOMIAL);

      ExodusII_IO exii(mesh);

      if (mesh.processor_id() == 0)
        exii.read("elemental_from_nodal.e");
      MeshCommunication().broadcast(mesh);
      mesh.prepare_for_use();

      es.init();

      for (auto i : index_range(file_var_names))
        exii.copy_elemental_solution
          (sys, sys.variable_name(i), file_var_names[i]);

      // Check that the values we read back in are as expected.
      for (const auto & elem : mesh.active_element_ptr_range())
        for (auto i : index_range(file_var_names))
          {
            Real read_val = sys.point_value(i, elem->centroid());
            LIBMESH_ASSERT_FP_EQUAL
              (expected_values[i], read_val, TOLERANCE*TOLERANCE);
          }
    } // end second scope
  } // end testExodusWriteElementDataFromDiscontinuousNodalData

#endif // !LIBMESH_USE_COMPLEX_NUMBERS
#endif // LIBMESH_HAVE_EXODUS_API


#if defined(LIBMESH_HAVE_EXODUS_API) && defined(LIBMESH_HAVE_NEMESIS_API)
  template <typename MeshType>
  void testNemesisReadImpl ()
  {
    // first scope: write file
    {
      MeshType mesh(*TestCommWorld);
      MeshTools::Generation::build_square (mesh, 3, 3, 0., 1., 0., 1.);
      mesh.write("test_nemesis_read.nem");
    }

    // Make sure that the writing is done before the reading starts.
    TestCommWorld->barrier();

    // second scope: read file
    {
      MeshType mesh(*TestCommWorld);
      Nemesis_IO nem(mesh);

      nem.read("test_nemesis_read.nem");
      mesh.prepare_for_use();
      CPPUNIT_ASSERT_EQUAL(mesh.n_elem(),  dof_id_type(9));
      CPPUNIT_ASSERT_EQUAL(mesh.n_nodes(), dof_id_type(16));
    }
  }

  void testNemesisReadReplicated ()
  { testNemesisReadImpl<ReplicatedMesh>(); }

  void testNemesisReadDistributed ()
  { testNemesisReadImpl<DistributedMesh>(); }
#endif


  void testMasterCenters (const MeshBase & mesh)
  {
    auto locator = mesh.sub_point_locator();

    const std::set<subdomain_id_type> manifold_subdomain { 0 };
    const std::set<subdomain_id_type> nodeelem_subdomain { 1 };

    for (auto & elem : mesh.element_ptr_range())
      {
        Point master_pt = {}; // center, for tensor product elements

        // But perturb it to try and trigger any mapping weirdness
        if (elem->dim() > 0)
          master_pt(0) = 0.25;

        if (elem->dim() > 1)
          master_pt(1) = -0.25;

        if (elem->dim() > 2)
          master_pt(2) = 0.75;

        FEMap fe_map;

        Point physical_pt = fe_map.map(elem->dim(), elem, master_pt);

        Point inverse_pt = fe_map.inverse_map(elem->dim(), elem,
                                              physical_pt);

        CPPUNIT_ASSERT((inverse_pt-master_pt).norm() < TOLERANCE);

        CPPUNIT_ASSERT(elem->contains_point(physical_pt));

        const std::set<subdomain_id_type> * sbd_set =
          (elem->type() == NODEELEM) ?
          &nodeelem_subdomain : &manifold_subdomain;

        const Elem * located_elem = (*locator)(physical_pt, sbd_set);

        CPPUNIT_ASSERT(located_elem == elem);
      }
  }



  void testDynaReadElem ()
  {
    Mesh mesh(*TestCommWorld);

    DynaIO dyna(mesh);

    // Make DynaIO::add_spline_constraints work on DistributedMesh
    mesh.allow_renumbering(false);
    mesh.allow_remote_element_removal(false);

    if (mesh.processor_id() == 0)
      dyna.read("meshes/1_quad.bxt.gz");
    MeshCommunication().broadcast (mesh);

    mesh.prepare_for_use();

    // We have 1 QUAD9 finite element, attached via a trivial map to 9
    // spline Node+NodeElem objects
    CPPUNIT_ASSERT_EQUAL(mesh.n_elem(), dof_id_type(10));
    CPPUNIT_ASSERT_EQUAL(mesh.n_nodes(), dof_id_type(18));

    CPPUNIT_ASSERT_EQUAL(mesh.default_mapping_type(),
                         RATIONAL_BERNSTEIN_MAP);

    unsigned char weight_index = mesh.default_mapping_data();

    for (auto & elem : mesh.element_ptr_range())
      {
        if (elem->type() == NODEELEM)
          continue;

        CPPUNIT_ASSERT_EQUAL(elem->type(), QUAD9);
        for (unsigned int n=0; n != 9; ++n)
          CPPUNIT_ASSERT_EQUAL
            (elem->node_ref(n).get_extra_datum<Real>(weight_index),
             Real(0.75));

        CPPUNIT_ASSERT_EQUAL(elem->point(0)(0), Real(0.5));
        CPPUNIT_ASSERT_EQUAL(elem->point(0)(1), Real(0.5));
        CPPUNIT_ASSERT_EQUAL(elem->point(1)(0), Real(1.5));
        CPPUNIT_ASSERT_EQUAL(elem->point(1)(1), Real(0.5));
        CPPUNIT_ASSERT_EQUAL(elem->point(2)(0), Real(1.5));
        CPPUNIT_ASSERT_EQUAL(elem->point(2)(1), Real(1.5));
        CPPUNIT_ASSERT_EQUAL(elem->point(3)(0), Real(0.5));
        CPPUNIT_ASSERT_EQUAL(elem->point(3)(1), Real(1.5));
        CPPUNIT_ASSERT(elem->has_affine_map());
#if LIBMESH_DIM > 2
        for (unsigned int v=0; v != 4; ++v)
          CPPUNIT_ASSERT_EQUAL(elem->point(v)(2), Real(0));
#endif
      }

    testMasterCenters(mesh);
  }


  void testDynaReadPatch ()
  {
    Mesh mesh(*TestCommWorld);

    // Make DynaIO::add_spline_constraints work on DistributedMesh
    mesh.allow_renumbering(false);
    mesh.allow_remote_element_removal(false);

    DynaIO dyna(mesh);
    if (mesh.processor_id() == 0)
      dyna.read("meshes/25_quad.bxt.gz");
    MeshCommunication().broadcast (mesh);

    mesh.prepare_for_use();

    // We have 5^2 QUAD9 elements, with 11^2 nodes,
    // tied to 49 Node/NodeElem spline nodes
    CPPUNIT_ASSERT_EQUAL(mesh.n_elem(), dof_id_type(25+49));
    CPPUNIT_ASSERT_EQUAL(mesh.n_nodes(), dof_id_type(121+49));

    CPPUNIT_ASSERT_EQUAL(mesh.default_mapping_type(),
                         RATIONAL_BERNSTEIN_MAP);

    unsigned char weight_index = mesh.default_mapping_data();

    for (const auto & elem : mesh.active_element_ptr_range())
      {
        if (elem->type() == NODEELEM)
          continue;
        LIBMESH_ASSERT_FP_EQUAL(libmesh_real(0.04), elem->volume(), TOLERANCE);

        for (unsigned int n=0; n != 9; ++n)
          CPPUNIT_ASSERT_EQUAL
            (elem->node_ref(n).get_extra_datum<Real>(weight_index),
             Real(1.0));

        unsigned int n_neighbors = 0, n_neighbors_expected = 2;
        for (unsigned int side=0; side != 4; ++side)
          if (elem->neighbor_ptr(side))
            n_neighbors++;
        Point c = elem->centroid();

        if (c(0) > 0.2 && c(0) < 0.8)
          n_neighbors_expected++;
        if (c(1) > 0.2 && c(1) < 0.8)
          n_neighbors_expected++;

        CPPUNIT_ASSERT_EQUAL(n_neighbors, n_neighbors_expected);
      }

    testMasterCenters(mesh);

#ifdef LIBMESH_HAVE_SOLVER
#ifdef LIBMESH_ENABLE_CONSTRAINTS
    // Now test whether we can assign the desired constraint equations
    EquationSystems es(mesh);
    System & sys = es.add_system<LinearImplicitSystem>("test");
    sys.add_variable("u", SECOND); // to match QUAD9
    es.init();
    dyna.add_spline_constraints(sys.get_dof_map(), 0, 0);

    // We should have a constraint on every FE dof
    CPPUNIT_ASSERT_EQUAL(sys.get_dof_map().n_constrained_dofs(), dof_id_type(121));
#endif // LIBMESH_ENABLE_CONSTRAINTS
#endif // LIBMESH_HAVE_SOLVER
  }

  void testProjectionRegression(MeshBase & mesh, DynaIO & dyna, std::array<Real, 4> expected_norms)
  {
    int order = 0;
    for (const auto elem : mesh.element_ptr_range())
      order = std::max(order, int(elem->default_order()));
    TestCommWorld->max(order);
    CPPUNIT_ASSERT (order > 0);

    EquationSystems es(mesh);
    System &sys = es.add_system<System> ("SimpleSystem");
    unsigned int n_var =
      sys.add_variable("n", Order(order), RATIONAL_BERNSTEIN);

    es.init();
    dyna.add_spline_constraints(sys.get_dof_map(), sys.number(), n_var);

    sys.project_solution(sin_x_plus_cos_y, nullptr, es.parameters);

    // Calculate some norms, skipping the spline points, and compare
    // to regression standard values
    std::set<unsigned int> skip_dimensions {0};
    const Real L2_norm =
      sys.calculate_norm(*sys.solution, 0, L2, &skip_dimensions);
//    std::cout.precision(16);
//    std::cout << "L2_norm = " << L2_norm << std::endl;
    LIBMESH_ASSERT_FP_EQUAL(L2_norm, expected_norms[0], TOLERANCE);
    const Real Linf_norm =
      sys.calculate_norm(*sys.solution, 0, L_INF, &skip_dimensions);
//    std::cout << "Linf_norm = " << Linf_norm << std::endl;
    LIBMESH_ASSERT_FP_EQUAL(Linf_norm, expected_norms[1], TOLERANCE);
    const Real H1_norm =
      sys.calculate_norm(*sys.solution, 0, H1_SEMINORM, &skip_dimensions);
//    std::cout << "H1_norm = " << H1_norm << std::endl;
    LIBMESH_ASSERT_FP_EQUAL(H1_norm, expected_norms[2], TOLERANCE);
    const Real W1inf_norm =
      sys.calculate_norm(*sys.solution, 0, W1_INF_SEMINORM, &skip_dimensions);
//    std::cout << "W1inf_norm = " << W1inf_norm << std::endl;
    // W1_inf seems more sensitive to FP error...
    LIBMESH_ASSERT_FP_EQUAL(W1inf_norm, expected_norms[3], 10*TOLERANCE);
  }

  void testDynaFileMappings (const std::string & filename, std::array<Real, 4> expected_norms)
  {
    Mesh mesh(*TestCommWorld);

    // Make DynaIO::add_spline_constraints work on DistributedMesh
    mesh.allow_renumbering(false);
    mesh.allow_remote_element_removal(false);

    DynaIO dyna(mesh);
    if (mesh.processor_id() == 0)
      dyna.read(filename);
    MeshCommunication().broadcast (mesh);

    mesh.prepare_for_use();

    CPPUNIT_ASSERT_EQUAL(mesh.default_mapping_type(),
                         RATIONAL_BERNSTEIN_MAP);

    // Useful when trying out different projection functions
    // std::cout << filename << ":" << std::endl;

    testMasterCenters(mesh);

    testProjectionRegression(mesh, dyna, expected_norms);
  }

  void testDynaFileMappingsFEMEx5 ()
  {
    testDynaFileMappings("meshes/PressurizedCyl_Patch6_256Elem.bxt.gz",
    // Regression values for sin_x_plus_cos_y
                         {0.9639857809698268, 1.839870171669186,
                          0.7089812562241862, 1.306121188539059});
  }

  void testDynaFileMappingsBlockWithHole ()
  {
    testDynaFileMappings("meshes/BlockWithHole_Patch9.bxt.gz",
    // Regression values for sin_x_plus_cos_y
                         {3.226125496262302, 1.97405596521291,
                          2.533759662135491, 1.413785069495184});
  }

  void testDynaFileMappingsPlateWithHole ()
  {
    testDynaFileMappings("meshes/PlateWithHole_Patch8.bxt.gz",
    // Regression values for sin_x_plus_cos_y
                         {2.2812154374012, 1.974049990211937,
                          1.791640772215248, 1.413679237529376});
  }

  void testDynaFileMappingsCyl3d ()
  {
    testDynaFileMappings("meshes/PressurizedCyl3d_Patch1_8Elem.bxt.gz",
    // Regression values for sin_x_plus_cos_y
                         {0.9636130896326653, 1.823294442918401,
                          0.7080084233124895, 1.314114853940283});
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION( MeshInputTest );
