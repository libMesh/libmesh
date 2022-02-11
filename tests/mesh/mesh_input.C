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
#include <libmesh/vtk_io.h>

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
#if LIBMESH_DIM > 2
  CPPUNIT_TEST( testExodusIGASidesets );
#endif
#ifndef LIBMESH_USE_COMPLEX_NUMBERS
  CPPUNIT_TEST( testExodusCopyElementVectorDistributed );
  CPPUNIT_TEST( testExodusCopyElementVectorReplicated );

  // Eventually this will support complex numbers.
  CPPUNIT_TEST( testExodusWriteElementDataFromDiscontinuousNodalData );
#endif // !LIBMESH_USE_COMPLEX_NUMBERS

  CPPUNIT_TEST( testExodusFileMappingsPlateWithHole);
  CPPUNIT_TEST( testExodusFileMappingsTwoBlocks);
  CPPUNIT_TEST( testExodusFileMappingsCyl3d);
#endif // LIBMESH_HAVE_EXODUS_API

#if defined(LIBMESH_HAVE_EXODUS_API) && defined(LIBMESH_HAVE_NEMESIS_API)
  CPPUNIT_TEST( testNemesisReadReplicated );
  CPPUNIT_TEST( testNemesisReadDistributed );

  CPPUNIT_TEST( testNemesisCopyNodalSolutionDistributed );
  CPPUNIT_TEST( testNemesisCopyNodalSolutionReplicated );
  CPPUNIT_TEST( testNemesisCopyElementSolutionDistributed );
  CPPUNIT_TEST( testNemesisCopyElementSolutionReplicated );

  CPPUNIT_TEST( testNemesisSingleElementDistributed );
  CPPUNIT_TEST( testNemesisSingleElementReplicated );
#ifndef LIBMESH_USE_COMPLEX_NUMBERS
  CPPUNIT_TEST( testNemesisCopyElementVectorDistributed );
  CPPUNIT_TEST( testNemesisCopyElementVectorReplicated );
#endif // !LIBMESH_USE_COMPLEX_NUMBERS
#endif // defined(LIBMESH_HAVE_EXODUS_API) && defined(LIBMESH_HAVE_NEMESIS_API)

#ifdef LIBMESH_HAVE_GZSTREAM
  CPPUNIT_TEST( testDynaReadElem );
  CPPUNIT_TEST( testDynaNoSplines );
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
      CPPUNIT_ASSERT_EQUAL(header_info.num_node_sets, 4);
      CPPUNIT_ASSERT_EQUAL(header_info.num_side_sets, 4);
      CPPUNIT_ASSERT_EQUAL(header_info.num_edge_blk, 0);
      CPPUNIT_ASSERT_EQUAL(header_info.num_edge, 0);
    }
  }


  void testExodusIGASidesets ()
  {
    Mesh mesh(*TestCommWorld);

    // Block here so we trigger exii destructor early; I thought I
    // might have had a bug in there at one point
    {
      ExodusII_IO exii(mesh);
      // IGA Exodus meshes require ExodusII 8 or higher
      if (exii.get_exodus_version() < 800)
        return;

      if (mesh.processor_id() == 0)
        exii.read("meshes/Cube_With_Sidesets.e");

    }

    MeshCommunication().broadcast(mesh);
    mesh.prepare_for_use();

    // 5^3 spline nodes + 7^3 Rational Bezier nodes
    CPPUNIT_ASSERT_EQUAL(mesh.n_nodes(), dof_id_type(468));
    // 5^3 spline elements + 3^3 Rational Bezier elements
    CPPUNIT_ASSERT_EQUAL(mesh.n_elem(),  dof_id_type(152));

    // Check that we see the boundary ids we expect
    BoundaryInfo & bi = mesh.get_boundary_info();

    // On a ReplicatedMesh, we should see all 6 boundary ids on each processor
    if (mesh.is_serial())
      CPPUNIT_ASSERT_EQUAL(static_cast<std::size_t>(6), bi.n_boundary_ids());

    // On any mesh, we should see each id on *some* processor
    {
      const std::set<boundary_id_type> & bc_ids = bi.get_boundary_ids();
      // CoreForm gave me a file with 1-based numbering! (faints)
      for (boundary_id_type i = 1 ; i != 7; ++i)
        {
          bool has_bcid = bc_ids.count(i);
          mesh.comm().max(has_bcid);
          CPPUNIT_ASSERT(has_bcid);
        }
    }

    // Indexed by bcid-1, because we count from 0, like God and
    // Dijkstra intended!
    std::vector<int> side_counts(6, 0);

    // Map from our side numbers to the file's BCIDs
    const boundary_id_type libmesh_side_to_bcid[] = {1, 4, 6, 3, 5, 2};

    for (const auto & elem : mesh.active_local_element_ptr_range())
      {
        if (elem->type() == NODEELEM)
          continue;

        for (unsigned short side=0; side<elem->n_sides(); side++)
          {
            if (elem->neighbor_ptr(side))
              CPPUNIT_ASSERT_EQUAL(bi.n_boundary_ids(elem, side), 0u);
            else
              {
                CPPUNIT_ASSERT_EQUAL(bi.n_boundary_ids(elem, side), 1u);
                std::vector<boundary_id_type> bids;
                bi.boundary_ids(elem, side, bids);
                side_counts[bids[0]-1]++;
                CPPUNIT_ASSERT_EQUAL(libmesh_side_to_bcid[side], bids[0]);
              }
          }
      }

    for (auto bc_count : side_counts)
      {
        // We should have 3^2 sides with each id
        mesh.comm().sum(bc_count);
        CPPUNIT_ASSERT_EQUAL(bc_count, 9);
      }

    // Test a write when we're done reading; I was getting weirdness
    // from NetCDF at this point in a Moose output test.
    {
      ExodusII_IO exii(mesh);

      exii.write("Cube_With_Sidesets_out.e");
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


  // This tests that a single-element mesh solution makes it through all of the API to write out a
  // set of Nemesis files. When this executes in parallel, there are ((number of processors) - 1)
  // subdomains with zero elements and therefore nothing to write, but the API should handle this.
  template <typename MeshType, typename IOType>
  void testSingleElementImpl(const std::string & filename)
  {
    {
      // Generate a single 1x1 square element mesh
      MeshType mesh(*TestCommWorld);
      MeshTools::Generation::build_square(mesh, 1, 1);

      EquationSystems es(mesh);
      auto & sys = es.add_system<System>("SimpleSystem");
      sys.add_variable("e", CONSTANT, MONOMIAL);

      // Set an arbitrary solution for the single element DOF
      es.init();
      sys.project_solution(six_x_plus_sixty_y, nullptr, es.parameters);

      // Write the solution to Nemesis file(s) - only proc 0 should have anything to write!
      Nemesis_IO nem_io(mesh);
      std::set<std::string> sys_list;
      nem_io.write_equation_systems(filename, es, &sys_list);
      nem_io.write_element_data(es);
    }

    // If we can read and copy the correct element value back into the mesh, we know the Nemesis
    // file(s) were written properly.
    {
      MeshType mesh(*TestCommWorld);
      EquationSystems es(mesh);
      auto & sys = es.add_system<System>("SimpleSystem");
      sys.add_variable("teste", CONSTANT, MONOMIAL);

      Nemesis_IO nem_io(mesh);
      if (mesh.processor_id() == 0 || nem_io.is_parallel_format())
        nem_io.read(filename);

      mesh.prepare_for_use();
      es.init();

#ifdef LIBMESH_USE_COMPLEX_NUMBERS
      nem_io.copy_elemental_solution(sys, "teste", "r_e");
#else
      nem_io.copy_elemental_solution(sys, "teste", "e");
#endif

      // The result should be '\frac{6 + 60}{2} = 33' at all points in the element domain
      CPPUNIT_ASSERT_EQUAL(int(sys.solution->size()), 1);
      CPPUNIT_ASSERT_EQUAL(libmesh_real(sys.point_value(0, Point(0.5, 0.5))), Real(33));
    }
  }


  void testNemesisSingleElementReplicated ()
  { testSingleElementImpl<ReplicatedMesh,Nemesis_IO>("repl_with_single_elem.nem"); }

  void testNemesisSingleElementDistributed ()
  { testSingleElementImpl<DistributedMesh,Nemesis_IO>("dist_with_single_elem.nem"); }
#endif //defined(LIBMESH_HAVE_NEMESIS_API)


#ifndef LIBMESH_USE_COMPLEX_NUMBERS
  // So this tester runs through pretty much the same process as 'testCopyElementSolutionImpl()'
  // except for a CONSTANT MONOMIAL_VEC variable. It mainly serves as a test for writing elemental
  // vector variables to elements in ExodusII files, and is not actually a test for reading and
  // copying an elemental solution.
  template <typename MeshType, typename IOType>
  void testCopyElementVectorImpl(const std::string & filename)
  {
    {
      MeshType mesh(*TestCommWorld);

      EquationSystems es(mesh);
      System & sys = es.add_system<System> ("SimpleSystem");
      auto e_var = sys.add_variable("e", CONSTANT, MONOMIAL_VEC);

      MeshTools::Generation::build_square (mesh,
                                           3, 3,
                                           0., 1., 0., 1.);

      es.init();

      // Here, we're going to manually set up the solution because the 'project_solution()' and
      // 'project_vector()' methods don't work so well with CONSTANT MONOMIAL_VEC variables. They
      // each lead to an error downstream asserting positive-definiteness when Cholesky decomposing.
      // Interestingly, the error is only invoked for CONSTANT MONOMIAL_VEC, and not, e.g.,
      // CONSTANT MONOMIAL nor FIRST LAGRANGE_VEC.
      //
      // Anyways, the important thing here is that we test the ExodusII and Nemesis writers, how
      // the solution is set is hardly important, and we're pretty much following the same
      // philosophy as the 'test2DProjectVectorFE()' unit tester in 'systems_test.C'
      Parameters params;
      for (const auto & elem : mesh.active_local_element_ptr_range())
      {
        const Point & p = elem->vertex_average();

        // Set the x-component with the value from 'six_x_plus_sixty_y()' and the y-component
        // with that from 'sin_x_plus_cos_y()' at the element centroid (vertex average)
        sys.current_local_solution->set(
          elem->dof_number(sys.number(), e_var, 0), six_x_plus_sixty_y(p, params, "", ""));
        sys.current_local_solution->set(
          elem->dof_number(sys.number(), e_var, 1), sin_x_plus_cos_y(p, params, "", ""));
      }

      // After setting values, we need to assemble
      sys.current_local_solution->close();

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
      System & sys = es.add_system<System> ("SimpleSystem");

      // We have to read the CONSTANT MONOMIAL_VEC var "e" into separate CONSTANT MONOMIAL vars
      // "e_x" and "e_y" because 'copy_elemental_solution()' currently doesn't support vectors.
      // Again, this isn't a test for reading/copying an elemental vector solution, only writing.
      sys.add_variable("teste_x", CONSTANT, MONOMIAL);
      sys.add_variable("teste_y", CONSTANT, MONOMIAL);

      if (mesh.processor_id() == 0 || meshinput.is_parallel_format())
        meshinput.read(filename);
      if (!meshinput.is_parallel_format())
        MeshCommunication().broadcast(mesh);
      mesh.prepare_for_use();

      es.init();

      // Read the solution e_x and e_y into variable teste_x and teste_y, respectively.
      meshinput.copy_elemental_solution(sys, "teste_x", "e_x");
      meshinput.copy_elemental_solution(sys, "teste_y", "e_y");

      // Exodus only handles double precision
      Real exotol = std::max(TOLERANCE*TOLERANCE, Real(1e-12));

      for (Real x = Real(1.L/6.L); x < 1; x += Real(1.L/3.L))
        for (Real y = Real(1.L/6.L); y < 1; y += Real(1.L/3.L))
          {
            Point p(x,y);
            LIBMESH_ASSERT_FP_EQUAL(libmesh_real(sys.point_value(0,p)),
                                    libmesh_real(6*x+60*y),
                                    exotol);
            LIBMESH_ASSERT_FP_EQUAL(libmesh_real(sys.point_value(1,p)),
                                    libmesh_real(sin(x)+cos(y)),
                                    exotol);
          }
    }
  }

  void testExodusCopyElementVectorReplicated ()
  { testCopyElementVectorImpl<ReplicatedMesh, ExodusII_IO>("repl_with_elem_vec.e"); }

  void testExodusCopyElementVectorDistributed ()
  { testCopyElementVectorImpl<DistributedMesh,ExodusII_IO>("dist_with_elem_vec.e"); }

#if defined(LIBMESH_HAVE_NEMESIS_API)
  void testNemesisCopyElementVectorReplicated ()
  { testCopyElementVectorImpl<ReplicatedMesh,Nemesis_IO>("repl_with_elem_vec.nem"); }

  void testNemesisCopyElementVectorDistributed ()
  { testCopyElementVectorImpl<DistributedMesh,Nemesis_IO>("dist_with_elem_vec.nem"); }
#endif


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
            Real read_val = sys.point_value(i, elem->vertex_average());
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

        std::set<subdomain_id_type> my_subdomain { elem->subdomain_id() };

        const Elem * located_elem = (*locator)(physical_pt, &my_subdomain);

        CPPUNIT_ASSERT(located_elem == elem);
      }
  }



  void helperTestingDynaQuad (const MeshBase & mesh)
  {
    CPPUNIT_ASSERT_EQUAL(mesh.default_mapping_type(),
                         RATIONAL_BERNSTEIN_MAP);

    unsigned char weight_index = mesh.default_mapping_data();

    bool found_the_quad = false;

    for (auto & elem : mesh.element_ptr_range())
      {
        if (elem->type() == NODEELEM)
          continue;

        CPPUNIT_ASSERT_EQUAL(elem->type(), QUAD9);
        found_the_quad = true;

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

    TestCommWorld->max(found_the_quad);
    CPPUNIT_ASSERT(found_the_quad);

    testMasterCenters(mesh);
  }


  void testDynaReadElem ()
  {
    Mesh mesh(*TestCommWorld);

    DynaIO dyna(mesh);

    if (mesh.processor_id() == 0)
      dyna.read("meshes/1_quad.bxt.gz");
    MeshCommunication().broadcast (mesh);

    mesh.prepare_for_use();

    // We have 1 QUAD9 finite element, attached via a trivial map to 9
    // spline Node+NodeElem objects
    CPPUNIT_ASSERT_EQUAL(mesh.n_elem(), dof_id_type(10));
    CPPUNIT_ASSERT_EQUAL(mesh.n_nodes(), dof_id_type(18));

    helperTestingDynaQuad(mesh);
  }


  void testDynaNoSplines ()
  {
    Mesh mesh(*TestCommWorld);

    DynaIO dyna(mesh, /* keep_spline_nodes = */ false);

    if (mesh.processor_id() == 0)
      dyna.read("meshes/1_quad.bxt.gz");
    MeshCommunication().broadcast (mesh);

    mesh.prepare_for_use();

    // We have 1 QUAD9 finite element
    CPPUNIT_ASSERT_EQUAL(mesh.n_elem(), dof_id_type(1));
    CPPUNIT_ASSERT_EQUAL(mesh.n_nodes(), dof_id_type(9));

    helperTestingDynaQuad(mesh);
  }


  void testDynaReadPatch ()
  {
    Mesh mesh(*TestCommWorld);

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
        Point c = elem->vertex_average();

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

    // We should have a constraint on every FE dof
    CPPUNIT_ASSERT_EQUAL(sys.get_dof_map().n_constrained_dofs(), dof_id_type(121));
#endif // LIBMESH_ENABLE_CONSTRAINTS
#endif // LIBMESH_HAVE_SOLVER
  }

  void testProjectionRegression(MeshBase & mesh, std::array<Real, 4> expected_norms)
  {
    int order = 0;
    for (const auto elem : mesh.element_ptr_range())
      order = std::max(order, int(elem->default_order()));
    TestCommWorld->max(order);
    CPPUNIT_ASSERT (order > 0);

    EquationSystems es(mesh);
    System &sys = es.add_system<System> ("SimpleSystem");
    sys.add_variable("n", Order(order), RATIONAL_BERNSTEIN);

    es.init();

    sys.project_solution(sin_x_plus_cos_y, nullptr, es.parameters);

    // Make this easy to tweak in the future
    const double my_tolerance = TOLERANCE;

    // Calculate some norms, skipping the spline points, and compare
    // to regression standard values
    std::set<unsigned int> skip_dimensions {0};
    const Real L2_norm =
      sys.calculate_norm(*sys.solution, 0, L2, &skip_dimensions);
//    std::cout.precision(16);
//    std::cout << "L2_norm = " << L2_norm << std::endl;
    LIBMESH_ASSERT_FP_EQUAL(L2_norm, expected_norms[0], my_tolerance);
    const Real Linf_norm =
      sys.calculate_norm(*sys.solution, 0, L_INF, &skip_dimensions);
//    std::cout << "Linf_norm = " << Linf_norm << std::endl;
    LIBMESH_ASSERT_FP_EQUAL(Linf_norm, expected_norms[1], my_tolerance);
    const Real H1_norm =
      sys.calculate_norm(*sys.solution, 0, H1_SEMINORM, &skip_dimensions);
//    std::cout << "H1_norm = " << H1_norm << std::endl;
    LIBMESH_ASSERT_FP_EQUAL(H1_norm, expected_norms[2], my_tolerance);
    const Real W1inf_norm =
      sys.calculate_norm(*sys.solution, 0, W1_INF_SEMINORM, &skip_dimensions);
//    std::cout << "W1inf_norm = " << W1inf_norm << std::endl;
    // W1_inf seems more sensitive to FP error...
    LIBMESH_ASSERT_FP_EQUAL(W1inf_norm, expected_norms[3], 10*my_tolerance);
  }

  void testDynaFileMappings (const std::string & filename, std::array<Real, 4> expected_norms)
  {
    Mesh mesh(*TestCommWorld);

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

    testProjectionRegression(mesh, expected_norms);
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
                         {3.22612556930183, 1.97405365384733,
                          2.53376235803176, 1.41374070517223});
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
                         {0.963612880188165, 1.82329452603503,
                          0.707998701597943, 1.31399222566683});
  }

  void testExodusFileMappings (const std::string & filename, std::array<Real, 4> expected_norms)
  {
    Mesh mesh(*TestCommWorld);

    ExodusII_IO exii(mesh);
    // IGA Exodus meshes require ExodusII 8 or higher
    if (exii.get_exodus_version() < 800)
      return;

    if (mesh.processor_id() == 0)
      exii.read(filename);
    MeshCommunication().broadcast (mesh);

    mesh.prepare_for_use();

    CPPUNIT_ASSERT_EQUAL(mesh.default_mapping_type(),
                         RATIONAL_BERNSTEIN_MAP);

    testMasterCenters(mesh);

    testProjectionRegression(mesh, expected_norms);

    // Test a write when we're done reading; I was getting weirdness
    // from NetCDF at this point in a Moose output test.
    {
      ExodusII_IO exii(mesh);

      exii.write("exodus_file_mapping_out.e");
    }

    {
      VTKIO vtkout(mesh);

      vtkout.write("vtk_file_mapping_out.pvtu");
    }
  }

  void testExodusFileMappingsPlateWithHole ()
  {
    testExodusFileMappings("meshes/PlateWithHole_Patch8.e",
    // Regression values for sin_x_plus_cos_y
                           {2.2812154374012, 1.974049990211937,
                            1.791640772215248, 1.413679237529376});
  }

  void testExodusFileMappingsTwoBlocks ()
  {
    testExodusFileMappings("meshes/two_quads_two_blocks.e",
    // Regression values for sin_x_plus_cos_y
                           {2.03496953073072, 1.97996853164955,
                            1.18462134113435, 1.03085301158959});
  }

  void testExodusFileMappingsCyl3d ()
  {
    testExodusFileMappings("meshes/PressurizedCyl3d_Patch1_8Elem.e",
                           {0.963612880188165, 1.82329452603503,
                            0.707998701597943, 1.31399222566683});
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION( MeshInputTest );
