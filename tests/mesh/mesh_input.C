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
#include <libmesh/enum_to_string.h>

#include <libmesh/abaqus_io.h>
#include <libmesh/dyna_io.h>
#include <libmesh/exodusII_io.h>
#include <libmesh/gmsh_io.h>
#include <libmesh/nemesis_io.h>
#include <libmesh/stl_io.h>
#include <libmesh/vtk_io.h>
#include <libmesh/tetgen_io.h>

#include "test_comm.h"
#include "libmesh_cppunit.h"

#include <regex>

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

constexpr int added_sides_nxyz[] = {2,2,2};

Number designed_for_side_elems (const Point & p,
                                const Parameters & param,
                                const std::string &,
                                const std::string &)
{
  const Real & x = p(0);
  const Real & y = p(1);
  const Real & z = p(2);

  short facedim = param.have_parameter<short>("face") ?
    param.get<short>("face") : -1;

  // What face are we on?
  auto is_on_face = [facedim](Real r, short rdim) {
    if (facedim == rdim)
      return true;
    if (facedim >= 0)
      return false;
    const Real numerator = r * added_sides_nxyz[rdim];
    return (std::abs(numerator - std::round(numerator)) <
            TOLERANCE*TOLERANCE);
  };



  // x/y/z components of a div-free flux,
  // curl([x^2yz, xy^2z, xyz])
  if (is_on_face(x, 0))
    {
      libmesh_assert(!is_on_face(y, 1));
      libmesh_assert(!is_on_face(z, 2));
      return (x*z-x*y*y);
    }
  if (is_on_face(y, 1))
    {
      libmesh_assert(!is_on_face(z, 2));
      return (x*x*y-y*z);
    }

  libmesh_assert(is_on_face(z, 2));
  return (y*y*z-x*x*z);
}


class MeshInputTest : public CppUnit::TestCase {
public:
  LIBMESH_CPPUNIT_TEST_SUITE( MeshInputTest );

#if LIBMESH_DIM > 1
#ifdef LIBMESH_HAVE_VTK
  CPPUNIT_TEST( testVTKPreserveElemIds );
  CPPUNIT_TEST( testVTKPreserveSubdomainIds );
#endif

#ifdef LIBMESH_HAVE_EXODUS_API
  CPPUNIT_TEST( testExodusCopyNodalSolutionDistributed );
  CPPUNIT_TEST( testExodusCopyElementSolutionDistributed );
  CPPUNIT_TEST( testExodusCopyNodalSolutionReplicated );
  CPPUNIT_TEST( testExodusCopyElementSolutionReplicated );
  CPPUNIT_TEST( testExodusReadHeader );
#if LIBMESH_DIM > 2
  CPPUNIT_TEST( testExodusIGASidesets );
  CPPUNIT_TEST( testLowOrderEdgeBlocks );
#endif
#ifndef LIBMESH_USE_COMPLEX_NUMBERS
  CPPUNIT_TEST( testExodusCopyElementVectorDistributed );
  CPPUNIT_TEST( testExodusCopyElementVectorReplicated );

  // Eventually this will support complex numbers.
  CPPUNIT_TEST( testExodusWriteElementDataFromDiscontinuousNodalData );
#endif // !LIBMESH_USE_COMPLEX_NUMBERS

  CPPUNIT_TEST( testExodusWriteAddedSidesEdgeC0 );
  CPPUNIT_TEST( testExodusDiscWriteAddedSidesEdgeC0 );
  CPPUNIT_TEST( testExodusWriteAddedSidesMixedEdgeC0 );
  CPPUNIT_TEST( testExodusDiscWriteAddedSidesMixedEdgeC0 );
  // CPPUNIT_TEST( testExodusDiscWriteAddedSidesEdgeDisc ); // need is_on_face fixes
  // CPPUNIT_TEST( testExodusDiscWriteAddedSidesEdgeDisc ); // need is_on_face fixes
  CPPUNIT_TEST( testExodusWriteAddedSidesTriC0 );
  CPPUNIT_TEST( testExodusDiscWriteAddedSidesTriC0 );
  CPPUNIT_TEST( testExodusWriteAddedSidesMixedTriC0 );
  CPPUNIT_TEST( testExodusDiscWriteAddedSidesMixedTriC0 );
  // CPPUNIT_TEST( testExodusWriteAddedSidesTriDisc ); // Need aligned faces
  // CPPUNIT_TEST( testExodusDiscWriteAddedSidesTriDisc ); // Need aligned faces
  CPPUNIT_TEST( testExodusWriteAddedSidesQuadC0 );
  CPPUNIT_TEST( testExodusDiscWriteAddedSidesQuadC0 );
  CPPUNIT_TEST( testExodusWriteAddedSidesMixedQuadC0 );
  CPPUNIT_TEST( testExodusDiscWriteAddedSidesMixedQuadC0 );
  // CPPUNIT_TEST( testExodusWriteAddedSidesQuadDisc ); // need is_on_face fixes
  // CPPUNIT_TEST( testExodusDiscWriteAddedSidesQuadDisc ); // need is_on_face fixes
  // CPPUNIT_TEST( testExodusWriteAddedSidesTetC0 ); // BROKEN!?!  WHY!?!
  // CPPUNIT_TEST( testExodusDiscWriteAddedSidesTetC0 ); // BROKEN!?!  WHY!?!
  // CPPUNIT_TEST( testExodusWriteAddedSidesTetDisc );
  // CPPUNIT_TEST( testExodusDiscWriteAddedSidesTetDisc );
  CPPUNIT_TEST( testExodusWriteAddedSidesHexC0 );
  CPPUNIT_TEST( testExodusDiscWriteAddedSidesHexC0 );
  CPPUNIT_TEST( testExodusWriteAddedSidesMixedHexC0 );
  CPPUNIT_TEST( testExodusDiscWriteAddedSidesMixedHexC0 );
  CPPUNIT_TEST( testExodusWriteAddedSidesHexDisc );
  CPPUNIT_TEST( testExodusDiscWriteAddedSidesHexDisc );

  CPPUNIT_TEST( testExodusFileMappingsPlateWithHole);
  CPPUNIT_TEST( testExodusFileMappingsTwoBlocks);
  CPPUNIT_TEST( testExodusFileMappingsTwoElemIGA);
  CPPUNIT_TEST( testExodusFileMappingsCyl3d);

  CPPUNIT_TEST( testExodusDiscPlateWithHole);
  CPPUNIT_TEST( testExodusDiscTwoBlocks);
  CPPUNIT_TEST( testExodusDiscTwoElemIGA);
  CPPUNIT_TEST( testExodusDiscCyl3d);
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
  CPPUNIT_TEST( testAbaqusReadFirst );
  CPPUNIT_TEST( testAbaqusReadSecond );
  CPPUNIT_TEST( testDynaReadElem );
  CPPUNIT_TEST( testDynaNoSplines );
  CPPUNIT_TEST( testDynaReadPatch );
  CPPUNIT_TEST( testDynaFileMappingsFEMEx5);
  CPPUNIT_TEST( testDynaFileMappingsBlockWithHole);
  CPPUNIT_TEST( testDynaFileMappingsPlateWithHole);
  CPPUNIT_TEST( testDynaFileMappingsCyl3d);
#endif // LIBMESH_HAVE_GZSTREAM
#endif // LIBMESH_DIM > 1
       //
#if LIBMESH_DIM > 1
  CPPUNIT_TEST( testBadGmsh );
  CPPUNIT_TEST( testGoodGmsh );
#endif

#if LIBMESH_DIM > 2
  CPPUNIT_TEST( testGoodSTL );
  CPPUNIT_TEST( testGoodSTLBinary );

#ifdef LIBMESH_HAVE_TETGEN
  CPPUNIT_TEST( testTetgenIO );
#endif
#endif

  CPPUNIT_TEST_SUITE_END();

private:

public:
  void setUp()
  {}

  void tearDown()
  {}

#ifdef LIBMESH_HAVE_VTK
  void testVTKPreserveElemIds ()
  {
    LOG_UNIT_TEST;

    // Come up with some crazy numbering.  Make all the new ids higher
    // than the existing ids so we don't have to worry about conflicts
    // while renumbering.
    dof_id_type start_id;

    // first scope: write file
    {
      Mesh mesh(*TestCommWorld);
      mesh.allow_renumbering(false);
      MeshTools::Generation::build_square (mesh, 3, 3, 0., 1., 0., 1.);

      start_id = mesh.max_elem_id();

      // Use a separate container that won't invalidate iterators when
      // we renumber
      std::set<Elem *> elements {mesh.elements_begin(), mesh.elements_end()};
      for (Elem * elem : elements)
      {
        const Point center = elem->vertex_average();
        const int xn = center(0)*3;
        const int yn = center(1)*3;
        const dof_id_type new_id = start_id + yn*5 + xn;
        mesh.renumber_elem(elem->id(), new_id);
      }

      // Explicit writer object here to be absolutely sure we get VTK
      VTKIO vtk(mesh);
      vtk.write("read_elem_ids_test.pvtu");
    }

    // Make sure that the writing is done before the reading starts.
    TestCommWorld->barrier();

    // second scope: read file
    {
      Mesh mesh(*TestCommWorld);
      mesh.allow_renumbering(false);

      mesh.read("read_elem_ids_test.pvtu");
      mesh.prepare_for_use();

      CPPUNIT_ASSERT_EQUAL(mesh.n_nodes(), dof_id_type(16));
      CPPUNIT_ASSERT_EQUAL(mesh.n_elem(), dof_id_type(9));

      for (const auto & elem : mesh.element_ptr_range())
      {
        const Point center = elem->vertex_average();
        const int xn = center(0)*3;
        const int yn = center(1)*3;
        const dof_id_type expected_id = start_id + yn*5 + xn;
        CPPUNIT_ASSERT_EQUAL(elem->id(), expected_id);
      }
    }
  }

  void testVTKPreserveSubdomainIds ()
  {
    LOG_UNIT_TEST;

    // first scope: write file
    {
      Mesh mesh(*TestCommWorld);
      mesh.allow_renumbering(false);
      MeshTools::Generation::build_square (mesh, 3, 3, 0., 1., 0., 1.);

      for (const auto & elem : mesh.element_ptr_range())
      {
        const Point center = elem->vertex_average();
        const int xn = center(0)*3;
        const int yn = center(1)*3;
        const subdomain_id_type new_id = yn*4 + xn;
        elem->subdomain_id() = new_id;
      }

      // Explicit writer object here to be absolutely sure we get VTK
      VTKIO vtk(mesh);
      vtk.write("read_sbd_ids_test.pvtu");
    }

    // Make sure that the writing is done before the reading starts.
    TestCommWorld->barrier();

    // second scope: read file
    {
      Mesh mesh(*TestCommWorld);
      mesh.allow_renumbering(false);

      mesh.read("read_sbd_ids_test.pvtu");
      mesh.prepare_for_use();

      CPPUNIT_ASSERT_EQUAL(mesh.n_nodes(), dof_id_type(16));
      CPPUNIT_ASSERT_EQUAL(mesh.n_elem(), dof_id_type(9));

      for (const auto & elem : mesh.element_ptr_range())
      {
        const Point center = elem->vertex_average();
        const int xn = center(0)*3;
        const int yn = center(1)*3;
        const subdomain_id_type expected_id = yn*4 + xn;
        CPPUNIT_ASSERT_EQUAL(elem->subdomain_id(), expected_id);
      }
    }
  }
#endif // LIBMESH_HAVE_VTK


#ifdef LIBMESH_HAVE_EXODUS_API
  void testExodusReadHeader ()
  {
    LOG_UNIT_TEST;

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

  void testLowOrderEdgeBlocks ()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
    ExodusII_IO exii(mesh);

    if (mesh.processor_id() == 0)
      exii.read("meshes/mesh_with_low_order_edge_blocks.e");

    MeshCommunication().broadcast(mesh);
    mesh.prepare_for_use();

    // Check that we see the boundary ids we expect
    BoundaryInfo & bi = mesh.get_boundary_info();

    // On a ReplicatedMesh, check that the number of edge boundary
    // conditions is as expected. The real test is that we can read
    // this file in at all. Prior to the changes in #3491, the Exodus
    // reader threw an exception while trying to read this mesh.
    if (mesh.is_serial())
    {
      // Mesh has 26 boundary ids total (including edge and side ids).
      // ss_prop1 = 200, 201 ;
      // ed_prop1 = 8000, 8001, 8002, 8003, 8004, 8005, 8006, 8007, 8008, 8009, 8010,
      //            8011, 9001, 9002, 9003, 9004, 9005, 9006, 9007, 9008, 9009, 9010,
      //            9011, 9012 ;
      CPPUNIT_ASSERT_EQUAL(static_cast<std::size_t>(26), bi.n_boundary_ids());

      // We can binary_search() the build_edge_list() which is sorted
      // in lexicographical order before it's returned.
      auto edge_list = bi.build_edge_list();

      // Search for some tuples we expect to be present
      CPPUNIT_ASSERT(std::binary_search(edge_list.begin(), edge_list.end(), std::make_tuple(4, 1, 8007)));
      CPPUNIT_ASSERT(std::binary_search(edge_list.begin(), edge_list.end(), std::make_tuple(10, 6, 8001)));

      // And make sure we don't have entries we shouldn't have
      CPPUNIT_ASSERT(!std::binary_search(edge_list.begin(), edge_list.end(), std::make_tuple(1, 8, 8009)));
      CPPUNIT_ASSERT(!std::binary_search(edge_list.begin(), edge_list.end(), std::make_tuple(2, 10, 9011)));
    }
  }

  void testExodusIGASidesets ()
  {
    LOG_UNIT_TEST;

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
    CPPUNIT_ASSERT_EQUAL(mesh.n_nodes(), static_cast<dof_id_type>(468));
    // 5^3 spline elements + 3^3 Rational Bezier elements
    CPPUNIT_ASSERT_EQUAL(mesh.n_elem(),  static_cast<dof_id_type>(152));

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
  { LOG_UNIT_TEST; testCopyNodalSolutionImpl<ReplicatedMesh,ExodusII_IO>("repl_with_nodal_soln.e"); }

  void testExodusCopyNodalSolutionDistributed ()
  { LOG_UNIT_TEST; testCopyNodalSolutionImpl<DistributedMesh,ExodusII_IO>("dist_with_nodal_soln.e"); }

#if defined(LIBMESH_HAVE_NEMESIS_API)
  void testNemesisCopyNodalSolutionReplicated ()
  { LOG_UNIT_TEST; testCopyNodalSolutionImpl<ReplicatedMesh,Nemesis_IO>("repl_with_nodal_soln.nem"); }

  void testNemesisCopyNodalSolutionDistributed ()
  { LOG_UNIT_TEST; testCopyNodalSolutionImpl<DistributedMesh,Nemesis_IO>("dist_with_nodal_soln.nem"); }
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
  { LOG_UNIT_TEST; testCopyElementSolutionImpl<ReplicatedMesh,ExodusII_IO>("repl_with_elem_soln.e"); }

  void testExodusCopyElementSolutionDistributed ()
  { LOG_UNIT_TEST; testCopyElementSolutionImpl<DistributedMesh,ExodusII_IO>("dist_with_elem_soln.e"); }

#if defined(LIBMESH_HAVE_NEMESIS_API)
  void testNemesisCopyElementSolutionReplicated ()
  { LOG_UNIT_TEST; testCopyElementSolutionImpl<ReplicatedMesh,Nemesis_IO>("repl_with_elem_soln.nem"); }

  void testNemesisCopyElementSolutionDistributed ()
  { LOG_UNIT_TEST; testCopyElementSolutionImpl<DistributedMesh,Nemesis_IO>("dist_with_elem_soln.nem"); }


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
  { LOG_UNIT_TEST; testSingleElementImpl<ReplicatedMesh,Nemesis_IO>("repl_with_single_elem.nem"); }

  void testNemesisSingleElementDistributed ()
  { LOG_UNIT_TEST; testSingleElementImpl<DistributedMesh,Nemesis_IO>("dist_with_single_elem.nem"); }
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
  { LOG_UNIT_TEST; testCopyElementVectorImpl<ReplicatedMesh, ExodusII_IO>("repl_with_elem_vec.e"); }

  void testExodusCopyElementVectorDistributed ()
  { LOG_UNIT_TEST; testCopyElementVectorImpl<DistributedMesh,ExodusII_IO>("dist_with_elem_vec.e"); }

#if defined(LIBMESH_HAVE_NEMESIS_API)
  void testNemesisCopyElementVectorReplicated ()
  { LOG_UNIT_TEST; testCopyElementVectorImpl<ReplicatedMesh,Nemesis_IO>("repl_with_elem_vec.nem"); }

  void testNemesisCopyElementVectorDistributed ()
  { LOG_UNIT_TEST; testCopyElementVectorImpl<DistributedMesh,Nemesis_IO>("dist_with_elem_vec.nem"); }
#endif


  void testExodusWriteElementDataFromDiscontinuousNodalData()
  {
    LOG_UNIT_TEST;

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


  void testExodusWriteAddedSides
    (Number (*exact_sol)(const Point &, const Parameters &, const
                         std::string &, const std::string &),
     const ElemType elem_type,
     const Order order,
     const bool write_discontinuous = false,
     const std::vector<FEType> earlier_vars = {},
     const std::vector<FEType> later_vars = {})
  {
    constexpr unsigned int nx = added_sides_nxyz[0],
                           ny = added_sides_nxyz[1],
                           nz = added_sides_nxyz[2];

    const unsigned int dim = Elem::build(elem_type)->dim();
    const bool is_tensor = (Elem::build(elem_type)->n_sides() == dim * 2);

    // Figure out how many fake and true elements to expect
    dof_id_type n_fake_elem = 0;
    dof_id_type n_true_elem = 0;

    dof_id_type n_fake_nodes = 0;
    dof_id_type n_true_nodes = 0;

    const std::string filename =
      "side_discontinuous_"+Utility::enum_to_string<ElemType>(elem_type)+(write_discontinuous?"_disc":"")+".e";

    // first scope: write file
    {
      Mesh mesh(*TestCommWorld);

      EquationSystems es(mesh);
      System & sys = es.add_system<System> ("SimpleSystem");
      int varnum = 1;
      for (auto vartype : earlier_vars)
        sys.add_variable("earlier_"+std::to_string(varnum++), vartype);

      sys.add_variable("u", order, SIDE_HIERARCHIC);

      varnum = 1;
      for (auto vartype : later_vars)
        sys.add_variable("later_"+std::to_string(varnum++), vartype);

      if (dim == 3)
        MeshTools::Generation::build_cube
          (mesh, nx, ny, nz, 0., 1., 0., 1., 0., 1., elem_type);
      else if (dim == 2)
        MeshTools::Generation::build_square
          (mesh, nx, ny, 0., 1., 0., 1., elem_type);
      else
        MeshTools::Generation::build_line
          (mesh, nx, 0., 1., elem_type);

      n_true_elem = mesh.n_elem();
      n_true_nodes = mesh.n_nodes();
      CPPUNIT_ASSERT_LESS(n_true_nodes, n_true_elem); // Ne < Nn

      const unsigned int our_ny = dim>1 ? ny : 1;
      const unsigned int our_nz = dim>2 ? nz : 1;

      dof_id_type min_n_elem = nx * our_ny * our_nz;
      CPPUNIT_ASSERT_LESSEQUAL(n_true_elem, min_n_elem); // "backwards" API...

      for (const auto & elem : mesh.active_local_element_ptr_range())
        {
          for (auto s : make_range(elem->n_sides()))
          if (!elem->neighbor_ptr(s) || elem->neighbor_ptr(s)->id() < elem->id())
            {
              ++n_fake_elem;
              auto side = elem->build_side_ptr(s);
              n_fake_nodes += side->n_nodes();
            }
        }
      mesh.comm().sum(n_fake_elem);
      mesh.comm().sum(n_fake_nodes);

      const dof_id_type expected_fakes = [elem_type]() {
        switch (elem_type)
        {
          case EDGE3:
            return nx+1;
          case TRI6:
            return 3*nx*ny + nx + ny;
          case QUAD8:
          case QUAD9:
            return 2*nx*ny + nx + ny;
          case TET14:
            return 48*nx*ny*nz + 4*(nx*ny+nx*nz+ny*nz);
          case HEX27:
            return 3*nx*ny*nz + nx*ny + nx*nz + ny*nz;
          default:
            libmesh_error();
        }
      } (); // Invoke anonymous lambda

      CPPUNIT_ASSERT_EQUAL(n_fake_elem, expected_fakes); // "backwards" API...

      es.init();
      sys.project_solution(exact_sol, nullptr,
                           es.parameters);

      // Set solution u^e_i = i, for the ith vertex of a given element e.

      // Now write to file.
      ExodusII_IO exii(mesh);
      exii.write_added_sides(true);

      if (write_discontinuous)
        exii.write_discontinuous_equation_systems(filename, es);
      else
        exii.write_equation_systems(filename, es);
    } // end first scope

    // second scope: read file, verify extra elements exist
    {
      Mesh mesh(*TestCommWorld);
      mesh.allow_renumbering(false);

      ExodusII_IO exii(mesh);

      if (mesh.processor_id() == 0)
        exii.read(filename);
      MeshCommunication().broadcast(mesh);
      mesh.prepare_for_use();

      CPPUNIT_ASSERT_EQUAL(mesh.n_elem(), n_true_elem + n_fake_elem);
      if (write_discontinuous)
        {
          const dof_id_type nodes_per_elem = Elem::build(elem_type)->n_nodes();
          CPPUNIT_ASSERT_EQUAL(mesh.n_nodes(),
                               n_true_elem*nodes_per_elem + n_fake_nodes);
        }
      else
        CPPUNIT_ASSERT_EQUAL(mesh.n_nodes(), n_true_nodes + n_fake_nodes);

      EquationSystems es(mesh);
      System & sys = es.add_system<System> ("SimpleSystem");
      // Read back into a LAGRANGE variable for testing; we still
      // can't use Exodus for a proper restart.
      sys.add_variable("ul", SECOND);
      es.init();

      const DofMap & dof_map = sys.get_dof_map();

#ifdef LIBMESH_USE_COMPLEX_NUMBERS
      exii.copy_nodal_solution(sys, "ul", "r_u");
#else
      exii.copy_nodal_solution(sys, "ul", "u");
#endif

      dof_id_type n_side_nodes = 0;
      const std::string nullstr;
      const std::string facestr = "face";

      // Debugging this in parallel is tricky.  Let's make sure that
      // if we have a failure on one rank we see it on all the others
      // and we can go on to other tests.
#ifdef LIBMESH_ENABLE_EXCEPTIONS
      bool threw_exception = false;
      try
#endif // LIBMESH_ENABLE_EXCEPTIONS
      {
      for (const auto & elem : mesh.active_local_element_ptr_range())
        {
          // Just look at side elements, not interiors
          if (elem->dim() == dim)
            continue;

          std::vector<dof_id_type> dof_indices;
          dof_map.dof_indices(elem, dof_indices, 0);

          // Find what face direction we're looking at, to
          // disambiguate when testing against a discontinuous
          // function, since we're evaluating on nodes that overlap
          // multiple faces
          const Point normal = [elem](){
            if (elem->dim() == 2)
              return Point((elem->point(1) - elem->point(0)).cross
                           (elem->point(2) - elem->point(0)));
            else if (elem->dim() == 1)
              return Point
                (elem->point(1)(1)-elem->point(0)(1),
                 elem->point(0)(0)-elem->point(1)(0));
            else
              return Point(1);
          } (); // Invoke anonymous lambda

          short faceval = -1;
          if (is_tensor)
            {
              if (std::abs(normal(0)) > TOLERANCE)
                {
                  faceval = 0;
                  libmesh_assert_less(std::abs(normal(1)), TOLERANCE);
                  libmesh_assert_less(std::abs(normal(2)), TOLERANCE);
                }
              else if (std::abs(normal(1)) > TOLERANCE)
                {
                  faceval = 1;
                  libmesh_assert_less(std::abs(normal(2)), TOLERANCE);
                }
              else
                {
                  faceval = 2;
                  libmesh_assert_greater(std::abs(normal(2)), TOLERANCE);
                }
              libmesh_assert_greater_equal(faceval, 0);
              es.parameters.set<short>(facestr) = faceval;
            }

          for (auto i : index_range(dof_indices))
            {
              const Point node_pt = elem->point(i);
              const Real nodal_coef =
                libmesh_real((*sys.current_local_solution)(dof_indices[i]));
              const Real exact_val =
                libmesh_real(exact_sol
                             (node_pt, es.parameters, nullstr,
                              nullstr));
              LIBMESH_ASSERT_FP_EQUAL
                (nodal_coef, exact_val,
                 std::max(Real(2),nodal_coef+exact_val)*
                 TOLERANCE*std::sqrt(TOLERANCE));
              ++n_side_nodes;
            }
        }
      }
#ifdef LIBMESH_ENABLE_EXCEPTIONS
      catch (...)
      {
        threw_exception = true;
        TestCommWorld->max(threw_exception);
        throw;
      }
      if (!threw_exception)
        TestCommWorld->max(threw_exception);
      CPPUNIT_ASSERT(!threw_exception);
#endif // LIBMESH_ENABLE_EXCEPTIONS

      TestCommWorld->sum(n_side_nodes);
      CPPUNIT_ASSERT_EQUAL(n_side_nodes, n_fake_nodes);
    } // end second scope
  } // end testExodusWriteAddedSides

  void testExodusWriteAddedSidesEdgeC0()
  {
    testExodusWriteAddedSides(six_x_plus_sixty_y, EDGE3, FIRST);
    testExodusWriteAddedSides(six_x_plus_sixty_y, EDGE3, SECOND);
  }

  void testExodusDiscWriteAddedSidesEdgeC0()
  {
    testExodusWriteAddedSides(six_x_plus_sixty_y, EDGE3, FIRST, true);
    testExodusWriteAddedSides(six_x_plus_sixty_y, EDGE3, SECOND, true);
  }

  void testExodusWriteAddedSidesMixedEdgeC0()
  {
    testExodusWriteAddedSides(six_x_plus_sixty_y, EDGE3, FIRST, false, {{FIRST, LAGRANGE}});
    testExodusWriteAddedSides(six_x_plus_sixty_y, EDGE3, SECOND, false, {}, {{FIRST, LAGRANGE}});
  }

  void testExodusDiscWriteAddedSidesMixedEdgeC0()
  {
    testExodusWriteAddedSides(six_x_plus_sixty_y, EDGE3, FIRST, true, {{FIRST, LAGRANGE}});
    testExodusWriteAddedSides(six_x_plus_sixty_y, EDGE3, SECOND, true, {}, {{FIRST, LAGRANGE}});
  }

  void testExodusWriteAddedSidesEdgeDisc()
  {
    testExodusWriteAddedSides(designed_for_side_elems, EDGE3, SECOND);
  }

  void testExodusDiscWriteAddedSidesEdgeDisc()
  {
    testExodusWriteAddedSides(designed_for_side_elems, EDGE3, SECOND, true);
  }

  void testExodusWriteAddedSidesTriC0()
  {
    testExodusWriteAddedSides(six_x_plus_sixty_y, TRI6, FIRST);
    testExodusWriteAddedSides(six_x_plus_sixty_y, TRI6, SECOND);
  }

  void testExodusDiscWriteAddedSidesTriC0()
  {
    testExodusWriteAddedSides(six_x_plus_sixty_y, TRI6, FIRST, true);
    testExodusWriteAddedSides(six_x_plus_sixty_y, TRI6, SECOND, true);
  }

  void testExodusWriteAddedSidesMixedTriC0()
  {
    testExodusWriteAddedSides(six_x_plus_sixty_y, TRI6, FIRST, false, {{SECOND, HIERARCHIC}});
    testExodusWriteAddedSides(six_x_plus_sixty_y, TRI6, SECOND, false, {}, {{SECOND, SZABAB}});
  }

  void testExodusDiscWriteAddedSidesMixedTriC0()
  {
    testExodusWriteAddedSides(six_x_plus_sixty_y, TRI6, FIRST, true, {{SECOND, HIERARCHIC}});
    testExodusWriteAddedSides(six_x_plus_sixty_y, TRI6, SECOND, true, {}, {{SECOND, SZABAB}});
  }

  void testExodusWriteAddedSidesTriDisc()
  {
    testExodusWriteAddedSides(designed_for_side_elems, TRI6, SECOND);
  }

  void testExodusDiscWriteAddedSidesTriDisc()
  {
    testExodusWriteAddedSides(designed_for_side_elems, TRI6, SECOND, true);
  }

  void testExodusWriteAddedSidesQuadC0()
  {
    testExodusWriteAddedSides(six_x_plus_sixty_y, QUAD9, FIRST);
    testExodusWriteAddedSides(six_x_plus_sixty_y, QUAD9, SECOND);
  }

  void testExodusDiscWriteAddedSidesQuadC0()
  {
    testExodusWriteAddedSides(six_x_plus_sixty_y, QUAD9, FIRST, true);
    testExodusWriteAddedSides(six_x_plus_sixty_y, QUAD9, SECOND, true);
  }

  void testExodusWriteAddedSidesMixedQuadC0()
  {
    testExodusWriteAddedSides(six_x_plus_sixty_y, QUAD9, FIRST, false, {{SECOND, LAGRANGE}});
    testExodusWriteAddedSides(six_x_plus_sixty_y, QUAD9, SECOND, false, {}, {{FIRST, LAGRANGE}});
  }

  void testExodusDiscWriteAddedSidesMixedQuadC0()
  {
    testExodusWriteAddedSides(six_x_plus_sixty_y, QUAD9, FIRST, true, {{SECOND, LAGRANGE}});
    testExodusWriteAddedSides(six_x_plus_sixty_y, QUAD9, SECOND, true, {}, {{FIRST, LAGRANGE}});
  }

  void testExodusWriteAddedSidesQuadDisc()
  {
    testExodusWriteAddedSides(designed_for_side_elems, QUAD9, SECOND);
  }

  void testExodusDiscWriteAddedSidesQuadDisc()
  {
    testExodusWriteAddedSides(designed_for_side_elems, QUAD9, SECOND, true);
  }

  void testExodusWriteAddedSidesTetC0()
  {
    testExodusWriteAddedSides(six_x_plus_sixty_y, TET14, FIRST);
    testExodusWriteAddedSides(six_x_plus_sixty_y, TET14, SECOND);
  }

  void testExodusDiscWriteAddedSidesTetC0()
  {
    testExodusWriteAddedSides(six_x_plus_sixty_y, TET14, FIRST, true);
    testExodusWriteAddedSides(six_x_plus_sixty_y, TET14, SECOND, true);
  }

  void testExodusWriteAddedSidesTetDisc()
  {
    testExodusWriteAddedSides(designed_for_side_elems, TET14, SECOND);
  }

  void testExodusDiscWriteAddedSidesTetDisc()
  {
    testExodusWriteAddedSides(designed_for_side_elems, TET14, SECOND, true);
  }

  void testExodusWriteAddedSidesHexC0()
  {
    testExodusWriteAddedSides(six_x_plus_sixty_y, HEX27, FIRST);
    testExodusWriteAddedSides(six_x_plus_sixty_y, HEX27, SECOND);
  }

  void testExodusDiscWriteAddedSidesHexC0()
  {
    testExodusWriteAddedSides(six_x_plus_sixty_y, HEX27, FIRST, true);
    testExodusWriteAddedSides(six_x_plus_sixty_y, HEX27, SECOND, true);
  }

  void testExodusWriteAddedSidesMixedHexC0()
  {
    testExodusWriteAddedSides(six_x_plus_sixty_y, HEX27, FIRST, false, {{FIRST, LAGRANGE}});
    testExodusWriteAddedSides(six_x_plus_sixty_y, HEX27, SECOND, false, {}, {{SECOND, HIERARCHIC}});
  }

  void testExodusDiscWriteAddedSidesMixedHexC0()
  {
    testExodusWriteAddedSides(six_x_plus_sixty_y, HEX27, FIRST, true, {{FIRST, LAGRANGE}});
    testExodusWriteAddedSides(six_x_plus_sixty_y, HEX27, SECOND, true, {}, {{SECOND, HIERARCHIC}});
  }

  void testExodusWriteAddedSidesHexDisc()
  {
    testExodusWriteAddedSides(designed_for_side_elems, HEX27, SECOND);
  }

  void testExodusDiscWriteAddedSidesHexDisc()
  {
    testExodusWriteAddedSides(designed_for_side_elems, HEX27, SECOND, true);
  }

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
      CPPUNIT_ASSERT_EQUAL(mesh.n_elem(),  static_cast<dof_id_type>(9));
      CPPUNIT_ASSERT_EQUAL(mesh.n_nodes(), static_cast<dof_id_type>(16));
    }
  }

  void testNemesisReadReplicated ()
  { LOG_UNIT_TEST; testNemesisReadImpl<ReplicatedMesh>(); }

  void testNemesisReadDistributed ()
  { LOG_UNIT_TEST; testNemesisReadImpl<DistributedMesh>(); }
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

        Point physical_pt = FEMap::map(elem->dim(), elem, master_pt);

        Point inverse_pt = FEMap::inverse_map(elem->dim(), elem,
                                              physical_pt);

        CPPUNIT_ASSERT((inverse_pt-master_pt).norm() < TOLERANCE);

        CPPUNIT_ASSERT(elem->contains_point(physical_pt));

        // We only want to find elements in the same block
        std::set<subdomain_id_type> my_subdomain { elem->subdomain_id() };

        // We can *still* have overlapping NodeElem from a slit mesh
        // input file; better check them all
        std::set<const Elem * > located_elems;
        (*locator)(physical_pt, located_elems, &my_subdomain);

        CPPUNIT_ASSERT(located_elems.count(elem));
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


  void testAbaqusRead (const std::string & fname,
                       dof_id_type n_elem,
                       dof_id_type n_nodes)
  {
    Mesh mesh(*TestCommWorld);

    AbaqusIO abaqus(mesh);

    if (mesh.processor_id() == 0)
      abaqus.read(fname);
    MeshCommunication().broadcast (mesh);

    mesh.prepare_for_use();

    CPPUNIT_ASSERT_EQUAL(mesh.n_elem(),  n_elem);
    CPPUNIT_ASSERT_EQUAL(mesh.n_nodes(), n_nodes);
  }


  void testAbaqusReadFirst()
  {
    LOG_UNIT_TEST;
    testAbaqusRead("meshes/tensile_sample_test1.inp.gz", 728, 1166);
  }


  void testAbaqusReadSecond()
  {
    LOG_UNIT_TEST;
    testAbaqusRead("meshes/poly_sample_test2.inp.gz", 1280, 1625);
  }


  void testDynaReadElem ()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);

    DynaIO dyna(mesh);

    if (mesh.processor_id() == 0)
      dyna.read("meshes/1_quad.bxt.gz");
    MeshCommunication().broadcast (mesh);

    mesh.prepare_for_use();

    // We have 1 QUAD9 finite element, attached via a trivial map to 9
    // spline Node+NodeElem objects
    CPPUNIT_ASSERT_EQUAL(mesh.n_elem(),  static_cast<dof_id_type>(10));
    CPPUNIT_ASSERT_EQUAL(mesh.n_nodes(), static_cast<dof_id_type>(18));

    helperTestingDynaQuad(mesh);
  }

  void testBadGmsh ()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);

    GmshIO gmsh_io(mesh);

#ifdef LIBMESH_ENABLE_EXCEPTIONS
    std::string what = "";
    try
    {
      if (mesh.processor_id() == 0)
        gmsh_io.read("meshes/block.msh");
    }
    catch (libMesh::LogicError & e)
    {
      what = e.what();
    }

    TestCommWorld->broadcast(what);
    std::regex msg_regex("outside entity physical bounding box");
    CPPUNIT_ASSERT(std::regex_search(what, msg_regex));
#endif
  }

  void testGoodGmsh ()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);

    GmshIO gmsh_io(mesh);

    if (mesh.processor_id() == 0)
      gmsh_io.read("meshes/circle.msh");
    MeshCommunication().broadcast(mesh);

    CPPUNIT_ASSERT_EQUAL(mesh.n_elem(), dof_id_type(14));
  }

  void testGoodSTL ()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);

    STLIO stl_io(mesh);

    if (mesh.processor_id() == 0)
      stl_io.read("meshes/Cluster_34.stl");
    MeshCommunication().broadcast(mesh);

    CPPUNIT_ASSERT_EQUAL(mesh.n_elem(), dof_id_type(40));
  }

  void testGoodSTLBinary ()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);

    STLIO stl_io(mesh);

    if (mesh.processor_id() == 0)
      stl_io.read("meshes/engraving.stl");
    MeshCommunication().broadcast(mesh);

    CPPUNIT_ASSERT_EQUAL(mesh.n_elem(), dof_id_type(426));
    CPPUNIT_ASSERT_EQUAL(mesh.n_nodes(), dof_id_type(215));
  }

  void testTetgenIO ()
  {
#ifdef LIBMESH_HAVE_TETGEN
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);

    TetGenIO tetgen_io(mesh);

    if (mesh.processor_id() == 0)
      tetgen_io.read("meshes/tetgen_one_tet10.ele");
    MeshCommunication().broadcast(mesh);

    mesh.prepare_for_use();

    // Mesh should contain 1 TET10 finite element
    CPPUNIT_ASSERT_EQUAL(mesh.n_elem(),  static_cast<dof_id_type>(1));
    CPPUNIT_ASSERT_EQUAL(mesh.n_nodes(), static_cast<dof_id_type>(10));

    // Element should have TET10 reference element volume
    const Elem * const elem = mesh.query_elem_ptr(0);

    // On a serial mesh we have every element everywhere
    if (mesh.is_serial())
      CPPUNIT_ASSERT(elem);
    else
      {
        bool have_elem = elem;
        mesh.comm().max(have_elem);
        CPPUNIT_ASSERT(have_elem);
      }

    if (elem)
      {
        const Real vol = elem->volume();
        CPPUNIT_ASSERT_DOUBLES_EQUAL(vol, 1./6, TOLERANCE*TOLERANCE);
      }
#endif
  }

  void testDynaNoSplines ()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);

    DynaIO dyna(mesh, /* keep_spline_nodes = */ false);

    if (mesh.processor_id() == 0)
      dyna.read("meshes/1_quad.bxt.gz");
    MeshCommunication().broadcast (mesh);

    mesh.prepare_for_use();

    // We have 1 QUAD9 finite element
    CPPUNIT_ASSERT_EQUAL(mesh.n_elem(),  static_cast<dof_id_type>(1));
    CPPUNIT_ASSERT_EQUAL(mesh.n_nodes(), static_cast<dof_id_type>(9));

    helperTestingDynaQuad(mesh);
  }


  void testDynaReadPatch ()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);

    DynaIO dyna(mesh);
    if (mesh.processor_id() == 0)
      dyna.read("meshes/25_quad.bxt.gz");
    MeshCommunication().broadcast (mesh);

    mesh.prepare_for_use();

    // We have 5^2 QUAD9 elements, with 11^2 nodes,
    // tied to 49 Node/NodeElem spline nodes
    CPPUNIT_ASSERT_EQUAL(mesh.n_elem(),  static_cast<dof_id_type>(25+49));
    CPPUNIT_ASSERT_EQUAL(mesh.n_nodes(), static_cast<dof_id_type>(121+49));

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
    CPPUNIT_ASSERT_EQUAL(sys.get_dof_map().n_constrained_dofs(), static_cast<dof_id_type>(121));
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

    // Let's test that IGA constraints are preserved (in a relative
    // sense) when we clone a mesh.
    std::unique_ptr<MeshBase> mesh_clone = mesh.clone();
    CPPUNIT_ASSERT(*mesh_clone == mesh);

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
    LOG_UNIT_TEST;

    testDynaFileMappings("meshes/PressurizedCyl_Patch6_256Elem.bxt.gz",
    // Regression values for sin_x_plus_cos_y
                         {{0.9639857809698268, 1.839870171669186,
                           0.7089812562241862, 1.306121188539059}});
  }

  void testDynaFileMappingsBlockWithHole ()
  {
    LOG_UNIT_TEST;

    testDynaFileMappings("meshes/BlockWithHole_Patch9.bxt.gz",
    // Regression values for sin_x_plus_cos_y
                         {{3.22612556930183, 1.97405365384733,
                           2.53376235803176, 1.41374070517223}});
  }

  void testDynaFileMappingsPlateWithHole ()
  {
    LOG_UNIT_TEST;

    testDynaFileMappings("meshes/PlateWithHole_Patch8.bxt.gz",
    // Regression values for sin_x_plus_cos_y
                         {{2.2812154374012, 1.974049990211937,
                           1.791640772215248, 1.413679237529376}});
  }

  void testDynaFileMappingsCyl3d ()
  {
    LOG_UNIT_TEST;

    testDynaFileMappings("meshes/PressurizedCyl3d_Patch1_8Elem.bxt.gz",
    // Regression values for sin_x_plus_cos_y
                         {{0.963612880188165, 1.82329452603503,
                           0.707998701597943, 1.31399222566683}});
  }

  void testExodusFileMappings (const std::string & filename,
                               std::array<Real, 4> expected_norms,
                               bool use_disc_bex = false)
  {
    Mesh mesh(*TestCommWorld);

    ExodusII_IO exii(mesh);
    // IGA Exodus meshes require ExodusII 8 or higher
    if (exii.get_exodus_version() < 800)
      return;

    // This should default to false
    if (use_disc_bex)
      exii.set_discontinuous_bex(true);

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

#ifdef LIBMESH_HAVE_VTK
    {
      VTKIO vtkout(mesh);

      vtkout.write("vtk_file_mapping_out.pvtu");
    }
#endif
  }

  void testExodusFileMappingsPlateWithHole ()
  {
    LOG_UNIT_TEST;

    testExodusFileMappings("meshes/PlateWithHole_Patch8.e",
    // Regression values for sin_x_plus_cos_y
                           {{2.2812154374012, 1.974049990211937,
                             1.791640772215248, 1.413679237529376}});
  }

  void testExodusFileMappingsTwoBlocks ()
  {
    LOG_UNIT_TEST;

    testExodusFileMappings("meshes/two_quads_two_blocks.e",
    // Regression values for sin_x_plus_cos_y
                           {{2.03496953073072, 1.97996853164955,
                             1.18462134113435, 1.03085301158959}});
  }
  void testExodusFileMappingsTwoElemIGA()
  {
    LOG_UNIT_TEST;

    testExodusFileMappings("meshes/two_element_iga_in.e",
    // Regression values for sin_x_plus_cos_y
                           {{1.26865962862531, 1.42562070158386,
                             1.54905363492342, 1.29782906548366}});
  }

  void testExodusFileMappingsCyl3d ()
  {
    LOG_UNIT_TEST;

    testExodusFileMappings("meshes/PressurizedCyl3d_Patch1_8Elem.e",
                           {{0.963612880188165, 1.82329452603503,
                             0.707998701597943, 1.31399222566683}});
  }

  void testExodusDiscPlateWithHole ()
  {
    LOG_UNIT_TEST;

    testExodusFileMappings("meshes/PlateWithHole_Patch8.e",
    // Regression values for sin_x_plus_cos_y
    //
    // These are *not* the same as for the continuous plate, because
    // we do the C^TKCx=C^Tf trick to pull back projections to the
    // spline nodes, and the pseudoinverse here minimizes the error in
    // a discretization-dependent norm, not in a Sobolev norm.  For
    // these coarse meshes, our Sobolev norms can end up being ~0.1%
    // different.
                           {{2.28234312456534, 1.97439548757586,
                             1.79290449809266, 1.41075128955985}},
                             true);
  }

  void testExodusDiscTwoBlocks ()
  {
    LOG_UNIT_TEST;

    testExodusFileMappings("meshes/two_quads_two_blocks.e",
    // Regression values for sin_x_plus_cos_y
                           {{2.03496953073072, 1.97996853164955,
                             1.18462134113435, 1.03085301158959}},
                             true);
  }
  void testExodusDiscTwoElemIGA()
  {
    LOG_UNIT_TEST;

    testExodusFileMappings("meshes/two_element_iga_in.e",
    // Regression values for sin_x_plus_cos_y
                           {{1.26877626663365, 1.42553698909339,
                             1.54810114917177, 1.29792704408979}},
                             true);
  }

  void testExodusDiscCyl3d ()
  {
    LOG_UNIT_TEST;

    testExodusFileMappings("meshes/PressurizedCyl3d_Patch1_8Elem.e",
                           {{0.963855209590556, 1.8234396424318,
                             0.708286572453382, 1.31468940958327}},
                             true);
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION( MeshInputTest );
