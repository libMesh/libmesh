#include <libmesh/equation_systems.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_communication.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/replicated_mesh.h>
#include <libmesh/dyna_io.h>
#include <libmesh/exodusII_io.h>
#include <libmesh/dof_map.h>

#include "test_comm.h"
#include "libmesh_cppunit.h"


using namespace libMesh;


Number x_plus_y (const Point& p,
                 const Parameters&,
                 const std::string&,
                 const std::string&)
{
  const Real & x = p(0);
  const Real & y = p(1);

  return x + y;
}


class MeshInputTest : public CppUnit::TestCase {
public:
  CPPUNIT_TEST_SUITE( MeshInputTest );

#if LIBMESH_DIM > 1
#ifdef LIBMESH_HAVE_EXODUS_API
  CPPUNIT_TEST( testExodusCopyElementSolution );
#ifndef LIBMESH_USE_COMPLEX_NUMBERS
  // Eventually this will support complex numbers.
  CPPUNIT_TEST( testExodusWriteElementDataFromDiscontinuousNodalData );
#endif
#endif
  CPPUNIT_TEST( testDynaReadElem );
  CPPUNIT_TEST( testDynaReadPatch );

  CPPUNIT_TEST( testMeshMoveConstructor );
#endif // LIBMESH_DIM > 1

  CPPUNIT_TEST_SUITE_END();

private:

public:
  void setUp()
  {}

  void tearDown()
  {}

#ifdef LIBMESH_HAVE_EXODUS_API
  void testExodusCopyElementSolution ()
  {
    {
      Mesh mesh(*TestCommWorld);

      EquationSystems es(mesh);
      System &sys = es.add_system<System> ("SimpleSystem");
      sys.add_variable("e", CONSTANT, MONOMIAL);

      MeshTools::Generation::build_square (mesh,
                                           3, 3,
                                           0., 1., 0., 1.);

      es.init();
      sys.project_solution(x_plus_y, nullptr, es.parameters);

      ExodusII_IO exii(mesh);

      // Don't try to write element data as nodal data
      std::set<std::string> sys_list;
      exii.write_equation_systems("mesh_with_soln.e", es, &sys_list);

      // Just write it as element data
      exii.write_element_data(es);
    }

    // copy_elemental_solution currently requires ReplicatedMesh
    {
      ReplicatedMesh mesh(*TestCommWorld);

      EquationSystems es(mesh);
      System &sys = es.add_system<System> ("SimpleSystem");
      sys.add_variable("teste", CONSTANT, MONOMIAL);

      ExodusII_IO exii(mesh);

      if (mesh.processor_id() == 0)
        exii.read("mesh_with_soln.e");
      MeshCommunication().broadcast(mesh);
      mesh.prepare_for_use();

      es.init();

      // Read the solution e into variable teste.
      //
      // With complex numbers, we'll only bother reading the real
      // part.
#ifdef LIBMESH_USE_COMPLEX_NUMBERS
      exii.copy_elemental_solution(sys, "teste", "r_e");
#else
      exii.copy_elemental_solution(sys, "teste", "e");
#endif

      // Exodus only handles double precision
      Real exotol = std::max(TOLERANCE*TOLERANCE, Real(1e-12));

      for (Real x = Real(1.L/6.L); x < 1; x += Real(1.L/3.L))
        for (Real y = Real(1.L/6.L); y < 1; y += Real(1.L/3.L))
          {
            Point p(x,y);
            LIBMESH_ASSERT_FP_EQUAL(libmesh_real(sys.point_value(0,p)),
                                    libmesh_real(x+y),
                                    exotol);
          }
    }
  }


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
              (read_val, expected_values[i], TOLERANCE*TOLERANCE);
          }
    } // end second scope
  } // end testExodusWriteElementDataFromDiscontinuousNodalData

#endif // !LIBMESH_USE_COMPLEX_NUMBERS
#endif // LIBMESH_HAVE_EXODUS_API

  void testDynaReadElem ()
  {
    Mesh mesh(*TestCommWorld);

    DynaIO dyna(mesh);
    if (mesh.processor_id() == 0)
      dyna.read("1_quad.dyn");
    MeshCommunication().broadcast (mesh);

    mesh.prepare_for_use();

    CPPUNIT_ASSERT_EQUAL(mesh.n_elem(), dof_id_type(1));
    CPPUNIT_ASSERT_EQUAL(mesh.n_nodes(), dof_id_type(9));

    CPPUNIT_ASSERT_EQUAL(mesh.default_mapping_type(),
                         RATIONAL_BERNSTEIN_MAP);

    unsigned char weight_index = mesh.default_mapping_data();

    if (mesh.query_elem_ptr(0))
      {
        const Elem & elem = mesh.elem_ref(0);

        CPPUNIT_ASSERT_EQUAL(elem.type(), QUAD9);
        for (unsigned int n=0; n != 9; ++n)
          CPPUNIT_ASSERT_EQUAL
            (elem.node_ref(n).get_extra_datum<Real>(weight_index),
             Real(0.75));

        CPPUNIT_ASSERT_EQUAL(elem.point(0)(0), Real(0.5));
        CPPUNIT_ASSERT_EQUAL(elem.point(0)(1), Real(0.5));
        CPPUNIT_ASSERT_EQUAL(elem.point(1)(0), Real(1.5));
        CPPUNIT_ASSERT_EQUAL(elem.point(1)(1), Real(0.5));
        CPPUNIT_ASSERT_EQUAL(elem.point(2)(0), Real(1.5));
        CPPUNIT_ASSERT_EQUAL(elem.point(2)(1), Real(1.5));
        CPPUNIT_ASSERT_EQUAL(elem.point(3)(0), Real(0.5));
        CPPUNIT_ASSERT_EQUAL(elem.point(3)(1), Real(1.5));
        CPPUNIT_ASSERT(elem.has_affine_map());
#if LIBMESH_DIM > 2
        for (unsigned int v=0; v != 4; ++v)
          CPPUNIT_ASSERT_EQUAL(elem.point(v)(2), Real(0));
#endif
      }
  }


  void testDynaReadPatch ()
  {
    Mesh mesh(*TestCommWorld);

    DynaIO dyna(mesh);
    if (mesh.processor_id() == 0)
      dyna.read("25_quad.bxt");
    MeshCommunication().broadcast (mesh);

    mesh.prepare_for_use();

    CPPUNIT_ASSERT_EQUAL(mesh.n_elem(), dof_id_type(25));
    CPPUNIT_ASSERT_EQUAL(mesh.n_nodes(), dof_id_type(121));

    CPPUNIT_ASSERT_EQUAL(mesh.default_mapping_type(),
                         RATIONAL_BERNSTEIN_MAP);

    unsigned char weight_index = mesh.default_mapping_data();

    for (const auto & elem : mesh.active_element_ptr_range())
      {
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
  }


  void testMeshMoveConstructor ()
  {
    Mesh mesh(*TestCommWorld);
    MeshTools::Generation::build_square (mesh,
                                         3, 3,
                                         0., 1., 0., 1.);

    // Construct mesh2, stealing the resources of the original.
    Mesh mesh2(std::move(mesh));

    // Make sure mesh2 now has the 9 elements.
    CPPUNIT_ASSERT_EQUAL(mesh2.n_elem(),
                         static_cast<dof_id_type>(9));

    // Verify that the moved-from mesh's Partitioner and BoundaryInfo
    // objects were successfully stolen.  Note: moved-from unique_ptrs
    // are guaranteed to compare equal to nullptr, see e.g. Section
    // 20.8.1/4 of the standard.
    // https://stackoverflow.com/questions/24061767/is-unique-ptr-guaranteed-to-store-nullptr-after-move
    CPPUNIT_ASSERT(!mesh.partitioner());
    CPPUNIT_ASSERT(!mesh.boundary_info);
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION( MeshInputTest );
