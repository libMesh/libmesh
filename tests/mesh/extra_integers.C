#include <libmesh/libmesh.h>
#include <libmesh/mesh.h>
#include <libmesh/elem.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_refinement.h>

#include <libmesh/exodusII_io.h>

#include "test_comm.h"
#include "libmesh_cppunit.h"

#include <array>
#include <cstddef> // std::ptrdiff_t

using namespace libMesh;

class ExtraIntegersTest : public CppUnit::TestCase
{
  /**
   * The goal of this test is to verify the ability to add extra
   * integer storage to a Mesh object, then set and query those extra
   * integers in the objects within the mesh.
   */
public:
  LIBMESH_CPPUNIT_TEST_SUITE( ExtraIntegersTest );

  CPPUNIT_TEST( testExtraIntegersEdge2 );
  CPPUNIT_TEST( testExtraIntegersTri6 );

#ifdef LIBMESH_HAVE_EXODUS_API
  CPPUNIT_TEST( testExtraIntegersExodusReading );
#endif
#if defined(LIBMESH_HAVE_EXODUS_API) && defined(LIBMESH_ENABLE_EXCEPTIONS)
  CPPUNIT_TEST( testBadExtraIntegersExodusReading );
#endif

#ifdef LIBMESH_HAVE_XDR
  CPPUNIT_TEST( testExtraIntegersCheckpointEdge3 );
  CPPUNIT_TEST( testExtraIntegersCheckpointHex8 );
#endif

  CPPUNIT_TEST_SUITE_END();

protected:
  // Helper functions called by the test implementations, saves a few lines of code.
  std::array<unsigned int, 6>
    build_mesh(Mesh & mesh, ElemType elem_type, unsigned int n_elem_per_side)
  {
    // Request some extra integers before building
    unsigned int i1 = mesh.add_elem_integer("i1");
    unsigned int r1 = mesh.add_elem_datum<Real>("r1");
    unsigned int ni1 = mesh.add_node_integer("ni1");
    unsigned int ni2 = mesh.add_node_integer("ni2");
    unsigned int nr1 = mesh.add_node_datum<Real>("nr1");
    unsigned int nr2 = mesh.add_node_datum<Real>("nr2");

    const std::unique_ptr<Elem> test_elem = Elem::build(elem_type);
    const unsigned int ymax = test_elem->dim() > 1;
    const unsigned int zmax = test_elem->dim() > 2;
    const unsigned int ny = ymax * n_elem_per_side;
    const unsigned int nz = zmax * n_elem_per_side;

    MeshTools::Generation::build_cube (mesh,
                                       n_elem_per_side,
                                       ny,
                                       nz,
                                       0., 1.,
                                       0., ymax,
                                       0., zmax,
                                       elem_type);
    return {{i1, r1, ni1, ni2, nr1, nr2}};
  }


  void test_and_set_initial_data
    (Mesh & mesh, std::array<unsigned int, 6> ini)
  {
    const unsigned int i1 = ini[0],
                       r1 = ini[1],
                       ni1 = ini[2],
                       ni2 = ini[3],
                       nr1 = ini[4];
    for (const auto & elem : mesh.element_ptr_range())
      {
        const unsigned int expected_extra_ints =
          2 + (sizeof(Real)-1)/sizeof(dof_id_type);
        CPPUNIT_ASSERT_EQUAL(elem->n_extra_integers(), expected_extra_ints);
        CPPUNIT_ASSERT_EQUAL(elem->get_extra_integer(i1), DofObject::invalid_id);
        elem->set_extra_integer(i1, dof_id_type(elem->point(0)(0)*100));
        CPPUNIT_ASSERT_EQUAL(elem->get_extra_integer(i1), dof_id_type(elem->point(0)(0)*100));
        elem->set_extra_datum<Real>(r1, elem->point(0)(0)*1000);
        CPPUNIT_ASSERT_EQUAL(elem->get_extra_datum<Real>(r1), elem->point(0)(0)*1000);
      }

    for (const auto & node : mesh.node_ptr_range())
      {
        const unsigned int expected_extra_ints =
          4 + (sizeof(Real)-1)/sizeof(dof_id_type)*2;
        CPPUNIT_ASSERT_EQUAL(node->n_extra_integers(), expected_extra_ints);
        CPPUNIT_ASSERT_EQUAL(node->get_extra_integer(ni1), DofObject::invalid_id);
        CPPUNIT_ASSERT_EQUAL(node->get_extra_integer(ni2), DofObject::invalid_id);
        node->set_extra_datum<Real>(nr1, (*node)(0)*1000);
        CPPUNIT_ASSERT_EQUAL(node->get_extra_datum<Real>(nr1), (*node)(0)*1000);
      }

  }


  void test_final_integers(Mesh & mesh, unsigned int i1)
  {
    // Make sure old (level 0) elements have the same integers and any new
    // elements have inherited their ancestors' settings
    for (const auto & elem : mesh.element_ptr_range())
      {
        const Elem * top_parent = elem;
#ifdef LIBMESH_ENABLE_AMR
        top_parent = elem->top_parent();

        if (!elem->level())
#endif
          for (auto & node : elem->node_ref_range())
            {
              const unsigned int expected_extra_ints =
                4 + (sizeof(Real)-1)/sizeof(dof_id_type)*2;
              CPPUNIT_ASSERT_EQUAL(node.n_extra_integers(), expected_extra_ints);
            }

        const unsigned int expected_extra_ints =
          2 + (sizeof(Real)-1)/sizeof(dof_id_type);
        CPPUNIT_ASSERT_EQUAL(elem->n_extra_integers(), expected_extra_ints);
        CPPUNIT_ASSERT_EQUAL(elem->get_extra_integer(i1), dof_id_type(top_parent->point(0)(0)*100));
      }
  }


  void test_helper(ElemType elem_type, unsigned int n_elem_per_side)
  {
    Mesh mesh(*TestCommWorld);

    std::array<unsigned int, 6> ini = build_mesh(mesh, elem_type, n_elem_per_side);
    const unsigned int i1 = ini[0],
                       ni1 = ini[2],
                       ni2 = ini[3];

    test_and_set_initial_data(mesh, ini);

    // Force (in parallel) a different partitioning - we'll simply put
    // everything on rank 0, which hopefully is not what our default
    // partitioner did!
    mesh.partition(1);

    CPPUNIT_ASSERT_EQUAL(i1, mesh.add_elem_integer("i1"));
    CPPUNIT_ASSERT_EQUAL(ni1, mesh.add_node_integer("ni1"));
    CPPUNIT_ASSERT_EQUAL(ni2, mesh.add_node_integer("ni2"));

    // Make sure we didn't screw up any extra integers thereby.
    test_final_integers(mesh, i1);

    // Try out switching to 2nd order and back
    mesh.all_second_order();

    test_final_integers(mesh, i1);

    mesh.all_first_order();

    test_final_integers(mesh, i1);

    // And try refining if we can
#ifdef LIBMESH_ENABLE_AMR
    MeshRefinement mr(mesh);
    mr.uniformly_refine(1);

    test_final_integers(mesh, i1);
#endif
  }

  void checkpoint_helper(ElemType elem_type, unsigned int n_elem_per_side, bool binary)
  {
    Mesh mesh(*TestCommWorld);

    std::array<unsigned int, 6> ini = build_mesh(mesh, elem_type, n_elem_per_side);
    const unsigned int i1 = ini[0], ni1 = ini[2], ni2 = ini[3];

    test_and_set_initial_data(mesh, ini);

    const std::string filename =
      std::string("extra_integers.cp") + (binary ? "r" : "a");

    mesh.write(filename);

    TestCommWorld->barrier();

    Mesh mesh2(*TestCommWorld);

    mesh2.read(filename);

    // Make sure the integers got transferred to the second mesh
    CPPUNIT_ASSERT_EQUAL(i1, mesh2.add_elem_integer("i1"));
    CPPUNIT_ASSERT_EQUAL(ni1, mesh2.add_node_integer("ni1"));
    CPPUNIT_ASSERT_EQUAL(ni2, mesh2.add_node_integer("ni2"));

    // Make sure we didn't screw up any extra integers thereby.
    test_final_integers(mesh2, i1);
  }

public:
  void setUp() {}

  void tearDown() {}

  void testExtraIntegersEdge2() { LOG_UNIT_TEST; test_helper(EDGE2, 5); }

  void testExtraIntegersTri6() { LOG_UNIT_TEST; test_helper(TRI6, 4); }

  void testExtraIntegersCheckpointEdge3() { LOG_UNIT_TEST; checkpoint_helper(EDGE3, 5, false); }

  void testExtraIntegersCheckpointHex8() { LOG_UNIT_TEST; checkpoint_helper(HEX8, 2, true); }

#ifdef LIBMESH_HAVE_EXODUS_API
  void testExtraIntegersExodusReading()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
    mesh.allow_renumbering(false);

    // This 3-by-3 mesh contains the following element integers:
    std::vector<std::ptrdiff_t> material_id = {0, -1, 2,
                                               1,  3, 450359962,
                                               2,  3, 450359963};
    const std::string filename = "meshes/good_32bit_elem_integers.e";
    ExodusII_IO exreader(mesh);
    exreader.set_extra_integer_vars({"material_id"});
    exreader.read(filename);

    // Test that the ExodusII_IO::get_{elem,node}_num_map() APIs give
    // us something sensible.  Note: this is unrelated to reading
    // extra integers, but, in my opinion, it does not warrant its own
    // standalone test either.
    const auto & elem_num_map = exreader.get_elem_num_map();
    const auto & node_num_map = exreader.get_node_num_map();

    // This mesh has trivial elem_num_map and node_num_map
    CPPUNIT_ASSERT_EQUAL(int(elem_num_map.size()), 9);
    CPPUNIT_ASSERT_EQUAL(int(node_num_map.size()), 16);
    for (int i=0; i != 9; ++i)
      CPPUNIT_ASSERT_EQUAL(elem_num_map[i], i+1);
    for (int i=0; i != 16; ++i)
      CPPUNIT_ASSERT_EQUAL(node_num_map[i], i+1);

    CPPUNIT_ASSERT(mesh.has_elem_integer("material_id"));
    unsigned int int_idx = mesh.get_elem_integer_index("material_id");
    for (dof_id_type i=0; i != 9; ++i)
      {
        Elem * elem = mesh.query_elem_ptr(i);
        if (!elem)
          continue;
        if (material_id[i] == -1)
          CPPUNIT_ASSERT_EQUAL(DofObject::invalid_id, elem->get_extra_integer(int_idx));
        else
          CPPUNIT_ASSERT_EQUAL(dof_id_type(material_id[i]), elem->get_extra_integer(int_idx));
      }
  }
#endif

#if defined(LIBMESH_HAVE_EXODUS_API) && defined(LIBMESH_ENABLE_EXCEPTIONS)
  void testBadExtraIntegersExodusReading()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
    /*
     This 3-by-3 mesh contains the following element integers:
     material_id = '0 -1    (unsigned long long)(-2) (bad)
                    1  3    4503599627370496
                    2  3    4503599627370497'
     Real(-2) != Real(-1) (used for invalid_id), so we can tell
     this is an error.
    */
    const std::string filename = "meshes/bad_64bit_elem_integers.e";
    ExodusII_IO exreader(mesh);
    exreader.set_extra_integer_vars({"material_id"});
    CPPUNIT_ASSERT_THROW_MESSAGE("Bad elem integer not detected",
                                 exreader.read(filename),
                                 libMesh::LogicError);
  }
#endif
};


CPPUNIT_TEST_SUITE_REGISTRATION( ExtraIntegersTest );
