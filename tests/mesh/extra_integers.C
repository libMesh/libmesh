// Ignore unused parameter warnings coming from cppunit headers
#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

#include <libmesh/libmesh.h>
#include <libmesh/mesh.h>
#include <libmesh/elem.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_refinement.h>

#include "test_comm.h"

// THE CPPUNIT_TEST_SUITE_END macro expands to code that involves
// std::auto_ptr, which in turn produces -Wdeprecated-declarations
// warnings.  These can be ignored in GCC as long as we wrap the
// offending code in appropriate pragmas.  We can't get away with a
// single ignore_warnings.h inclusion at the beginning of this file,
// since the libmesh headers pull in a restore_warnings.h at some
// point.  We also don't bother restoring warnings at the end of this
// file since it's not a header.
#include <libmesh/ignore_warnings.h>

using namespace libMesh;

class ExtraIntegersTest : public CppUnit::TestCase
{
  /**
   * The goal of this test is to verify the ability to add extra
   * integer storage to a Mesh object, then set and query those extra
   * integers in the objects within the mesh.
   */
public:
  CPPUNIT_TEST_SUITE( ExtraIntegersTest );

  CPPUNIT_TEST( testExtraIntegersEdge2 );

  CPPUNIT_TEST( testExtraIntegersTri6 );

  CPPUNIT_TEST( testExtraIntegersCheckpointEdge3 );

  CPPUNIT_TEST( testExtraIntegersCheckpointHex8 );

  CPPUNIT_TEST_SUITE_END();

protected:
  // Helper functions called by the test implementations, saves a few lines of code.
  std::array<unsigned int, 3>
    build_mesh(Mesh & mesh, ElemType elem_type, unsigned int n_elem_per_side)
  {
    // Request some extra integers before building
    unsigned int i1 = mesh.add_elem_integer("i1");
    unsigned int ni1 = mesh.add_node_integer("ni1");
    unsigned int ni2 = mesh.add_node_integer("ni2");

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
    return {i1, ni1, ni2};
  }


  void test_and_set_initial_integers
    (Mesh & mesh, unsigned int i1, unsigned int ni1, unsigned int ni2)
  {
    for (const auto & elem : mesh.element_ptr_range())
      {
        CPPUNIT_ASSERT_EQUAL(elem->n_extra_integers(), 1u);
        CPPUNIT_ASSERT_EQUAL(elem->get_extra_integer(i1), DofObject::invalid_id);
        elem->set_extra_integer(i1, dof_id_type(elem->point(0)(0)*100));
        CPPUNIT_ASSERT_EQUAL(elem->get_extra_integer(i1), dof_id_type(elem->point(0)(0)*100));
      }

    for (const auto & node : mesh.node_ptr_range())
      {
        CPPUNIT_ASSERT_EQUAL(node->n_extra_integers(), 2u);
        CPPUNIT_ASSERT_EQUAL(node->get_extra_integer(ni1), DofObject::invalid_id);
        CPPUNIT_ASSERT_EQUAL(node->get_extra_integer(ni2), DofObject::invalid_id);
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
            CPPUNIT_ASSERT_EQUAL(node.n_extra_integers(), 2u);

        CPPUNIT_ASSERT_EQUAL(elem->n_extra_integers(), 1u);
        CPPUNIT_ASSERT_EQUAL(elem->get_extra_integer(i1), dof_id_type(top_parent->point(0)(0)*100));
      }
  }


  void test_helper(ElemType elem_type, unsigned int n_elem_per_side)
  {
    Mesh mesh(*TestCommWorld);

    std::array<unsigned int, 3> ini = build_mesh(mesh, elem_type, n_elem_per_side);
    const unsigned int i1 = ini[0], ni1 = ini[1], ni2 = ini[2];

    test_and_set_initial_integers(mesh, i1, ni1, ni2);

    // Force (in parallel) a different partitioning - we'll simply put
    // everything on rank 0, which hopefully is not what our default
    // partitioner did!
    mesh.partition(1);

    CPPUNIT_ASSERT_EQUAL(i1, mesh.add_elem_integer("i1"));
    CPPUNIT_ASSERT_EQUAL(ni1, mesh.add_node_integer("ni1"));
    CPPUNIT_ASSERT_EQUAL(ni2, mesh.add_node_integer("ni2"));

    // Make sure we didn't screw up any extra integers thereby.
    test_final_integers(mesh, i1);

#ifdef LIBMESH_ENABLE_AMR
    MeshRefinement mr(mesh);
    mr.uniformly_refine(1);

    test_final_integers(mesh, i1);
#endif
  }

  void checkpoint_helper(ElemType elem_type, unsigned int n_elem_per_side, bool binary)
  {
    Mesh mesh(*TestCommWorld);

    std::array<unsigned int, 3> ini = build_mesh(mesh, elem_type, n_elem_per_side);
    const unsigned int i1 = ini[0], ni1 = ini[1], ni2 = ini[2];

    test_and_set_initial_integers(mesh, i1, ni1, ni2);

    const std::string filename =
      std::string("extra_integers.cp") + (binary ? "r" : "a");

    mesh.write(filename);

    TestCommWorld->barrier();

    Mesh mesh2(*TestCommWorld);

    mesh2.read(filename);

    // Make sure everything got transferred to the second mesh
    CPPUNIT_ASSERT_EQUAL(i1, mesh2.add_elem_integer("i1"));
    CPPUNIT_ASSERT_EQUAL(ni1, mesh2.add_node_integer("ni1"));
    CPPUNIT_ASSERT_EQUAL(ni2, mesh2.add_node_integer("ni2"));

    // Make sure we didn't screw up any extra integers thereby.
    test_final_integers(mesh2, i1);
  }

public:
  void setUp() {}

  void tearDown() {}

  void testExtraIntegersEdge2() { test_helper(EDGE2, 5); }

  void testExtraIntegersTri6() { test_helper(TRI6, 4); }

  void testExtraIntegersCheckpointEdge3() { checkpoint_helper(EDGE3, 5, false); }

  void testExtraIntegersCheckpointHex8() { checkpoint_helper(HEX8, 2, true); }

};


CPPUNIT_TEST_SUITE_REGISTRATION( ExtraIntegersTest );
