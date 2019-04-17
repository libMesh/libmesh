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

  CPPUNIT_TEST_SUITE_END();

protected:
  // Helper function called by the test implementations, saves a few lines of code.
  void test_helper_1D(ElemType elem_type)
  {
    Mesh mesh(*TestCommWorld, /*dim=*/2);

    // Request some extra integers before building
    unsigned int i1 = mesh.add_elem_integer("i1");
    unsigned int ni1 = mesh.add_node_integer("ni1");
    unsigned int ni2 = mesh.add_node_integer("ni2");

    MeshTools::Generation::build_line(mesh,
                                      /*nx=*/10,
                                      /*xmin=*/0., /*xmax=*/1.,
                                      elem_type);

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

    // Force (in parallel) a different partitioning - we'll simply put
    // everything on rank 0, which hopefully is not what our default
    // partitioner did!
    mesh.partition(1);

    CPPUNIT_ASSERT_EQUAL(i1, mesh.add_elem_integer("i1"));
    CPPUNIT_ASSERT_EQUAL(ni1, mesh.add_node_integer("ni1"));
    CPPUNIT_ASSERT_EQUAL(ni2, mesh.add_node_integer("ni2"));

    // Make sure we didn't screw up any extra integers thereby.
    for (const auto & elem : mesh.element_ptr_range())
      {
        CPPUNIT_ASSERT_EQUAL(elem->n_extra_integers(), 1u);
        CPPUNIT_ASSERT_EQUAL(elem->get_extra_integer(i1), dof_id_type(elem->point(0)(0)*100));
      }

    for (const auto & node : mesh.node_ptr_range())
      {
        CPPUNIT_ASSERT_EQUAL(node->n_extra_integers(), 2u);
        CPPUNIT_ASSERT_EQUAL(node->get_extra_integer(ni1), DofObject::invalid_id);
        CPPUNIT_ASSERT_EQUAL(node->get_extra_integer(ni2), DofObject::invalid_id);
      }

#ifdef LIBMESH_ENABLE_AMR
    MeshRefinement mr(mesh);
    mr.uniformly_refine(1);

    // Make sure the old elements have the same integers and the new
    // elements have be initialized as undefined.
    for (const auto & elem : mesh.element_ptr_range())
      {
        CPPUNIT_ASSERT_EQUAL(elem->n_extra_integers(), 1u);
        CPPUNIT_ASSERT_EQUAL(elem->get_extra_integer(i1), dof_id_type(elem->top_parent()->point(0)(0)*100));
        if (!elem->level())
          for (auto & node : elem->node_ref_range())
            CPPUNIT_ASSERT_EQUAL(node.n_extra_integers(), 2u);
      }
#endif
  }

public:
  void setUp() {}

  void tearDown() {}

  void testExtraIntegersEdge2() { test_helper_1D(EDGE2); }
};


CPPUNIT_TEST_SUITE_REGISTRATION( ExtraIntegersTest );
