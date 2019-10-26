#include <libmesh/elem.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/node.h>
#include <libmesh/replicated_mesh.h>

#include "test_comm.h"
#include "libmesh_cppunit.h"


using namespace libMesh;


class MeshStitchTest : public CppUnit::TestCase {
public:
  CPPUNIT_TEST_SUITE( MeshStitchTest );

#if LIBMESH_DIM > 2
  CPPUNIT_TEST( testMeshStitch );
#endif // LIBMESH_DIM > 2

  CPPUNIT_TEST_SUITE_END();

private:

public:
  void setUp()
  {}

  void tearDown()
  {}

  void testMeshStitch ()
  {
    // Generate four meshes to be stitched together
    ReplicatedMesh mesh0(*TestCommWorld),
                   mesh1(*TestCommWorld),
                   mesh2(*TestCommWorld),
                   mesh3(*TestCommWorld);

    // Give the meshes different extra integers to make sure those
    // merge.  Reuse names between nodes and elements to make sure
    // those don't mix.  Add some integers before and others after
    // generation to test flexibility there.

    std::vector<std::string> names2 {"bar", "baz"};
    mesh2.add_elem_integers(names2);

    std::vector<std::string> names3 {"bar", "foo"};
    mesh3.add_elem_integers(names3);

    int ps = 2;
    MeshTools::Generation::build_cube (mesh0, ps, ps, ps, -1,    0,    0,  1,  0, 1, HEX27);
    MeshTools::Generation::build_cube (mesh1, ps, ps, ps,    0,  1,    0,  1,  0, 1, HEX27);
    MeshTools::Generation::build_cube (mesh2, ps, ps, ps, -1,    0, -1,    0,  0, 1, HEX27);
    MeshTools::Generation::build_cube (mesh3, ps, ps, ps,    0,  1, -1,    0,  0, 1, HEX27);

    struct trivially_copyable_pair // std::pair triggers -Wclass-memaccess
    {
      dof_id_type first, second;
    };

    mesh0.add_node_integer("baz");
    unsigned int foo1e_idx = mesh1.add_elem_integer("foo");
    mesh2.add_elem_datum<trivially_copyable_pair>("qux");
    unsigned int qux2n_idx = mesh2.add_node_datum<trivially_copyable_pair>("qux");
    mesh3.add_node_integers(names3);

    for (const auto & elem : mesh1.element_ptr_range())
      elem->set_extra_integer(foo1e_idx, 2);

    for (const auto & node : mesh2.node_ptr_range())
      node->set_extra_datum<trivially_copyable_pair>
        (qux2n_idx, {3, 4});

    // We stitch the meshes in a hierarchical way.
    mesh0.stitch_meshes(mesh1, 2, 4, TOLERANCE, true, true, false, false);
    mesh2.stitch_meshes(mesh3, 2, 4, TOLERANCE, true, true, false, false);
    mesh0.stitch_meshes(mesh2, 1, 3, TOLERANCE, true, true, false, false);

    CPPUNIT_ASSERT_EQUAL(mesh0.n_elem(), dof_id_type(32));
    CPPUNIT_ASSERT_EQUAL(mesh0.n_nodes(), dof_id_type(405));
    CPPUNIT_ASSERT_EQUAL(mesh0.n_elem_integers(), 5u); // that pair counts 2x
    CPPUNIT_ASSERT_EQUAL(mesh0.n_node_integers(), 5u);
    std::vector<std::string> all_names {"foo", "bar", "baz", "qux"};
    std::vector<unsigned int> node_name_indices {4, 3, 0, 1};
    for (unsigned int i=0; i != 4; ++i)
      {
        CPPUNIT_ASSERT(mesh0.has_elem_integer(all_names[i]));
        CPPUNIT_ASSERT_EQUAL(mesh0.get_elem_integer_index(all_names[i]), i);
        CPPUNIT_ASSERT(mesh0.has_node_integer(all_names[i]));
        CPPUNIT_ASSERT_EQUAL(mesh0.get_node_integer_index(all_names[i]), node_name_indices[i]);
      }

    unsigned int foo0e_idx = mesh0.get_elem_integer_index("foo");
    for (const auto & elem : mesh0.element_ptr_range())
      {
        CPPUNIT_ASSERT_EQUAL(elem->n_extra_integers(), 5u);
        const Point c = elem->centroid();
        if (c(0) > 0 && c(1) > 0) // this came from mesh1
          CPPUNIT_ASSERT_EQUAL(elem->get_extra_integer(foo0e_idx), dof_id_type(2));
        else
          CPPUNIT_ASSERT_EQUAL(elem->get_extra_integer(foo0e_idx), DofObject::invalid_id);
      }

    unsigned int qux0n_idx = mesh0.get_node_integer_index("qux");
    for (const auto & node : mesh0.node_ptr_range())
      {
        CPPUNIT_ASSERT_EQUAL(node->n_extra_integers(), 5u);
        trivially_copyable_pair datum =
          node->get_extra_datum<trivially_copyable_pair>(qux0n_idx);
        if ((*node)(0) <= 0 && (*node)(1) < 0) // this came from mesh2
          {
            CPPUNIT_ASSERT_EQUAL(datum.first, dof_id_type(3));
            CPPUNIT_ASSERT_EQUAL(datum.second, dof_id_type(4));
          }
        else
          {
            CPPUNIT_ASSERT_EQUAL(datum.first, DofObject::invalid_id);
            CPPUNIT_ASSERT_EQUAL(datum.second, DofObject::invalid_id);
          }
      }
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION( MeshStitchTest );
