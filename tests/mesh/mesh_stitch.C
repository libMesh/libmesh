#include <libmesh/elem.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/node.h>
#include <libmesh/replicated_mesh.h>
#include <libmesh/mesh_modification.h>
#include <libmesh/boundary_info.h>

#include "test_comm.h"
#include "libmesh_cppunit.h"

#include <algorithm>

using namespace libMesh;


class MeshStitchTest : public CppUnit::TestCase {
public:
  LIBMESH_CPPUNIT_TEST_SUITE( MeshStitchTest );

#if LIBMESH_DIM > 2
  CPPUNIT_TEST( testMeshStitch );
  CPPUNIT_TEST( testBoundaryInfo );
#endif // LIBMESH_DIM > 2

  CPPUNIT_TEST_SUITE_END();

private:

public:
  void setUp()
  {}

  void tearDown()
  {}

  void renameAndShift(ReplicatedMesh & mesh,
                      const boundary_id_type boundary_id_offset,
                      const std::string & boundary_name_prefix)
  {
    BoundaryInfo & boundary_info = mesh.get_boundary_info();
    const auto & mesh_boundary_ids = boundary_info.get_boundary_ids();
    for (auto rit = mesh_boundary_ids.rbegin(); rit != mesh_boundary_ids.rend(); ++rit)
    {
      boundary_info.sideset_name(*rit + boundary_id_offset) =
          boundary_name_prefix + boundary_info.sideset_name(*rit);
      boundary_info.nodeset_name(*rit + boundary_id_offset) =
          boundary_name_prefix + boundary_info.nodeset_name(*rit);

      MeshTools::Modification::change_boundary_id(mesh, *rit, *rit + boundary_id_offset);
    }
  }

  void testBoundaryInfo()
  {
    LOG_UNIT_TEST;

    ReplicatedMesh mesh0(*TestCommWorld), mesh1(*TestCommWorld);

    int ps = 2;
    MeshTools::Generation::build_cube(mesh0, ps, ps, ps, -1, 0, 0, 1, 0, 1, HEX8);
    MeshTools::Generation::build_cube(mesh1, ps, ps, ps, 0, 1, 0, 1, 0, 1, HEX8);

    // rename and shift boundaries
    renameAndShift(mesh0, 0, "zero_");
    renameAndShift(mesh1, 6, "one_");

    mesh0.stitch_meshes(mesh1, 2, 10, TOLERANCE, true, true, false, false);

    CPPUNIT_ASSERT_EQUAL(mesh0.n_elem(), dof_id_type(16));
    CPPUNIT_ASSERT_EQUAL(mesh0.n_nodes(), dof_id_type(45));

    const BoundaryInfo & bi = mesh0.get_boundary_info();
    const auto & sbi = bi.get_side_boundary_ids();
    typename std::decay<decltype(sbi.size())>::type expected_size = 10;
    CPPUNIT_ASSERT_EQUAL(expected_size, sbi.size());

    const auto & nbi = bi.get_node_boundary_ids();
    CPPUNIT_ASSERT_EQUAL(expected_size, nbi.size());

    const auto & ss_id_to_name = bi.get_sideset_name_map();
    std::set<std::string> ss_names;
    std::for_each(ss_id_to_name.begin(),
                  ss_id_to_name.end(),
                  [&ss_names](const std::pair<boundary_id_type, std::string> & map_pr) {
                    ss_names.insert(map_pr.second);
                  });
    std::set<std::string> expected_names = {{"zero_left",
                                             "zero_top",
                                             "zero_front",
                                             "zero_back",
                                             "zero_bottom",
                                             "one_right",
                                             "one_top",
                                             "one_front",
                                             "one_back",
                                             "one_bottom"}};
    CPPUNIT_ASSERT(ss_names == expected_names);

    const auto & ns_id_to_name = bi.get_nodeset_name_map();
    std::set<std::string> ns_names;
    std::for_each(ns_id_to_name.begin(),
                  ns_id_to_name.end(),
                  [&ns_names](const std::pair<boundary_id_type, std::string> & map_pr) {
                    ns_names.insert(map_pr.second);
                  });
    CPPUNIT_ASSERT(ns_names == expected_names);
  }

  void testMeshStitch ()
  {
    LOG_UNIT_TEST;

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
        const Point c = elem->vertex_average();
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
