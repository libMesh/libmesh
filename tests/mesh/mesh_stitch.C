// libmesh includes
#include <libmesh/boundary_info.h>
#include <libmesh/distributed_mesh.h>
#include <libmesh/elem.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_modification.h>
#include <libmesh/node.h>
#include <libmesh/replicated_mesh.h>
#include <libmesh/utility.h>

// cppunit includes
#include "test_comm.h"
#include "libmesh_cppunit.h"

#include <algorithm>

using namespace libMesh;


class MeshStitchTest : public CppUnit::TestCase {
public:
  LIBMESH_CPPUNIT_TEST_SUITE( MeshStitchTest );

#if LIBMESH_DIM > 2
  CPPUNIT_TEST( testReplicatedMeshStitch );
  CPPUNIT_TEST( testDistributedMeshStitch );
  CPPUNIT_TEST( testReplicatedBoundaryInfo );
  CPPUNIT_TEST( testDistributedBoundaryInfo );
  CPPUNIT_TEST( testReplicatedMeshStitchElemsets );
#endif // LIBMESH_DIM > 2

  CPPUNIT_TEST_SUITE_END();

private:

public:
  void setUp()
  {}

  void tearDown()
  {}

  void renameAndShift(UnstructuredMesh & mesh,
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


  template <typename MeshType>
  void testBoundaryInfo()
  {
    LOG_UNIT_TEST;

    MeshType mesh0(*TestCommWorld), mesh1(*TestCommWorld);

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
    std::set<std::string> ss_names;
    for (const auto & pr : bi.get_sideset_name_map())
      ss_names.insert(pr.second);
    CPPUNIT_ASSERT(ss_names == expected_names);

    std::set<std::string> ns_names;
    for (const auto & pr : bi.get_nodeset_name_map())
      ns_names.insert(pr.second);
    CPPUNIT_ASSERT(ns_names == expected_names);
  }


  void testReplicatedBoundaryInfo()
  {
    testBoundaryInfo<ReplicatedMesh>();
  }


  void testDistributedBoundaryInfo()
  {
    testBoundaryInfo<DistributedMesh>();
  }


  template <typename MeshType>
  void testMeshStitch ()
  {
    LOG_UNIT_TEST;

    // Generate four meshes to be stitched together
    MeshType mesh0(*TestCommWorld),
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
    mesh2.template add_elem_datum<trivially_copyable_pair>("qux");
    unsigned int qux2n_idx = mesh2.template add_node_datum<trivially_copyable_pair>("qux");
    mesh3.add_node_integers(names3);

    for (const auto & elem : mesh1.element_ptr_range())
      elem->set_extra_integer(foo1e_idx, 2);

    for (const auto & node : mesh2.node_ptr_range())
      node->template set_extra_datum<trivially_copyable_pair>
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
          node->template get_extra_datum<trivially_copyable_pair>(qux0n_idx);
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

  void testReplicatedMeshStitch ()
  {
    testMeshStitch<ReplicatedMesh>();
  }

  void testDistributedMeshStitch ()
  {
    testMeshStitch<DistributedMesh>();
  }

  template <typename MeshType>
  void testMeshStitchElemsets (unsigned int ps)
  {
    LOG_UNIT_TEST;

    // Generate meshes to be stitched together. We are going to clone
    // these so work with unique_ptrs directly.
    auto mesh0 = std::make_unique<MeshType>(*TestCommWorld);

    // If the user tries to stitch meshes with overlapping codes, we
    // allow this as long as the codes refer to the same underlying
    // set ids.

    // Build a mesh on the unit cube
    MeshTools::Generation::build_cube (*mesh0, ps, ps, ps,
                                       /*xmin=*/0., /*xmax=*/1.,
                                       /*ymin=*/0., /*ymax=*/1.,
                                       /*zmin=*/0., /*zmax=*/1.,
                                       HEX27);

    // Make a copy
    auto mesh1 = mesh0->clone();

    // Shift copy one unit to the right
    MeshTools::Modification::translate(*mesh1, /*x-dir*/1.0);

    // For both meshes:
    // .) Put odd-numbered Elems in elmset 1
    // .) Put even-numbered Elems in elemset 2
    // We use the trivial encoding: elemset id == elemset code for simplicity
    auto place_elems = [](MeshBase & mesh)
      {
        unsigned int elemset_index =
          mesh.add_elem_integer("elemset_code", /*allocate_data=*/true);

        mesh.add_elemset_code(/*code=*/1, /*set_ids*/{1});
        mesh.add_elemset_code(/*code=*/2, /*set_ids*/{2});

        for (const auto & elem : mesh.element_ptr_range())
          {
            if (elem->id() % 2) // id odd
              elem->set_extra_integer(elemset_index, 1);
            else // id even
              elem->set_extra_integer(elemset_index, 2);
          }
      };

    place_elems(*mesh0);
    place_elems(*mesh1);

    // Stitch the meshes together at the indicated boundary ids
    mesh0->stitch_meshes(dynamic_cast<UnstructuredMesh &>(*mesh1),
                         /*this boundary=*/2,
                         /*other boundary=*/4,
                         TOLERANCE,
                         /*clear_stitched_boundary_ids=*/true,
                         /*verbose=*/true,
                         /*use_binary_search=*/false,
                         /*enforce_all_nodes_match_on_boundaries=*/false);

    // Number of elements in each mesh pre-stitch
    dof_id_type n_elem_prestitch = Utility::pow<3>(ps);

    // mesh0 should contain 2 * ps**3 total elements after stitching
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(2 * n_elem_prestitch), mesh0->n_elem());

    // Check that the stitched mesh still stores "elemset_code" in the
    // same index (0) as it was before the meshes were stitched.
    unsigned int elemset_index = mesh0->get_elem_integer_index("elemset_code");
    CPPUNIT_ASSERT_EQUAL(0u, elemset_index);

    // There should still be 2 elemset ids, 1 and 2, on the stitched mesh
    CPPUNIT_ASSERT_EQUAL(2u, mesh0->n_elemsets());

    MeshBase::elemset_type id_set_to_fill;

    // Make sure that elemset_code 1 == {1} and 2 == {2} on the new mesh
    for (int i=1; i<3; ++i)
      {
        mesh0->get_elemsets(static_cast<dof_id_type>(i), id_set_to_fill);
        CPPUNIT_ASSERT_EQUAL(static_cast<std::size_t>(1), id_set_to_fill.size());
        CPPUNIT_ASSERT_EQUAL(static_cast<elemset_id_type>(i), *id_set_to_fill.begin());
      }

    // Debugging
    // libMesh::out << "In stitched mesh, elemset_index = " << elemset_index << std::endl;

    // The elems from mesh1 will all have n_elem_prestitch added to their ids.
    for (const auto & elem : mesh0->element_ptr_range())
      {
        dof_id_type elemset_code = elem->get_extra_integer(elemset_index);

        // Debugging
        // libMesh::out << "Elem " << elem->id() << " in stitched mesh has elemset_code = " << elemset_code << std::endl;

        if (ps % 2) // ps == odd
          {
            // i.) If ps == odd, then n_elem_prestitch == odd, and even mesh1
            // elem ids will become odd, and odd mesh1 elem ids will become
            // even..
            if (elem->id() < n_elem_prestitch) // lower half id
              {
                if (elem->id() % 2) // id odd
                  CPPUNIT_ASSERT_EQUAL(1u, elemset_code);
                else // id even
                  CPPUNIT_ASSERT_EQUAL(2u, elemset_code);
              }
            else // upper half id (sets reversed)
              {
                if (elem->id() % 2) // id odd
                  CPPUNIT_ASSERT_EQUAL(2u, elemset_code);
                else // id even
                  CPPUNIT_ASSERT_EQUAL(1u, elemset_code);
              }
          }
        else // ps == even
          {
            // ii.) If ps == even, then n_elem_prestitch == even, and even mesh1
            // elem ids will remain even, odd mesh1 elem ids will remain odd.
            if (elem->id() % 2) // id odd
              CPPUNIT_ASSERT_EQUAL(1u, elemset_code);
            else // id even
              CPPUNIT_ASSERT_EQUAL(2u, elemset_code);
          }
      }
  }

  void testReplicatedMeshStitchElemsets ()
  {
    testMeshStitchElemsets<ReplicatedMesh>(/*ps=*/2);
    testMeshStitchElemsets<ReplicatedMesh>(/*ps=*/3);
  }


};

CPPUNIT_TEST_SUITE_REGISTRATION( MeshStitchTest );
