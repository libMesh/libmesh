// libmesh includes
#include <libmesh/boundary_info.h>
#include <libmesh/distributed_mesh.h>
#include <libmesh/elem.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_modification.h>
#include <libmesh/parallel_implementation.h>
#include <libmesh/node.h>
#include <libmesh/replicated_mesh.h>
#include <libmesh/utility.h>

// cppunit includes
#include "test_comm.h"
#include "libmesh_cppunit.h"

#include <algorithm>
#include <regex>

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
  CPPUNIT_TEST( testRemappingStitch );
  CPPUNIT_TEST( testAmbiguousRemappingStitch );
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
    const auto mesh_boundary_ids = boundary_info.get_global_boundary_ids();
    for (auto rit = mesh_boundary_ids.rbegin(); rit != mesh_boundary_ids.rend(); ++rit)
    {
      const auto old_sideset_name = boundary_info.sideset_name(*rit);
      const auto old_nodeset_name = boundary_info.nodeset_name(*rit);

      MeshTools::Modification::change_boundary_id(mesh, *rit, *rit + boundary_id_offset);

      boundary_info.sideset_name(*rit + boundary_id_offset) =
        boundary_name_prefix + old_sideset_name;
      boundary_info.nodeset_name(*rit + boundary_id_offset) =
        boundary_name_prefix + old_nodeset_name;
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

    CPPUNIT_ASSERT_EQUAL(mesh0.n_elem(),  static_cast<dof_id_type>(16));
    CPPUNIT_ASSERT_EQUAL(mesh0.n_nodes(), static_cast<dof_id_type>(45));

    const BoundaryInfo & bi = mesh0.get_boundary_info();
    std::set<boundary_id_type> sbi = bi.get_side_boundary_ids();
    TestCommWorld->set_union(sbi);

    typename std::decay<decltype(sbi.size())>::type expected_size = 10;
    CPPUNIT_ASSERT_EQUAL(expected_size, sbi.size());

    std::set<boundary_id_type> nbi = bi.get_node_boundary_ids();
    TestCommWorld->set_union(nbi);
    CPPUNIT_ASSERT_EQUAL(expected_size, nbi.size());

    // We expect that the "zero_right" and "one_left" boundaries have
    // disappeared after being stitched together.
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

    CPPUNIT_ASSERT_EQUAL(mesh0.n_elem(),  static_cast<dof_id_type>(32));
    CPPUNIT_ASSERT_EQUAL(mesh0.n_nodes(), static_cast<dof_id_type>(405));
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
          CPPUNIT_ASSERT_EQUAL(elem->get_extra_integer(foo0e_idx), static_cast<dof_id_type>(2));
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
            CPPUNIT_ASSERT_EQUAL(datum.first,  static_cast<dof_id_type>(3));
            CPPUNIT_ASSERT_EQUAL(datum.second, static_cast<dof_id_type>(4));
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

    // Before stitching, change the elemset codes on mesh1 so they
    // don't overlap with the codes on mesh0.
    mesh1->change_elemset_code(/*old*/1, /*new*/3); // 1 -> 3
    mesh1->change_elemset_code(/*old*/2, /*new*/4); // 2 -> 4

    // Before stitching, change the elemset ids on mesh1 so they
    // don't overlap with the elemset ids on mesh0.
    mesh1->change_elemset_id(/*old*/1, /*new*/100); // 1 -> 100
    mesh1->change_elemset_id(/*old*/2, /*new*/200); // 2 -> 200

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

    // Check that the stitched mesh has merged elemset codes and ids as expected
    MeshBase::elemset_type id_set_to_fill;
    const elemset_id_type code_to_type[] = {0,1,2,100,200};
    for (dof_id_type elemset_code=1; elemset_code<5; ++elemset_code)
      {
        mesh0->get_elemsets(elemset_code, id_set_to_fill);

        // Assert one elemset id in each set, and that set contains the correct id
        CPPUNIT_ASSERT_EQUAL(static_cast<std::size_t>(1), id_set_to_fill.size());
        CPPUNIT_ASSERT(id_set_to_fill.count(code_to_type[elemset_code]));
      }

    bool ps_odd = ps % 2;

    for (const auto & elem : mesh0->element_ptr_range())
      {
        dof_id_type elemset_code = elem->get_extra_integer(elemset_index);
        bool elem_id_odd = elem->id() % 2;

        // Debugging
        // libMesh::out << "Elem " << elem->id() << " in stitched mesh has elemset_code = " << elemset_code << std::endl;

        // Verify that the stitched mesh elemset codes match their pre-stitched values
        if (elem->id() < n_elem_prestitch) // lower half id
          {
            if (elem_id_odd)
              CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(1), elemset_code);
            else
              CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(2), elemset_code);
          }
        else // upper half id
          {
            // i.) If ps == odd, then n_elem_prestitch == odd, and even mesh1
            // elem ids will become odd, and odd mesh1 elem ids will become
            // even..
            // ii.) If ps == even, then n_elem_prestitch == even, and even mesh1
            // elem ids will remain even, odd mesh1 elem ids will remain odd.
            if (elem_id_odd)
              CPPUNIT_ASSERT_EQUAL(ps_odd ? static_cast<dof_id_type>(4) : static_cast<dof_id_type>(3), elemset_code);
            else
              CPPUNIT_ASSERT_EQUAL(ps_odd ? static_cast<dof_id_type>(3) : static_cast<dof_id_type>(4), elemset_code);
          }
      }
  }

  void testReplicatedMeshStitchElemsets ()
  {
    testMeshStitchElemsets<ReplicatedMesh>(/*ps=*/2);
    testMeshStitchElemsets<ReplicatedMesh>(/*ps=*/3);
  }


  void testRemappingStitch()
  {
    LOG_UNIT_TEST;

    Mesh mesh0(*TestCommWorld), mesh1(*TestCommWorld);

    int ps = 2;
    MeshTools::Generation::build_cube(mesh0, ps, ps, ps, -1, 0, 0, 1, 0, 1, HEX8);
    MeshTools::Generation::build_cube(mesh1, ps, ps, ps, 0, 1, 0, 1, 0, 1, HEX8);

    // rename and shift boundaries
    renameAndShift(mesh0, 0, "zero_");
    renameAndShift(mesh1, 6, "one_");

    // Create "auto" generated subdomain ids
    for (const auto & elem : mesh0.element_ptr_range())
      elem->subdomain_id() = 123;

    for (const auto & elem : mesh1.element_ptr_range())
      elem->subdomain_id() = 456;

    // Resolve them to the same name
    mesh0.subdomain_name(123) = "OneTwoThree";
    mesh1.subdomain_name(456) = "OneTwoThree"; // silly autogen

    mesh0.stitch_meshes(mesh1, 2, 10, TOLERANCE, true, false, false,
                        false, false, /* remap_subdomain_ids = */ true);

    CPPUNIT_ASSERT_EQUAL(mesh0.n_elem(),  static_cast<dof_id_type>(16));
    CPPUNIT_ASSERT_EQUAL(mesh0.n_nodes(), static_cast<dof_id_type>(45));

    // Ensure they still map to the same name but now with the same id
    for (const auto & elem : mesh0.element_ptr_range())
      CPPUNIT_ASSERT_EQUAL(elem->subdomain_id(), subdomain_id_type(123));
  }


  void testAmbiguousRemappingStitch()
  {
    LOG_UNIT_TEST;

    Mesh mesh0(*TestCommWorld), mesh1(*TestCommWorld);

    int ps = 2;
    MeshTools::Generation::build_cube(mesh0, ps, ps, ps, -1, 0, 0, 1, 0, 1, HEX8);
    MeshTools::Generation::build_cube(mesh1, ps, ps, ps, 0, 1, 0, 1, 0, 1, HEX8);

    // rename and shift boundaries
    renameAndShift(mesh0, 0, "zero_");
    renameAndShift(mesh1, 6, "one_");

    // Create matching subdomain ids
    for (const auto & elem : mesh0.element_ptr_range())
      elem->subdomain_id() = 123;

    for (const auto & elem : mesh1.element_ptr_range())
      elem->subdomain_id() = 123;

    // Create a conflict when only one is named
    mesh1.subdomain_name(123) = "OneTwoThree";

#ifdef LIBMESH_ENABLE_EXCEPTIONS
    bool threw_error = false;
    try
    {
      mesh0.stitch_meshes(mesh1, 2, 10, TOLERANCE, true, false, false,
                          false, false, /* remap_subdomain_ids = */ true);
    }
    catch (libMesh::LogicError & e)
    {
      std::regex msg_regex("safely stitch with a mesh");
      CPPUNIT_ASSERT(std::regex_search(e.what(), msg_regex));
      threw_error = true;
    }

    CPPUNIT_ASSERT(threw_error);
#endif // LIBMESH_ENABLE_EXCEPTIONS
  }

};

CPPUNIT_TEST_SUITE_REGISTRATION( MeshStitchTest );
