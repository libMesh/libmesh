// libMesh includes
#include <libmesh/libmesh.h>
#include <libmesh/distributed_mesh.h>
#include <libmesh/replicated_mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh.h>
#include <libmesh/point.h>
#include <libmesh/elem.h>

// cppunit includes
#include "test_comm.h"
#include "libmesh_cppunit.h"

using namespace libMesh;

class AllSecondOrderTest : public CppUnit::TestCase
{
public:
  LIBMESH_CPPUNIT_TEST_SUITE( AllSecondOrderTest );

  CPPUNIT_TEST( allSecondOrder );
  CPPUNIT_TEST( allSecondOrderRange );
  CPPUNIT_TEST( allSecondOrderDoNothing );
  CPPUNIT_TEST( allSecondOrderMixed );
  CPPUNIT_TEST( allSecondOrderMixedFixing );
  CPPUNIT_TEST( allSecondOrderMixedFixing3D );

  CPPUNIT_TEST( allCompleteOrder );
  CPPUNIT_TEST( allCompleteOrderRange );
  CPPUNIT_TEST( allCompleteOrderDoNothing );
  CPPUNIT_TEST( allCompleteOrderMixed );
  CPPUNIT_TEST( allCompleteOrderMixedFixing );
  CPPUNIT_TEST( allCompleteOrderMixedFixing3D );

  CPPUNIT_TEST_SUITE_END();

public:
  void setUp() {}

  void tearDown() {}

  void allSecondOrder()
  {
    LOG_UNIT_TEST;

    DistributedMesh mesh(*TestCommWorld, /*dim=*/2);

    mesh.allow_remote_element_removal(false);

    MeshTools::Generation::build_square(mesh, 2, 2);

    mesh.all_second_order();
  }

  void allSecondOrderRange()
  {
    ReplicatedMesh mesh(*TestCommWorld);

    // Construct a single-element Quad4 mesh
    MeshTools::Generation::build_square(mesh, /*nx=*/1, /*ny=*/1);

    // Add extra nodes above the four original nodes
    Point z_trans(0,0,1);
    for (dof_id_type n=0; n<4; ++n)
      mesh.add_point(mesh.point(n) + z_trans, /*id=*/4+n, /*proc_id=*/0);

    // Construct Edge2 elements attached to Quad4 at individual points
    // and in subdomain 1.
    for (dof_id_type n=0; n<4; ++n)
      {
        Elem * elem = mesh.add_elem(Elem::build(EDGE2));
        elem->set_node(0) = mesh.node_ptr(n);
        elem->set_node(1) = mesh.node_ptr(n+4);
        elem->subdomain_id() = 1;
      }

    // Convert only Edge2 elements (all elements in subdomain 1) to
    // SECOND-order.
    mesh.all_second_order_range(mesh.active_subdomain_elements_ptr_range(1),
                                /*full_ordered=*/true);

    // Make sure we still have the expected total number of elements
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(5), mesh.n_elem());

    // Make sure we added the right number of nodes
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(12), mesh.n_nodes());

    // Make sure that the type of the Quad4 is unchanged
    CPPUNIT_ASSERT_EQUAL(QUAD4, mesh.elem_ptr(0)->type());

    // Make sure that the other elements are now Edge3s
    for (dof_id_type e=1; e<5; ++e)
      CPPUNIT_ASSERT_EQUAL(EDGE3, mesh.elem_ptr(e)->type());
  }

  void allSecondOrderDoNothing()
  {
    DistributedMesh mesh(*TestCommWorld, /*dim=*/2);

    mesh.allow_remote_element_removal(false);

    MeshTools::Generation::build_square(mesh,
                                        /*nx=*/2, /*ny=*/2,
                                        /*xmin=*/0., /*xmax=*/1.,
                                        /*ymin=*/0., /*ymax=*/1.,
                                        QUAD9);

    // Test that all_second_order_range() correctly "does nothing" in
    // parallel when passed ranges of local elements. Here we should
    // hit one of the "early return" cases for this function.
    mesh.all_second_order_range(mesh.active_local_element_ptr_range(),
                                /*full_ordered=*/true);

    // Make sure we still have the same number of elements
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(4), mesh.n_elem());

    // Make sure we still have the same number of nodes
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(25), mesh.n_nodes());
  }

  void allSecondOrderMixed()
  {
    ReplicatedMesh mesh(*TestCommWorld);

    // Construct a single-element Quad9 mesh
    MeshTools::Generation::build_square(mesh,
                                        /*nx=*/1, /*ny=*/1,
                                        /*xmin=*/0., /*xmax=*/1.,
                                        /*ymin=*/0., /*ymax=*/1.,
                                        /*elem_type=*/QUAD9);

    // Pointer to the single Elem
    const Elem * elem = mesh.elem_ptr(0);

    // Amount to offset the indices of newly-added nodes by
    const unsigned int node_offset = elem->n_nodes();

    // Add an extra node above each of the original vertices
    for (dof_id_type n=0; n<elem->n_vertices(); ++n)
      mesh.add_point(mesh.point(n) + Point(0,0,1), /*id=*/node_offset + n, /*proc_id=*/0);

    // Construct Edge2 elements and attach to each vertex of Quad9
    for (dof_id_type n=0; n<elem->n_vertices(); ++n)
      {
        Elem * elem = mesh.add_elem(Elem::build(EDGE2));
        elem->set_node(0) = mesh.node_ptr(n);
        elem->set_node(1) = mesh.node_ptr(node_offset + n);
      }

    // Convert all elements to SECOND-order, automatically skipping those that are
    // already SECOND-order.
    mesh.all_second_order_range(mesh.element_ptr_range(), /*full_ordered=*/true);

    // Make sure we still have the expected total number of elements
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(5), mesh.n_elem());

    // Make sure we have the correct number of nodes, 9+4+4=17
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(17), mesh.n_nodes());

    // Make sure that the type of the original Elem is unchanged
    CPPUNIT_ASSERT_EQUAL(QUAD9, elem->type());

    // Make sure that the other elements are now Edge3s
    for (dof_id_type e=1; e<5; ++e)
      CPPUNIT_ASSERT_EQUAL(EDGE3, mesh.elem_ptr(e)->type());
  }

  void allSecondOrderMixedFixing()
  {
    Mesh mesh(*TestCommWorld);

    // Construct a multi-element Quad4 mesh
    MeshTools::Generation::build_square(mesh,
                                        /*nx=*/3, /*ny=*/3,
                                        /*xmin=*/0., /*xmax=*/1.,
                                        /*ymin=*/0., /*ymax=*/1.,
                                        /*elem_type=*/QUAD4);

    // Convert elements to SECOND-order, but do so a few at a time,
    // leaving a broken mesh after initial conversions to see if it
    // will be fixed by subsequent conversions.
    int flipcnt = 0;
    for (auto it = mesh.elements_begin(), end = mesh.elements_end(); it != end; ++it)
      {
        auto range_start = it, range_end = it;
        ++range_end;
        if (flipcnt++%2 && range_end != end)
          ++range_end;

        mesh.all_second_order_range({range_start, range_end},
                                    /*full_ordered=*/true);
      }

    // Make sure we still have the expected total number of elements
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(9), mesh.n_elem());

    // Make sure we have the correct number of nodes, 7*7
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(49), mesh.n_nodes());

    // Make sure that the elements are now upgraded
    for (const auto & elem : mesh.element_ptr_range())
      CPPUNIT_ASSERT_EQUAL(QUAD9, elem->type());
  }

  void allSecondOrderMixedFixing3D()
  {
    Mesh mesh(*TestCommWorld);

    // Construct a multi-element Hex8 mesh
    MeshTools::Generation::build_cube(mesh,
                                      /*nx=*/3, /*ny=*/3, /*nz=*/3,
                                      /*xmin=*/0., /*xmax=*/1.,
                                      /*ymin=*/0., /*ymax=*/1.,
                                      /*zmin=*/0., /*zmax=*/1.,
                                      /*elem_type=*/HEX8);

    // Convert elements to SECOND-order, but do so a few at a time,
    // leaving a broken mesh after initial conversions to see if it
    // will be fixed by subsequent conversions.
    int flipcnt = 0;
    for (auto it = mesh.elements_begin(), end = mesh.elements_end(); it != end; ++it)
      {
        auto range_start = it, range_end = it;
        ++range_end;
        if (flipcnt++%2 && range_end != end)
          ++range_end;

        mesh.all_second_order_range({range_start, range_end},
                                    /*full_ordered=*/true);
      }

    // Make sure we still have the expected total number of elements
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(27), mesh.n_elem());

    // Make sure we have the correct number of nodes, 7*7*7
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(343), mesh.n_nodes());

    // Make sure that the elements are now upgraded
    for (const auto & elem : mesh.element_ptr_range())
      CPPUNIT_ASSERT_EQUAL(HEX27, elem->type());
  }


  void allCompleteOrder()
  {
    LOG_UNIT_TEST;

    DistributedMesh mesh(*TestCommWorld, /*dim=*/2);

    mesh.allow_remote_element_removal(false);

    MeshTools::Generation::build_square(mesh, 2, 2);

    mesh.all_complete_order();
  }

  void allCompleteOrderRange()
  {
    ReplicatedMesh mesh(*TestCommWorld);

    // Construct a single-element Quad4 mesh
    MeshTools::Generation::build_square(mesh, /*nx=*/1, /*ny=*/1);

    // Add extra nodes above the four original nodes
    Point z_trans(0,0,1);
    for (dof_id_type n=0; n<4; ++n)
      mesh.add_point(mesh.point(n) + z_trans, /*id=*/4+n, /*proc_id=*/0);

    // Construct Edge2 elements attached to Quad4 at individual points
    // and in subdomain 1.
    for (dof_id_type n=0; n<4; ++n)
      {
        Elem * elem = mesh.add_elem(Elem::build(EDGE2));
        elem->set_node(0) = mesh.node_ptr(n);
        elem->set_node(1) = mesh.node_ptr(n+4);
        elem->subdomain_id() = 1;
      }

    // Convert only Edge2 elements (all elements in subdomain 1) to
    // higher-order.
    mesh.all_complete_order_range(mesh.active_subdomain_elements_ptr_range(1));

    // Make sure we still have the expected total number of elements
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(5), mesh.n_elem());

    // Make sure we added the right number of nodes
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(12), mesh.n_nodes());

    // Make sure that the type of the Quad4 is unchanged
    CPPUNIT_ASSERT_EQUAL(QUAD4, mesh.elem_ptr(0)->type());

    // Make sure that the other elements are now Edge3s
    for (dof_id_type e=1; e<5; ++e)
      CPPUNIT_ASSERT_EQUAL(EDGE3, mesh.elem_ptr(e)->type());
  }

  void allCompleteOrderDoNothing()
  {
    DistributedMesh mesh(*TestCommWorld, /*dim=*/2);

    mesh.allow_remote_element_removal(false);

    MeshTools::Generation::build_square(mesh,
                                        /*nx=*/2, /*ny=*/2,
                                        /*xmin=*/0., /*xmax=*/1.,
                                        /*ymin=*/0., /*ymax=*/1.,
                                        QUAD9);

    // Test that all_complete_order_range() correctly "does nothing" in
    // parallel when passed ranges of local elements. Here we should
    // hit one of the "early return" cases for this function.
    mesh.all_complete_order_range(mesh.active_local_element_ptr_range());

    // Make sure we still have the same number of elements
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(4), mesh.n_elem());

    // Make sure we still have the same number of nodes
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(25), mesh.n_nodes());
  }

  void allCompleteOrderMixed()
  {
    ReplicatedMesh mesh(*TestCommWorld);

    // Construct a single-element Quad9 mesh
    MeshTools::Generation::build_square(mesh,
                                        /*nx=*/1, /*ny=*/1,
                                        /*xmin=*/0., /*xmax=*/1.,
                                        /*ymin=*/0., /*ymax=*/1.,
                                        /*elem_type=*/QUAD9);

    // Pointer to the single Elem
    const Elem * elem = mesh.elem_ptr(0);

    // Amount to offset the indices of newly-added nodes by
    const unsigned int node_offset = elem->n_nodes();

    // Add an extra node above each of the original vertices
    for (dof_id_type n=0; n<elem->n_vertices(); ++n)
      mesh.add_point(mesh.point(n) + Point(0,0,1), /*id=*/node_offset + n, /*proc_id=*/0);

    // Construct Edge2 elements and attach to each vertex of Quad9
    for (dof_id_type n=0; n<elem->n_vertices(); ++n)
      {
        Elem * elem = mesh.add_elem(Elem::build(EDGE2));
        elem->set_node(0) = mesh.node_ptr(n);
        elem->set_node(1) = mesh.node_ptr(node_offset + n);
      }

    // Convert all elements, automatically skipping those that are
    // already complete-order.
    mesh.all_complete_order_range(mesh.element_ptr_range());

    // Make sure we still have the expected total number of elements
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(5), mesh.n_elem());

    // Make sure we have the correct number of nodes, 9+4+4=17
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(17), mesh.n_nodes());

    // Make sure that the type of the original Elem is unchanged
    CPPUNIT_ASSERT_EQUAL(QUAD9, elem->type());

    // Make sure that the other elements are now Edge3s
    for (dof_id_type e=1; e<5; ++e)
      CPPUNIT_ASSERT_EQUAL(EDGE3, mesh.elem_ptr(e)->type());
  }


  void allCompleteOrderMixedFixing()
  {
    Mesh mesh(*TestCommWorld);

    // Construct a multi-element Tri3 mesh
    MeshTools::Generation::build_square(mesh,
                                        /*nx=*/3, /*ny=*/3,
                                        /*xmin=*/0., /*xmax=*/1.,
                                        /*ymin=*/0., /*ymax=*/1.,
                                        /*elem_type=*/TRI3);

    // Convert elements to "complete"-order, but do so a few at a time,
    // leaving a broken mesh after initial conversions to see if it
    // will be fixed by subsequent conversions.
    int flipcnt = 0;
    for (auto it = mesh.elements_begin(), end = mesh.elements_end(); it != end; ++it)
      {
        auto range_start = it, range_end = it;
        ++range_end;
        if (flipcnt++%2 && range_end != end)
          ++range_end;

        mesh.all_complete_order_range({range_start, range_end});
      }

    // Make sure we still have the expected total number of elements
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(18), mesh.n_elem());

    // Make sure we have the correct number of nodes, 7*7+3*3*2
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(67), mesh.n_nodes());

    // Make sure that the elements are now upgraded
    for (const auto & elem : mesh.element_ptr_range())
      CPPUNIT_ASSERT_EQUAL(TRI7, elem->type());
  }

  void allCompleteOrderMixedFixing3D()
  {
    Mesh mesh(*TestCommWorld);

    // Construct a multi-element Prism6 mesh
    MeshTools::Generation::build_cube(mesh,
                                      /*nx=*/3, /*ny=*/3, /*nz=*/3,
                                      /*xmin=*/0., /*xmax=*/1.,
                                      /*ymin=*/0., /*ymax=*/1.,
                                      /*zmin=*/0., /*zmax=*/1.,
                                      /*elem_type=*/PRISM6);

    // Convert elements to "complete"-order, but do so a few at a time,
    // leaving a broken mesh after initial conversions to see if it
    // will be fixed by subsequent conversions.
    int flipcnt = 0;
    for (auto it = mesh.elements_begin(), end = mesh.elements_end(); it != end; ++it)
      {
        auto range_start = it, range_end = it;
        ++range_end;
        if (flipcnt++%2 && range_end != end)
          ++range_end;

        mesh.all_complete_order_range({range_start, range_end});
      }

    // Make sure we still have the expected total number of elements
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(54), mesh.n_elem());

    // Make sure we have the correct number of nodes, 7*7*7 + 3*3*2*7
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(469), mesh.n_nodes());

    // Make sure that the elements are now upgraded
    for (const auto & elem : mesh.element_ptr_range())
      CPPUNIT_ASSERT_EQUAL(PRISM21, elem->type());
  }


};


CPPUNIT_TEST_SUITE_REGISTRATION( AllSecondOrderTest );
