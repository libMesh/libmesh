// Ignore unused parameter warnings coming from cppunit headers
#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

#include <libmesh/equation_systems.h>
#include <libmesh/replicated_mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/edge_edge2.h>
#include <libmesh/face_quad4.h>
#include <libmesh/face_tri3.h>
#include <libmesh/cell_hex8.h>
#include <libmesh/dof_map.h>
#include <libmesh/linear_implicit_system.h>
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

class MixedDimensionMeshTest : public CppUnit::TestCase {
  /**
   * The goal of this test is to ensure that a 2D mesh with 1D elements overlapping
   * on the edge is consistent. That is, they share the same global node numbers and
   * the same dof numbers for a variable.
   */
public:
  CPPUNIT_TEST_SUITE( MixedDimensionMeshTest );

  CPPUNIT_TEST( testMesh );
  CPPUNIT_TEST( testDofOrdering );
  CPPUNIT_TEST( testPointLocatorTree );

  CPPUNIT_TEST_SUITE_END();

protected:

  ReplicatedMesh* _mesh;

  void build_mesh()
  {
    _mesh = new ReplicatedMesh(*TestCommWorld);

    // (0,1)           (1,1)
    // x---------------x
    // |               |
    // |               |
    // |               |
    // |               |
    // |               |
    // x---------------x
    // (0,0)           (1,0)
    // |               |
    // |               |
    // |               |
    // |               |
    // x---------------x
    // (0,-1)          (1,-1)

    _mesh->set_mesh_dimension(2);

    _mesh->add_point( Point(0.0,-1.0), 4 );
    _mesh->add_point( Point(1.0,-1.0), 5 );
    _mesh->add_point( Point(1.0, 0.0), 1 );
    _mesh->add_point( Point(1.0, 1.0), 2 );
    _mesh->add_point( Point(0.0, 1.0), 3 );
    _mesh->add_point( Point(0.0, 0.0), 0 );

    {
      Elem* elem_top = _mesh->add_elem( new Quad4 );
      elem_top->set_node(0) = _mesh->node_ptr(0);
      elem_top->set_node(1) = _mesh->node_ptr(1);
      elem_top->set_node(2) = _mesh->node_ptr(2);
      elem_top->set_node(3) = _mesh->node_ptr(3);

      Elem* elem_bottom = _mesh->add_elem( new Quad4 );
      elem_bottom->set_node(0) = _mesh->node_ptr(4);
      elem_bottom->set_node(1) = _mesh->node_ptr(5);
      elem_bottom->set_node(2) = _mesh->node_ptr(1);
      elem_bottom->set_node(3) = _mesh->node_ptr(0);

      Elem* edge = _mesh->add_elem( new Edge2 );
      edge->set_node(0) = _mesh->node_ptr(0);
      edge->set_node(1) = _mesh->node_ptr(1);

      // 2D elements will have subdomain id 0, this one will have 1
      edge->subdomain_id() = 1;
    }

    // libMesh will renumber, but we numbered according to its scheme
    // anyway. We do this because when we call uniformly_refine subsequenly,
    // it's going use skip_renumber=false.
    _mesh->prepare_for_use(false /*skip_renumber*/);
  }

public:
  void setUp()
  {
    this->build_mesh();
  }

  void tearDown()
  {
    delete _mesh;
  }

  void testMesh()
  {
    // There'd better be 3 elements
    CPPUNIT_ASSERT_EQUAL( (dof_id_type)3, _mesh->n_elem() );

    // There'd better be only 6 nodes
    CPPUNIT_ASSERT_EQUAL( (dof_id_type)6, _mesh->n_nodes() );

    /* The nodes for the EDGE2 element should have the same global ids
       as the bottom edge of the top QUAD4 element */
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(2).node_id(0), _mesh->elem_ref(0).node_id(0) );
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(2).node_id(1), _mesh->elem_ref(0).node_id(1) );

    /* The nodes for the EDGE2 element should have the same global ids
       as the top edge of the bottom QUAD4 element */
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(2).node_id(0), _mesh->elem_ref(1).node_id(3) );
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(2).node_id(1), _mesh->elem_ref(1).node_id(2) );

    /* The nodes for the bottom edge of the top QUAD4 element should have
       the same global ids as the top edge of the bottom QUAD4 element */
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(0).node_id(0), _mesh->elem_ref(1).node_id(3) );
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(0).node_id(1), _mesh->elem_ref(1).node_id(2) );

    // We didn't set an interior_parent on the edge element, so it
    // should default to NULL
    CPPUNIT_ASSERT( _mesh->elem_ref(2).interior_parent() );
  }

  void testDofOrdering()
  {
    EquationSystems es(*_mesh);
    es.add_system<LinearImplicitSystem>("TestDofSystem");
    es.get_system("TestDofSystem").add_variable("u",FIRST);
    es.init();

    DofMap& dof_map = es.get_system("TestDofSystem").get_dof_map();

    std::vector<dof_id_type> top_quad_dof_indices, bottom_quad_dof_indices, edge_dof_indices;

    dof_map.dof_indices( _mesh->elem_ptr(0), top_quad_dof_indices );
    dof_map.dof_indices( _mesh->elem_ptr(1), bottom_quad_dof_indices );
    dof_map.dof_indices( _mesh->elem_ptr(2), edge_dof_indices );

    /* The dofs for the EDGE2 element should be the same
       as the bottom edge of the top QUAD4 dofs */
    CPPUNIT_ASSERT_EQUAL( edge_dof_indices[0], top_quad_dof_indices[0] );
    CPPUNIT_ASSERT_EQUAL( edge_dof_indices[1], top_quad_dof_indices[1] );

    /* The dofs for the EDGE2 element should be the same
       as the top edge of the bottom QUAD4 dofs */
    CPPUNIT_ASSERT_EQUAL( edge_dof_indices[0], bottom_quad_dof_indices[3] );
    CPPUNIT_ASSERT_EQUAL( edge_dof_indices[1], bottom_quad_dof_indices[2] );

    /* The nodes for the bottom edge of the top QUAD4 element should have
       the same global ids as the top edge of the bottom QUAD4 element */
    CPPUNIT_ASSERT_EQUAL( top_quad_dof_indices[0], bottom_quad_dof_indices[3] );
    CPPUNIT_ASSERT_EQUAL( top_quad_dof_indices[1], bottom_quad_dof_indices[2] );
  }

  void testPointLocatorTree()
  {
    UniquePtr<PointLocatorBase> locator = _mesh->sub_point_locator();

    Point top_point(0.5, 0.5);
    const Elem* top_elem = (*locator)(top_point);
    CPPUNIT_ASSERT(top_elem);

    // We should have gotten back the top quad
    CPPUNIT_ASSERT_EQUAL( (dof_id_type)0, top_elem->id() );

    Point bottom_point(0.5, -0.5);
    const Elem* bottom_elem = (*locator)(bottom_point);
    CPPUNIT_ASSERT(bottom_elem);

    // We should have gotten back the bottom quad
    CPPUNIT_ASSERT_EQUAL( (dof_id_type)1, bottom_elem->id() );

    // Test getting back the edge
    {
      std::set<subdomain_id_type> subdomain_id; subdomain_id.insert(1);
      Point interface_point( 0.5, 0.0 );
      const Elem* interface_elem = (*locator)(interface_point, &subdomain_id);
      CPPUNIT_ASSERT(interface_elem);

      // We should have gotten back the overlapping edge element
      CPPUNIT_ASSERT_EQUAL( (dof_id_type)2, interface_elem->id() );
    }
  }

};

class MixedDimensionRefinedMeshTest : public MixedDimensionMeshTest {
  /**
   * The goal of this test is the same as the previous, but now we do a
   * uniform refinement and make sure the result mesh is consistent. i.e.
   * the new node shared between the 1D elements is the same as the node
   * shared on the underlying quads, and so on.
   */
public:
  CPPUNIT_TEST_SUITE( MixedDimensionRefinedMeshTest );

  CPPUNIT_TEST( testMesh );
  CPPUNIT_TEST( testDofOrdering );

  CPPUNIT_TEST_SUITE_END();

  // Yes, this is necessary. Somewhere in those macros is a protected/private
public:

  void setUp()
  {
    // 3-------10------2
    // |       |       |
    // |   5   |   6   |
    // 8-------7-------9
    // |       |       |
    // |   3   |   4   |
    // 0-------6-------1
    // |       |       |
    // |   9   |  10   |
    // 13------12-------14
    // |       |       |
    // |   7   |   8   |
    // 4-------11------5
    this->build_mesh();

#ifdef LIBMESH_ENABLE_AMR
    MeshRefinement(*_mesh).uniformly_refine(1);
#endif
  }

  void testMesh()
  {
#ifdef LIBMESH_ENABLE_AMR
    // We should have 13 total and 10 active elements.
    CPPUNIT_ASSERT_EQUAL( (dof_id_type)13, _mesh->n_elem() );
    CPPUNIT_ASSERT_EQUAL( (dof_id_type)10, _mesh->n_active_elem() );

    // We should have 15 nodes
    CPPUNIT_ASSERT_EQUAL( (dof_id_type)15, _mesh->n_nodes() );

    // EDGE2,id=11 should have same nodes of bottom of QUAD4, id=3
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(11).node_id(0),
                          _mesh->elem_ref(3).node_id(0) );
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(11).node_id(1),
                          _mesh->elem_ref(3).node_id(1) );

    // EDGE2,id=12 should have same nodes of bottom of QUAD4, id=4
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(12).node_id(0),
                          _mesh->elem_ref(4).node_id(0) );
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(12).node_id(1),
                          _mesh->elem_ref(4).node_id(1) );

    // EDGE2,id=11 should have same nodes of top of QUAD4, id=9
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(11).node_id(0),
                          _mesh->elem_ref(9).node_id(3) );
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(11).node_id(1),
                          _mesh->elem_ref(9).node_id(2) );

    // EDGE2,id=12 should have same nodes of top of QUAD4, id=10
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(12).node_id(0),
                          _mesh->elem_ref(10).node_id(3) );
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(12).node_id(1),
                          _mesh->elem_ref(10).node_id(2) );

    // Shared node between the EDGE2 elements should have the same global id
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(11).node_id(1),
                          _mesh->elem_ref(12).node_id(0) );

    // EDGE2 child elements should have the correct parent
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(11).parent(),
                          _mesh->elem_ptr(2) );
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(12).parent(),
                          _mesh->elem_ptr(2) );

    // EDGE2 child elements should have the correct interior_parent
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(11).interior_parent(),
                          _mesh->elem_ptr(3) );
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(12).interior_parent(),
                          _mesh->elem_ptr(4) );
#endif
  }

  void testDofOrdering()
  {
#ifdef LIBMESH_ENABLE_AMR
    EquationSystems es(*_mesh);
    es.add_system<LinearImplicitSystem>("TestDofSystem");
    es.get_system("TestDofSystem").add_variable("u",FIRST);
    es.init();

    DofMap& dof_map = es.get_system("TestDofSystem").get_dof_map();

    std::vector<dof_id_type> top_quad3_dof_indices, top_quad4_dof_indices;
    std::vector<dof_id_type> bottom_quad9_dof_indices, bottom_quad10_dof_indices;
    std::vector<dof_id_type> edge11_dof_indices, edge12_dof_indices;

    dof_map.dof_indices( _mesh->elem_ptr(3), top_quad3_dof_indices );
    dof_map.dof_indices( _mesh->elem_ptr(4), top_quad4_dof_indices );
    dof_map.dof_indices( _mesh->elem_ptr(9), bottom_quad9_dof_indices );
    dof_map.dof_indices( _mesh->elem_ptr(10), bottom_quad10_dof_indices );
    dof_map.dof_indices( _mesh->elem_ptr(11), edge11_dof_indices );
    dof_map.dof_indices( _mesh->elem_ptr(12), edge12_dof_indices );

    // EDGE2,id=11 should have same dofs as of bottom of QUAD4, id=3
    CPPUNIT_ASSERT_EQUAL( edge11_dof_indices[0], top_quad3_dof_indices[0] );
    CPPUNIT_ASSERT_EQUAL( edge11_dof_indices[1], top_quad3_dof_indices[1] );

    // EDGE2,id=12 should have same dofs of bottom of QUAD4, id=4
    CPPUNIT_ASSERT_EQUAL( edge12_dof_indices[0], top_quad4_dof_indices[0] );
    CPPUNIT_ASSERT_EQUAL( edge12_dof_indices[1], top_quad4_dof_indices[1] );

    // EDGE2,id=11 should have same dofs of top of QUAD4, id=9
    CPPUNIT_ASSERT_EQUAL( edge11_dof_indices[0], bottom_quad9_dof_indices[3] );
    CPPUNIT_ASSERT_EQUAL( edge11_dof_indices[1], bottom_quad9_dof_indices[2] );

    // EDGE2,id=12 should have same dofs of top of QUAD4, id=10
    CPPUNIT_ASSERT_EQUAL( edge12_dof_indices[0], bottom_quad10_dof_indices[3] );
    CPPUNIT_ASSERT_EQUAL( edge12_dof_indices[1], bottom_quad10_dof_indices[2] );

    //EDGE2 elements should have same shared dof number
    CPPUNIT_ASSERT_EQUAL( edge11_dof_indices[1], edge12_dof_indices[0] );
#endif
  }

};

class MixedDimensionNonUniformRefinement : public CppUnit::TestCase {
  /**
   * Given a mesh with four 2D elements and an overlapping 1D element this test ensures that when
   * a single 2D element that is flagged for refinement, which is a neighbor of the 1D element, that
   * the 1D element will is also be flagged for refinement due to an underrefined_boundary_limit of
   * 0 (default) and the neighboring 2D element will also refined due to an overrefined_boundary_limit
   * of 0 (default).
   */
public:
  CPPUNIT_TEST_SUITE( MixedDimensionNonUniformRefinement );

  CPPUNIT_TEST( testMesh );
  CPPUNIT_TEST( testDofOrdering );

  CPPUNIT_TEST_SUITE_END();

  // Yes, this is necessary. Somewhere in those macros is a protected/private
protected:

  ReplicatedMesh* _mesh;

  void build_mesh()
  {
    _mesh = new ReplicatedMesh(*TestCommWorld);
    // We start with this
    //
    //
    // (0,2)           (1,2)
    // 4---------------5
    // |               |
    // |               |
    // |       1       |
    // |               |
    // |               |
    // (0,1)           (1,1)
    // 3---------------2
    // |               |
    // |               |
    // |       0       |
    // |               |
    // |               |
    // 0---------------1
    // (0,0)           (1,0)
    // |               |
    // |       2       |
    // |               |
    // |               |
    // 6---------------7
    // (0,-1)          (1,-1)
    // |               |
    // |       3       |
    // |               |
    // |               |
    // 9---------------8
    // (0,-2)          (1,-2)
    //
    // But the single element refinement should result
    // with this for the default max_mismatch = 0 case
    //
    // 4---------------5
    // |               |
    // |               |
    // |       1       |
    // |               |
    // |               |
    // 3------14-------2
    // |       |       |
    // |   7   |   8   |
    // 12------11-------13
    // |       |       |
    // |   5   |   6   |
    // 0------10-------1
    // |       |       |
    // |   11  |   12  |
    // 17------16-------18
    // |       |       |
    // |   9   |   10  |
    // 6------15-------7
    // |               |
    // |               |
    // |       3       |
    // |               |
    // |               |
    // 9---------------8

    _mesh->set_mesh_dimension(2);

    _mesh->add_point( Point(0.0,0.0), 0 );
    _mesh->add_point( Point(1.0,0.0), 1 );
    _mesh->add_point( Point(1.0,1.0), 2 );
    _mesh->add_point( Point(0.0,1.0), 3 );
    _mesh->add_point( Point(0.0,2.0), 4 );
    _mesh->add_point( Point(1.0,2.0), 5 );
    _mesh->add_point( Point(0.0,-1.0), 6 );
    _mesh->add_point( Point(1.0,-1.0), 7 );
    _mesh->add_point( Point(1.0,-2.0), 8 );
    _mesh->add_point( Point(0.0,-2.0), 9 );


    {
      Elem* quad0 = _mesh->add_elem( new Quad4 );
      quad0->set_node(0) = _mesh->node_ptr(0);
      quad0->set_node(1) = _mesh->node_ptr(1);
      quad0->set_node(2) = _mesh->node_ptr(2);
      quad0->set_node(3) = _mesh->node_ptr(3);

      Elem* quad1 = _mesh->add_elem( new Quad4 );
      quad1->set_node(0) = _mesh->node_ptr(3);
      quad1->set_node(1) = _mesh->node_ptr(2);
      quad1->set_node(2) = _mesh->node_ptr(5);
      quad1->set_node(3) = _mesh->node_ptr(4);

      Elem* quad2 = _mesh->add_elem( new Quad4 );
      quad2->set_node(0) = _mesh->node_ptr(6);
      quad2->set_node(1) = _mesh->node_ptr(7);
      quad2->set_node(2) = _mesh->node_ptr(1);
      quad2->set_node(3) = _mesh->node_ptr(0);

      Elem* quad3 = _mesh->add_elem( new Quad4 );
      quad3->set_node(0) = _mesh->node_ptr(9);
      quad3->set_node(1) = _mesh->node_ptr(8);
      quad3->set_node(2) = _mesh->node_ptr(7);
      quad3->set_node(3) = _mesh->node_ptr(6);

      Elem* edge = _mesh->add_elem( new Edge2 );
      edge->set_node(0) = _mesh->node_ptr(0);
      edge->set_node(1) = _mesh->node_ptr(1);

      // 2D elements will have subdomain id 0, this one will have 1
      edge->subdomain_id() = 1;
    }

    // libMesh will renumber, but we numbered according to its scheme
    // anyway. We do this because when we call uniformly_refine subsequenly,
    // it's going use skip_renumber=false.
    _mesh->prepare_for_use(false /*skip_renumber*/);


#ifdef LIBMESH_ENABLE_AMR
    //Flag the bottom element for refinement
    _mesh->elem_ref(0).set_refinement_flag(Elem::REFINE);
    MeshRefinement(*_mesh).refine_and_coarsen_elements();
#endif

  }

public:
  void setUp()
  {
    this->build_mesh();
  }

  void tearDown()
  {
    delete _mesh;
  }

  void testMesh()
  {
#ifdef LIBMESH_ENABLE_AMR
    // We should have 13 total and 10 active elements.
    CPPUNIT_ASSERT_EQUAL( (dof_id_type)15, _mesh->n_elem() );
    CPPUNIT_ASSERT_EQUAL( (dof_id_type)12, _mesh->n_active_elem() );

    // We should have 15 nodes
    CPPUNIT_ASSERT_EQUAL( (dof_id_type)19, _mesh->n_nodes() );

    // EDGE2,id=13 should have same nodes of bottom of QUAD4, id=5
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(13).node_id(0),
                          _mesh->elem_ref(5).node_id(0) );
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(13).node_id(1),
                          _mesh->elem_ref(5).node_id(1) );

    // EDGE2,id=14 should have same nodes of bottom of QUAD4, id=6
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(14).node_id(0),
                          _mesh->elem_ref(6).node_id(0) );
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(14).node_id(1),
                          _mesh->elem_ref(6).node_id(1) );

    // EDGE2,id=13 should have same nodes of top of QUAD4, id=11
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(13).node_id(0),
                          _mesh->elem_ref(11).node_id(3) );
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(13).node_id(1),
                          _mesh->elem_ref(11).node_id(2) );

    // EDGE2,id=14 should have same nodes of top of QUAD4, id=12
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(14).node_id(0),
                          _mesh->elem_ref(12).node_id(3) );
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(14).node_id(1),
                          _mesh->elem_ref(12).node_id(2) );

    // Shared node between the EDGE2 elements should have the same global id
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(13).node_id(1),
                          _mesh->elem_ref(14).node_id(0) );

    // EDGE2 child elements should have the correct parent
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(13).parent(),
                          _mesh->elem_ptr(4) );
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(14).parent(),
                          _mesh->elem_ptr(4) );

    // EDGE2 child elements should have the correct interior_parent
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(13).interior_parent(),
                          _mesh->elem_ptr(5) );
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(14).interior_parent(),
                          _mesh->elem_ptr(6) );
#endif
  }

  void testDofOrdering()
  {
#ifdef LIBMESH_ENABLE_AMR
    EquationSystems es(*_mesh);
    es.add_system<LinearImplicitSystem>("TestDofSystem");
    es.get_system("TestDofSystem").add_variable("u",FIRST);
    es.init();

    DofMap& dof_map = es.get_system("TestDofSystem").get_dof_map();

    std::vector<dof_id_type> top_quad5_dof_indices, top_quad6_dof_indices;
    std::vector<dof_id_type> bottom_quad11_dof_indices, bottom_quad12_dof_indices;
    std::vector<dof_id_type> edge13_dof_indices, edge14_dof_indices;

    dof_map.dof_indices( _mesh->elem_ptr(5),  top_quad5_dof_indices );
    dof_map.dof_indices( _mesh->elem_ptr(6),  top_quad6_dof_indices );
    dof_map.dof_indices( _mesh->elem_ptr(11), bottom_quad11_dof_indices );
    dof_map.dof_indices( _mesh->elem_ptr(12), bottom_quad12_dof_indices );
    dof_map.dof_indices( _mesh->elem_ptr(13), edge13_dof_indices );
    dof_map.dof_indices( _mesh->elem_ptr(14), edge14_dof_indices );

    // EDGE2,id=13 should have same dofs as of bottom of QUAD4, id=5
    CPPUNIT_ASSERT_EQUAL( edge13_dof_indices[0], top_quad5_dof_indices[0] );
    CPPUNIT_ASSERT_EQUAL( edge13_dof_indices[1], top_quad5_dof_indices[1] );

    // EDGE2,id=14 should have same dofs of bottom of QUAD4, id=6
    CPPUNIT_ASSERT_EQUAL( edge14_dof_indices[0], top_quad6_dof_indices[0] );
    CPPUNIT_ASSERT_EQUAL( edge14_dof_indices[1], top_quad6_dof_indices[1] );

    // EDGE2,id=13 should have same dofs of top of QUAD4, id=11
    CPPUNIT_ASSERT_EQUAL( edge13_dof_indices[0], bottom_quad11_dof_indices[3] );
    CPPUNIT_ASSERT_EQUAL( edge13_dof_indices[1], bottom_quad11_dof_indices[2] );

    // EDGE2,id=14 should have same dofs of top of QUAD4, id=12
    CPPUNIT_ASSERT_EQUAL( edge14_dof_indices[0], bottom_quad12_dof_indices[3] );
    CPPUNIT_ASSERT_EQUAL( edge14_dof_indices[1], bottom_quad12_dof_indices[2] );

    //EDGE2 elements should have same shared dof number
    CPPUNIT_ASSERT_EQUAL( edge13_dof_indices[1], edge14_dof_indices[0] );
#endif
  }
};

class MixedDimensionNonUniformRefinementTriangle : public CppUnit::TestCase {
  /**
   * Given a mesh with four TRI3 elements and an overlapping EDG2 element this test ensures that when
   * a single TRI3 element that is flagged for refinement, which is a neighbor of the EDGE2 element, that
   * the EDGE2 element will is also be flagged for refinement due to an underrefined_boundary_limit of
   * 0 (default) and the neighboring TRI3 element will also refined due to an overrefined_boundary_limit
   * of 0 (default).
   */
public:
  CPPUNIT_TEST_SUITE( MixedDimensionNonUniformRefinementTriangle );

  CPPUNIT_TEST( testMesh );
  CPPUNIT_TEST( testDofOrdering );

  CPPUNIT_TEST_SUITE_END();

protected:

  ReplicatedMesh* _mesh;

  void build_mesh()
  {
    _mesh = new ReplicatedMesh(*TestCommWorld);

    /**
     * We start with this
     *
     *  (0,1)             (1,1)
     *        3---------2
     *        |        -|
     *        |      -  |
     *        |    -    |
     *        |  -      |
     *        |-        |
     *  (0,0) 0---------1 (1,0)
     *        |        -|
     *        |      -  |
     *        |    -    |
     *        |  -      |
     *        |-        |
     *        4---------5
     *  (0,-1)            (1,-1)
     *
     *
     * But the single element refinement should result
     * with this for the default max_mismatch = 0 case
     *
     * (0,1)             (1,1)
     *       3---------2
     *       |        -|
     *       |      -  |
     *       |    7----8
     *       |  - |  - |
     *       |-   |-   |
     * (0,0) 0----6----1 (1,0)
     *       |  - |   -|
     *       |-   | -  |
     *      10----9    |
     *       |  -      |
     *       |-        |
     *       4---------5
     * (0,-1)           (1,-1)
     */

    _mesh->set_mesh_dimension(2);

    _mesh->add_point( Point(0.0,-1.0), 4 );
    _mesh->add_point( Point(1.0,-1.0), 5 );
    _mesh->add_point( Point(1.0, 0.0), 1 );
    _mesh->add_point( Point(1.0, 1.0), 2 );
    _mesh->add_point( Point(0.0, 1.0), 3 );
    _mesh->add_point( Point(0.0, 0.0), 0 );

    {
      Elem* elem0 = _mesh->add_elem( new Tri3 );
      elem0->set_node(0) = _mesh->node_ptr(0);
      elem0->set_node(1) = _mesh->node_ptr(1);
      elem0->set_node(2) = _mesh->node_ptr(2);

      Elem* elem1 = _mesh->add_elem( new Tri3 );
      elem1->set_node(0) = _mesh->node_ptr(2);
      elem1->set_node(1) = _mesh->node_ptr(3);
      elem1->set_node(2) = _mesh->node_ptr(0);

      Elem* elem2 = _mesh->add_elem( new Tri3 );
      elem2->set_node(0) = _mesh->node_ptr(1);
      elem2->set_node(1) = _mesh->node_ptr(0);
      elem2->set_node(2) = _mesh->node_ptr(4);

      Elem* elem3 = _mesh->add_elem( new Tri3 );
      elem3->set_node(0) = _mesh->node_ptr(4);
      elem3->set_node(1) = _mesh->node_ptr(5);
      elem3->set_node(2) = _mesh->node_ptr(1);

      Elem* edge = _mesh->add_elem( new Edge2 );
      edge->set_node(0) = _mesh->node_ptr(0);
      edge->set_node(1) = _mesh->node_ptr(1);

      // 2D elements will have subdomain id 0, this one will have 1
      edge->subdomain_id() = 1;

    }

    // libMesh will renumber, but we numbered according to its scheme
    // anyway. We do this because when we call uniformly_refine subsequently,
    // it's going use skip_renumber=false.
    _mesh->prepare_for_use(false /*skip_renumber*/);

#ifdef LIBMESH_ENABLE_AMR
    //Flag the bottom element for refinement
    _mesh->elem_ref(4).set_refinement_flag(Elem::REFINE);
    MeshRefinement(*_mesh).refine_and_coarsen_elements();
#endif
  }

public:
  void setUp()
  {
    this->build_mesh();
  }

  void tearDown()
  {
    delete _mesh;
  }

  void testMesh()
  {
#ifdef LIBMESH_ENABLE_AMR
    // We should have 15 total and 12 active elements.
    CPPUNIT_ASSERT_EQUAL( (dof_id_type)15, _mesh->n_elem() );
    CPPUNIT_ASSERT_EQUAL( (dof_id_type)12, _mesh->n_active_elem() );

    // We should have 15 nodes
    CPPUNIT_ASSERT_EQUAL( (dof_id_type)11, _mesh->n_nodes() );

    // EDGE2,id=13 should have same nodes of the base of TRI3, id=5
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(13).node_id(0),
                          _mesh->elem_ref(5).node_id(0) );
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(13).node_id(1),
                          _mesh->elem_ref(5).node_id(1) );

    // EDGE2,id=13 should have same nodes of the base of TRI3, id=10
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(13).node_id(0),
                          _mesh->elem_ref(10).node_id(1) );
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(13).node_id(1),
                          _mesh->elem_ref(10).node_id(0) );

    // EDGE2,id=13 should have same node as the tip of TRI3, id=8 and id=12
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(13).node_id(1),
                          _mesh->elem_ref(8).node_id(0) );
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(13).node_id(1),
                          _mesh->elem_ref(12).node_id(0) );

    // EDGE2,id=14 should have same nodes of the base of TRI3, id=6
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(14).node_id(0),
                          _mesh->elem_ref(6).node_id(0) );
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(14).node_id(1),
                          _mesh->elem_ref(6).node_id(1) );

    // EDGE2,id=14 should have same nodes of the base of TRI3, id=9
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(14).node_id(0),
                          _mesh->elem_ref(9).node_id(1) );
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(14).node_id(1),
                          _mesh->elem_ref(9).node_id(0) );

    // EDGE2,id=14 should have same node as the tip of TRI3, id=8 and id=12
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(14).node_id(0),
                          _mesh->elem_ref(8).node_id(0) );
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(14).node_id(0),
                          _mesh->elem_ref(12).node_id(0) );

    // Shared node between the EDGE2 elements should have the same global id
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(13).node_id(1),
                          _mesh->elem_ref(14).node_id(0) );

    // EDGE2 child elements should have the correct parent
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(13).parent(),
                          _mesh->elem_ptr(4) );
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(14).parent(),
                          _mesh->elem_ptr(4) );

    // EDGE2 child elements should have the correct interior_parent
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(13).interior_parent(),
                          _mesh->elem_ptr(5) );
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(14).interior_parent(),
                          _mesh->elem_ptr(6) );

#endif // LIBMESH_ENABLE_AMR
  }

  void testDofOrdering()
  {
#ifdef LIBMESH_ENABLE_AMR
    EquationSystems es(*_mesh);
    es.add_system<LinearImplicitSystem>("TestDofSystem");
    es.get_system("TestDofSystem").add_variable("u",FIRST);
    es.init();

    DofMap& dof_map = es.get_system("TestDofSystem").get_dof_map();

    //Elements above the EDGE2 elements
    std::vector<dof_id_type> elem5_dof_indices, elem6_dof_indices, elem8_dof_indices;

    //Elements below the EDGE2 elements
    std::vector<dof_id_type> elem9_dof_indices, elem10_dof_indices, elem12_dof_indices;

    //EDGE2 Elements
    std::vector<dof_id_type> elem13_dof_indices, elem14_dof_indices;

    dof_map.dof_indices( _mesh->elem_ptr(5), elem5_dof_indices );
    dof_map.dof_indices( _mesh->elem_ptr(6), elem6_dof_indices );
    dof_map.dof_indices( _mesh->elem_ptr(8), elem8_dof_indices );
    dof_map.dof_indices( _mesh->elem_ptr(9), elem9_dof_indices );
    dof_map.dof_indices( _mesh->elem_ptr(10), elem10_dof_indices );
    dof_map.dof_indices( _mesh->elem_ptr(12), elem12_dof_indices );
    dof_map.dof_indices( _mesh->elem_ptr(13), elem13_dof_indices );
    dof_map.dof_indices( _mesh->elem_ptr(14), elem14_dof_indices );

    /* The dofs for the EDGE2 (id = 13 and id =14) element should be the same
       as the bottom edge of the top TRI3 (id=5 and id=6) and the tip of
       TRI3 id = 8 dofs */
    CPPUNIT_ASSERT_EQUAL( elem13_dof_indices[0], elem5_dof_indices[0] );
    CPPUNIT_ASSERT_EQUAL( elem13_dof_indices[1], elem5_dof_indices[1] );
    CPPUNIT_ASSERT_EQUAL( elem13_dof_indices[1], elem6_dof_indices[0] );
    CPPUNIT_ASSERT_EQUAL( elem13_dof_indices[1], elem8_dof_indices[0] );
    CPPUNIT_ASSERT_EQUAL( elem14_dof_indices[0], elem6_dof_indices[0] );
    CPPUNIT_ASSERT_EQUAL( elem14_dof_indices[1], elem6_dof_indices[1] );
    CPPUNIT_ASSERT_EQUAL( elem14_dof_indices[0], elem5_dof_indices[1] );
    CPPUNIT_ASSERT_EQUAL( elem14_dof_indices[0], elem8_dof_indices[0] );

    /* The dofs for the EDGE2 (id = 13 and id =14) element should be the same
       as the top edge of the bottom TRI3 (id=9 and id=10) and the tip of
       TRI3 id = 12 dofs */
    CPPUNIT_ASSERT_EQUAL( elem13_dof_indices[0], elem10_dof_indices[1] );
    CPPUNIT_ASSERT_EQUAL( elem13_dof_indices[1], elem10_dof_indices[0] );
    CPPUNIT_ASSERT_EQUAL( elem13_dof_indices[1], elem9_dof_indices[1] );
    CPPUNIT_ASSERT_EQUAL( elem13_dof_indices[1], elem12_dof_indices[0] );
    CPPUNIT_ASSERT_EQUAL( elem14_dof_indices[0], elem9_dof_indices[1] );
    CPPUNIT_ASSERT_EQUAL( elem14_dof_indices[1], elem9_dof_indices[0] );
    CPPUNIT_ASSERT_EQUAL( elem14_dof_indices[0], elem10_dof_indices[0] );
    CPPUNIT_ASSERT_EQUAL( elem14_dof_indices[0], elem12_dof_indices[0] );

    /* The nodes for the bottom edge of the top TRI3 elements should have
       the same global ids as the top edge of the bottom TRI3 elements. */
    CPPUNIT_ASSERT_EQUAL( elem5_dof_indices[0], elem10_dof_indices[1] );
    CPPUNIT_ASSERT_EQUAL( elem5_dof_indices[1], elem10_dof_indices[0] );
    CPPUNIT_ASSERT_EQUAL( elem6_dof_indices[0], elem9_dof_indices[1] );
    CPPUNIT_ASSERT_EQUAL( elem6_dof_indices[1], elem9_dof_indices[0] );
    CPPUNIT_ASSERT_EQUAL( elem8_dof_indices[0], elem12_dof_indices[0] );
#endif // LIBMESH_ENABLE_AMR
  }

};

class MixedDimensionNonUniformRefinement3D : public CppUnit::TestCase {
  /**
   * Given a mesh with 36 HEX8 elements and an overlapping QUAD4 element located in the center of the HEX8
   * elements, this test ensures that when a single HEX8 element is flagged for refinement, which
   * is a neighbor of the QUAD4 element, that the QUAD4 element will is also be flagged for refinement
   * due to an underrefined_boundary_limit of 0 (default) and the neighboring HEX8 element (with respect
   * to the QUAD4 element) will also refined due to an overrefined_boundary_limit of 0 (default).
   */
public:
  CPPUNIT_TEST_SUITE( MixedDimensionNonUniformRefinement3D );

  CPPUNIT_TEST( testMesh );
  CPPUNIT_TEST( testDofOrdering );

  CPPUNIT_TEST_SUITE_END();

  // Yes, this is necessary. Somewhere in those macros is a protected/private
protected:

  ReplicatedMesh* _mesh;

  void build_mesh()
  {
    _mesh = new ReplicatedMesh(*TestCommWorld);

    _mesh->set_mesh_dimension(3);

    //Add the nodes
    for (unsigned int z = 0; z < 5; z++)
      {
        for (unsigned int y = 0; y < 4; y++)
          {
            for (unsigned int x = 0; x < 4; x++)
              {
                _mesh->add_point( Point(Real(x),Real(y),Real(z)), 16*z+4*y+x);
              }
          }
      }

    {
      //Add the HEX8 elements
      for (unsigned int z = 0; z < 4; z++)
        {
          for (unsigned int y = 0; y < 3; y++)
            {
              for (unsigned int x = 0; x < 3; x++)
                {
                  Elem* hex = _mesh->add_elem( new Hex8 );
                  hex->set_node(0) = _mesh->node_ptr(x+4*y    +16*z        );
                  hex->set_node(1) = _mesh->node_ptr(x+4*y    +16*z     + 1);
                  hex->set_node(2) = _mesh->node_ptr(x+4*(y+1)+16*z     + 1);
                  hex->set_node(3) = _mesh->node_ptr(x+4*(y+1)+16*z        );
                  hex->set_node(4) = _mesh->node_ptr(x+4*y    +16*(z+1)    );
                  hex->set_node(5) = _mesh->node_ptr(x+4*y    +16*(z+1) + 1);
                  hex->set_node(6) = _mesh->node_ptr(x+4*(y+1)+16*(z+1) + 1);
                  hex->set_node(7) = _mesh->node_ptr(x+4*(y+1)+16*(z+1)    );
                }
            }
        }
      Elem* quad = _mesh->add_elem( new Quad4 );
      unsigned int x=1,y=1,z=2;
      quad->set_node(0) = _mesh->node_ptr(x+4*y    +16*z    );
      quad->set_node(1) = _mesh->node_ptr(x+4*y    +16*z + 1);
      quad->set_node(2) = _mesh->node_ptr(x+4*(y+1)+16*z + 1);
      quad->set_node(3) = _mesh->node_ptr(x+4*(y+1)+16*z    );

      // 2D elements will have subdomain id 0, this one will have 1
      quad->subdomain_id() = 1;
    }

    // libMesh will renumber, but we numbered according to its scheme
    // anyway. We do this because when we call uniformly_refine subsequently,
    // it's going use skip_renumber=false.
    _mesh->prepare_for_use(false /*skip_renumber*/);

#ifdef LIBMESH_ENABLE_AMR
    //Flag the bottom element for refinement
    _mesh->elem_ref(13).set_refinement_flag(Elem::REFINE);
    MeshRefinement(*_mesh).refine_and_coarsen_elements();
#endif
  }

public:
  void setUp()
  {
    this->build_mesh();
  }

  void tearDown()
  {
    delete _mesh;
  }

  void testMesh()
  {
#ifdef LIBMESH_ENABLE_AMR
    // We should have 57 total and 54 active elements.
    CPPUNIT_ASSERT_EQUAL( (dof_id_type)57, _mesh->n_elem() );
    CPPUNIT_ASSERT_EQUAL( (dof_id_type)54, _mesh->n_active_elem() );

    // We should have 113 nodes
    CPPUNIT_ASSERT_EQUAL( (dof_id_type)113, _mesh->n_nodes() );

    // QUAD4,id=53 should have same nodes as a face in HEX8, id=39
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(53).node_id(0),
                          _mesh->elem_ref(41).node_id(4) );
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(53).node_id(1),
                          _mesh->elem_ref(41).node_id(5) );
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(53).node_id(2),
                          _mesh->elem_ref(41).node_id(6) );
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(53).node_id(3),
                          _mesh->elem_ref(41).node_id(7) );

    // QUAD4,id=53 should have same nodes as a face in HEX8, id=45
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(53).node_id(0),
                          _mesh->elem_ref(45).node_id(0) );
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(53).node_id(1),
                          _mesh->elem_ref(45).node_id(1) );
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(53).node_id(2),
                          _mesh->elem_ref(45).node_id(2) );
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(53).node_id(3),
                          _mesh->elem_ref(45).node_id(3) );

    // QUAD4,id=54 should have same nodes as a face in HEX8, id=42
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(54).node_id(0),
                          _mesh->elem_ref(42).node_id(4) );
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(54).node_id(1),
                          _mesh->elem_ref(42).node_id(5) );
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(54).node_id(2),
                          _mesh->elem_ref(42).node_id(6) );
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(54).node_id(3),
                          _mesh->elem_ref(42).node_id(7) );

    // QUAD4,id=54 should have same nodes as a face in HEX8, id=46
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(54).node_id(0),
                          _mesh->elem_ref(46).node_id(0) );
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(54).node_id(1),
                          _mesh->elem_ref(46).node_id(1) );
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(54).node_id(2),
                          _mesh->elem_ref(46).node_id(2) );
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(54).node_id(3),
                          _mesh->elem_ref(46).node_id(3) );

    // QUAD4,id=55 should have same nodes as a face in HEX8, id=43
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(55).node_id(0),
                          _mesh->elem_ref(43).node_id(4) );
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(55).node_id(1),
                          _mesh->elem_ref(43).node_id(5) );
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(55).node_id(2),
                          _mesh->elem_ref(43).node_id(6) );
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(55).node_id(3),
                          _mesh->elem_ref(43).node_id(7) );

    // QUAD4,id=55 should have same nodes as a face in HEX8, id=47
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(55).node_id(0),
                          _mesh->elem_ref(47).node_id(0) );
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(55).node_id(1),
                          _mesh->elem_ref(47).node_id(1) );
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(55).node_id(2),
                          _mesh->elem_ref(47).node_id(2) );
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(55).node_id(3),
                          _mesh->elem_ref(47).node_id(3) );

    // QUAD4,id=56 should have same nodes as a face in HEX8, id=44
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(56).node_id(0),
                          _mesh->elem_ref(44).node_id(4) );
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(56).node_id(1),
                          _mesh->elem_ref(44).node_id(5) );
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(56).node_id(2),
                          _mesh->elem_ref(44).node_id(6) );
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(56).node_id(3),
                          _mesh->elem_ref(44).node_id(7) );

    // QUAD4,id=56 should have same nodes as a face in HEX8, id=48
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(56).node_id(0),
                          _mesh->elem_ref(48).node_id(0) );
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(56).node_id(1),
                          _mesh->elem_ref(48).node_id(1) );
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(56).node_id(2),
                          _mesh->elem_ref(48).node_id(2) );
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(56).node_id(3),
                          _mesh->elem_ref(48).node_id(3) );

    // Shared node between the QUAD4 elements should have the same global id
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(53).node_id(1),
                          _mesh->elem_ref(54).node_id(0) );
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(53).node_id(2),
                          _mesh->elem_ref(54).node_id(3) );
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(53).node_id(3),
                          _mesh->elem_ref(55).node_id(0) );
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(53).node_id(2),
                          _mesh->elem_ref(55).node_id(1) );
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(54).node_id(3),
                          _mesh->elem_ref(56).node_id(0) );
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(54).node_id(2),
                          _mesh->elem_ref(56).node_id(1) );
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(55).node_id(1),
                          _mesh->elem_ref(56).node_id(0) );
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(55).node_id(2),
                          _mesh->elem_ref(56).node_id(3) );

    // QUAD4 child elements should have the correct parent
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(53).parent(),
                          _mesh->elem_ptr(36) );
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(54).parent(),
                          _mesh->elem_ptr(36) );
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(55).parent(),
                          _mesh->elem_ptr(36) );
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(56).parent(),
                          _mesh->elem_ptr(36) );

    // QUAD4 child elements should have the correct interior_parent
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(53).interior_parent(),
                          _mesh->elem_ptr(41) );
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(54).interior_parent(),
                          _mesh->elem_ptr(42) );
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(55).interior_parent(),
                          _mesh->elem_ptr(43) );
    CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(56).interior_parent(),
                          _mesh->elem_ptr(44) );

#endif
  }

  void testDofOrdering()
  {
#ifdef LIBMESH_ENABLE_AMR
    EquationSystems es(*_mesh);
    es.add_system<LinearImplicitSystem>("TestDofSystem");
    es.get_system("TestDofSystem").add_variable("u",FIRST);
    es.init();

    DofMap& dof_map = es.get_system("TestDofSystem").get_dof_map();

    //Elements to the left of the QUAD4 elements
    std::vector<dof_id_type> elem41_dof_indices, elem42_dof_indices, elem43_dof_indices, elem44_dof_indices;
    //Elements to the right of the QUAD4 elements
    std::vector<dof_id_type> elem45_dof_indices, elem46_dof_indices, elem47_dof_indices, elem48_dof_indices;
    //QUAD4 elements
    std::vector<dof_id_type> elem53_dof_indices, elem54_dof_indices, elem55_dof_indices, elem56_dof_indices;

    dof_map.dof_indices( _mesh->elem_ptr(41), elem41_dof_indices );
    dof_map.dof_indices( _mesh->elem_ptr(42), elem42_dof_indices );
    dof_map.dof_indices( _mesh->elem_ptr(43), elem43_dof_indices );
    dof_map.dof_indices( _mesh->elem_ptr(44), elem44_dof_indices );
    dof_map.dof_indices( _mesh->elem_ptr(45), elem45_dof_indices );
    dof_map.dof_indices( _mesh->elem_ptr(46), elem46_dof_indices );
    dof_map.dof_indices( _mesh->elem_ptr(47), elem47_dof_indices );
    dof_map.dof_indices( _mesh->elem_ptr(48), elem48_dof_indices );
    dof_map.dof_indices( _mesh->elem_ptr(53), elem53_dof_indices );
    dof_map.dof_indices( _mesh->elem_ptr(54), elem54_dof_indices );
    dof_map.dof_indices( _mesh->elem_ptr(55), elem55_dof_indices );
    dof_map.dof_indices( _mesh->elem_ptr(56), elem56_dof_indices );

    /* The dofs for the QUAD4 (ids = 53, 54, 55, and 56) element should be the same
       as the face of the HEX8 elements HEX8 (id=41, 42, 43, and 44) left of the
       QUAD4 elements. */
    CPPUNIT_ASSERT_EQUAL( elem53_dof_indices[0], elem41_dof_indices[4] );
    CPPUNIT_ASSERT_EQUAL( elem53_dof_indices[1], elem41_dof_indices[5] );
    CPPUNIT_ASSERT_EQUAL( elem53_dof_indices[2], elem41_dof_indices[6] );
    CPPUNIT_ASSERT_EQUAL( elem53_dof_indices[3], elem41_dof_indices[7] );

    CPPUNIT_ASSERT_EQUAL( elem54_dof_indices[0], elem42_dof_indices[4] );
    CPPUNIT_ASSERT_EQUAL( elem54_dof_indices[1], elem42_dof_indices[5] );
    CPPUNIT_ASSERT_EQUAL( elem54_dof_indices[2], elem42_dof_indices[6] );
    CPPUNIT_ASSERT_EQUAL( elem54_dof_indices[3], elem42_dof_indices[7] );

    CPPUNIT_ASSERT_EQUAL( elem55_dof_indices[0], elem43_dof_indices[4] );
    CPPUNIT_ASSERT_EQUAL( elem55_dof_indices[1], elem43_dof_indices[5] );
    CPPUNIT_ASSERT_EQUAL( elem55_dof_indices[2], elem43_dof_indices[6] );
    CPPUNIT_ASSERT_EQUAL( elem55_dof_indices[3], elem43_dof_indices[7] );

    CPPUNIT_ASSERT_EQUAL( elem56_dof_indices[0], elem44_dof_indices[4] );
    CPPUNIT_ASSERT_EQUAL( elem56_dof_indices[1], elem44_dof_indices[5] );
    CPPUNIT_ASSERT_EQUAL( elem56_dof_indices[2], elem44_dof_indices[6] );
    CPPUNIT_ASSERT_EQUAL( elem56_dof_indices[3], elem44_dof_indices[7] );

    /* The dofs for the QUAD4 (ids = 53, 54, 55, and 56) element should be the same
       as the face of the HEX8 elements HEX8 (id=45, 46, 47, and 49) left of the
       QUAD4 elements. */
    CPPUNIT_ASSERT_EQUAL( elem53_dof_indices[0], elem45_dof_indices[0] );
    CPPUNIT_ASSERT_EQUAL( elem53_dof_indices[1], elem45_dof_indices[1] );
    CPPUNIT_ASSERT_EQUAL( elem53_dof_indices[2], elem45_dof_indices[2] );
    CPPUNIT_ASSERT_EQUAL( elem53_dof_indices[3], elem45_dof_indices[3] );

    CPPUNIT_ASSERT_EQUAL( elem54_dof_indices[0], elem46_dof_indices[0] );
    CPPUNIT_ASSERT_EQUAL( elem54_dof_indices[1], elem46_dof_indices[1] );
    CPPUNIT_ASSERT_EQUAL( elem54_dof_indices[2], elem46_dof_indices[2] );
    CPPUNIT_ASSERT_EQUAL( elem54_dof_indices[3], elem46_dof_indices[3] );

    CPPUNIT_ASSERT_EQUAL( elem55_dof_indices[0], elem47_dof_indices[0] );
    CPPUNIT_ASSERT_EQUAL( elem55_dof_indices[1], elem47_dof_indices[1] );
    CPPUNIT_ASSERT_EQUAL( elem55_dof_indices[2], elem47_dof_indices[2] );
    CPPUNIT_ASSERT_EQUAL( elem55_dof_indices[3], elem47_dof_indices[3] );

    CPPUNIT_ASSERT_EQUAL( elem56_dof_indices[0], elem48_dof_indices[0] );
    CPPUNIT_ASSERT_EQUAL( elem56_dof_indices[1], elem48_dof_indices[1] );
    CPPUNIT_ASSERT_EQUAL( elem56_dof_indices[2], elem48_dof_indices[2] );
    CPPUNIT_ASSERT_EQUAL( elem56_dof_indices[3], elem48_dof_indices[3] );

    /* The dofs for the HEX8 elements (id=41, 42, 43, and 44) should be the same
       on the left side of the QUAD4 elements as the HEX8 elements (id=45, 46, 47, and 48)
       on the right as QUAD4 elements. */
    CPPUNIT_ASSERT_EQUAL( elem41_dof_indices[4], elem45_dof_indices[0] );
    CPPUNIT_ASSERT_EQUAL( elem41_dof_indices[5], elem45_dof_indices[1] );
    CPPUNIT_ASSERT_EQUAL( elem41_dof_indices[6], elem45_dof_indices[2] );
    CPPUNIT_ASSERT_EQUAL( elem41_dof_indices[7], elem45_dof_indices[3] );

    CPPUNIT_ASSERT_EQUAL( elem42_dof_indices[4], elem46_dof_indices[0] );
    CPPUNIT_ASSERT_EQUAL( elem42_dof_indices[5], elem46_dof_indices[1] );
    CPPUNIT_ASSERT_EQUAL( elem42_dof_indices[6], elem46_dof_indices[2] );
    CPPUNIT_ASSERT_EQUAL( elem42_dof_indices[7], elem46_dof_indices[3] );

    CPPUNIT_ASSERT_EQUAL( elem43_dof_indices[4], elem47_dof_indices[0] );
    CPPUNIT_ASSERT_EQUAL( elem43_dof_indices[5], elem47_dof_indices[1] );
    CPPUNIT_ASSERT_EQUAL( elem43_dof_indices[6], elem47_dof_indices[2] );
    CPPUNIT_ASSERT_EQUAL( elem43_dof_indices[7], elem47_dof_indices[3] );

    CPPUNIT_ASSERT_EQUAL( elem44_dof_indices[4], elem48_dof_indices[0] );
    CPPUNIT_ASSERT_EQUAL( elem44_dof_indices[5], elem48_dof_indices[1] );
    CPPUNIT_ASSERT_EQUAL( elem44_dof_indices[6], elem48_dof_indices[2] );
    CPPUNIT_ASSERT_EQUAL( elem44_dof_indices[7], elem48_dof_indices[3] );

    /* The dofs for the QUAD4 elements (ids = 53, 54, 55, and 56) should be the
       same for shared nodes. */
    CPPUNIT_ASSERT_EQUAL( elem53_dof_indices[1], elem54_dof_indices[0] );
    CPPUNIT_ASSERT_EQUAL( elem53_dof_indices[2], elem54_dof_indices[3] );
    CPPUNIT_ASSERT_EQUAL( elem53_dof_indices[3], elem55_dof_indices[0] );
    CPPUNIT_ASSERT_EQUAL( elem53_dof_indices[2], elem55_dof_indices[1] );
    CPPUNIT_ASSERT_EQUAL( elem53_dof_indices[2], elem56_dof_indices[0] );
    CPPUNIT_ASSERT_EQUAL( elem54_dof_indices[3], elem55_dof_indices[1] );
    CPPUNIT_ASSERT_EQUAL( elem54_dof_indices[3], elem56_dof_indices[0] );
    CPPUNIT_ASSERT_EQUAL( elem54_dof_indices[2], elem56_dof_indices[1] );
    CPPUNIT_ASSERT_EQUAL( elem55_dof_indices[1], elem56_dof_indices[0] );
    CPPUNIT_ASSERT_EQUAL( elem55_dof_indices[2], elem56_dof_indices[3] );


#endif
  }
};


CPPUNIT_TEST_SUITE_REGISTRATION( MixedDimensionMeshTest );
CPPUNIT_TEST_SUITE_REGISTRATION( MixedDimensionRefinedMeshTest );
CPPUNIT_TEST_SUITE_REGISTRATION( MixedDimensionNonUniformRefinement );
CPPUNIT_TEST_SUITE_REGISTRATION( MixedDimensionNonUniformRefinementTriangle );
CPPUNIT_TEST_SUITE_REGISTRATION( MixedDimensionNonUniformRefinement3D );
