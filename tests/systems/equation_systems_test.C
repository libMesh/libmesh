#include <libmesh/dof_map.h>
#include <libmesh/elem.h>
#include <libmesh/equation_systems.h>
#include <libmesh/ghost_point_neighbors.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_refinement.h>
#include <libmesh/remote_elem.h>
#include <libmesh/replicated_mesh.h>
#include <libmesh/node_elem.h>

#include "test_comm.h"
#include "libmesh_cppunit.h"


using namespace libMesh;

// Anonymous namespace to avoid linker conflicts
namespace {

Number bilinear_test (const Point& p,
                      const Parameters&,
                      const std::string&,
                      const std::string&)
{
  const Real & x = p(0);
  const Real & y = p(1);

  return 4*x*y - 3*x + 2*y - 1;
}

}

class EquationSystemsTest : public CppUnit::TestCase {
public:
  CPPUNIT_TEST_SUITE( EquationSystemsTest );

  CPPUNIT_TEST( testConstruction );
  CPPUNIT_TEST( testAddSystem );
  CPPUNIT_TEST( testInit );
  CPPUNIT_TEST( testPostInitAddSystem );
  CPPUNIT_TEST( testPostInitAddElem );
  CPPUNIT_TEST( testReinitWithNodeElem );
#if LIBMESH_DIM > 1
  CPPUNIT_TEST( testRefineThenReinitPreserveFlags );
#ifdef LIBMESH_ENABLE_AMR // needs project_solution, even for reordering
  CPPUNIT_TEST( testRepartitionThenReinit );
#endif
#endif
  CPPUNIT_TEST( testDisableDefaultGhosting );

  CPPUNIT_TEST_SUITE_END();

private:

public:
  void setUp()
  {}

  void tearDown()
  {}



  void testConstruction()
  {
    Mesh mesh(*TestCommWorld);
    EquationSystems es(mesh);
  }

  void testAddSystem()
  {
    Mesh mesh(*TestCommWorld);
    EquationSystems es(mesh);
    /*System &sys = */es.add_system<System> ("SimpleSystem");
  }

  void testInit()
  {
    Mesh mesh(*TestCommWorld);
    EquationSystems es(mesh);
    /*System &sys = */es.add_system<System> ("SimpleSystem");
    MeshTools::Generation::build_point(mesh);
    es.init();
  }

  void testPostInitAddSystem()
  {
    Mesh mesh(*TestCommWorld);
    MeshTools::Generation::build_point(mesh);
    EquationSystems es(mesh);
    /*System &sys1 = */es.add_system<System> ("SimpleSystem");
    es.init();
    /*System &sys2 = */es.add_system<System> ("SecondSystem");
    es.reinit();
  }

  void testPostInitAddRealSystem()
  {
    Mesh mesh(*TestCommWorld);
    MeshTools::Generation::build_point(mesh);
    EquationSystems es(mesh);
    System &sys1 = es.add_system<System> ("SimpleSystem");
    sys1.add_variable("u1", FIRST);
    es.init();
    System &sys2 = es.add_system<System> ("SecondSystem");
    sys2.add_variable("u2", FIRST);
    es.reinit();
  }

  void testPostInitAddElem()
  {
    ReplicatedMesh mesh(*TestCommWorld);

    EquationSystems es(mesh);
    System &sys = es.add_system<System> ("SimpleSystem");
    sys.add_variable("u", FIRST);

    MeshTools::Generation::build_line (mesh, 10, 0., 1., EDGE2);
    es.init();

    Elem* e = Elem::build(EDGE2).release();
    e->set_id(mesh.max_elem_id());
    e->processor_id() = 0;
    e->set_node(0) = mesh.node_ptr(2);
    e->set_node(1) = mesh.node_ptr(8);
    mesh.add_elem(e);
    mesh.prepare_for_use();

    es.reinit();
  }

  void testReinitWithNodeElem()
  {
    ReplicatedMesh mesh(*TestCommWorld);

    MeshTools::Generation::build_line (mesh, 10, 0., 1., EDGE2);
    Elem* node_elem = mesh.add_elem (new NodeElem);
    node_elem->set_node(0) = mesh.node_ptr(0);
    mesh.prepare_for_use();

    EquationSystems es(mesh);
    System &sys = es.add_system<System> ("SimpleSystem");
    sys.add_variable("u", CONSTANT, MONOMIAL);
    es.init();
    es.reinit();
  }

  void testRefineThenReinitPreserveFlags()
  {
    // This test requires AMR support since it sets refinement flags.
#ifdef LIBMESH_ENABLE_AMR
    Mesh mesh(*TestCommWorld);
    mesh.allow_renumbering(false);
    EquationSystems es(mesh);
    System & sys = es.add_system<System> ("SimpleSystem");
    sys.add_variable("u", FIRST);
    MeshTools::Generation::build_square(mesh,2,1);
    es.init();

    Elem * to_refine = mesh.query_elem_ptr(0);
    if (to_refine)
      to_refine->set_refinement_flag(Elem::REFINE);

    MeshRefinement mr(mesh);
    mr.refine_elements();
    es.disable_refine_in_reinit();
    es.reinit();

    if (mesh.query_elem_ptr(1))
      CPPUNIT_ASSERT( mesh.elem_ptr(1)->active() );

    const Elem * elem = mesh.query_elem_ptr(0);
    if (elem)
      {
        CPPUNIT_ASSERT_EQUAL( Elem::INACTIVE,elem->refinement_flag() );

        for (unsigned int c=0; c<elem->n_children(); c++)
          if (elem->child_ptr(c) != RemoteElem::get_instance())
            CPPUNIT_ASSERT_EQUAL(Elem::JUST_REFINED,
                                 elem->child_ptr(c)->refinement_flag());
      }
#endif
  }



  void testRepartitionThenReinit()
  {
    Mesh mesh(*TestCommWorld);
    mesh.allow_renumbering(false);
    EquationSystems es(mesh);
    System & sys = es.add_system<System> ("SimpleSystem");
    sys.add_variable("u", FIRST);
    MeshTools::Generation::build_square(mesh,5,5);
    es.init();
    sys.project_solution(bilinear_test, NULL, es.parameters);

    // Force (in parallel) a different partitioning - we'll simply put
    // everything on rank 0, which hopefully is not what our default
    // partitioner did!
    mesh.partition(1);

    // Make sure the solution is still intact after reinit
    es.reinit();

    for (Real x = 0.1; x < 1; x += 0.2)
      for (Real y = 0.1; y < 1; y += 0.2)
        {
          Point p(x,y);
          LIBMESH_ASSERT_FP_EQUAL(libmesh_real(sys.point_value(0,p)),
                                  libmesh_real(bilinear_test(p,es.parameters,"","")),
                                  TOLERANCE*TOLERANCE);
        }
  }

  void testDisableDefaultGhosting()
  {
    Mesh mesh(*TestCommWorld);
    EquationSystems es(mesh);

    auto n_ghosts = [&mesh]() {
      return int(std::distance(mesh.ghosting_functors_begin(),
                               mesh.ghosting_functors_end()));
    };

    auto n_evaluables = [](System &sys) {
      return int(std::distance(sys.get_dof_map().algebraic_ghosting_functors_begin(),
                               sys.get_dof_map().algebraic_ghosting_functors_end()));
    };

    auto n_couplings = [](System &sys) {
      return int(std::distance(sys.get_dof_map().coupling_functors_begin(),
                               sys.get_dof_map().coupling_functors_end()));
    };

    // One default ghosting functor on the mesh
    CPPUNIT_ASSERT_EQUAL(n_ghosts(), 1);

    // Add another functor, making two
    auto gpn = std::make_shared<GhostPointNeighbors>(mesh);
    mesh.add_ghosting_functor(gpn);
    CPPUNIT_ASSERT_EQUAL(n_ghosts(), 2);

    // Remove the default, leaving just the user functor
    es.enable_default_ghosting(false);
    CPPUNIT_ASSERT_EQUAL(n_ghosts(), 1);

    // Which can be removed too
    mesh.remove_ghosting_functor(*gpn);
    CPPUNIT_ASSERT_EQUAL(n_ghosts(), 0);

    // Adding a new system shouldn't add any default ghosting if the
    // EquationSystems disabled it
    System & sys1 = es.add_system<System>("System1");
    CPPUNIT_ASSERT_EQUAL(n_ghosts(), 0);
    CPPUNIT_ASSERT_EQUAL(n_evaluables(sys1), 0);
    CPPUNIT_ASSERT_EQUAL(n_couplings(sys1), 0);

    // But if we reenable it then now we should have three functors,
    // with the default algebraic and coupling functors from sys1.
    //
    // We currently iterate over coupling functors manually even when
    // using them for evaluability... this test will need to change
    // eventually, when that does.
    es.enable_default_ghosting(true);
    CPPUNIT_ASSERT_EQUAL(n_ghosts(), 3);
    CPPUNIT_ASSERT_EQUAL(n_evaluables(sys1), 1);
    CPPUNIT_ASSERT_EQUAL(n_couplings(sys1), 1);

    // Adding a second system with default ghosting reenabled should
    // give us 2 more functors.
    System & sys2 = es.add_system<System>("System2");
    CPPUNIT_ASSERT_EQUAL(n_ghosts(), 5);
    CPPUNIT_ASSERT_EQUAL(n_evaluables(sys2), 1);
    CPPUNIT_ASSERT_EQUAL(n_couplings(sys2), 1);

    // Adding a user functor to evaluables and couplings should add it
    // to the mesh
    GhostPointNeighbors gpn2(mesh);
    sys1.get_dof_map().add_algebraic_ghosting_functor(gpn2);
    CPPUNIT_ASSERT_EQUAL(n_ghosts(), 6);
    CPPUNIT_ASSERT_EQUAL(n_evaluables(sys1), 2);
    CPPUNIT_ASSERT_EQUAL(n_couplings(sys1), 1);

    // Unless we say not to.
    auto gpn3 = std::make_shared<GhostPointNeighbors>(mesh);
    sys1.get_dof_map().add_coupling_functor(gpn3, /*to_mesh=*/false);
    CPPUNIT_ASSERT_EQUAL(n_ghosts(), 6);
    CPPUNIT_ASSERT_EQUAL(n_evaluables(sys1), 2);
    CPPUNIT_ASSERT_EQUAL(n_couplings(sys1), 2);

    // Turning off default coupling again should get rid of everything
    // except the user functors.
    es.enable_default_ghosting(false);
    CPPUNIT_ASSERT_EQUAL(n_ghosts(), 1);
    CPPUNIT_ASSERT_EQUAL(n_evaluables(sys1), 1);
    CPPUNIT_ASSERT_EQUAL(n_couplings(sys1), 1);
    CPPUNIT_ASSERT_EQUAL(n_evaluables(sys2), 0);
    CPPUNIT_ASSERT_EQUAL(n_couplings(sys2), 0);
  }





};

CPPUNIT_TEST_SUITE_REGISTRATION( EquationSystemsTest );
