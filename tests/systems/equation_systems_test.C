// Ignore unused parameter warnings coming from cppunit headers
#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

#include <libmesh/elem.h>
#include <libmesh/equation_systems.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_refinement.h>
#include <libmesh/remote_elem.h>
#include <libmesh/replicated_mesh.h>

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

class EquationSystemsTest : public CppUnit::TestCase {
public:
  CPPUNIT_TEST_SUITE( EquationSystemsTest );

  CPPUNIT_TEST( testConstruction );
  CPPUNIT_TEST( testAddSystem );
  CPPUNIT_TEST( testInit );
  CPPUNIT_TEST( testPostInitAddSystem );
  CPPUNIT_TEST( testPostInitAddElem );
  CPPUNIT_TEST( testRefineThenReinitPreserveFlags );

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

  void testRefineThenReinitPreserveFlags()
  {
    Mesh mesh(*TestCommWorld);
    EquationSystems es(mesh);
    System & sys = es.add_system<System> ("SimpleSystem");
    sys.add_variable("u", FIRST);
    MeshTools::Generation::build_square(mesh,2,1);
    es.init();

    Elem * to_refine = mesh.query_elem_ptr(0);
    if (to_refine)
      to_refine->set_refinement_flag(Elem::RefinementState::REFINE);

    MeshRefinement mr(mesh);
    mr.refine_elements();
    es.disable_refine_in_reinit();
    es.reinit();

    if (mesh.query_elem_ptr(1))
    CPPUNIT_ASSERT( mesh.elem(1)->active() );

    const Elem * elem = mesh.query_elem_ptr(0);
    if (elem)
      {
        CPPUNIT_ASSERT_EQUAL( Elem::RefinementState::INACTIVE,elem->refinement_flag() );

        for (unsigned int c=0; c<elem->n_children(); c++)
          if (elem->child(c) != remote_elem)
            CPPUNIT_ASSERT_EQUAL(Elem::RefinementState::JUST_REFINED,
                                 elem->child(c)->refinement_flag());
      }
  }





};

CPPUNIT_TEST_SUITE_REGISTRATION( EquationSystemsTest );
