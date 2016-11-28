// Ignore unused parameter warnings coming from cppunit headers
#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

#include <libmesh/equation_systems.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/replicated_mesh.h>
#include <libmesh/elem.h>

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





};

CPPUNIT_TEST_SUITE_REGISTRATION( EquationSystemsTest );
