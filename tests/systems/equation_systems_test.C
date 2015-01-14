// Ignore unused parameter warnings coming from cppuint headers
#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

#include <libmesh/equation_systems.h>
#include <libmesh/mesh.h>

#include "test_comm.h"

using namespace libMesh;

class EquationSystemsTest : public CppUnit::TestCase {
public:
  CPPUNIT_TEST_SUITE( EquationSystemsTest );

  CPPUNIT_TEST( testConstruction );

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

};

CPPUNIT_TEST_SUITE_REGISTRATION( EquationSystemsTest );
