
#include "libmesh_cppunit.h"

#ifdef LIBMESH_HAVE_NETGEN
namespace nglib {
#include "netgen/nglib/nglib.h"
}
#endif

#include <numeric>

class LibMeshNetgenTest : public CppUnit::TestCase {
public:
  LIBMESH_CPPUNIT_TEST_SUITE( LibMeshNetgenTest );

#ifdef LIBMESH_HAVE_NETGEN
  CPPUNIT_TEST( testLibMeshNetgen );
#endif

  CPPUNIT_TEST_SUITE_END();

private:

public:
  void setUp()
  {}

  void tearDown()
  {}

#ifdef LIBMESH_HAVE_NETGEN
  void testLibMeshNetgen ()
  {
    LOG_UNIT_TEST;

    auto ngmesh = nglib::Ng_NewMesh();
    nglib::Ng_DeleteMesh(ngmesh);
  }
#endif
};

CPPUNIT_TEST_SUITE_REGISTRATION( LibMeshNetgenTest );
