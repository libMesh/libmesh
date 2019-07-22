// Ignore unused parameter warnings coming from cppunit headers
#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

#include <libmesh/node.h>

#include "../geom/point_test.h"
#include "../base/dof_object_test.h"

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

class NodeTest : public PointTestBase<Node>, public DofObjectTest<Node> {
public:
  CPPUNIT_TEST_SUITE( NodeTest );

  // These tests currently use the Node copy constructor, which is marked
  // as deprecated, so only instantiate them if deprecated code is allowed.
#ifdef LIBMESH_ENABLE_DEPRECATED
  POINTTEST
  DOFOBJECTTEST
#endif

  CPPUNIT_TEST_SUITE_END();

private:

  Node * dof_object_instance;

public:

  virtual void setUp()
  {
    PointTestBase<Node>::setUp();

    dof_object_instance = new Node(1,1,1);
    DofObjectTest<Node>::setUp(dof_object_instance);
  }

  virtual void tearDown()
  {
    PointTestBase<Node>::tearDown();

    delete dof_object_instance;
  }

};

CPPUNIT_TEST_SUITE_REGISTRATION( NodeTest );
