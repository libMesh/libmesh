// Ignore unused parameter warnings coming from cppunit headers
#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

#include <libmesh/node.h>

#include "../geom/point_test.h"
#include "../base/dof_object_test.h"

using namespace libMesh;

class NodeTest : public PointTestBase<Node>, public DofObjectTest<Node> {
public:
  CPPUNIT_TEST_SUITE( NodeTest );

  POINTTEST

  DOFOBJECTTEST

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
