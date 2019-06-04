#include <libmesh/node.h>

#include "../geom/point_test.h"
#include "../base/dof_object_test.h"

#include "libmesh_cppunit.h"

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
