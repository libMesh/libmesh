#include <libmesh/node.h>

#include "../geom/point_test.h"
#include "../base/dof_object_test.h"

#include "libmesh_cppunit.h"

using namespace libMesh;

class NodeTest : public PointTestBase<Node>, public DofObjectTest<Node> {
public:
  NodeTest() :
    PointTestBase<Node>(), DofObjectTest<Node>() {
    this->PointTestBase<Node>::libmesh_suite_name = "NodeTest";
    this->DofObjectTest<Node>::libmesh_suite_name = "NodeTest";
  }
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
