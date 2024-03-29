#include <libmesh/distributed_vector.h>

#include "numeric_vector_test.h"


using namespace libMesh;

class DistributedVectorTest : public NumericVectorTest<DistributedVector<Number>> {
public:
  DistributedVectorTest() :
    NumericVectorTest<DistributedVector<Number>>() {
    if (unitlog->summarized_logs_enabled())
      this->libmesh_suite_name = "NumericVectorTest";
    else
      this->libmesh_suite_name = "DistributedVectorTest";
  }
  CPPUNIT_TEST_SUITE( DistributedVectorTest );

  NUMERICVECTORTEST

  CPPUNIT_TEST_SUITE_END();
};

CPPUNIT_TEST_SUITE_REGISTRATION( DistributedVectorTest );
