#include <libmesh/laspack_vector.h>

#ifdef LIBMESH_HAVE_LASPACK

#include "numeric_vector_test.h"


using namespace libMesh;

class LaspackVectorTest : public NumericVectorTest<LaspackVector<Number>>
{
private:
  // This class manages the memory for the communicator it uses,
  // providing a dumb pointer to the managed resource for the base
  // class.
  std::unique_ptr<Parallel::Communicator> _managed_comm;

public:
  void setUp()
  {
    // Laspack doesn't support distributed parallel vectors, but we
    // can build a serial vector on each processor
    _managed_comm = std::make_unique<Parallel::Communicator>();

    // Base class communicator points to our managed communicator
    my_comm = _managed_comm.get();

    this->NumericVectorTest<LaspackVector<Number>>::setUp();
  }

  void tearDown() {}

  LaspackVectorTest() :
    NumericVectorTest<LaspackVector<Number>>() {
    if (unitlog->summarized_logs_enabled())
      this->libmesh_suite_name = "NumericVectorTest";
    else
      this->libmesh_suite_name = "LaspackVectorTest";
  }

  CPPUNIT_TEST_SUITE( LaspackVectorTest );

  NUMERICVECTORTEST

  CPPUNIT_TEST_SUITE_END();
};

CPPUNIT_TEST_SUITE_REGISTRATION( LaspackVectorTest );

#endif // #ifdef LIBMESH_HAVE_LASPACK

