#include <libmesh/laspack_vector.h>

#ifdef LIBMESH_HAVE_LASPACK

#include "numeric_vector_test.h"
#include "test_comm.h"

#include <regex>


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

  void testDistributedInit()
  {
    LOG_UNIT_TEST;

    const numeric_index_type n = 10;

#ifdef LIBMESH_ENABLE_EXCEPTIONS
    // Rank 0 owning every entry, as it does the DoFs of a mesh with a single element, still leaves
    // the vector distributed, and every rank has to report that
    if (TestCommWorld->size() > 1)
      {
        const std::string expected = "LaspackVectors can only be used in serial";
        bool threw = false;
        try
          {
            LaspackVector<Number> distributed(*TestCommWorld, n, TestCommWorld->rank() ? 0 : n);
          }
        catch (libMesh::LogicError & e)
          {
            CPPUNIT_ASSERT_MESSAGE(e.what(), std::regex_search(e.what(), std::regex(expected)));
            threw = true;
          }
        CPPUNIT_ASSERT_MESSAGE("Expected an error containing \"" + expected + "\"", threw);
      }
#endif

    // A vector that every rank holds whole is serial in effect
    LaspackVector<Number> replicated(*TestCommWorld, n, n);
    CPPUNIT_ASSERT_EQUAL(n, replicated.local_size());
  }

  LaspackVectorTest() :
    NumericVectorTest<LaspackVector<Number>>() {
    if (unitlog->summarized_logs_enabled())
      this->libmesh_suite_name = "NumericVectorTest";
    else
      this->libmesh_suite_name = "LaspackVectorTest";
  }

  CPPUNIT_TEST_SUITE( LaspackVectorTest );

  NUMERICVECTORTEST
  CPPUNIT_TEST( testDistributedInit );

  CPPUNIT_TEST_SUITE_END();
};

CPPUNIT_TEST_SUITE_REGISTRATION( LaspackVectorTest );

#endif // #ifdef LIBMESH_HAVE_LASPACK

