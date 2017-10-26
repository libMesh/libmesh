#include <libmesh/laspack_vector.h>

#ifdef LIBMESH_HAVE_LASPACK

#include "numeric_vector_test.h"

#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

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

class LaspackVectorTest : public NumericVectorTest<LaspackVector<Number>>
{
public:
  void setUp()
  {
    // Laspack doesn't support distributed parallel vectors, but we
    // can build a serial vector on each processor
    my_comm = new Parallel::Communicator();
  }

  void tearDown()
  {
    delete my_comm;
  }

  CPPUNIT_TEST_SUITE( LaspackVectorTest );

  NUMERICVECTORTEST

  CPPUNIT_TEST_SUITE_END();
};

CPPUNIT_TEST_SUITE_REGISTRATION( LaspackVectorTest );

#endif // #ifdef LIBMESH_HAVE_LASPACK

