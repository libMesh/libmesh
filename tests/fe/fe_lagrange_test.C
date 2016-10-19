
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>

#include "fe_test.h"

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

#define INSTANTIATE_FETEST(order, family, elemtype) \
class FETest_##order##_##family##_##elemtype : public FETest<order, family, elemtype> { \
public: \
  CPPUNIT_TEST_SUITE( FETest_##order##_##family##_##elemtype ); \
  FETEST \
  CPPUNIT_TEST_SUITE_END(); \
}; \
 \
CPPUNIT_TEST_SUITE_REGISTRATION( FETest_##order##_##family##_##elemtype );

INSTANTIATE_FETEST(FIRST, LAGRANGE, EDGE2);
INSTANTIATE_FETEST(FIRST, LAGRANGE, TRI3);
INSTANTIATE_FETEST(FIRST, LAGRANGE, QUAD4);
INSTANTIATE_FETEST(FIRST, LAGRANGE, TET4);
INSTANTIATE_FETEST(FIRST, LAGRANGE, HEX8);
INSTANTIATE_FETEST(FIRST, LAGRANGE, PRISM6);
