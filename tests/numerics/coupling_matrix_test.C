// Ignore unused parameter warnings coming from cppuint headers
#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

// libmesh includes
#include <libmesh/coupling_matrix.h>

using namespace libMesh;

class CouplingMatrixTest : public CppUnit::TestCase
{
public:
  void setUp() {}

  void tearDown() {}

  CPPUNIT_TEST_SUITE(CouplingMatrixTest);

  CPPUNIT_TEST(testSimpleAPI);

  CPPUNIT_TEST_SUITE_END();


private:
  void testSimpleAPI()
  {
    CouplingMatrix cm(2);

    cm(0,1) = 1;

    bool cm01 = cm(0,1);
    CPPUNIT_ASSERT_EQUAL(cm01, true);

    cm(1,0) = 1;

    for (unsigned i=0; i<2; ++i)
      for (unsigned j=0; j<2; ++j)
        {
          bool cmij = cm(i,j);
          CPPUNIT_ASSERT_EQUAL(cmij, (i != j));
        }

    cm.resize(8);

    for (unsigned i=0; i<8; ++i)
      for (unsigned j=0; j<8; ++j)
        {
          bool cmij = cm(i,j);
          CPPUNIT_ASSERT_EQUAL(cmij, false);
        }

    // Set some elements true, in a weird order.
    for (unsigned i=6; i>0; --i)
      {
        const unsigned int pi = i + (i > 4);
        for (unsigned j=0; j<6; ++j)
          {
            const unsigned int pj = j + (j > 3);
            cm(pi, pj) = true;
          }
      }

    // Now the tensor product of {1,2,3,4,6,7} with {0,1,2,3,5,6}
    // should be 1.
    for (unsigned i=0; i<8; ++i)
      for (unsigned j=0; j<8; ++j)
        {
          bool cmij = cm(i,j);
          if ((i != 0) && (i != 5) && (j != 4) && (j != 7))
            {
              CPPUNIT_ASSERT_EQUAL(cmij, true);
            }
          else
            {
              CPPUNIT_ASSERT_EQUAL(cmij, false);
            }
        }

    // Set some elements to false.
    for (unsigned k=0; k<8; ++k)
      {
        cm(3, k) = false;
        cm(k, 0) = false;
      }

    // Now the tensor product of {1,2,4,6,7} with {1,2,3,5,6}
    // should be 1.
    for (unsigned i=0; i<8; ++i)
      for (unsigned j=0; j<8; ++j)
        {
          bool cmij = cm(i,j);
          if ((i != 0) && (i != 3) && (i != 5) &&
              (j != 0) && (j != 4) && (j != 7))
            {
              CPPUNIT_ASSERT_EQUAL(cmij, true);
            }
          else
            {
              CPPUNIT_ASSERT_EQUAL(cmij, false);
            }
        }
  }

};

CPPUNIT_TEST_SUITE_REGISTRATION(CouplingMatrixTest);
