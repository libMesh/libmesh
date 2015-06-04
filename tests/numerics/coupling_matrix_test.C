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

  CPPUNIT_TEST(testIteratorAPI);

  CPPUNIT_TEST_SUITE_END();


private:
  void testSimpleAPI()
  {
    CouplingMatrix cm(2);

    // Use a constant reference to make sure we test both const and
    // non-const operator() implementations
    const CouplingMatrix& cmr = cm;

    cm(0,1) = 1;

    bool cm01 = cm(0,1);
    CPPUNIT_ASSERT_EQUAL(cm01, true);

    cm(1,0) = 1;

    for (unsigned i=0; i<2; ++i)
      for (unsigned j=0; j<2; ++j)
        {
          bool cmij = cm(i,j);
          bool cmrij = cmr(i,j);
          CPPUNIT_ASSERT_EQUAL(cmij, cmrij);
          CPPUNIT_ASSERT_EQUAL(cmij, (i != j));
        }

    cm.resize(8);

    for (unsigned i=0; i<8; ++i)
      for (unsigned j=0; j<8; ++j)
        {
          bool cmij = cm(i,j);
          bool cmrij = cmr(i,j);
          CPPUNIT_ASSERT_EQUAL(cmij, cmrij);
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
          bool cmrij = cmr(i,j);
          CPPUNIT_ASSERT_EQUAL(cmij, cmrij);
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
          bool cmrij = cmr(i,j);
          CPPUNIT_ASSERT_EQUAL(cmij, cmrij);
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

  void testIteratorAPI()
  {
    CouplingMatrix cm(8);

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

    // Set some elements to false.
    for (unsigned k=0; k<8; ++k)
      {
        cm(3, k) = false;
        cm(k, 0) = false;
      }

    // Now the tensor product of {1,2,4,6,7} with {1,2,3,5,6}
    // should be 1.
    const unsigned int ivals[] = {1,2,4,6,7};
    const unsigned int non_ivals[] = {0,3,5};
    const unsigned int jvals[] = {1,2,3,5,6};
    // const unsigned int non_jvals[] = {0,4,7};

    const unsigned int isize = sizeof(unsigned int);

    for (unsigned int pi = 0; pi != sizeof(non_ivals)/isize; ++pi)
      {
        unsigned int i = non_ivals[pi];
        ConstCouplingRow ccr(i,cm);
        CPPUNIT_ASSERT(ccr.begin() == ccr.end());
      }

    for (unsigned int pi = 0; pi != sizeof(ivals)/isize; ++pi)
      {
        unsigned int i = ivals[pi];
        ConstCouplingRow ccr(i,cm);

        ConstCouplingRow::const_iterator ccr_it = ccr.begin();

        for (unsigned int pj = 0; pj != sizeof(jvals)/isize; ++pj)
          {
            CPPUNIT_ASSERT(ccr_it != ccr.end());
            CPPUNIT_ASSERT_EQUAL(*ccr_it, jvals[pj]);
            ++ccr_it;
          }

        CPPUNIT_ASSERT(ccr_it == ccr.end());
      }
  }


};

CPPUNIT_TEST_SUITE_REGISTRATION(CouplingMatrixTest);
