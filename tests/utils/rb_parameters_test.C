// libMesh includes
#include "libmesh/rb_parameters.h"

// CPPUnit includes
#include "libmesh_cppunit.h"

using namespace libMesh;

class RBParametersTest : public CppUnit::TestCase
{
public:
  LIBMESH_CPPUNIT_TEST_SUITE ( RBParametersTest );
  CPPUNIT_TEST( testScalar );
  CPPUNIT_TEST( testOldConstructor );
  CPPUNIT_TEST( testIterators );
  CPPUNIT_TEST( testAppend );
  CPPUNIT_TEST( testNSteps );
  CPPUNIT_TEST_SUITE_END();

public:

  // virtual void setUp()
  // {}

  // virtual void tearDown()
  // {}

  void testScalar()
  {
    LOG_UNIT_TEST;

    // Test scalar-valued interfaces
    RBParameters params;
    params.set_value("a", 1.);
    params.set_value("b", 2.);
    params.set_value("c", 3.);

    // Expected result
    // a: 1.000000e+00
    // b: 2.000000e+00
    // c: 3.000000e+00
    CPPUNIT_ASSERT(params.has_value("a"));
    CPPUNIT_ASSERT(params.has_value("b"));
    CPPUNIT_ASSERT(params.has_value("c"));
    CPPUNIT_ASSERT_EQUAL(params.get_value("a"), 1.);
    CPPUNIT_ASSERT_EQUAL(params.get_value("b"), 2.);
    CPPUNIT_ASSERT_EQUAL(params.get_value("c"), 3.);
  }

  void testOldConstructor()
  {
    LOG_UNIT_TEST;

    // Test constructing an RBParameters object from a std::map
    std::map<std::string, Real> in = {{"a", 1.}, {"b", 2.}, {"c", 3.}};
    RBParameters params(in);

    // Expected result
    // a: 1.000000e+00
    // b: 2.000000e+00
    // c: 3.000000e+00
    CPPUNIT_ASSERT(params.has_value("a"));
    CPPUNIT_ASSERT(params.has_value("b"));
    CPPUNIT_ASSERT(params.has_value("c"));
    CPPUNIT_ASSERT_EQUAL(params.get_value("a"), 1.);
    CPPUNIT_ASSERT_EQUAL(params.get_value("b"), 2.);
    CPPUNIT_ASSERT_EQUAL(params.get_value("c"), 3.);

    // Test that RBParameters objects constructed with the old
    // constructor have the correct number of samples.
    CPPUNIT_ASSERT_EQUAL(params.n_samples(), 1u);
  }

  void testIterators()
  {
    LOG_UNIT_TEST;

    // Test creating a std::map using RBParameters iterators
    RBParameters params;
    params.set_value("a", 1.);
    params.set_value("b", 2.);
    params.set_value("c", 3.);

    std::map<std::string, Real> m;
    m.insert(params.begin(), params.end());

    // Expected result
    // a: 1.000000e+00
    // b: 2.000000e+00
    // c: 3.000000e+00
    CPPUNIT_ASSERT(m.count("a"));
    CPPUNIT_ASSERT(m.count("b"));
    CPPUNIT_ASSERT(m.count("c"));
    CPPUNIT_ASSERT_EQUAL(m["a"], 1.);
    CPPUNIT_ASSERT_EQUAL(m["b"], 2.);
    CPPUNIT_ASSERT_EQUAL(m["c"], 3.);
  }

  void testAppend()
  {
    LOG_UNIT_TEST;

    // Create first multi-sample RBParameters object
    RBParameters params1;
    for (int i=0; i<3; ++i)
      params1.push_back_value("a", Real(i));

    // Create second multi-sample RBParameters object
    // (must have same number of samples)
    RBParameters params2;
    for (int i=0; i<3; ++i)
      {
        params2.push_back_value("b", Real(i+3));
        params2.push_back_extra_value("c", Real(i*i));
      }

    // Append second onto first
    params1 += params2;

    // Print result
    // params1.print();

    // Check that the desired appending happened
    CPPUNIT_ASSERT(params1.has_value("b"));
    CPPUNIT_ASSERT(params1.has_extra_value("c"));
    for (int i=0; i<3; ++i)
      {
        CPPUNIT_ASSERT_EQUAL(params1.get_sample_value("b", i), static_cast<Real>(i+3));
        CPPUNIT_ASSERT_EQUAL(params1.get_extra_sample_value("c", i), static_cast<Real>(i*i));
      }
  }

  void testNSteps()
  {
    LOG_UNIT_TEST;

    // A default-constructed RBparameters object has 1 sample by definition
    RBParameters params;
    CPPUNIT_ASSERT_EQUAL(params.n_samples(), static_cast<unsigned int>(1));

    // Set the number of samples to use in the no-parameters case
    params.set_n_samples(10);
    CPPUNIT_ASSERT_EQUAL(params.n_samples(), static_cast<unsigned int>(10));

    // Define multiple samples for a single parameter. Now we no longer
    // use the set_n_samples() value, since we have actual samples.
    params.push_back_value("a", 1.);
    params.push_back_value("a", 2.);
    CPPUNIT_ASSERT_EQUAL(params.n_samples(), static_cast<unsigned int>(2));
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION ( RBParametersTest );
