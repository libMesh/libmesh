// libMesh includes
#include "libmesh/libmesh_exceptions.h"
#include "libmesh/rb_parameters.h"
#include "libmesh/rb_parametrized.h"
#include "libmesh/simple_range.h"
#include "libmesh/int_range.h"

// CPPUnit includes
#include "libmesh_cppunit.h"
#include <cppunit/TestAssert.h>

using namespace libMesh;

class RBParametersTest : public CppUnit::TestCase
{
public:
  LIBMESH_CPPUNIT_TEST_SUITE ( RBParametersTest );
  CPPUNIT_TEST( testScalar );
  CPPUNIT_TEST( testOldConstructor );
  CPPUNIT_TEST( testIterators );
  CPPUNIT_TEST( testIteratorsWithSamples );
  CPPUNIT_TEST( testAppend );
  CPPUNIT_TEST( testNSamples );
  CPPUNIT_TEST( testMultiValued );
  CPPUNIT_TEST( testRBParametrized );
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
    m.insert(params.begin_serialized(), params.end_serialized());

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

  void testIteratorsWithSamples()
  {
    LOG_UNIT_TEST;

    RBParameters params;
    params.set_value("a", 2, 1.);   // "a" : [0,0,1]
    params.set_value("b", 2, 2.);   // "b" : [0,0,2]
    params.set_value("c", 2, 3.);   // "c" : [0,0,3]

    // The iterators work on a serialized version of the map, meaning we see an
    // iterator for each sample of each parameter, e.g. {a,0}, {a,0}, {a,1}.
    // map.insert(begin(),end()) says "it is unspecified which element is inserted" in this case,
    // so instead we manually iterate and overwrite with the value from the latest sample.
    std::map<std::string, Real> m;
    // m.insert(params.begin(), params.end());  // unspecified behavior
    for (const auto & [key,val] : as_range(params.begin_serialized(), params.end_serialized()))
      m[key] = val;

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

  void testNSamples()
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
    params.push_back_value("a", 0.);
    params.push_back_value("a", 1.);
    params.push_back_value("a", 2.);
    params.push_back_value("a", 3.);
    CPPUNIT_ASSERT_EQUAL(params.n_samples(), static_cast<unsigned int>(4));
    CPPUNIT_ASSERT_EQUAL(params.get_sample_value("a", 2), static_cast<Real>(2.));

    // Test set_value() with an index.
    params.set_value("b", 3, 0.);
    params.set_value("b", 2, 1.);
    params.set_value("b", 1, 2.);
    params.set_value("b", 0, 3.);
    CPPUNIT_ASSERT_EQUAL(params.n_samples(), static_cast<unsigned int>(4));
    CPPUNIT_ASSERT_EQUAL(params.n_parameters(), static_cast<unsigned int>(2));
    CPPUNIT_ASSERT_EQUAL(params.get_sample_value("a", 2), static_cast<Real>(2.));
    CPPUNIT_ASSERT_EQUAL(params.get_sample_value("b", 2), static_cast<Real>(1.));

    // Test the default getters.
    CPPUNIT_ASSERT_EQUAL(params.get_sample_value("a", 5, 100.), static_cast<Real>(100.));
    CPPUNIT_ASSERT_EQUAL(params.get_sample_value("b", 5, 200.), static_cast<Real>(200.));
    CPPUNIT_ASSERT_EQUAL(params.get_sample_value("c", 5, 300.), static_cast<Real>(300.));

    // Test some errors.
#ifdef LIBMESH_ENABLE_EXCEPTIONS
    CPPUNIT_ASSERT_THROW(params.get_sample_value("b", 5), libMesh::LogicError); // sample_idx 5 does not exist.
    CPPUNIT_ASSERT_THROW(params.get_sample_value("c", 0), libMesh::LogicError); // parameter "c" does not exist.
    CPPUNIT_ASSERT_THROW(params.get_value("a"), libMesh::LogicError);           // a has multiple samples.
    CPPUNIT_ASSERT_THROW(params.get_value("a", 1.0), libMesh::LogicError);      // a has multiple samples.
#endif
  }

  void testMultiValued()
  {
    LOG_UNIT_TEST;

    RBParameters params;
    std::vector<Real> test_vec1 = {2.1, 2.2, 2.3};
    params.set_value("a", 0, {1.1, 1.2, 1.3});
    params.set_value("a", 1, test_vec1);
    params.set_value("a", 2, {3.1, 3.2, 3.3});
    params.set_value("b", 2, {3.1, 3.2, 3.3});

#ifdef LIBMESH_ENABLE_EXCEPTIONS
    // Test some errors.
    CPPUNIT_ASSERT_THROW(params.get_sample_value("b", 0), libMesh::LogicError); // single-value requested, but multi-valued.
    CPPUNIT_ASSERT_THROW(params.get_sample_value("b", 5), libMesh::LogicError); // sample_idx 5 does not exist.
    CPPUNIT_ASSERT_THROW(params.get_sample_value("c", 0), libMesh::LogicError); // parameter "c" does not exist.
    CPPUNIT_ASSERT_THROW(params.get_value("a"), libMesh::LogicError);           // a has multiple samples.
    CPPUNIT_ASSERT_THROW(params.get_value("a", 1.0), libMesh::LogicError);      // a has multiple samples.
    CPPUNIT_ASSERT_THROW(params.get_sample_vector_value("a", 3), libMesh::LogicError);  // sample_idx 3 does not exist.
    CPPUNIT_ASSERT_THROW(params.get_sample_vector_value("c", 0), libMesh::LogicError);  // parameter "c" does not exist.
#endif

    const auto & vector_value_1 = params.get_sample_vector_value("a", 1);
    for (const auto index : index_range(vector_value_1))
      CPPUNIT_ASSERT_EQUAL(vector_value_1[index], test_vec1[index]);
    const auto & vector_value_2 = params.get_sample_vector_value("a", 3, test_vec1);
    for (const auto index : index_range(vector_value_2))
      CPPUNIT_ASSERT_EQUAL(vector_value_2[index], test_vec1[index]);

    std::vector<Real> test_vec2 = {5.1, 5.2, 5.3};
    params.push_back_value("a", test_vec2);
    params.push_back_value("b", test_vec2);
    const auto & vector_value_3 = params.get_sample_vector_value("a", 3);
    for (const auto index : index_range(vector_value_3))
      CPPUNIT_ASSERT_EQUAL(vector_value_3[index], test_vec2[index]);
  }

  void testRBParametrized()
  {
    LOG_UNIT_TEST;

    RBParameters mu_min, mu_max;

    mu_min.set_value("a", -10.);
    mu_max.set_value("a", -11.);

    RBParametrized rb_parametrized;
    // rb_parametrized.verbose_mode = true; // Enable for more printed details.

#ifdef LIBMESH_ENABLE_EXCEPTIONS
    // Throw an error due to invalid min/max.
    CPPUNIT_ASSERT_THROW(rb_parametrized.initialize_parameters(mu_min, mu_max, {}), libMesh::LogicError);
#endif

    // Define an invalid max RBParameter with multiple steps.
    mu_max.set_value("a", 2, -40.);

#ifdef LIBMESH_ENABLE_EXCEPTIONS
    // Throw an error due to invalid max value.
    CPPUNIT_ASSERT_THROW(rb_parametrized.initialize_parameters(mu_min,mu_max, {}), libMesh::LogicError);
#endif

    // Set the max value correctly and initialize the RBParametrized object
    mu_max = RBParameters();
    mu_max.set_value("a",  10.);
    rb_parametrized.initialize_parameters(mu_min, mu_max, {});

    RBParameters params;

    // Value == min --> OK.
    params.set_value("a", 0,  0.);
    params.set_value("a", 1, -10.);
    CPPUNIT_ASSERT(rb_parametrized.set_parameters(params));

    // Value == max --> OK.
    params.set_value("a", 2,  10.);
    CPPUNIT_ASSERT(rb_parametrized.set_parameters(params));

    // Value < min --> not OK.
    params.set_value("a", 3, -40.);
    CPPUNIT_ASSERT(!rb_parametrized.set_parameters(params));

#ifdef LIBMESH_ENABLE_EXCEPTIONS
    // Throw an error due to different number of parameters.
    mu_max.set_value("b", 3,  40.);
    CPPUNIT_ASSERT_THROW(rb_parametrized.initialize_parameters(mu_min, mu_max, {}), libMesh::LogicError);
#endif
  }

};

CPPUNIT_TEST_SUITE_REGISTRATION ( RBParametersTest );
