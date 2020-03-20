#include <libmesh/parameters.h>

#include "libmesh_cppunit.h"

#include <map>
#include <string>

using namespace libMesh;

class ParametersTest : public CppUnit::TestCase {
public:
  CPPUNIT_TEST_SUITE( ParametersTest );

  CPPUNIT_TEST( testInt );
  CPPUNIT_TEST( testFloat );
  CPPUNIT_TEST( testDouble );

  CPPUNIT_TEST( testMap );

  CPPUNIT_TEST_SUITE_END();

private:

public:
  void setUp()
  {}

  void tearDown()
  {}

  template <typename T>
  void testScalar ()
  {
    Parameters param;

    T t = 10, t_orig = 10;

    param.set<T>("first") = t;

    ++t;

    param.set<T>("second") = t;

    CPPUNIT_ASSERT_EQUAL(t_orig, param.get<T>("first"));
    CPPUNIT_ASSERT_EQUAL(t, param.get<T>("second"));
  }

  void testMap ()
  {
    Parameters param;

    std::map<int, std::string> m;
    m[2] = "two";
    m[4] = "four";
    param.set<std::map<int, std::string>>("mymap") = m;
    const std::map<int, std::string> & gotten =
      param.get<std::map<int, std::string>>("mymap");

    CPPUNIT_ASSERT_EQUAL(gotten.size(), std::size_t(2));
    CPPUNIT_ASSERT_EQUAL(gotten.at(2), std::string("two"));
    CPPUNIT_ASSERT_EQUAL(gotten.at(4), std::string("four"));
  }

  void testInt () { testScalar<int>(); }
  void testFloat () { testScalar<float>(); }
  void testDouble () { testScalar<double>(); }
};

CPPUNIT_TEST_SUITE_REGISTRATION( ParametersTest );
