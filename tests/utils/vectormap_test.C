#include "libmesh/vectormap.h"
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>

#define VECTORMAPOBJECTTEST \
  CPPUNIT_TEST( testCreate );			\

using namespace libMesh;

class VectormapTest : public CppUnit::TestCase
{
public:
  CPPUNIT_TEST_SUITE ( VectormapTest );

  CPPUNIT_TEST( testCreate );
  CPPUNIT_TEST( testInsert );
  CPPUNIT_TEST( testIterate );

  CPPUNIT_TEST_SUITE_END();

private:

  template <typename Key, typename Val>
  void create()
  {
    vectormap<Key,Val> vm;
  }

  template <typename Key, typename Val>
  void insert()
  {
    vectormap<Key,Val> vm;

    Val val(0); // requires default constructor for val type.

    for (Key key=1; key<32; key*=2)
      vm.insert (std::make_pair(key,val));

    vm.sort();
  }

  template <typename Key, typename Val>
  void iterate(const Val &default_value=0)
  {
    vectormap<Key,Val> vm;

    Val val(default_value); // requires default constructor for val type.

    for (Key key=1; key<32; key*=2)
      vm.insert (std::make_pair(key,val));

    vm.sort();

    for (typename vectormap<Key,Val>::const_iterator it=vm.begin();
	 it != vm.end(); ++it)
      {
	const Key &ikey = it->first;
	const Val &ival = it->second;

	CPPUNIT_ASSERT       ( vm.count(ikey) == 1 );
	CPPUNIT_ASSERT_EQUAL (vm[ikey], ival);
	CPPUNIT_ASSERT_EQUAL (ival, val);
      }
  }

public:

  // virtual void setUp()
  // {}

  // virtual void tearDown()
  // {}


  void testCreate()
  {
    create<int, int> ();
    create<int*,int> ();
    create<int*,int*>();
    create<int, std::vector<int> >();
  }

  void testInsert()
  {
    insert<int, int> ();
    insert<char,int> ();
    insert<long,int*>();
    insert<int, std::vector<int> >();
  }

  void testIterate()
  {
    iterate<int, int> ();
    iterate<char,int> ();
    iterate<long,int*>();
    iterate<int, std::string>("test_string");
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION ( VectormapTest );
