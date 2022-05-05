#include "libmesh/vectormap.h"

#include "libmesh_cppunit.h"


#define VECTORMAPOBJECTTEST                     \
  CPPUNIT_TEST( testCreate );                   \

using namespace libMesh;

class VectormapTest : public CppUnit::TestCase
{
public:
  LIBMESH_CPPUNIT_TEST_SUITE ( VectormapTest );

  CPPUNIT_TEST( testCreate );
  CPPUNIT_TEST( testInsert );
  CPPUNIT_TEST( testEmplace );
  CPPUNIT_TEST( testIterate );
  CPPUNIT_TEST( testFind );

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
  void emplace()
  {
    vectormap<Key,Val> vm;

    Val val(0); // requires default constructor for val type.

    for (Key key=1; key<32; key*=2)
      vm.emplace(key, val);

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
    LOG_UNIT_TEST;

    create<int, int> ();
    create<int*,int> ();
    create<int*,int*>();
    create<int, std::vector<int>>();
  }

  void testInsert()
  {
    LOG_UNIT_TEST;

    insert<int, int> ();
    insert<char,int> ();
    insert<long,int*>();
    insert<int, std::vector<int>>();
  }

  void testEmplace()
  {
    LOG_UNIT_TEST;

    emplace<int, int> ();
    emplace<char,int> ();
    emplace<long,int*>();
    emplace<int, std::vector<int>>();
  }

  void testIterate()
  {
    LOG_UNIT_TEST;

    iterate<int, int> ();
    iterate<char,int> ();
    iterate<long,int*>();
    iterate<int, std::string>("test_string");
  }

  void testFind()
  {
    LOG_UNIT_TEST;

    vectormap<int, int> vm;
    for (int i=16; i<32; ++i)
      vm.insert(std::make_pair(i,i));

    vectormap<int, int>::iterator
      it1 = vm.find(24),
      it2 = vm.find(4);

    CPPUNIT_ASSERT(it1 != vm.end());
    CPPUNIT_ASSERT(it2 == vm.end());
    CPPUNIT_ASSERT(vm.count(24) == 1);
    CPPUNIT_ASSERT(vm.count(4) == 0);
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION ( VectormapTest );
