#include <libmesh/elem_range.h>
#include <libmesh/replicated_mesh.h>
#include <libmesh/mesh_generation.h>

#include "test_comm.h"
#include "libmesh_cppunit.h"

using namespace libMesh;

class StoredRangeTest : public CppUnit::TestCase {

public:
  LIBMESH_CPPUNIT_TEST_SUITE( StoredRangeTest );
#if LIBMESH_DIM > 1
  CPPUNIT_TEST( testEquality );
#endif
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp() {}

  void tearDown() {}

  /**
   * Ranges are equal when they hold the same objects in the same order. A range compares equal to
   * a range over a copy of its objects, wherever the two are stored, including when one range
   * begins where the other ends.
   */
  void testEquality()
  {
    LOG_UNIT_TEST;

    // Every processor needs every element
    ReplicatedMesh mesh(*TestCommWorld);
    MeshTools::Generation::build_square(mesh, 2, 2);

    std::vector<const Elem *> all(mesh.elements_begin(), mesh.elements_end());
    std::vector<const Elem *> copy(all);
    std::vector<const Elem *> fewer(all.begin(), all.end() - 1);
    std::vector<const Elem *> reversed(all.rbegin(), all.rend());

    const ConstElemRange all_range(&all), copy_range(&copy), fewer_range(&fewer),
      reversed_range(&reversed);

    CPPUNIT_ASSERT(all_range == copy_range);
    CPPUNIT_ASSERT(!(all_range == fewer_range));
    CPPUNIT_ASSERT(!(fewer_range == all_range));
    CPPUNIT_ASSERT(!(all_range == reversed_range));

    // Two copies of the same objects, stored back to back
    std::vector<const Elem *> twice(all);
    twice.insert(twice.end(), all.begin(), all.end());
    const ConstElemRange twice_range(&twice);
    const auto middle = twice_range.begin() + all.size();
    const ConstElemRange first_half(twice_range, twice_range.begin(), middle),
      second_half(twice_range, middle, twice_range.end());

    CPPUNIT_ASSERT(first_half == second_half);
    CPPUNIT_ASSERT(second_half == first_half);
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION( StoredRangeTest );
