#include "test_comm.h"

#include <libmesh/elem.h>
#include <libmesh/enum_elem_type.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>

// THE CPPUNIT_TEST_SUITE_END macro expands to code that involves
// std::auto_ptr, which in turn produces -Wdeprecated-declarations
// warnings.  These can be ignored in GCC as long as we wrap the
// offending code in appropriate pragmas.  We can't get away with a
// single ignore_warnings.h inclusion at the beginning of this file,
// since the libmesh headers pull in a restore_warnings.h at some
// point.  We also don't bother restoring warnings at the end of this
// file since it's not a header.
#include <libmesh/ignore_warnings.h>

#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>

using namespace libMesh;

template <ElemType elem_type>
class ElemTest : public CppUnit::TestCase {

private:
  Mesh * _mesh;

public:
  void setUp()
  {
    const Real minpos = 1.5, maxpos = 5.5;
    const unsigned int N = 2;

    _mesh = new Mesh(*TestCommWorld);
    const UniquePtr<Elem> test_elem = Elem::build(elem_type);
    const unsigned int dim = test_elem->dim();
    const unsigned int use_y = dim > 1;
    const unsigned int use_z = dim > 2;

    MeshTools::Generation::build_cube (*_mesh,
                                       N, N*use_y, N*use_z,
                                       minpos, maxpos,
                                       minpos, use_y*maxpos,
                                       minpos, use_z*maxpos,
                                       elem_type);
  }

  void tearDown()
  {
    delete _mesh;
  }

  void test_bounding_box()
  {
    MeshBase::const_element_iterator
      elem_it  = _mesh->active_local_elements_begin(),
      elem_end = _mesh->active_local_elements_end();
    for (; elem_it != elem_end; ++elem_it)
      {
        const Elem & elem = **elem_it;

        const BoundingBox bbox = elem.loose_bounding_box();

        const Point centroid = elem.centroid();

        // The "loose" bounding box should actually be pretty tight
        // in most of these cases, but for weirdly aligned triangles
        // (such as occur in pyramid elements) it won't be, so we'll
        // just test against a widened bounding box.
        BoundingBox wide_bbox(elem.point(0), elem.point(0));

        for (unsigned int n = 0; n != elem.n_nodes(); ++n)
          {
            const Point & p = elem.point(n);

            CPPUNIT_ASSERT(bbox.contains_point(p));

            wide_bbox.union_with
              (BoundingBox(elem.point(n), elem.point(n)));
          }

        for (unsigned int d=0; d != LIBMESH_DIM; ++d)
          {
            const Real widening =
              (wide_bbox.max()(d) - wide_bbox.min()(d)) / 3;
            wide_bbox.min()(d) -= widening;
            wide_bbox.max()(d) += widening;
          }

        CPPUNIT_ASSERT(!bbox.contains_point(wide_bbox.min()));
        CPPUNIT_ASSERT(!bbox.contains_point(wide_bbox.max()));
      }
  }
};

#define ELEMTEST                                \
  CPPUNIT_TEST( test_bounding_box )

#define INSTANTIATE_ELEMTEST(elemtype)                          \
  class ElemTest_##elemtype : public ElemTest<elemtype> {       \
  public:                                                       \
  CPPUNIT_TEST_SUITE( ElemTest_##elemtype );                    \
  ELEMTEST;                                                     \
  CPPUNIT_TEST_SUITE_END();                                     \
  };                                                            \
                                                                \
  CPPUNIT_TEST_SUITE_REGISTRATION( ElemTest_##elemtype )

INSTANTIATE_ELEMTEST(EDGE2);
INSTANTIATE_ELEMTEST(EDGE3);
INSTANTIATE_ELEMTEST(EDGE4);

INSTANTIATE_ELEMTEST(TRI3);
INSTANTIATE_ELEMTEST(TRI6);

INSTANTIATE_ELEMTEST(QUAD4);
INSTANTIATE_ELEMTEST(QUAD8);
INSTANTIATE_ELEMTEST(QUAD9);

INSTANTIATE_ELEMTEST(TET4);
INSTANTIATE_ELEMTEST(TET10);

INSTANTIATE_ELEMTEST(HEX8);
INSTANTIATE_ELEMTEST(HEX20);
INSTANTIATE_ELEMTEST(HEX27);

INSTANTIATE_ELEMTEST(PRISM6);
INSTANTIATE_ELEMTEST(PRISM15);
INSTANTIATE_ELEMTEST(PRISM18);

INSTANTIATE_ELEMTEST(PYRAMID5);
INSTANTIATE_ELEMTEST(PYRAMID13);
INSTANTIATE_ELEMTEST(PYRAMID14);
