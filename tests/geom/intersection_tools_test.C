#include "libmesh_cppunit.h"

#include "test_comm.h"

#include <libmesh/intersection_tools.h>
#include <libmesh/point.h>
#include <libmesh/int_range.h>
#include <libmesh/enum_elem_type.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/elem.h>
#include <libmesh/elem_corner.h>

using namespace libMesh;

class IntersectionToolsTest : public CppUnit::TestCase
{

public:
  LIBMESH_CPPUNIT_TEST_SUITE( IntersectionToolsTest );
  CPPUNIT_TEST( within_segment );
  CPPUNIT_TEST( collinear );
  CPPUNIT_TEST_SUITE_END();

public:

  void within_segment()
  {
    LOG_UNIT_TEST;

    const Point s1(1.0, 2.0, 3.0);
    const Point s2(2.0, 3.0, 4.0);
    const auto length_vec = s2 - s1;
    const auto length = length_vec.norm();
    const auto s1_to_s2 = length_vec / length;

    int segments = 10;
    Real dx = (Real)1 / segments * length;
    for (const auto i : make_range(-1, segments + 1))
    {
      const auto p = s1 + Real(i) * dx * s1_to_s2;
      IntersectionTools::WithinSegmentResult within_result = IntersectionTools::WithinSegmentResult::NOT_WITHIN;
      if (i == 0)
        within_result = IntersectionTools::WithinSegmentResult::AT_BEGINNING;
      else if (i > 0 && i < segments)
        within_result = IntersectionTools::WithinSegmentResult::BETWEEN;
      else if (i == segments)
        within_result = IntersectionTools::WithinSegmentResult::AT_END;

      CPPUNIT_ASSERT_EQUAL(IntersectionTools::within_segment(s1, s2, p), within_result);
    }

    CPPUNIT_ASSERT_EQUAL(IntersectionTools::within_segment(s1, s2, Point(9.9, 5, 3)),
                         IntersectionTools::WithinSegmentResult::NOT_WITHIN);
  }

  void collinear()
  {
    LOG_UNIT_TEST;

    const auto do_assert = [](const Point & p1,
                              const Point & p2,
                              const Point & p3,
                              const bool collinear)
    {
      CPPUNIT_ASSERT_EQUAL(IntersectionTools::collinear(p1, p2, p3), collinear);
      CPPUNIT_ASSERT_EQUAL(IntersectionTools::collinear(p1, p3, p2), collinear);
      CPPUNIT_ASSERT_EQUAL(IntersectionTools::collinear(p2, p1, p3), collinear);
      CPPUNIT_ASSERT_EQUAL(IntersectionTools::collinear(p2, p3, p1), collinear);
      CPPUNIT_ASSERT_EQUAL(IntersectionTools::collinear(p3, p1, p2), collinear);
      CPPUNIT_ASSERT_EQUAL(IntersectionTools::collinear(p3, p2, p1), collinear);
    };

    // two of the same points
    do_assert(Point(1, 2, 3), Point(1, 2, 3), Point(4, 5, 6), true);

    // three of the same points
    do_assert(Point(1, 2, 3), Point(1, 2, 3), Point(1, 2, 3), true);

    // all in a line
    do_assert(Point(1, 1, 2), Point(1, 1, 3), Point(1, 1, 0), true);

    // not collinear
    do_assert(Point(0, 1, 2), Point(0, 1, 3), Point(1, 5, 10), false);
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION( IntersectionToolsTest );

template <ElemType elem_type>
class MeshedIntersectionToolsTest : public CppUnit::TestCase {

private:
  std::unique_ptr<Mesh> _mesh;

protected:
  std::string libmesh_suite_name;

public:
  void setUp()
  {
    const Real minpos = 1.5, maxpos = 4.86;
    const unsigned int N = 3;

    _mesh = std::make_unique<Mesh>(*TestCommWorld);
    auto test_elem = Elem::build(elem_type);

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

  void test_within_edge_on_side()
  {
    LOG_UNIT_TEST;

    if (_mesh->mesh_dimension() != 3)
      return;

    // check locations at every node
    for (const auto & elem : _mesh->active_local_element_ptr_range())
      for (const auto s : elem->side_index_range())
        for (const auto e : elem->edge_index_range())
          for (const auto n : elem->nodes_on_edge(e))
          {
            ElemCorner corner;
            const auto within = IntersectionTools::within_edge_on_side(*elem,
                                                                        elem->point(n),
                                                                        s,
                                                                        corner);

            CPPUNIT_ASSERT_EQUAL(within, elem->is_node_on_side(n, s));
            if (elem->is_node_on_side(n, s))
            {
              CPPUNIT_ASSERT_EQUAL(elem->is_vertex(n), corner.at_vertex(n));
              CPPUNIT_ASSERT_EQUAL(elem->is_vertex(n), !corner.at_edge(*elem, e));
            }
          }

    // cut edges into segments
    for (const auto & elem : _mesh->active_local_element_ptr_range())
      for (const auto e : elem->edge_index_range())
        for (const auto s : elem->side_index_range())
          if (elem->is_edge_on_side(e, s))
          {
            const auto nodes_on_edge = elem->nodes_on_edge(e);
            const auto & p1 = elem->point(nodes_on_edge[0]);
            const auto & p2 = elem->point(nodes_on_edge[1]);
            const auto length_vec = p2 - p1;
            const auto length = length_vec.norm();
            const auto p1_to_p2 = length_vec / length;

            int segments = 5;
            Real dx = (Real)1 / segments * length;
            for (const auto i : make_range(-1, segments + 1))
            {
              const auto p = p1 + Real(i) * dx * p1_to_p2;
              ElemCorner corner;
              const auto within = IntersectionTools::within_edge_on_side(*elem,
                                                                         p,
                                                                         s,
                                                                         corner);

              CPPUNIT_ASSERT_EQUAL(within, i >= 0 && i <= segments);
              CPPUNIT_ASSERT_EQUAL(corner.at_vertex(nodes_on_edge[0]), i == 0);
              CPPUNIT_ASSERT_EQUAL(corner.at_vertex(nodes_on_edge[1]), i == segments);
              CPPUNIT_ASSERT_EQUAL(corner.at_edge(*elem, e), i > 0 && i < segments);
            }
          }

    // check elem centroids
    for (const auto & elem : _mesh->active_local_element_ptr_range())
      for (const auto s : elem->side_index_range())
      {
        ElemCorner corner;
        CPPUNIT_ASSERT(!IntersectionTools::within_edge_on_side(*elem,
                                                               elem->vertex_average(),
                                                               s,
                                                               corner));
      }
  }

};

#define MESHEDINTERSECTIONTOOLSTEST               \
  CPPUNIT_TEST( test_within_edge_on_side );

#define INSTANTIATE_MESHEDINTERSECTIONTOOLSTEST(elemtype)                                       \
  class MeshedIntersectionToolsTest_##elemtype : public MeshedIntersectionToolsTest<elemtype> { \
  public:                                                                                       \
  MeshedIntersectionToolsTest_##elemtype() :                                                    \
    MeshedIntersectionToolsTest<elemtype>() {                                                   \
    if (unitlog->summarized_logs_enabled())                                                     \
      this->libmesh_suite_name = "MeshedIntersectionToolsTest";                                 \
    else                                                                                        \
      this->libmesh_suite_name = "MeshedIntersectionToolsTest_" #elemtype;                      \
  }                                                                                             \
  CPPUNIT_TEST_SUITE( MeshedIntersectionToolsTest_##elemtype );                                 \
  MESHEDINTERSECTIONTOOLSTEST;                                                                  \
  CPPUNIT_TEST_SUITE_END();                                                                     \
  };                                                                                            \
                                                                                                \
  CPPUNIT_TEST_SUITE_REGISTRATION( MeshedIntersectionToolsTest_##elemtype )


#if LIBMESH_DIM > 2
INSTANTIATE_MESHEDINTERSECTIONTOOLSTEST(TET4);
INSTANTIATE_MESHEDINTERSECTIONTOOLSTEST(TET10);
INSTANTIATE_MESHEDINTERSECTIONTOOLSTEST(TET14);

INSTANTIATE_MESHEDINTERSECTIONTOOLSTEST(HEX8);
INSTANTIATE_MESHEDINTERSECTIONTOOLSTEST(HEX20);
INSTANTIATE_MESHEDINTERSECTIONTOOLSTEST(HEX27);

INSTANTIATE_MESHEDINTERSECTIONTOOLSTEST(PRISM6);
INSTANTIATE_MESHEDINTERSECTIONTOOLSTEST(PRISM15);
INSTANTIATE_MESHEDINTERSECTIONTOOLSTEST(PRISM18);
INSTANTIATE_MESHEDINTERSECTIONTOOLSTEST(PRISM20);
INSTANTIATE_MESHEDINTERSECTIONTOOLSTEST(PRISM21);

INSTANTIATE_MESHEDINTERSECTIONTOOLSTEST(PYRAMID5);
INSTANTIATE_MESHEDINTERSECTIONTOOLSTEST(PYRAMID13);
INSTANTIATE_MESHEDINTERSECTIONTOOLSTEST(PYRAMID14);
#endif // LIBMESH_DIM > 2
