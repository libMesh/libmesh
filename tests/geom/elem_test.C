#include "test_comm.h"

#include <libmesh/elem.h>
#include <libmesh/enum_elem_type.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/elem_side_builder.h>
#include <libmesh/auto_ptr.h>

#include "libmesh_cppunit.h"

using namespace libMesh;

template <ElemType elem_type>
class ElemTest : public CppUnit::TestCase {

private:
  std::unique_ptr<Mesh> _mesh;

protected:
  std::string libmesh_suite_name;

public:
  void setUp()
  {
    const Real minpos = 1.5, maxpos = 5.5;
    const unsigned int N = 2;

    _mesh = libmesh_make_unique<Mesh>(*TestCommWorld);
    std::unique_ptr<Elem> test_elem = Elem::build(elem_type);

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
#if LIBMESH_DIM > 1
    if (test_elem->infinite())
      {
        Elem * elem = _mesh->add_elem(std::move(test_elem));

        const auto add_point =
          [this, elem](const unsigned int i,
                       const Real x,
                       const Real y,
                       const Real
#if LIBMESH_DIM == 3
                                  z
#endif
                      )
          {
#if LIBMESH_DIM == 2
            auto node = _mesh->add_point(Point(x, y), i);
#else
            auto node = _mesh->add_point(Point(x, y, z), i);
#endif
            elem->set_node(i) = node;
          };

        const Real halfpos = (minpos + maxpos) / 2.;

        if (elem_type == INFQUAD4 || elem_type == INFQUAD6 ||
            elem_type == INFHEX8 || elem_type == INFHEX16 || elem_type == INFHEX18)
          {
            add_point(0, minpos, minpos, minpos);
            add_point(1, maxpos, minpos, minpos);
            add_point(2, minpos, maxpos, minpos);
            add_point(3, maxpos, maxpos, minpos);

            if (elem_type == INFQUAD6)
              {
                add_point(4, halfpos, minpos, minpos);
                add_point(5, halfpos, maxpos, minpos);
              }
          }
        if (elem_type == INFHEX8 || elem_type == INFHEX16 || elem_type == INFHEX18)
          {
            add_point(4, minpos, minpos, maxpos);
            add_point(5, maxpos, minpos, maxpos);
            add_point(6, minpos, maxpos, maxpos);
            add_point(7, maxpos, maxpos, maxpos);

            if (elem_type == INFHEX16 || elem_type == INFHEX18)
              {
                add_point(8, halfpos, minpos, minpos);
                add_point(9, maxpos, halfpos, minpos);
                add_point(10, halfpos, maxpos, minpos);
                add_point(11, minpos, halfpos, minpos);
                add_point(12, halfpos, minpos, maxpos);
                add_point(13, maxpos, halfpos, maxpos);
                add_point(14, halfpos, maxpos, maxpos);
                add_point(15, minpos, halfpos, maxpos);
              }
            if (elem_type == INFHEX18)
              {
                add_point(16, halfpos, halfpos, minpos);
                add_point(17, halfpos, halfpos, maxpos);
              }
          }
        if (elem_type == INFPRISM6 || elem_type == INFPRISM12)
          {
            add_point(0, minpos, minpos, minpos);
            add_point(1, maxpos, minpos, minpos);
            add_point(2, halfpos, maxpos, minpos);
            add_point(3, minpos, minpos, maxpos);
            add_point(4, maxpos, minpos, maxpos);
            add_point(5, halfpos, maxpos, maxpos);

            if (elem_type == INFPRISM12)
              {
                add_point(6, halfpos, minpos, minpos);
                add_point(7, (halfpos + maxpos) / 2., halfpos, minpos);
                add_point(8, (halfpos + minpos) / 2., halfpos, minpos);
                add_point(9, halfpos, minpos, maxpos);
                add_point(10, (halfpos + maxpos) / 2., halfpos, maxpos);
                add_point(11, (halfpos + minpos) / 2., halfpos, maxpos);
              }
          }

        _mesh->prepare_for_use();
      }
    else
#endif
#endif
      {
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
  }

  void test_bounding_box()
  {
    LOG_UNIT_TEST;

    for (const auto & elem : _mesh->active_local_element_ptr_range())
      {
        const BoundingBox bbox = elem->loose_bounding_box();

        // The "loose" bounding box should actually be pretty tight
        // in most of these cases, but for weirdly aligned triangles
        // (such as occur in pyramid elements) it won't be, so we'll
        // just test against a widened bounding box.
        BoundingBox wide_bbox(elem->point(0), elem->point(0));

        for (unsigned int n = 0; n != elem->n_nodes(); ++n)
          {
            const Point & p = elem->point(n);

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
            if (!elem->infinite())
#endif
              {
                CPPUNIT_ASSERT(bbox.contains_point(p));
              }

            wide_bbox.union_with
              (BoundingBox(elem->point(n), elem->point(n)));
          }

        wide_bbox.scale(1. / 3.);

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
        if (!elem->infinite())
#endif
          {
            CPPUNIT_ASSERT(!bbox.contains_point(wide_bbox.min()));
            CPPUNIT_ASSERT(!bbox.contains_point(wide_bbox.max()));
          }
      }
  }

  void test_maps()
  {
    LOG_UNIT_TEST;

    for (const auto & elem : _mesh->active_local_element_ptr_range())
      {
        for (const auto edge : elem->edge_index_range())
          for (const auto side_on_edge : elem->sides_on_edge(edge))
            for (const auto node_on_edge : elem->nodes_on_edge(edge))
              CPPUNIT_ASSERT(elem->is_node_on_side(node_on_edge, side_on_edge));

        for (const auto side : elem->side_index_range())
          for (const auto node_on_side : elem->nodes_on_side(side))
            CPPUNIT_ASSERT(elem->is_node_on_side(node_on_side, side));

        for (const auto edge : elem->edge_index_range())
          for (const auto node_on_edge : elem->nodes_on_edge(edge))
            CPPUNIT_ASSERT(elem->is_node_on_edge(node_on_edge, edge));

        for (const auto edge : elem->edge_index_range())
          for (const auto side_on_edge : elem->sides_on_edge(edge))
            CPPUNIT_ASSERT(elem->is_edge_on_side(edge, side_on_edge));
      }
  }

  void test_contains_point_node()
  {
    LOG_UNIT_TEST;

    for (const auto & elem : _mesh->active_local_element_ptr_range())
      {
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
        if (elem->infinite())
          continue;
#endif

        for (const auto n : elem->node_index_range())
#ifndef LIBMESH_ENABLE_EXCEPTIONS
          // If this node has a singular Jacobian, we need exceptions in order
          // to catch the failed inverse_map solve and return the singular
          // master point. Therefore, if we don't have exceptions and we're
          // at a singular node, we can't test this. As of the writing of
          // this comment, this issue exists for only Pyramid elements at
          // the apex.
          if (elem->local_singular_node(elem->point(n), TOLERANCE*TOLERANCE) == invalid_uint)
#endif
            CPPUNIT_ASSERT(elem->contains_point(elem->point(n)));
      }
  }

  void test_permute()
  {
    LOG_UNIT_TEST;

    for (const auto & elem : _mesh->active_local_element_ptr_range())
      {
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
        if (elem->infinite())
          continue;
#endif
        const Point centroid = elem->true_centroid();
        const Point vertex_avg = elem->vertex_average();
        Point quasicc;
        if (elem->dim() < 3)
          quasicc = elem->quasicircumcenter();

        for (const auto p : IntRange<unsigned int>(0, elem->n_permutations()))
          {
            elem->permute(p);
            CPPUNIT_ASSERT(elem->has_invertible_map());
            const Point new_centroid = elem->true_centroid();
            const Point new_vertex_avg = elem->vertex_average();
            Point new_quasicc;
            if (elem->dim() < 3)
              new_quasicc = elem->quasicircumcenter();
            for (const auto d : make_range(LIBMESH_DIM))
              {
                LIBMESH_ASSERT_FP_EQUAL(centroid(d), new_centroid(d),
                                        TOLERANCE*TOLERANCE);
                LIBMESH_ASSERT_FP_EQUAL(vertex_avg(d), new_vertex_avg(d),
                                        TOLERANCE*TOLERANCE);
                LIBMESH_ASSERT_FP_EQUAL(quasicc(d), new_quasicc(d),
                                        TOLERANCE*TOLERANCE);
              }
          }
      }
  }

  void test_center_node_on_side()
  {
    LOG_UNIT_TEST;

    for (const auto & elem : _mesh->active_local_element_ptr_range())
      for (const auto s : elem->side_index_range())
        {
          if (elem->type() == EDGE2 || elem->type() == EDGE3 || elem->type() == EDGE4)
            CPPUNIT_ASSERT_EQUAL((unsigned int)s, elem->center_node_on_side(s));
          else if (elem->type() == TRI6 || elem->type() == TRI7)
            CPPUNIT_ASSERT_EQUAL((unsigned int)(s + 3), elem->center_node_on_side(s));
          else if (elem->type() == QUAD8 || elem->type() == QUAD9 || elem->type() == QUADSHELL8)
            CPPUNIT_ASSERT_EQUAL((unsigned int)(s + 4), elem->center_node_on_side(s));
          else if (elem->type() == HEX27)
            CPPUNIT_ASSERT_EQUAL((unsigned int)(s + 20), elem->center_node_on_side(s));
          else if (elem->type() == PRISM18 && s >= 1 && s <= 3)
            CPPUNIT_ASSERT_EQUAL((unsigned int)(s + 14), elem->center_node_on_side(s));
          else if (elem->type() == PYRAMID14 && s == 4)
            CPPUNIT_ASSERT_EQUAL((unsigned int)(13), elem->center_node_on_side(s));
          else
            CPPUNIT_ASSERT_EQUAL(invalid_uint, elem->center_node_on_side(s));
        }
  }

  void test_side_type()
  {
    LOG_UNIT_TEST;

    for (const auto & elem : _mesh->active_local_element_ptr_range())
      for (const auto s : elem->side_index_range())
        CPPUNIT_ASSERT_EQUAL(elem->build_side_ptr(s)->type(), elem->side_type(s));
  }

  void test_elem_side_builder()
  {
    LOG_UNIT_TEST;

    ElemSideBuilder cache;
    for (auto & elem : _mesh->active_local_element_ptr_range())
      for (const auto s : elem->side_index_range())
      {
        const auto side = elem->build_side_ptr(s);

        auto & cached_side = cache(*elem, s);
        CPPUNIT_ASSERT_EQUAL(side->type(), cached_side.type());
        for (const auto n : side->node_index_range())
          CPPUNIT_ASSERT_EQUAL(side->node_ref(n), cached_side.node_ref(n));

        const auto & const_cached_side = cache(const_cast<const Elem &>(*elem), s);
        CPPUNIT_ASSERT_EQUAL(side->type(), const_cached_side.type());
        for (const auto n : side->node_index_range())
          CPPUNIT_ASSERT_EQUAL(side->node_ref(n), const_cached_side.node_ref(n));
      }
  }
};

#define ELEMTEST                                \
  CPPUNIT_TEST( test_bounding_box );            \
  CPPUNIT_TEST( test_maps );                    \
  CPPUNIT_TEST( test_permute );                 \
  CPPUNIT_TEST( test_contains_point_node );     \
  CPPUNIT_TEST( test_center_node_on_side );     \
  CPPUNIT_TEST( test_side_type );               \
  CPPUNIT_TEST( test_elem_side_builder );

#define INSTANTIATE_ELEMTEST(elemtype)                          \
  class ElemTest_##elemtype : public ElemTest<elemtype> {       \
  public:                                                       \
  ElemTest_##elemtype() :                                       \
    ElemTest<elemtype>() {                                      \
    this->libmesh_suite_name = "ElemTest_" #elemtype;           \
  }                                                             \
  CPPUNIT_TEST_SUITE( ElemTest_##elemtype );                    \
  ELEMTEST;                                                     \
  CPPUNIT_TEST_SUITE_END();                                     \
  };                                                            \
                                                                \
  CPPUNIT_TEST_SUITE_REGISTRATION( ElemTest_##elemtype )

INSTANTIATE_ELEMTEST(EDGE2);
INSTANTIATE_ELEMTEST(EDGE3);
INSTANTIATE_ELEMTEST(EDGE4);

#if LIBMESH_DIM > 1
INSTANTIATE_ELEMTEST(TRI3);
INSTANTIATE_ELEMTEST(TRI6);
INSTANTIATE_ELEMTEST(TRI7);

INSTANTIATE_ELEMTEST(QUAD4);
INSTANTIATE_ELEMTEST(QUAD8);
INSTANTIATE_ELEMTEST(QUAD9);

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
INSTANTIATE_ELEMTEST(INFQUAD4);
INSTANTIATE_ELEMTEST(INFQUAD6);
#endif
#endif // LIBMESH_DIM > 1

#if LIBMESH_DIM > 2
INSTANTIATE_ELEMTEST(TET4);
INSTANTIATE_ELEMTEST(TET10);
INSTANTIATE_ELEMTEST(TET14);

INSTANTIATE_ELEMTEST(HEX8);
INSTANTIATE_ELEMTEST(HEX20);
INSTANTIATE_ELEMTEST(HEX27);

INSTANTIATE_ELEMTEST(PRISM6);
INSTANTIATE_ELEMTEST(PRISM15);
INSTANTIATE_ELEMTEST(PRISM18);

INSTANTIATE_ELEMTEST(PYRAMID5);
INSTANTIATE_ELEMTEST(PYRAMID13);
INSTANTIATE_ELEMTEST(PYRAMID14);

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
INSTANTIATE_ELEMTEST(INFHEX8);
INSTANTIATE_ELEMTEST(INFHEX16);
INSTANTIATE_ELEMTEST(INFHEX18);

INSTANTIATE_ELEMTEST(INFPRISM6);
INSTANTIATE_ELEMTEST(INFPRISM12);
#endif
#endif // LIBMESH_DIM > 2
