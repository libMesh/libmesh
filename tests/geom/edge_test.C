#include <libmesh/elem.h>

#include <libmesh/cell_hex20.h>
#include <libmesh/cell_hex27.h>
#include <libmesh/cell_hex8.h>
#include <libmesh/cell_inf_hex16.h>
#include <libmesh/cell_inf_hex18.h>
#include <libmesh/cell_inf_hex8.h>
#include <libmesh/cell_inf_prism12.h>
#include <libmesh/cell_inf_prism6.h>
#include <libmesh/cell_prism15.h>
#include <libmesh/cell_prism18.h>
#include <libmesh/cell_prism6.h>
#include <libmesh/cell_pyramid13.h>
#include <libmesh/cell_pyramid14.h>
#include <libmesh/cell_pyramid5.h>
#include <libmesh/cell_tet10.h>
#include <libmesh/cell_tet14.h>
#include <libmesh/cell_tet4.h>
#include <libmesh/edge_edge2.h>
#include <libmesh/edge_edge3.h>
#include <libmesh/edge_edge4.h>
#include <libmesh/edge_inf_edge2.h>
#include <libmesh/face_inf_quad4.h>
#include <libmesh/face_inf_quad6.h>
#include <libmesh/face_quad4.h>
#include <libmesh/face_quad8.h>
#include <libmesh/face_quad9.h>
#include <libmesh/face_tri3.h>
#include <libmesh/face_tri6.h>
#include <libmesh/face_tri7.h>

#include <vector>

#include "libmesh_cppunit.h"

#define EDGETEST                                \
  CPPUNIT_TEST( testIsNodeOnEdge );             \
  CPPUNIT_TEST( testNodesOnEdge );              \
  CPPUNIT_TEST( testBuildEdgePtr );             \

using namespace libMesh;


template <typename ElemClass, ElemType edge_type,
          unsigned short indexbegin, unsigned short indexend>
class EdgeTest : public CppUnit::TestCase {

private:
  ElemClass elem;
  std::vector<std::unique_ptr<Node>> nodes;

protected:
  std::string libmesh_suite_name;

public:
  void setUp() {
    elem.set_id() = 0;
#ifdef LIBMESH_ENABLE_AMR
    // Do tests with an Elem having a non-default p_level to ensure
    // that edges which are built have a matching p_level. p-refinement
    // is only avaiable if LIBMESH_ENABLE_AMR is defined.
    elem.set_p_level(1);
#endif
    Point dummy;
    for (auto i : elem.node_index_range())
      {
        nodes.push_back(libmesh_make_unique<Node>(dummy, /*id=*/i));
        elem.set_node(i) = nodes[i].get();
      }
  }

  void tearDown() {}

  void testIsNodeOnEdge()
  {
    LOG_UNIT_TEST;

    for (auto e : make_range(indexbegin, indexend))
      {
        std::unique_ptr<Elem> edge = elem.build_edge_ptr(e);
        for (auto n : elem.node_index_range())
          {
            const Node * node = elem.node_ptr(n);
            bool found_node = false;
            for (auto en : edge->node_index_range())
              if (node == edge->node_ptr(en))
                {
                  found_node = true;
                  break;
                }

            if (elem.is_node_on_edge(n, e))
              {
                CPPUNIT_ASSERT(found_node);
              }
            else
              {
                CPPUNIT_ASSERT(!found_node);
              }
          }
      }
  }

  void testNodesOnEdge()
  {
    LOG_UNIT_TEST;

    for (auto e : make_range(indexbegin, indexend))
      {
        std::unique_ptr<Elem> edge = elem.build_edge_ptr(e);
        std::vector<unsigned int> edge_nodes = elem.nodes_on_edge(e);

        CPPUNIT_ASSERT_EQUAL(edge_nodes.size(), std::size_t(edge->n_nodes()));

        for (auto en : edge->node_index_range())
          {
            const Node * node = edge->node_ptr(en);
            bool found_node = false;
            for (auto ei : edge_nodes)
              if (node == elem.node_ptr(ei))
                {
                  found_node = true;
                  break;
                }
            CPPUNIT_ASSERT(found_node);
          }
      }
  }

  void testBuildEdgePtr()
  {
    LOG_UNIT_TEST;

    for (auto e : make_range(indexbegin, indexend))
      {
        std::unique_ptr<Elem> edge = elem.build_edge_ptr(e);

        CPPUNIT_ASSERT(edge->type() == edge_type);
        CPPUNIT_ASSERT(edge->subdomain_id() == elem.subdomain_id());

#ifdef LIBMESH_ENABLE_AMR
        // p-refinement is only avaiable if LIBMESH_ENABLE_AMR is defined.
        CPPUNIT_ASSERT(edge->p_level() == elem.p_level());
#endif
      }
  }
};


#define INSTANTIATE_EDGETEST(elemclass, edgetype, indexbegin, indexend)                \
  class EdgeTest_##elemclass##_##edgetype##_##indexbegin##_##indexend :                \
    public EdgeTest<elemclass, edgetype, indexbegin, indexend> {                       \
  public:                                                                              \
  EdgeTest_##elemclass##_##edgetype##_##indexbegin##_##indexend() :                    \
    EdgeTest<elemclass,edgetype,indexbegin,indexend>() {                               \
    if (unitlog->summarized_logs_enabled())                                            \
      this->libmesh_suite_name = "EdgeTest";                                           \
    else                                                                               \
      this->libmesh_suite_name = "EdgeTest_" #elemclass"_" #edgetype "_" #indexbegin "_" #indexend; \
  }                                                                                    \
  CPPUNIT_TEST_SUITE( EdgeTest_##elemclass##_##edgetype##_##indexbegin##_##indexend ); \
  EDGETEST                                                                             \
  CPPUNIT_TEST_SUITE_END();                                                            \
  };                                                                                   \
                                                                                       \
  CPPUNIT_TEST_SUITE_REGISTRATION( EdgeTest_##elemclass##_##edgetype##_##indexbegin##_##indexend );

INSTANTIATE_EDGETEST(Hex20,     EDGE3, 0, 12);
INSTANTIATE_EDGETEST(Hex27,     EDGE3, 0, 12);
INSTANTIATE_EDGETEST(Hex8,      EDGE2, 0, 12);
INSTANTIATE_EDGETEST(Prism15,   EDGE3, 0, 9);
INSTANTIATE_EDGETEST(Prism18,   EDGE3, 0, 9);
INSTANTIATE_EDGETEST(Prism6,    EDGE2, 0, 9);
INSTANTIATE_EDGETEST(Pyramid13, EDGE3, 0, 8);
INSTANTIATE_EDGETEST(Pyramid14, EDGE3, 0, 8);
INSTANTIATE_EDGETEST(Pyramid5,  EDGE2, 0, 8);
INSTANTIATE_EDGETEST(Tet10,     EDGE3, 0, 6);
INSTANTIATE_EDGETEST(Tet14,     EDGE3, 0, 6);
INSTANTIATE_EDGETEST(Tet4,      EDGE2, 0, 6);
INSTANTIATE_EDGETEST(Quad4,     EDGE2, 0, 4);
INSTANTIATE_EDGETEST(Quad8,     EDGE3, 0, 4);
INSTANTIATE_EDGETEST(Quad9,     EDGE3, 0, 4);
INSTANTIATE_EDGETEST(Tri3,      EDGE2, 0, 3);
INSTANTIATE_EDGETEST(Tri6,      EDGE3, 0, 3);
INSTANTIATE_EDGETEST(Tri7,      EDGE3, 0, 3);

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
INSTANTIATE_EDGETEST(InfHex16,   EDGE3,    0, 4);
INSTANTIATE_EDGETEST(InfHex16,   INFEDGE2, 4, 8);
INSTANTIATE_EDGETEST(InfHex18,   EDGE3,    0, 4);
INSTANTIATE_EDGETEST(InfHex18,   INFEDGE2, 4, 8);
INSTANTIATE_EDGETEST(InfHex8,    EDGE2,    0, 4);
INSTANTIATE_EDGETEST(InfHex8,    INFEDGE2, 4, 8);
INSTANTIATE_EDGETEST(InfPrism12, EDGE3,    0, 3);
INSTANTIATE_EDGETEST(InfPrism12, INFEDGE2, 3, 6);
INSTANTIATE_EDGETEST(InfPrism6,  EDGE2,    0, 3);
INSTANTIATE_EDGETEST(InfPrism6,  INFEDGE2, 3, 6);
INSTANTIATE_EDGETEST(InfQuad4,   EDGE2,    0, 1);
INSTANTIATE_EDGETEST(InfQuad4,   INFEDGE2, 1, 3);
INSTANTIATE_EDGETEST(InfQuad6,   EDGE3,    0, 1);
INSTANTIATE_EDGETEST(InfQuad6,   INFEDGE2, 1, 3);
#endif // LIBMESH_ENABLE_INFINITE_ELEMENTS
