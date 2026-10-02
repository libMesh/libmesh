
#include "fe_test.h"

#include <libmesh/fe_type.h>
#include <libmesh/node.h>

INSTANTIATE_FETEST(FIRST, HIERARCHIC, EDGE2);
INSTANTIATE_FETEST(SECOND, HIERARCHIC, EDGE3);
INSTANTIATE_FETEST(THIRD, HIERARCHIC, EDGE3);
INSTANTIATE_FETEST(FOURTH, HIERARCHIC, EDGE3);

#if LIBMESH_DIM > 1
INSTANTIATE_FETEST(FIRST, HIERARCHIC, TRI3);
INSTANTIATE_FETEST(SECOND, HIERARCHIC, TRI6);
INSTANTIATE_FETEST(THIRD, HIERARCHIC, TRI6);
INSTANTIATE_FETEST(FOURTH, HIERARCHIC, TRI6);

INSTANTIATE_FETEST(SECOND, HIERARCHIC, TRI7);
INSTANTIATE_FETEST(THIRD, HIERARCHIC, TRI7);
INSTANTIATE_FETEST(FOURTH, HIERARCHIC, TRI7);

INSTANTIATE_FETEST(FIRST, HIERARCHIC, QUAD4);
INSTANTIATE_FETEST(SECOND, HIERARCHIC, QUAD9);
INSTANTIATE_FETEST(THIRD, HIERARCHIC, QUAD9);
INSTANTIATE_FETEST(FOURTH, HIERARCHIC, QUAD9);
#endif

#if LIBMESH_DIM > 2
INSTANTIATE_FETEST(FIRST, HIERARCHIC, HEX8);
INSTANTIATE_FETEST(SECOND, HIERARCHIC, HEX27);
INSTANTIATE_FETEST(THIRD, HIERARCHIC, HEX27);
INSTANTIATE_FETEST(FOURTH, HIERARCHIC, HEX27);

INSTANTIATE_FETEST(FIRST, HIERARCHIC, PRISM6);
INSTANTIATE_FETEST(FIRST, HIERARCHIC, PRISM15);
INSTANTIATE_FETEST(SECOND, HIERARCHIC, PRISM18);
INSTANTIATE_FETEST(SECOND, HIERARCHIC, PRISM20);
INSTANTIATE_FETEST(THIRD, HIERARCHIC, PRISM20);
INSTANTIATE_FETEST(FOURTH, HIERARCHIC, PRISM20);
INSTANTIATE_FETEST(SECOND, HIERARCHIC, PRISM21);
INSTANTIATE_FETEST(THIRD, HIERARCHIC, PRISM21);
INSTANTIATE_FETEST(FOURTH, HIERARCHIC, PRISM21);

INSTANTIATE_FETEST(FIRST, HIERARCHIC, TET4);
INSTANTIATE_FETEST(SECOND, HIERARCHIC, TET10);
INSTANTIATE_FETEST(THIRD, HIERARCHIC, TET14);
INSTANTIATE_FETEST(FOURTH, HIERARCHIC, TET14);
#endif

#if LIBMESH_DIM > 1

/**
 * Checks that the HIERARCHIC simplex edge functions are continuous where the barycentric
 * coordinates of an edge's two vertices sum to zero. There the general expression is singular
 * and the shape function is evaluated from its limit, which must carry the same orientation sign
 * as the general expression.
 */
class HierarchicSimplexEdgeTest : public CppUnit::TestCase
{
public:
  LIBMESH_CPPUNIT_TEST_SUITE( HierarchicSimplexEdgeTest );
  CPPUNIT_TEST( testTriEdgeLimit );
#if LIBMESH_DIM > 2
  CPPUNIT_TEST( testTetEdgeLimit );
#endif
  CPPUNIT_TEST_SUITE_END();

  void testTriEdgeLimit()
  {
    LOG_UNIT_TEST;
    checkEdgeLimit(TRI6);
  }

  void testTetEdgeLimit()
  {
    LOG_UNIT_TEST;
    checkEdgeLimit(TET14);
  }

private:
  void checkEdgeLimit(const ElemType elem_type)
  {
    const auto elem = Elem::build(elem_type);
    std::vector<std::unique_ptr<Node>> nodes;
    for (const auto n : elem->node_index_range())
      {
        nodes.push_back(std::make_unique<Node>(elem->master_point(n), n));
        elem->set_node(n, nodes.back().get());
      }

    // Total order four covers edge functions of both parities
    const unsigned int totalorder = 4;
    const FEType fe_type(totalorder, HIERARCHIC);
    const unsigned int n_vertices = elem->n_vertices();

    // Distance from the singular locus at which the general expression is evaluated, and a
    // tolerance comfortably above the O(delta) difference between it and the limit
    const Real delta = 1e-6;
    const Real tol = 1e-4;

    bool checked_flipped_edge = false;

    for (const auto e : make_range(elem->n_edges()))
      {
        const unsigned int v0 = elem->local_edge_node(e, 0),
                           v1 = elem->local_edge_node(e, 1);
        checked_flipped_edge |= elem->positive_edge_orientation(e);

        // Barycentric coordinates with zeta[v0] + zeta[v1] equal to c and
        // zeta[v1] - zeta[v0] equal to one; the other vertices share the remainder
        auto edge_point = [&](const Real c)
          {
            std::vector<Real> zeta(n_vertices, (1. - c) / (n_vertices - 2));
            zeta[v0] = (c - 1.) / 2.;
            zeta[v1] = (c + 1.) / 2.;

            Point p;
            for (const auto v : make_range(n_vertices))
              p += zeta[v] * elem->master_point(v);
            return p;
          };

        const Point p_limit = edge_point(0.), p_near = edge_point(delta);

        for (const auto basisorder : make_range(2u, totalorder + 1))
          {
            const unsigned int i = n_vertices + e * (totalorder - 1) + (basisorder - 2);
            const Real limit = FEInterface::shape(fe_type, elem.get(), i, p_limit);
            const Real near = FEInterface::shape(fe_type, elem.get(), i, p_near);
            LIBMESH_ASSERT_FP_EQUAL(near, limit, tol);
          }
      }

    // The reference vertex placement orients some edges positively, so the odd edge functions
    // there carry a negative sign
    CPPUNIT_ASSERT(checked_flipped_edge);
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION( HierarchicSimplexEdgeTest );

#endif // LIBMESH_DIM > 1
