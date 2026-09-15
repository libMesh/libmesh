// At orders one and two on the tensor product element types, LAGRANGE and the Gauss-Lobatto
// nodal families have the same shape functions, because the Gauss-Lobatto points of those
// orders, {-1,1} and {-1,0,1}, are the element nodes. The implementations are independent:
// LAGRANGE evaluates a closed form for each order, and the Gauss-Lobatto basis evaluates
// the barycentric formula on its point table. This test pairs each Gauss-Lobatto shape
// function with the LAGRANGE shape function of the node at the same reference point and
// checks that their values agree.
//
// It catches errors in the point table, the barycentric evaluation, and the tensor index
// map.

#include <libmesh/elem.h>
#include <libmesh/enum_fe_family.h>
#include <libmesh/enum_order.h>
#include <libmesh/fe_interface.h>
#include <libmesh/fe_type.h>
#include <libmesh/int_range.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/parallel_implementation.h>
#include <libmesh/quadrature_gauss_lobatto.h>

#include "test_comm.h"
#include "libmesh_cppunit.h"

#include <memory>
#include <vector>

using namespace libMesh;

template <ElemType elem_type, Order order>
class FEGLLLagrangeAgreementTest : public CppUnit::TestCase
{
protected:
  std::unique_ptr<Mesh> _mesh;
  std::string libmesh_suite_name;

public:
  void setUp()
  {
    _mesh = std::make_unique<Mesh>(*TestCommWorld);

    const unsigned int dim = Elem::type_to_dim_map[elem_type];

    if (dim == 1)
      MeshTools::Generation::build_line(*_mesh, 1, 0., 1., elem_type);
    else if (dim == 2)
      MeshTools::Generation::build_square(*_mesh, 1, 1, 0., 1., 0., 1., elem_type);
    else
      MeshTools::Generation::build_cube(*_mesh, 1, 1, 1, 0., 1., 0., 1., 0., 1., elem_type);
  }

  void tearDown() { _mesh.reset(); }

  /// The interpolation point of Gauss-Lobatto shape function \p i
  Point gll_point(const unsigned int i, const unsigned int dim)
  {
    const unsigned int n = static_cast<unsigned int>(order) + 1;

    // QGaussLobatto chooses n points for order 2n-3
    QGaussLobatto rule(1, static_cast<Order>(2*n - 3));
    rule.init(EDGE2);

    Point q;
    unsigned int stride = 1;

    for (const auto d : make_range(dim))
      {
        q(d) = rule.qp((i / stride) % n)(0);
        stride *= n;
      }

    return q;
  }

  /// The node of \p elem whose reference position is \p q, or invalid_uint if it has none
  unsigned int node_at(const Elem & elem, const Point & q)
  {
    unsigned int found = invalid_uint;

    for (const auto n : elem.node_index_range())
      if (elem.master_point(n).relative_fuzzy_equals(q, TOLERANCE))
        {
          // Each shape function must pair with exactly one node
          CPPUNIT_ASSERT_EQUAL(invalid_uint, found);
          found = n;
        }

    return found;
  }

  void testAgreesWithLagrange()
  {
    LOG_UNIT_TEST;

    const unsigned int dim = Elem::type_to_dim_map[elem_type];

    const FEType gll_type(order, L2_LAGRANGE_GLL);
    const FEType lagrange_type(order, LAGRANGE);

    // Sample on a grid spread over the reference element, so that most points lie between
    // the interpolation points
    std::vector<Point> samples;
    for (const auto a : make_range(-3, 4))
      for (const auto b : make_range(-3, 4))
        for (const auto c : make_range(-3, 4))
          samples.emplace_back(a/7., (dim > 1) ? b/7. : 0., (dim > 2) ? c/7. : 0.);

    unsigned int n_paired = 0;

    for (const auto * elem : _mesh->active_local_element_ptr_range())
      {
        const unsigned int n_gll = FEInterface::n_shape_functions(gll_type, elem);
        const unsigned int n_lag = FEInterface::n_shape_functions(lagrange_type, elem);

        // The two bases span the same space here, so their sizes match
        CPPUNIT_ASSERT_EQUAL(n_lag, n_gll);

        for (const auto i : make_range(n_gll))
          {
            const unsigned int n = node_at(*elem, gll_point(i, dim));
            CPPUNIT_ASSERT(n != invalid_uint);

            for (const Point & q : samples)
              LIBMESH_ASSERT_FP_EQUAL(FEInterface::shape(lagrange_type, elem, n, q),
                                      FEInterface::shape(gll_type, elem, i, q),
                                      TOLERANCE*TOLERANCE);

            ++n_paired;
          }
      }

    _mesh->comm().sum(n_paired);
    CPPUNIT_ASSERT(n_paired > 0);
  }
};

#define INSTANTIATE_FEGLLAGREEMENTTEST(elemtype, order)                         \
  class FEGLLLagrangeAgreementTest_##order##_##elemtype :                       \
    public FEGLLLagrangeAgreementTest<elemtype, order> {                        \
  public:                                                                       \
  FEGLLLagrangeAgreementTest_##order##_##elemtype() :                           \
    FEGLLLagrangeAgreementTest<elemtype, order>() {                             \
    if (unitlog->summarized_logs_enabled())                                     \
      this->libmesh_suite_name = "FEGLLLagrangeAgreementTest";                  \
    else                                                                        \
      this->libmesh_suite_name =                                                \
        "FEGLLLagrangeAgreementTest_" #order "_" #elemtype;                     \
  }                                                                             \
  CPPUNIT_TEST_SUITE( FEGLLLagrangeAgreementTest_##order##_##elemtype );        \
  CPPUNIT_TEST( testAgreesWithLagrange );                                       \
  CPPUNIT_TEST_SUITE_END();                                                     \
  };                                                                            \
                                                                                \
  CPPUNIT_TEST_SUITE_REGISTRATION( FEGLLLagrangeAgreementTest_##order##_##elemtype )

// Element types whose nodes are the Gauss-Lobatto points of the given order. QUAD8 and
// HEX20 are excluded because LAGRANGE on those serendipity types spans a smaller space than
// the tensor product basis.
INSTANTIATE_FEGLLAGREEMENTTEST(EDGE2, FIRST);
INSTANTIATE_FEGLLAGREEMENTTEST(EDGE3, SECOND);

#if LIBMESH_DIM > 1
INSTANTIATE_FEGLLAGREEMENTTEST(QUAD4, FIRST);
INSTANTIATE_FEGLLAGREEMENTTEST(QUAD9, SECOND);
#endif

#if LIBMESH_DIM > 2
INSTANTIATE_FEGLLAGREEMENTTEST(HEX8, FIRST);
INSTANTIATE_FEGLLAGREEMENTTEST(HEX27, SECOND);
#endif
