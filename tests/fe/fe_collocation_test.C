// A nodal basis whose interpolation points are the points of a quadrature rule is
// collocated with that rule. Each shape function is one at its own point and zero at the
// others. The matrix of shape function values at the quadrature points is thus a
// permutation matrix. The mass matrix assembled with the rule is then diagonal. This holds for any mapping
// and element shape, because the Jacobian only scales the diagonal entries.
//
// The Gauss-Lobatto families are collocated with their default Gauss-Lobatto rule.
// L2_LAGRANGE_GLL numbers its shape functions in the order of the rule points, so its
// permutation is the identity. LAGRANGE_GLL numbers them by the entity that owns them, so its
// permutation is general. Both families read their points from QGaussLobatto, so the
// property holds exactly. These cases compare against 1 and 0 without a tolerance.
//
// LAGRANGE and L2_LAGRANGE are collocated with QNODAL. The QNODAL points are the element
// nodes, in an order that can differ from the node numbering. The shape functions of these
// families are evaluated in floating point at the nodes, so these cases use a tolerance.

#include <libmesh/elem.h>
#include <libmesh/enum_fe_family.h>
#include <libmesh/enum_order.h>
#include <libmesh/fe_base.h>
#include <libmesh/fe_type.h>
#include <libmesh/int_range.h>
#include <libmesh/libmesh_common.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/node.h>
#include <libmesh/parallel_implementation.h>
#include <libmesh/quadrature.h>
#include <libmesh/quadrature_nodal.h>

#include "test_comm.h"
#include "libmesh_cppunit.h"

#include <cmath>
#include <memory>
#include <vector>

using namespace libMesh;

template <FEFamily family, ElemType elem_type, Order order>
class FECollocationTest : public CppUnit::TestCase
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
      MeshTools::Generation::build_line(*_mesh, 2, 0., 1., elem_type);
    else if (dim == 2)
      MeshTools::Generation::build_square(*_mesh, 2, 2, 0., 1., 0., 1., elem_type);
    else
      MeshTools::Generation::build_cube(*_mesh, 2, 2, 2, 0., 1., 0., 1., 0., 1., elem_type);

    // Shear the mesh so that the mapping to each element has off-diagonal Jacobian terms
    for (auto * node : _mesh->node_ptr_range())
      {
        Node & p = *node;
        const Real y = (LIBMESH_DIM > 1) ? p(1) : 0.;
        const Real z = (LIBMESH_DIM > 2) ? p(2) : 0.;

        p(0) += 0.3 * y + 0.17 * z;
        if (LIBMESH_DIM > 1)
          p(1) += 0.23 * z;
      }
  }

  void tearDown() { _mesh.reset(); }

  /// The quadrature rule collocated with the basis
  std::unique_ptr<QBase> collocated_rule(const FEType & fe_type, const unsigned int dim)
  {
    if (gll)
      return fe_type.default_quadrature_rule(dim);

    return std::make_unique<QNodal>(dim);
  }

  /// Whether the family is one of the Gauss-Lobatto nodal families
  static constexpr bool gll = (family == L2_LAGRANGE_GLL || family == LAGRANGE_GLL);

  /// The tolerance on shape function values at the collocation points
  static constexpr Real tol = gll ? Real(0) : TOLERANCE*TOLERANCE;

  /**
   * \returns Entry i is the index of the quadrature point at which shape function i is one,
   * after checking that the shape function values at the quadrature points form a
   * permutation matrix
   */
  std::vector<unsigned int>
  collocation_points(const std::vector<std::vector<Real>> & phi)
  {
    std::vector<unsigned int> shape_qp(phi.size(), invalid_uint);

    // Whether some shape function is already one at each quadrature point. A second shape
    // function that is one at the same point would leave another point without one.
    std::vector<bool> qp_claimed(phi.size(), false);

    for (const auto i : index_range(phi))
      {
        // The collocated rule has one point per degree of freedom
        CPPUNIT_ASSERT_EQUAL(phi.size(), phi[i].size());

        for (const auto q : index_range(phi[i]))
          if (std::abs(phi[i][q] - 1) <= tol)
            {
              CPPUNIT_ASSERT_EQUAL(invalid_uint, shape_qp[i]);
              shape_qp[i] = q;
            }
          else
            LIBMESH_ASSERT_FP_EQUAL(Real(0), phi[i][q], tol);

        CPPUNIT_ASSERT(shape_qp[i] != invalid_uint);
        CPPUNIT_ASSERT(!qp_claimed[shape_qp[i]]);
        qp_claimed[shape_qp[i]] = true;
      }

    return shape_qp;
  }

  void testShapesArePermutation()
  {
    LOG_UNIT_TEST;

    const FEType fe_type(order, family);
    const unsigned int dim = Elem::type_to_dim_map[elem_type];

    std::unique_ptr<FEBase> fe = FEBase::build(dim, fe_type);
    std::unique_ptr<QBase> qrule = collocated_rule(fe_type, dim);
    fe->attach_quadrature_rule(qrule.get());

    // Request phi before reinit so that reinit computes it
    const std::vector<std::vector<Real>> & phi = fe->get_phi();

    unsigned int n_elem = 0;

    for (const auto * elem : _mesh->active_local_element_ptr_range())
      {
        fe->reinit(elem);

        const std::vector<unsigned int> shape_qp = collocation_points(phi);

        if (family == L2_LAGRANGE_GLL)
          for (const auto i : index_range(shape_qp))
            CPPUNIT_ASSERT_EQUAL(cast_int<unsigned int>(i), shape_qp[i]);

        ++n_elem;
      }

    _mesh->comm().sum(n_elem);
    CPPUNIT_ASSERT(n_elem > 0);
  }

  void testMassMatrixIsDiagonal()
  {
    LOG_UNIT_TEST;

    const FEType fe_type(order, family);
    const unsigned int dim = Elem::type_to_dim_map[elem_type];

    std::unique_ptr<FEBase> fe = FEBase::build(dim, fe_type);
    std::unique_ptr<QBase> qrule = collocated_rule(fe_type, dim);
    fe->attach_quadrature_rule(qrule.get());

    const std::vector<std::vector<Real>> & phi = fe->get_phi();
    const std::vector<Real> & JxW = fe->get_JxW();

    for (const auto * elem : _mesh->active_local_element_ptr_range())
      {
        fe->reinit(elem);

        const std::size_t n_dofs = phi.size();
        const std::vector<unsigned int> shape_qp = collocation_points(phi);

        for (const auto i : make_range(n_dofs))
          for (const auto j : make_range(n_dofs))
            {
              Real m = 0.;
              for (const auto q : index_range(JxW))
                m += JxW[q] * phi[i][q] * phi[j][q];

              if (i == j)
                {
                  // Each diagonal entry equals JxW at its own point, which is positive on
                  // a valid element with positive quadrature weights
                  CPPUNIT_ASSERT(m > 0.);
                  LIBMESH_ASSERT_FP_EQUAL(JxW[shape_qp[i]], m, TOLERANCE*TOLERANCE);
                }
              else
                LIBMESH_ASSERT_FP_EQUAL(Real(0), m, tol);
            }
      }
  }
};

#define INSTANTIATE_FECOLLOCATIONTEST(family, elemtype, order)                 \
  class FECollocationTest_##family##_##order##_##elemtype :                     \
    public FECollocationTest<family, elemtype, order> {                         \
  public:                                                                       \
  FECollocationTest_##family##_##order##_##elemtype() :                         \
    FECollocationTest<family, elemtype, order>() {                              \
    if (unitlog->summarized_logs_enabled())                                     \
      this->libmesh_suite_name = "FECollocationTest";                           \
    else                                                                        \
      this->libmesh_suite_name =                                                \
        "FECollocationTest_" #family "_" #order "_" #elemtype;                  \
  }                                                                             \
  CPPUNIT_TEST_SUITE( FECollocationTest_##family##_##order##_##elemtype );      \
  CPPUNIT_TEST( testShapesArePermutation );                                     \
  CPPUNIT_TEST( testMassMatrixIsDiagonal );                                     \
  CPPUNIT_TEST_SUITE_END();                                                     \
  };                                                                            \
                                                                                \
  CPPUNIT_TEST_SUITE_REGISTRATION( FECollocationTest_##family##_##order##_##elemtype )

INSTANTIATE_FECOLLOCATIONTEST(L2_LAGRANGE_GLL, EDGE2, FIRST);
INSTANTIATE_FECOLLOCATIONTEST(L2_LAGRANGE_GLL, EDGE3, THIRD);
INSTANTIATE_FECOLLOCATIONTEST(L2_LAGRANGE_GLL, EDGE3, EIGHTH);
INSTANTIATE_FECOLLOCATIONTEST(LAGRANGE_GLL, EDGE2, FIRST);
INSTANTIATE_FECOLLOCATIONTEST(LAGRANGE_GLL, EDGE3, THIRD);
INSTANTIATE_FECOLLOCATIONTEST(LAGRANGE_GLL, EDGE3, EIGHTH);
INSTANTIATE_FECOLLOCATIONTEST(LAGRANGE, EDGE2, FIRST);
INSTANTIATE_FECOLLOCATIONTEST(LAGRANGE, EDGE3, SECOND);
INSTANTIATE_FECOLLOCATIONTEST(LAGRANGE, EDGE4, THIRD);
INSTANTIATE_FECOLLOCATIONTEST(L2_LAGRANGE, EDGE3, SECOND);

#if LIBMESH_DIM > 1
INSTANTIATE_FECOLLOCATIONTEST(L2_LAGRANGE_GLL, QUAD4, FIRST);
INSTANTIATE_FECOLLOCATIONTEST(L2_LAGRANGE_GLL, QUAD9, SECOND);
INSTANTIATE_FECOLLOCATIONTEST(L2_LAGRANGE_GLL, QUAD9, FOURTH);
INSTANTIATE_FECOLLOCATIONTEST(L2_LAGRANGE_GLL, QUAD9, EIGHTH);
INSTANTIATE_FECOLLOCATIONTEST(LAGRANGE_GLL, QUAD4, FIRST);
INSTANTIATE_FECOLLOCATIONTEST(LAGRANGE_GLL, QUAD9, SECOND);
INSTANTIATE_FECOLLOCATIONTEST(LAGRANGE_GLL, QUAD9, FOURTH);
INSTANTIATE_FECOLLOCATIONTEST(LAGRANGE_GLL, QUAD9, EIGHTH);
INSTANTIATE_FECOLLOCATIONTEST(LAGRANGE, TRI3, FIRST);
INSTANTIATE_FECOLLOCATIONTEST(LAGRANGE, TRI6, SECOND);
INSTANTIATE_FECOLLOCATIONTEST(LAGRANGE, TRI7, THIRD);
INSTANTIATE_FECOLLOCATIONTEST(LAGRANGE, QUAD4, FIRST);
INSTANTIATE_FECOLLOCATIONTEST(LAGRANGE, QUAD8, SECOND);
INSTANTIATE_FECOLLOCATIONTEST(LAGRANGE, QUAD9, SECOND);
INSTANTIATE_FECOLLOCATIONTEST(L2_LAGRANGE, TRI6, SECOND);
INSTANTIATE_FECOLLOCATIONTEST(L2_LAGRANGE, QUAD9, SECOND);
#endif

#if LIBMESH_DIM > 2
INSTANTIATE_FECOLLOCATIONTEST(L2_LAGRANGE_GLL, HEX8, FIRST);
INSTANTIATE_FECOLLOCATIONTEST(L2_LAGRANGE_GLL, HEX27, SECOND);
INSTANTIATE_FECOLLOCATIONTEST(L2_LAGRANGE_GLL, HEX27, FOURTH);
INSTANTIATE_FECOLLOCATIONTEST(LAGRANGE_GLL, HEX8, FIRST);
INSTANTIATE_FECOLLOCATIONTEST(LAGRANGE_GLL, HEX27, SECOND);
INSTANTIATE_FECOLLOCATIONTEST(LAGRANGE_GLL, HEX27, FOURTH);
INSTANTIATE_FECOLLOCATIONTEST(LAGRANGE, TET4, FIRST);
INSTANTIATE_FECOLLOCATIONTEST(LAGRANGE, TET10, SECOND);
INSTANTIATE_FECOLLOCATIONTEST(LAGRANGE, TET14, THIRD);
INSTANTIATE_FECOLLOCATIONTEST(LAGRANGE, HEX8, FIRST);
INSTANTIATE_FECOLLOCATIONTEST(LAGRANGE, HEX20, SECOND);
INSTANTIATE_FECOLLOCATIONTEST(LAGRANGE, HEX27, SECOND);
INSTANTIATE_FECOLLOCATIONTEST(LAGRANGE, PRISM6, FIRST);
INSTANTIATE_FECOLLOCATIONTEST(LAGRANGE, PRISM15, SECOND);
INSTANTIATE_FECOLLOCATIONTEST(LAGRANGE, PRISM18, SECOND);
INSTANTIATE_FECOLLOCATIONTEST(LAGRANGE, PRISM20, THIRD);
INSTANTIATE_FECOLLOCATIONTEST(LAGRANGE, PRISM21, THIRD);
INSTANTIATE_FECOLLOCATIONTEST(L2_LAGRANGE, HEX27, SECOND);
#endif
