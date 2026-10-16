// A nodal basis whose interpolation points are the points of a quadrature rule is
// collocated with that rule. Each shape function is one at its own point and zero at the
// others, so the matrix of shape function values at the quadrature points is the identity.
// The mass matrix assembled with the rule is then diagonal. This holds for any mapping and
// element shape, because the Jacobian only scales the diagonal entries.
//
// The basis reads its points from QGaussLobatto, so the property holds exactly. These tests
// compare shape function values and off-diagonal mass matrix entries against 1 and 0
// without a tolerance.

#include <libmesh/elem.h>
#include <libmesh/enum_fe_family.h>
#include <libmesh/enum_order.h>
#include <libmesh/fe_base.h>
#include <libmesh/fe_type.h>
#include <libmesh/int_range.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/node.h>
#include <libmesh/parallel_implementation.h>
#include <libmesh/quadrature.h>

#include "test_comm.h"
#include "libmesh_cppunit.h"

#include <cmath>
#include <memory>
#include <vector>

using namespace libMesh;

template <ElemType elem_type, Order order>
class FEGLLCollocationTest : public CppUnit::TestCase
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

  /// The family's default quadrature rule, which is collocated with the basis
  std::unique_ptr<QBase> collocated_rule(const FEType & fe_type, const unsigned int dim)
  {
    return fe_type.default_quadrature_rule(dim);
  }

  void testShapesAreIdentity()
  {
    LOG_UNIT_TEST;

    const FEType fe_type(order, L2_LAGRANGE_GLL);
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

        // The collocated rule has one point per degree of freedom
        CPPUNIT_ASSERT_EQUAL(phi.size(), std::size_t(qrule->n_points()));

        for (const auto i : index_range(phi))
          for (const auto q : index_range(phi[i]))
            CPPUNIT_ASSERT_EQUAL((i == q) ? Real(1) : Real(0), phi[i][q]);

        ++n_elem;
      }

    _mesh->comm().sum(n_elem);
    CPPUNIT_ASSERT(n_elem > 0);
  }

  void testMassMatrixIsDiagonal()
  {
    LOG_UNIT_TEST;

    const FEType fe_type(order, L2_LAGRANGE_GLL);
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

        for (const auto i : make_range(n_dofs))
          for (const auto j : make_range(n_dofs))
            {
              Real m = 0.;
              for (const auto q : index_range(JxW))
                m += JxW[q] * phi[i][q] * phi[j][q];

              if (i == j)
                {
                  // Each diagonal entry equals JxW at its own point, which is positive on
                  // a valid element
                  CPPUNIT_ASSERT(m > 0.);
                  LIBMESH_ASSERT_FP_EQUAL(JxW[i], m, TOLERANCE*TOLERANCE);
                }
              else
                CPPUNIT_ASSERT_EQUAL(Real(0), m);
            }
      }
  }
};

#define INSTANTIATE_FEGLLCOLLOCATIONTEST(elemtype, order)                       \
  class FEGLLCollocationTest_##order##_##elemtype :                             \
    public FEGLLCollocationTest<elemtype, order> {                              \
  public:                                                                       \
  FEGLLCollocationTest_##order##_##elemtype() :                                 \
    FEGLLCollocationTest<elemtype, order>() {                                   \
    if (unitlog->summarized_logs_enabled())                                     \
      this->libmesh_suite_name = "FEGLLCollocationTest";                        \
    else                                                                        \
      this->libmesh_suite_name =                                                \
        "FEGLLCollocationTest_" #order "_" #elemtype;                           \
  }                                                                             \
  CPPUNIT_TEST_SUITE( FEGLLCollocationTest_##order##_##elemtype );              \
  CPPUNIT_TEST( testShapesAreIdentity );                                        \
  CPPUNIT_TEST( testMassMatrixIsDiagonal );                                     \
  CPPUNIT_TEST_SUITE_END();                                                     \
  };                                                                            \
                                                                                \
  CPPUNIT_TEST_SUITE_REGISTRATION( FEGLLCollocationTest_##order##_##elemtype )

INSTANTIATE_FEGLLCOLLOCATIONTEST(EDGE2, FIRST);
INSTANTIATE_FEGLLCOLLOCATIONTEST(EDGE3, THIRD);
INSTANTIATE_FEGLLCOLLOCATIONTEST(EDGE3, EIGHTH);

#if LIBMESH_DIM > 1
INSTANTIATE_FEGLLCOLLOCATIONTEST(QUAD4, FIRST);
INSTANTIATE_FEGLLCOLLOCATIONTEST(QUAD9, SECOND);
INSTANTIATE_FEGLLCOLLOCATIONTEST(QUAD9, FOURTH);
INSTANTIATE_FEGLLCOLLOCATIONTEST(QUAD9, EIGHTH);
#endif

#if LIBMESH_DIM > 2
INSTANTIATE_FEGLLCOLLOCATIONTEST(HEX8, FIRST);
INSTANTIATE_FEGLLCOLLOCATIONTEST(HEX27, SECOND);
INSTANTIATE_FEGLLCOLLOCATIONTEST(HEX27, FOURTH);
#endif
