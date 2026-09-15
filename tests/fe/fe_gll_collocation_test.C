// A nodal basis whose interpolation points are the points of a quadrature rule is collocated
// with it: each shape function is one at its own point and zero at every other, so the mass
// matrix that rule assembles is diagonal. That holds for any mapping and any element shape,
// since the Jacobian only scales a diagonal entry, and it holds for both Gauss-Lobatto nodal
// families, which differ in the order they number their degrees of freedom rather than in
// where those degrees of freedom sit.
//
// The L2 family numbers them in the order the rule lays out its points, so its matrix of shape
// function values is the identity. The C0 family numbers them by the entity that owns them, so
// the matrix is a permutation of the identity. What both hold to is that the map from a degree
// of freedom to its point is one to one, which is what the diagonal rests on.
//
// The property is exact rather than approximate, so the tests compare against 1 and 0 with no
// tolerance. That is what the basis reading its points from QGaussLobatto buys: a basis that
// computed its own would agree only to roundoff.

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

template <ElemType elem_type, Order order, FEFamily family>
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

    // Skew the mesh so that no element is a scaled copy of the reference element, which
    // holds the diagonal mass matrix to being a property of the basis rather than of a
    // Jacobian that happens to be constant.
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

  /**
   * The point of the collocated rule that each degree of freedom interpolates at, taken from
   * the shape function values and holding them to being one there and zero everywhere else.
   *
   * The L2 family numbers its degrees of freedom in the order the rule lays out its points, so
   * this map is the identity for it. The C0 family numbers them by the entity that owns them,
   * so the map is a permutation. Either way it is one to one, which is what makes the mass
   * matrix diagonal.
   */
  std::vector<unsigned int> collocation_map(const std::vector<std::vector<Real>> & phi)
  {
    std::vector<unsigned int> point_of_dof(phi.size(), invalid_uint);
    std::vector<bool> taken(phi.size(), false);

    for (const auto i : index_range(phi))
      {
        CPPUNIT_ASSERT_EQUAL(phi.size(), phi[i].size());

        for (const auto q : index_range(phi[i]))
          {
            if (phi[i][q] == 1.)
              {
                CPPUNIT_ASSERT_EQUAL(invalid_uint, point_of_dof[i]);
                point_of_dof[i] = cast_int<unsigned int>(q);
              }
            else
              // Exactly zero, not nearly: the basis and the rule carry the same points
              CPPUNIT_ASSERT_EQUAL(Real(0), phi[i][q]);
          }

        CPPUNIT_ASSERT(point_of_dof[i] != invalid_uint);

        // One degree of freedom per point, so that no point is claimed twice
        CPPUNIT_ASSERT(!taken[point_of_dof[i]]);
        taken[point_of_dof[i]] = true;
      }

    return point_of_dof;
  }

  void testShapesAreCollocated()
  {
    LOG_UNIT_TEST;

    const FEType fe_type(order, family);
    const unsigned int dim = Elem::type_to_dim_map[elem_type];

    std::unique_ptr<FEBase> fe = FEBase::build(dim, fe_type);
    std::unique_ptr<QBase> qrule = fe_type.default_quadrature_rule(dim);
    fe->attach_quadrature_rule(qrule.get());

    // The accessors have to be reached before reinit, which is what tells the element what to
    // compute
    const std::vector<std::vector<Real>> & phi = fe->get_phi();

    unsigned int n_elem = 0;

    for (const auto * elem : _mesh->active_local_element_ptr_range())
      {
        fe->reinit(elem);

        // The rule carries one point per degree of freedom
        CPPUNIT_ASSERT_EQUAL(phi.size(), std::size_t(qrule->n_points()));

        const std::vector<unsigned int> point_of_dof = collocation_map(phi);

        // The L2 family also agrees with the rule on the order of the points
        if (family == L2_LAGRANGE_GLL)
          for (const auto i : index_range(point_of_dof))
            CPPUNIT_ASSERT_EQUAL(cast_int<unsigned int>(i), point_of_dof[i]);

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
    std::unique_ptr<QBase> qrule = fe_type.default_quadrature_rule(dim);
    fe->attach_quadrature_rule(qrule.get());

    const std::vector<std::vector<Real>> & phi = fe->get_phi();
    const std::vector<Real> & JxW = fe->get_JxW();

    for (const auto * elem : _mesh->active_local_element_ptr_range())
      {
        fe->reinit(elem);

        const std::vector<unsigned int> point_of_dof = collocation_map(phi);
        const std::size_t n_dofs = phi.size();

        // Summed the way an assembly loop would sum it, so that the diagonal is a property of
        // the arithmetic that assembles a mass matrix rather than one read back from the shape
        // function values the assertions above already hold
        for (const auto i : make_range(n_dofs))
          for (const auto j : make_range(n_dofs))
            {
              Real m = 0.;
              for (const auto q : index_range(JxW))
                m += JxW[q] * phi[i][q] * phi[j][q];

              if (i == j)
                {
                  // A diagonal entry is the weight of its own point times the Jacobian there,
                  // which is positive on a valid element, so the mass matrix is positive
                  // definite and its inverse costs a reciprocal per entry
                  CPPUNIT_ASSERT(m > 0.);
                  LIBMESH_ASSERT_FP_EQUAL(JxW[point_of_dof[i]], m, TOLERANCE*TOLERANCE);
                }
              else
                CPPUNIT_ASSERT_EQUAL(Real(0), m);
            }
      }
  }
};

#define INSTANTIATE_FEGLLCOLLOCATIONTEST(elemtype, order, family)                  \
  class FEGLLCollocationTest_##family##_##order##_##elemtype :                     \
    public FEGLLCollocationTest<elemtype, order, family> {                         \
  public:                                                                          \
  FEGLLCollocationTest_##family##_##order##_##elemtype() :                         \
    FEGLLCollocationTest<elemtype, order, family>() {                              \
    if (unitlog->summarized_logs_enabled())                                        \
      this->libmesh_suite_name = "FEGLLCollocationTest";                           \
    else                                                                           \
      this->libmesh_suite_name =                                                   \
        "FEGLLCollocationTest_" #family "_" #order "_" #elemtype;                   \
  }                                                                                \
  CPPUNIT_TEST_SUITE( FEGLLCollocationTest_##family##_##order##_##elemtype );       \
  CPPUNIT_TEST( testShapesAreCollocated );                                         \
  CPPUNIT_TEST( testMassMatrixIsDiagonal );                                        \
  CPPUNIT_TEST_SUITE_END();                                                        \
  };                                                                               \
                                                                                   \
  CPPUNIT_TEST_SUITE_REGISTRATION( FEGLLCollocationTest_##family##_##order##_##elemtype )

INSTANTIATE_FEGLLCOLLOCATIONTEST(EDGE2, FIRST,  L2_LAGRANGE_GLL);
INSTANTIATE_FEGLLCOLLOCATIONTEST(EDGE3, THIRD,  L2_LAGRANGE_GLL);
INSTANTIATE_FEGLLCOLLOCATIONTEST(EDGE3, EIGHTH, L2_LAGRANGE_GLL);
INSTANTIATE_FEGLLCOLLOCATIONTEST(EDGE2, FIRST,  LAGRANGE_GLL);
INSTANTIATE_FEGLLCOLLOCATIONTEST(EDGE3, THIRD,  LAGRANGE_GLL);
INSTANTIATE_FEGLLCOLLOCATIONTEST(EDGE3, EIGHTH, LAGRANGE_GLL);

#if LIBMESH_DIM > 1
INSTANTIATE_FEGLLCOLLOCATIONTEST(QUAD4, FIRST,  L2_LAGRANGE_GLL);
INSTANTIATE_FEGLLCOLLOCATIONTEST(QUAD9, SECOND, L2_LAGRANGE_GLL);
INSTANTIATE_FEGLLCOLLOCATIONTEST(QUAD9, FOURTH, L2_LAGRANGE_GLL);
INSTANTIATE_FEGLLCOLLOCATIONTEST(QUAD9, EIGHTH, L2_LAGRANGE_GLL);
INSTANTIATE_FEGLLCOLLOCATIONTEST(QUAD4, FIRST,  LAGRANGE_GLL);
INSTANTIATE_FEGLLCOLLOCATIONTEST(QUAD9, SECOND, LAGRANGE_GLL);
INSTANTIATE_FEGLLCOLLOCATIONTEST(QUAD9, FOURTH, LAGRANGE_GLL);
INSTANTIATE_FEGLLCOLLOCATIONTEST(QUAD9, EIGHTH, LAGRANGE_GLL);
#endif

#if LIBMESH_DIM > 2
INSTANTIATE_FEGLLCOLLOCATIONTEST(HEX8,  FIRST,  L2_LAGRANGE_GLL);
INSTANTIATE_FEGLLCOLLOCATIONTEST(HEX27, SECOND, L2_LAGRANGE_GLL);
INSTANTIATE_FEGLLCOLLOCATIONTEST(HEX27, FOURTH, L2_LAGRANGE_GLL);
INSTANTIATE_FEGLLCOLLOCATIONTEST(HEX8,  FIRST,  LAGRANGE_GLL);
INSTANTIATE_FEGLLCOLLOCATIONTEST(HEX27, SECOND, LAGRANGE_GLL);
INSTANTIATE_FEGLLCOLLOCATIONTEST(HEX27, FOURTH, LAGRANGE_GLL);
#endif
