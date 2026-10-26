// On an adaptively refined mesh, a conforming basis constrains the degrees of freedom on a
// hanging node to the coarse side that contains the node. A conforming mesh does not exercise
// these constraints. This test checks the property that the constraints must preserve. The
// constrained space contains the polynomials of its order, so the projection of such a
// polynomial reproduces it exactly on both sides of the interface.
//
// A wrong constraint moves the hanging degrees of freedom away from the values that
// interpolate the polynomial. The solution then differs from the polynomial.

#include <libmesh/dof_map.h>
#include <libmesh/elem.h>
#include <libmesh/enum_fe_family.h>
#include <libmesh/enum_order.h>
#include <libmesh/equation_systems.h>
#include <libmesh/fe_type.h>
#include <libmesh/int_range.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_refinement.h>
#include <libmesh/system.h>
#include <libmesh/utility.h>

#include "test_comm.h"
#include "libmesh_cppunit.h"

#include <cmath>
#include <memory>
#include <vector>

using namespace libMesh;

namespace {

/**
 * A polynomial of total degree \p D. The tensor product space of order \p D contains it, so
 * a projection onto a family of that order reproduces it.
 */
template <unsigned int D>
Number amr_test_poly(const Point & p,
                     const Parameters &,
                     const std::string &,
                     const std::string &)
{
  return Utility::pow<D>(1. + p(0) + 2.*p(1) + 3.*p(2));
}

}

#ifdef LIBMESH_ENABLE_AMR

template <ElemType elem_type, Order order, FEFamily family>
class FEAMRConstraintTest : public CppUnit::TestCase
{
protected:
  std::unique_ptr<Mesh> _mesh;
  std::string libmesh_suite_name;

public:
  void setUp()
  {
    _mesh = std::make_unique<Mesh>(*TestCommWorld);

    const unsigned int dim = Elem::type_to_dim_map[elem_type];

    // For a second-order type, build the first-order mesh and then convert it. Families above
    // order one need the mid-edge and face nodes of the second-order type.
    ElemType base = elem_type;
    if (elem_type == QUAD9)  base = QUAD4;
    if (elem_type == HEX27)  base = HEX8;
    if (elem_type == EDGE3)  base = EDGE2;

    if (dim == 1)
      MeshTools::Generation::build_line(*_mesh, 4, 0., 1., base);
    else if (dim == 2)
      MeshTools::Generation::build_square(*_mesh, 2, 2, 0., 1., 0., 1., base);
    else
      MeshTools::Generation::build_cube(*_mesh, 2, 2, 2, 0., 1., 0., 1., 0., 1., base);

    if (base != elem_type)
      _mesh->all_second_order();

    // Refine part of the mesh. This puts hanging nodes on the interface between the refined
    // elements and their unrefined neighbors.
    MeshRefinement refine(*_mesh);

    for (auto * elem : _mesh->active_element_ptr_range())
      if (elem->id() % 4 == 0)
        elem->set_refinement_flag(Elem::REFINE);

    refine.refine_and_coarsen_elements();
  }

  void tearDown() { _mesh.reset(); }

  void testConstrainedSpaceHoldsPolynomials()
  {
    LOG_UNIT_TEST;

    EquationSystems es(*_mesh);
    System & sys = es.add_system<System>("amr");
    sys.add_variable("u", FEType(order, family));
    es.init();

    sys.project_solution(amr_test_poly<static_cast<unsigned int>(order)>,
                         nullptr, es.parameters);

    // The projection enforced the constraints. Check that the solution still equals the
    // polynomial at points in the refined and the unrefined parts of the mesh.
    Parameters dummy;
    unsigned int n_checked = 0;

    for (const auto a : make_range(1, 8))
      for (const auto b : make_range(1, 8))
        {
          const unsigned int dim = Elem::type_to_dim_map[elem_type];
          const Point q(a/8., (dim > 1) ? b/8. : 0., (dim > 2) ? (a + b)/16. : 0.);

          const Number expected =
            amr_test_poly<static_cast<unsigned int>(order)>(q, dummy, "", "");
          const Number actual = sys.point_value(0, q);

          LIBMESH_ASSERT_NUMBERS_EQUAL(expected, actual, TOLERANCE*std::sqrt(TOLERANCE));
          ++n_checked;
        }

    CPPUNIT_ASSERT(n_checked > 0);
  }
};

#define INSTANTIATE_FEAMRCONSTRAINTTEST(elemtype, order, family)                    \
  class FEAMRConstraintTest_##family##_##order##_##elemtype :                       \
    public FEAMRConstraintTest<elemtype, order, family> {                           \
  public:                                                                           \
  FEAMRConstraintTest_##family##_##order##_##elemtype() :                           \
    FEAMRConstraintTest<elemtype, order, family>() {                                \
    if (unitlog->summarized_logs_enabled())                                         \
      this->libmesh_suite_name = "FEAMRConstraintTest";                             \
    else                                                                            \
      this->libmesh_suite_name =                                                    \
        "FEAMRConstraintTest_" #family "_" #order "_" #elemtype;                     \
  }                                                                                 \
  CPPUNIT_TEST_SUITE( FEAMRConstraintTest_##family##_##order##_##elemtype );         \
  CPPUNIT_TEST( testConstrainedSpaceHoldsPolynomials );                             \
  CPPUNIT_TEST_SUITE_END();                                                         \
  };                                                                                \
                                                                                    \
  CPPUNIT_TEST_SUITE_REGISTRATION( FEAMRConstraintTest_##family##_##order##_##elemtype )

// LAGRANGE_GLL is the family under test. LAGRANGE is the control, because its constraints
// are well established. A LAGRANGE failure thus indicates a defect in the test. HIERARCHIC is
// the other family that adds a separate degree of freedom at a hanging node.
INSTANTIATE_FEAMRCONSTRAINTTEST(QUAD9, SECOND, LAGRANGE);
INSTANTIATE_FEAMRCONSTRAINTTEST(QUAD9, SECOND, LAGRANGE_GLL);
INSTANTIATE_FEAMRCONSTRAINTTEST(QUAD9, THIRD,  LAGRANGE_GLL);
INSTANTIATE_FEAMRCONSTRAINTTEST(QUAD9, FOURTH, LAGRANGE_GLL);
INSTANTIATE_FEAMRCONSTRAINTTEST(QUAD9, SECOND, HIERARCHIC);
INSTANTIATE_FEAMRCONSTRAINTTEST(QUAD9, THIRD,  HIERARCHIC);
INSTANTIATE_FEAMRCONSTRAINTTEST(QUAD9, FOURTH, HIERARCHIC);

#if LIBMESH_DIM > 2
INSTANTIATE_FEAMRCONSTRAINTTEST(HEX27, SECOND, LAGRANGE);
INSTANTIATE_FEAMRCONSTRAINTTEST(HEX27, SECOND, LAGRANGE_GLL);
INSTANTIATE_FEAMRCONSTRAINTTEST(HEX27, THIRD,  LAGRANGE_GLL);
INSTANTIATE_FEAMRCONSTRAINTTEST(HEX27, FOURTH, LAGRANGE_GLL);
INSTANTIATE_FEAMRCONSTRAINTTEST(HEX27, SECOND, HIERARCHIC);
INSTANTIATE_FEAMRCONSTRAINTTEST(HEX27, THIRD,  HIERARCHIC);
INSTANTIATE_FEAMRCONSTRAINTTEST(HEX27, FOURTH, HIERARCHIC);
#endif

#endif // LIBMESH_ENABLE_AMR
