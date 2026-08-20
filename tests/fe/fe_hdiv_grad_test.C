// libmesh includes
#include "libmesh/elem.h"
#include "libmesh/fe.h"
#include "libmesh/fe_map.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/node.h"
#include "libmesh/quadrature_gauss.h"

// unit test includes
#include "test_comm.h"

#include "libmesh_cppunit.h"

using namespace libMesh;

/**
 * Verifies HDivFETransformation::map_dphi (the physical-space gradients of
 * H(div)-conforming Raviart-Thomas shape functions) by comparing
 * fe->get_dphi() against a central-difference approximation built from
 * map_phi's own output (fe->get_phi()), evaluated at physical points
 * perturbed via FEMap::inverse_map.
 */
class HDivGradTest : public CppUnit::TestCase
{
public:
  LIBMESH_CPPUNIT_TEST_SUITE(HDivGradTest);

#if LIBMESH_DIM > 1
  CPPUNIT_TEST(testQuad8Affine);
  CPPUNIT_TEST(testQuad8NonAffine);
  CPPUNIT_TEST(testTri6Affine);
  CPPUNIT_TEST(testTri6NonAffine);
#endif
#if LIBMESH_DIM > 2
  CPPUNIT_TEST(testHex27Affine);
  CPPUNIT_TEST(testHex27NonAffine);
  CPPUNIT_TEST(testTet14Affine);
  CPPUNIT_TEST(testTet14NonAffine);
#endif

  CPPUNIT_TEST_SUITE_END();

public:
  void setUp() {}
  void tearDown() {}

#if LIBMESH_DIM > 1
  // RAVIART_THOMAS dof placement piggybacks on the second-order edge/face
  // nodes, so the geometric element types must be the second-order ones
  // (QUAD8, TRI6, ...) even though the RT approximation order used here
  // is FIRST.
  void testQuad8Affine() { testGradByFiniteDifference(QUAD8, false); }
  void testQuad8NonAffine() { testGradByFiniteDifference(QUAD8, true); }
  void testTri6Affine() { testGradByFiniteDifference(TRI6, false); }
  void testTri6NonAffine() { testGradByFiniteDifference(TRI6, true); }
#endif

#if LIBMESH_DIM > 2
  void testHex27Affine() { testGradByFiniteDifference(HEX27, false); }
  void testHex27NonAffine() { testGradByFiniteDifference(HEX27, true); }
  void testTet14Affine() { testGradByFiniteDifference(TET14, false); }
  void testTet14NonAffine() { testGradByFiniteDifference(TET14, true); }
#endif

private:
  // Builds a single-element mesh of the given type, optionally distorting
  // one node so the element's map is non-affine, then checks the HDiv
  // shape function gradients returned by fe->get_dphi() at every
  // quadrature point and coordinate direction.
  void testGradByFiniteDifference(const ElemType elem_type, const bool distort)
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
    unsigned int dim = 0;

    switch (elem_type)
      {
      case QUAD8:
      case TRI6:
        dim = 2;
        MeshTools::Generation::build_square
          (mesh, 1, 1, 0., 1., 0., 1., elem_type);
        break;
      case HEX27:
      case TET14:
        dim = 3;
        MeshTools::Generation::build_cube
          (mesh, 1, 1, 1, 0., 1., 0., 1., 0., 1., elem_type);
        break;
      default:
        libmesh_error_msg("Unsupported element type " << elem_type);
      }

    if (distort)
      {
        // Move a single vertex node without touching the associated
        // edge/face nodes, turning the element's boundary curved and
        // giving the map nonzero second derivatives -- this exercises
        // the terms in map_dphi (A and B in the derivation) that vanish
        // identically on an affine element.
        Node & node = mesh.node_ref(0);
        for (auto d : make_range(dim))
          node(d) += Real(0.1) * (d + 1);
      }

    auto elem_range = mesh.active_local_element_ptr_range();
    if (elem_range.begin() == elem_range.end())
      return;
    const Elem * elem = *(elem_range.begin());

    const FEType fe_type(FIRST, RAVIART_THOMAS);
    std::unique_ptr<FEVectorBase> fe(FEVectorBase::build(dim, fe_type));

    const std::vector<std::vector<RealTensor>> & dphi = fe->get_dphi();
    const std::vector<std::vector<RealGradient>> & phi = fe->get_phi();
    const std::vector<Point> & xyz = fe->get_xyz();

    QGauss qrule(dim, fe_type.default_quadrature_order());
    fe->attach_quadrature_rule(&qrule);

    fe->reinit(elem);

    // Copy out what we need before the perturbed reinit() calls below
    // overwrite fe's internal state.
    const std::vector<std::vector<RealTensor>> dphi_base = dphi;
    const std::vector<Point> qpoints = xyz;

    const Real h = 1e-6;
    const Real tol = 5e-4;

    for (auto qp : index_range(qpoints))
      {
        const Point x0 = qpoints[qp];

        for (auto l : make_range(dim))
          {
            Point x_plus = x0;
            x_plus(l) += h;
            Point x_minus = x0;
            x_minus(l) -= h;

            std::vector<Point> pts_plus
              (1, FEMap::inverse_map(dim, elem, x_plus));
            fe->reinit(elem, &pts_plus);
            const std::vector<std::vector<RealGradient>> phi_plus = phi;

            std::vector<Point> pts_minus
              (1, FEMap::inverse_map(dim, elem, x_minus));
            fe->reinit(elem, &pts_minus);
            const std::vector<std::vector<RealGradient>> phi_minus = phi;

            for (auto i : index_range(dphi_base))
              for (auto k : make_range(LIBMESH_DIM))
                {
                  const Real fd =
                    (phi_plus[i][0](k) - phi_minus[i][0](k)) / (2*h);

                  LIBMESH_ASSERT_FP_EQUAL(fd, dphi_base[i][qp](k, l), tol);
                }
          }
      }
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION(HDivGradTest);
