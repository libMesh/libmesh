// libmesh includes
#include <libmesh/elem.h>
#include <libmesh/enum_elem_type.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_modification.h>
#include <libmesh/mesh.h>
#include <libmesh/reference_elem.h>
#include <libmesh/replicated_mesh.h>
#include <libmesh/node.h>
#include <libmesh/enum_to_string.h>
#include <libmesh/tensor_value.h>
#include <libmesh/enum_elem_quality.h>

// unit test includes
#include "test_comm.h"
#include "libmesh_cppunit.h"

// C++ includes
#include <iomanip>

using namespace libMesh;

class VolumeTest : public CppUnit::TestCase
{

public:
  LIBMESH_CPPUNIT_TEST_SUITE( VolumeTest );
  CPPUNIT_TEST( testTwistedVolume );
  CPPUNIT_TEST( testEdge3Volume );
  CPPUNIT_TEST( testEdge3Invertible );
  CPPUNIT_TEST( testEdge4Invertible );
  CPPUNIT_TEST( testQuad4Invertible );
  CPPUNIT_TEST( testTri3TrueCentroid );
  CPPUNIT_TEST( testQuad4TrueCentroid );
  CPPUNIT_TEST( testPyramid5TrueCentroid );
  CPPUNIT_TEST( testHex8TrueCentroid );
  CPPUNIT_TEST( testPrism6TrueCentroid );
  CPPUNIT_TEST( testHex20PLevelTrueCentroid );
  CPPUNIT_TEST( testQuad4AspectRatio );
  CPPUNIT_TEST( testQuad4Warpage );
  CPPUNIT_TEST( testQuad4MinMaxAngle );
  CPPUNIT_TEST( testQuad4Jacobian );
  CPPUNIT_TEST( testTri3AspectRatio );
  CPPUNIT_TEST( testTet4DihedralAngle );
  CPPUNIT_TEST( testTet4Jacobian );
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp()
  {
  }

  void tearDown()
  {
  }

  void testTri3TrueCentroid()
  {
    LOG_UNIT_TEST;

    // The true_centroid() == vertex_average() == (1/3, 1/3) for reference Tri3
    {
      const Elem & tri3 = ReferenceElem::get(TRI3);
      Point true_centroid = tri3.true_centroid();
      LIBMESH_ASSERT_FP_EQUAL(Real(1)/3, true_centroid(0), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(Real(1)/3, true_centroid(1), TOLERANCE*TOLERANCE);
    }
  }

  void testQuad4TrueCentroid()
  {
    LOG_UNIT_TEST;

    // Test Quad4::true_centroid() override
    {
      const Elem & quad4 = ReferenceElem::get(QUAD4);
      Point true_centroid = quad4.true_centroid();
      LIBMESH_ASSERT_FP_EQUAL(0, true_centroid(0), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(0, true_centroid(1), TOLERANCE*TOLERANCE);

      // Compare to centroid computed via generic base class implementation
      Point base_class_centroid = quad4.Elem::true_centroid();
      CPPUNIT_ASSERT(true_centroid.absolute_fuzzy_equals(base_class_centroid, TOLERANCE*TOLERANCE));
    }

    // Check that "optimized" Quad4::true_centroid() gives same result
    // as "generic" Elem::true_centroid() on a mesh of 10% distorted elements.
    {
      ReplicatedMesh mesh(*TestCommWorld);

      MeshTools::Generation::build_square(mesh,
                                          /*nx=*/3, /*ny=*/3,
                                          /*xmin=*/0., /*xmax=*/1.,
                                          /*ymin=*/0., /*ymax=*/1.,
                                          QUAD4);

      MeshTools::Modification::distort(mesh,
                                       /*factor=*/0.1,
                                       /*perturb_boundary=*/false);

      for (const auto & elem : mesh.element_ptr_range())
        {
          Point derived_centroid = elem->true_centroid();
          Point base_centroid = elem->Elem::true_centroid();

          // Debugging: check results in detail
          // auto flags = libMesh::out.flags();
          // libMesh::out << std::scientific << std::setprecision(16);
          // libMesh::out << "derived_centroid = " << derived_centroid << std::endl;
          // libMesh::out << "base_centroid = " << base_centroid << std::endl;
          // libMesh::out.flags(flags);

          CPPUNIT_ASSERT(derived_centroid.absolute_fuzzy_equals(base_centroid, TOLERANCE*TOLERANCE));

          Real derived_volume = elem->volume();
          Real base_volume = elem->Elem::volume();
          LIBMESH_ASSERT_FP_EQUAL(base_volume, derived_volume, TOLERANCE*TOLERANCE);
        }
    }
  }

  void testPyramid5TrueCentroid()
  {
    LOG_UNIT_TEST;

    // Test Pyramid5::true_centroid() gives the correct result for a reference element
    {
      const Elem & pyr5 = ReferenceElem::get(PYRAMID5);
      Point true_centroid = pyr5.true_centroid();
      LIBMESH_ASSERT_FP_EQUAL(0, true_centroid(0), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(0, true_centroid(1), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(0.25, true_centroid(2), TOLERANCE*TOLERANCE);
    }

    // Currently there is not an optimized Pyramid5::true_centroid() to compare against
    test_true_centroid_and_volume(PYRAMID5);
  }

  void testHex8TrueCentroid() { LOG_UNIT_TEST; test_true_centroid_and_volume(HEX8); }
  void testPrism6TrueCentroid() { LOG_UNIT_TEST; test_true_centroid_and_volume(PRISM6); }

  void testHex20PLevelTrueCentroid()
  {
    LOG_UNIT_TEST;

    // Test that Elem base class true_centroid() implementation works
    // for an elevated p_level HEX20
    {
#ifdef LIBMESH_ENABLE_AMR
      ReplicatedMesh mesh(*TestCommWorld);
      MeshTools::Generation::build_cube(mesh,
                                        /*nelem=*/1, /*nelem=*/1, /*nelem=*/1,
                                        /*xmin=*/-1, /*xmax=*/1,
                                        /*ymin=*/-1, /*ymax=*/1,
                                        /*zmin=*/-1, /*zmax=*/1,
                                        HEX20);
      Elem * hex20 = mesh.elem_ptr(0);
      hex20->set_p_level(1);
      Point true_centroid = hex20->true_centroid();
      LIBMESH_ASSERT_FP_EQUAL(0, true_centroid(0), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(0, true_centroid(1), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(0, true_centroid(2), TOLERANCE*TOLERANCE);
#endif // LIBMESH_ENABLE_AMR
    }
  }

  void testTwistedVolume()
  {
    LOG_UNIT_TEST;

    ReplicatedMesh mesh(*TestCommWorld);

    // Build an element type that will fall back on our generic
    // quadrature-based Elem::volume()
    MeshTools::Generation::build_cube(mesh,
                                      /*nelem=*/1, /*nelem=*/1, /*nelem=*/1,
                                      /*xmin=*/-1, /*xmax=*/1,
                                      /*ymin=*/-1, /*ymax=*/1,
                                      /*zmin=*/-1, /*zmax=*/1,
                                      PRISM21);

    // Pick an element and twist it
    Elem * prism6 = mesh.elem_ptr(0);
    prism6->point(1) *= -1;
    prism6->point(1) += 2*prism6->point(0);

    // The real test here is that volume() doesn't throw
    const Real vol = prism6->volume();
    const Real gold_vol = 3+Real(5)/9;
    CPPUNIT_ASSERT_LESS(TOLERANCE, std::abs(vol-gold_vol));
  }

  void testEdge3Volume()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
    MeshTools::Generation::build_line (mesh, /*nelem=*/1, 0., 1., EDGE3);
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(3), mesh.n_nodes());

    auto edge3 = mesh.query_elem_ptr(0);
    if (!edge3) // We may be on a distributed mesh
      return;

    // Check unperturbed, straight edge case
    LIBMESH_ASSERT_FP_EQUAL(1.0, edge3->volume(), TOLERANCE*TOLERANCE);

    // Get references to the individual Edge3 nodes
    auto & middle_node = edge3->node_ref(2);
    auto & right_node = edge3->node_ref(1);

    // Check middle node perturbed in +x direction case. This should
    // not change the volume because it's still a straight line
    // element.
    middle_node = Point(0.5 + 1.e-3, 0., 0.);
    LIBMESH_ASSERT_FP_EQUAL(1.0, edge3->volume(), TOLERANCE*TOLERANCE);

    // Check middle node perturbed in -x direction case. This should
    // not change the volume because it's still a straight line
    // element.
    middle_node = Point(0.5 - 1.e-3, 0., 0.);
    LIBMESH_ASSERT_FP_EQUAL(1.0, edge3->volume(), TOLERANCE*TOLERANCE);

    // Check volume of actual curved element against pre-computed value.
    middle_node = Point(0.5, 0.25, 0.);
    right_node = Point(1., 1., 0.);
    LIBMESH_ASSERT_FP_EQUAL(1.4789428575446, edge3->volume(), TOLERANCE);

    // Compare with volume computed by base class Elem::volume() call
    // which uses quadrature.  We don't expect this to have full
    // floating point accuracy.
    middle_node = Point(0.5, 0.1, 0.);
    right_node = Point(1., 0., 0.);
    LIBMESH_ASSERT_FP_EQUAL(edge3->Elem::volume(), edge3->volume(), std::sqrt(TOLERANCE));
  }

  void testEdge3Invertible()
  {
    LOG_UNIT_TEST;

    // 1.) This is the original test which started the investigation
    // of determining invertibility.  In this test, the actual
    // midpoint of nodes 0 and 1 is 0.5*(1.100328e2 + 1.176528e2) =
    // 113.8428, so we can see that the middle node is closer to the
    // left endpoint. In this case, it is too close and the element is
    // not invertible.
    bool invertible = test_elem_invertible({
      Point(-3.566160e1, -6.690970e-1, 1.100328e2),
      Point(-3.566160e1, -6.690970e-1, 1.176528e2),
      Point(-3.566160e1, -6.690970e-1, 1.115568e2)}, EDGE3);
    CPPUNIT_ASSERT(!invertible);

    // 2.) Just like case 1, but now node 2 is at the midpoint, so
    // this case is invertible.
    invertible = test_elem_invertible({
      Point(-3.566160e1, -6.690970e-1, 1.100328e2),
      Point(-3.566160e1, -6.690970e-1, 1.176528e2),
      Point(-3.566160e1, -6.690970e-1, 113.8428)}, EDGE3);
    CPPUNIT_ASSERT(invertible);

    // 3.) Non-collinear case where the mid-edge node is "above" and "way
    // past" the right endpoint. This case is not invertible
    invertible = test_elem_invertible({Point(0, 0, 0), Point(1, 0, 0), Point(3.5, 1.5, 0)}, EDGE3);
    CPPUNIT_ASSERT(!invertible);
  }

  void testEdge4Invertible()
  {
    LOG_UNIT_TEST;

    // Reference Elem should be invertible
    {
      const Elem & edge4 = ReferenceElem::get(EDGE4);
      CPPUNIT_ASSERT(edge4.has_invertible_map());
    }

    // If node 2 goes to the left past -5/9 = -.555, the element becomes non-invertible
    {
      // x2 > -5/9, the map is still invertible
      bool invertible =
        test_elem_invertible({Point(-1, 0, 0), Point(1, 0, 0), Point(-0.5, 0, 0), Point(Real(1)/3, 0, 0)},
                  EDGE4);
      CPPUNIT_ASSERT(invertible);

      // x2 < -5/9, it is too close to x0 now
      invertible =
        test_elem_invertible({Point(-1, 0, 0), Point(1, 0, 0), Point(-0.57, 0, 0), Point(Real(1)/3, 0, 0)},
                  EDGE4);
      CPPUNIT_ASSERT(!invertible);
    }

    // If node 2 goes to the right past 5/21 ~ 0.2381, the element becomes non-invertible
    {
      // x2 < 5/21, the map should still be invertible
      bool invertible =
        test_elem_invertible({Point(-1, 0, 0), Point(1, 0, 0), Point(Real(3)/21, 0, 0), Point(Real(1)/3, 0, 0)},
                  EDGE4);
      CPPUNIT_ASSERT(invertible);

      // x2 > 5/21, x2 is too close to x3 now
      invertible =
        test_elem_invertible({Point(-1, 0, 0), Point(1, 0, 0), Point(Real(6)/21, 0, 0), Point(Real(1)/3, 0, 0)},
                  EDGE4);
      CPPUNIT_ASSERT(!invertible);
    }
  }

  void testQuad4Invertible()
  {
    LOG_UNIT_TEST;

    // Case 1: Test that rigid body rotations have no effect on the
    // invertibility of the reference element
    {
      // 1a) The reference element rotated into various different different planes.
      std::vector<Point> pts = {Point(0, 0, 0), Point(1, 0, 0), Point(1, 1, 0), Point(0, 1, 0)};
      bool invertible = test_elem_invertible(pts, QUAD4);
      CPPUNIT_ASSERT(invertible);

      // 1b) Rotate all points about x-axis by 90 degrees
      Real cost = std::cos(.5*libMesh::pi);
      Real sint = std::sin(.5*libMesh::pi);
      RealTensorValue Rx(1, 0, 0,
                         0, cost, sint,
                         0, sint, cost);

      for (auto & pt : pts)
        pt = Rx * pt;

      invertible = test_elem_invertible(pts, QUAD4);
      CPPUNIT_ASSERT(invertible);

      // 1c) Rotate all points about z-axis by 90 degrees
      RealTensorValue Rz(cost, -sint, 0,
                         sint,  cost, 0,
                         0,        0, 1);

      for (auto & pt : pts)
        pt = Rz * pt;

      invertible = test_elem_invertible(pts, QUAD4);
      CPPUNIT_ASSERT(invertible);

      // 1d) Rotate all points about y-axis by 270 degrees
      RealTensorValue Ry(cost,  0, sint,
                         0,     1, 0,
                         -sint, 0, cost);

      for (int cnt=0; cnt<3; ++cnt)
        for (auto & pt : pts)
          pt = Ry * pt;

      invertible = test_elem_invertible(pts, QUAD4);
      CPPUNIT_ASSERT(invertible);
    }

    // Case 2: Planar quad with top right vertex displaced to the position
    // (alpha, alpha). Some different cases are described below.
    // .) alpha==1: affine case, always invertible
    // .) 1/2 < alpha < 1: planar case, invertible
    // .) alpha<=1/2: planar case but node is now at center of the
    //    element, should give a zero/negative Jacobian on the displaced
    //    Node -> not invertible.
    {
      const Real alpha = .5;

      bool invertible =
        test_elem_invertible({Point(0, 0, 0), Point(1, 0, 0), Point(alpha, alpha, 0), Point(0, 1, 0)}, QUAD4);

      CPPUNIT_ASSERT(!invertible);
    }

    // Case 3) Top right corner is moved to (alpha, 1, 0). Element
    // becomes non-invertible when alpha < 0.
    {
      const Real alpha = -0.25;

      bool invertible =
        test_elem_invertible({Point(0, 0, 0), Point(1, 0, 0), Point(alpha, 1, 0), Point(0, 1, 0)}, QUAD4);

      CPPUNIT_ASSERT(!invertible);
    }

    // Case 4) Degenerate case - all 4 points at same location. This
    // zero-volume element does not have an invertible map.
    {
      const Real alpha = std::log(2);

      bool invertible =
        test_elem_invertible({Point(alpha, alpha, alpha),
                   Point(alpha, alpha, alpha),
                   Point(alpha, alpha, alpha),
                   Point(alpha, alpha, alpha)}, QUAD4);

      CPPUNIT_ASSERT(!invertible);
    }
  }

  void testQuad4AspectRatio()
  {
    LOG_UNIT_TEST;

    // Case 1: Test that rigid body rotations of a unit square
    // quadrilateral that have no effect on the quality of the
    // element.
    {
      // Construct unit square QUAD4
      std::vector<Point> pts = {Point(0, 0, 0), Point(1, 0, 0), Point(1, 1, 0), Point(0, 1, 0)};
      auto [elem, nodes] = this->construct_elem(pts, QUAD4);
      libmesh_ignore(nodes);

      // 1a) Unit square aspect ratio should be == 1
      Real aspect_ratio = elem->quality(ASPECT_RATIO);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(/*expected=*/1.0, /*actual=*/aspect_ratio, TOLERANCE);

      // 1b) Rotate all points about x-axis by 90 degrees
      Real cost = std::cos(.5*libMesh::pi);
      Real sint = std::sin(.5*libMesh::pi);
      RealTensorValue Rx(1, 0, 0,
                         0, cost, sint,
                         0, sint, cost);

      for (auto & pt : pts)
        pt = Rx * pt;

      aspect_ratio = elem->quality(ASPECT_RATIO);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(/*expected=*/1.0, /*actual=*/aspect_ratio, TOLERANCE);

      // 1c) Rotate all points about z-axis by 90 degrees
      RealTensorValue Rz(cost, -sint, 0,
                         sint,  cost, 0,
                         0,        0, 1);

      for (auto & pt : pts)
        pt = Rz * pt;

      aspect_ratio = elem->quality(ASPECT_RATIO);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(/*expected=*/1.0, /*actual=*/aspect_ratio, TOLERANCE);

      // 1d) Rotate all points about y-axis by 270 degrees
      RealTensorValue Ry(cost,  0, sint,
                         0,     1, 0,
                         -sint, 0, cost);

      for (int cnt=0; cnt<3; ++cnt)
        for (auto & pt : pts)
          pt = Ry * pt;

      aspect_ratio = elem->quality(ASPECT_RATIO);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(/*expected=*/1.0, /*actual=*/aspect_ratio, TOLERANCE);
    }

    // Case 2: Rhombus QUAD4. This case should have an aspect ratio of
    // 1/sin(theta), where theta is the acute interior angle of the
    // rhombus.
    {
      // Helper lambda function that constructs a rhombus quad with
      // interior acute angle theta.
      auto test_rhombus_quad = [this](Real theta)
      {
        Real ct = std::cos(theta);
        Real st = std::sin(theta);
        std::vector<Point> pts = {
          Point(0, 0, 0),
          Point(1, 0, 0),
          Point(1. + ct, st, 0),
          Point(     ct, st, 0)};
        auto [elem, nodes] = this->construct_elem(pts, QUAD4);
        libmesh_ignore(nodes);

        // The expected aspect ratio for the rhombus is 1/sin(theta)
        Real aspect_ratio = elem->quality(ASPECT_RATIO);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(/*expected=*/1.0/st, /*actual=*/aspect_ratio, TOLERANCE);
      };

      // 2a) Rhombus with interior angle theta=pi/6. The expected
      // aspect ratio in this case is 1/std::sin(pi/6) = 2
      test_rhombus_quad(libMesh::pi / 6);

      // 2b) Rhombus with interior angle theta=pi/3. The expected
      // aspect ratio in this case is 1/std::sin(pi/3) = 2/sqrt(3) = 1.155
      test_rhombus_quad(libMesh::pi / 3);
    }

    // Case 3) Rectangle QUAD4. The "old" and "new" aspect ratio metrics
    // return the same result for this case, whih is simply the ratio of
    // the longest to shortest side.
    {
      auto test_rectangle_quad = [this](Real a)
      {
        std::vector<Point> pts = {Point(0, 0, 0), Point(a, 0, 0), Point(a, 1, 0), Point(0, 1, 0)};
        auto [elem, nodes] = this->construct_elem(pts, QUAD4);
        libmesh_ignore(nodes);

        // The expected aspect ratio for the rectangle is "a"
        Real aspect_ratio = elem->quality(ASPECT_RATIO);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(/*expected=*/a, /*actual=*/aspect_ratio, TOLERANCE);
      };

      // 3a) Test case twice as long as it is tall
      test_rectangle_quad(2.0);

      // 3b) Test no funny numerical stuff with extremely stretched case
      test_rectangle_quad(1e6);
    }

    // Case 4) Degenerate QUAD with zero length side. In the "old"
    // aspect ratio metric we'd just get zero (as a stand-in for
    // infinity) for this case since the minimum side length is zero,
    // but using the new metric we get a non-zero value of 2.5. Thus,
    // in the new metric it is possible to compare the aspect ratios
    // of different degenerate QUAD elements rather than simply
    // assigning all such elements a quality value of 0. Degenerate
    // quadrilaterals are sometimes used as a "hack" to avoid having a
    // separate subdomain of TRIs, and (surprisingly) many finite
    // element operations just "work" on such elements.
    {
      std::vector<Point> pts = {Point(0, 0, 0), Point(1, 0, 0), Point(1, 0, 0), Point(0, 1, 0)};
      auto [elem, nodes] = this->construct_elem(pts, QUAD4);
      libmesh_ignore(nodes);

      Real aspect_ratio = elem->quality(ASPECT_RATIO);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(/*expected=*/2.5, /*actual=*/aspect_ratio, TOLERANCE);
    }

    // Case 5) Trapezoid QUAD
    {
      auto test_trapezoid_quad = [this](Real a)
      {
        // 0 <= a <= 1/2
        std::vector<Point> pts = {Point(0, 0, 0), Point(1, 0, 0), Point(1-a, 1, 0), Point(a, 1, 0)};
        auto [elem, nodes] = this->construct_elem(pts, QUAD4);
        libmesh_ignore(nodes);

        Real aspect_ratio = elem->quality(ASPECT_RATIO);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(/*expected=*/1./(1. - a), /*actual=*/aspect_ratio, TOLERANCE);
      };

      // 3a) Test "simple" trapezoid with expected aspect ratio 1.5
      test_trapezoid_quad(1./3);

      // 3b) Test "degenerate" trapezoid with one zero length base
      test_trapezoid_quad(0.5);
    }
  }

  void testTri3AspectRatio()
  {
    LOG_UNIT_TEST;

    // Case 1) Reference TRI3
    {
      std::vector<Point> pts = {Point(0, 0, 0), Point(1, 0, 0), Point(0, 1, 0)};
      auto [elem, nodes] = this->construct_elem(pts, TRI3);
      libmesh_ignore(nodes);

      // Compute the aspect ratio for the reference Tri
      // The expected value is ~ 1.44338
      Real aspect_ratio = elem->quality(ASPECT_RATIO);
      // libMesh::out << "Unit TRI3 aspect ratio = " << aspect_ratio << std::endl;
      CPPUNIT_ASSERT_DOUBLES_EQUAL(/*expected=*/Real(5)/2/std::sqrt(Real(3)), /*actual=*/aspect_ratio, TOLERANCE);
    }

    // Case 2) Equilateral TRI3
    {
      std::vector<Point> pts = {Point(0, 0, 0), Point(1, 0, 0), Point(0.5, 0.5*std::sqrt(3), 0)};
      auto [elem, nodes] = this->construct_elem(pts, TRI3);
      libmesh_ignore(nodes);

      // Compute the aspect ratio for the reference Tri
      // The expected value is 1.0
      Real aspect_ratio = elem->quality(ASPECT_RATIO);
      // libMesh::out << "Equilateral TRI3 aspect ratio = " << aspect_ratio << std::endl;
      CPPUNIT_ASSERT_DOUBLES_EQUAL(/*expected=*/Real(1), /*actual=*/aspect_ratio, TOLERANCE);
    }

    // Case 3) Reference TRI3 with one leg length = L >> 1
    {
      Real L = 10.;
      std::vector<Point> pts = {Point(0, 0, 0), Point(L, 0, 0), Point(0, 1, 0)};
      auto [elem, nodes] = this->construct_elem(pts, TRI3);
      libmesh_ignore(nodes);

      // Compute the aspect ratio for the reference Tri
      // The expected value is ~ 11.5759
      Real aspect_ratio = elem->quality(ASPECT_RATIO);
      // libMesh::out << "TRI3 with leg length L = " << L << ", aspect ratio = " << aspect_ratio << std::endl;
      CPPUNIT_ASSERT_DOUBLES_EQUAL(/*expected=*/Real(1)/std::sqrt(Real(3)) * (Real(2)*L + Real(1)/(2*L)), /*actual=*/aspect_ratio, TOLERANCE);
    }
  }

  void testQuad4Warpage()
  {
    LOG_UNIT_TEST;

    // Case 1: Test that rigid body rotations of a unit square
    // quadrilateral that have no effect on the quality of the
    // element.
    {
      // Construct unit square QUAD4
      std::vector<Point> pts = {Point(0, 0, 0), Point(1, 0, 0), Point(1, 1, 0), Point(0, 1, 0)};
      auto [elem, nodes] = this->construct_elem(pts, QUAD4);
      libmesh_ignore(nodes);

      // 1a) Any flat element should have warp == 1
      Real warpage = elem->quality(WARP);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(/*expected=*/1.0, /*actual=*/warpage, TOLERANCE);

      // 1b) Rotate all points about x-axis by 90 degrees
      Real cost = std::cos(.5*libMesh::pi);
      Real sint = std::sin(.5*libMesh::pi);
      RealTensorValue Rx(1, 0, 0,
                         0, cost, sint,
                         0, sint, cost);

      for (auto & pt : pts)
        pt = Rx * pt;

      warpage = elem->quality(WARP);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(/*expected=*/1.0, /*actual=*/warpage, TOLERANCE);

      // 1c) Rotate all points about z-axis by 90 degrees
      RealTensorValue Rz(cost, -sint, 0,
                         sint,  cost, 0,
                         0,        0, 1);

      for (auto & pt : pts)
        pt = Rz * pt;

      warpage = elem->quality(WARP);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(/*expected=*/1.0, /*actual=*/warpage, TOLERANCE);

      // 1d) Rotate all points about y-axis by 270 degrees
      RealTensorValue Ry(cost,  0, sint,
                         0,     1, 0,
                         -sint, 0, cost);

      for (int cnt=0; cnt<3; ++cnt)
        for (auto & pt : pts)
          pt = Ry * pt;

      warpage = elem->quality(WARP);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(/*expected=*/1.0, /*actual=*/warpage, TOLERANCE);
    }

    // Case 2: Unit square quadrilateral with Node 2 displaced by a distance h in the z-direction.
    {
      auto test_warped_quad = [this](Real h)
      {
        std::vector<Point> pts = {Point(0, 0, 0), Point(1, 0, 0), Point(1, 1, h), Point(0, 1, 0)};
        auto [elem, nodes] = this->construct_elem(pts, QUAD4);
        libmesh_ignore(nodes);

        Real warpage = elem->quality(WARP);
        // libMesh::out << "QUAD with node 3 displaced by h = "
        //              << h << " in the z-direction , warpage = " << warpage
        //              << std::endl;
        CPPUNIT_ASSERT_DOUBLES_EQUAL(/*expected=*/Real(1) / (h*h + 1), /*actual=*/warpage, TOLERANCE);
      };

      // For h = 0.1, the expected warpage value is ~ 0.991, so not
      // much different than the value of 1.0 for a flat element.
      test_warped_quad(0.1);

      // For h = 0.3, the expected warpage value is ~ 0.917431, which
      // is pretty close to the lower bound (0.9) of "acceptable"
      // warpage, as suggested in the Cubit manual.
      test_warped_quad(0.3);
    }
  }

  void testQuad4MinMaxAngle()
  {
    LOG_UNIT_TEST;

    // Case 1: Test that rigid body rotations of a unit square
    // quadrilateral that have no effect on the quality of the
    // element.
    {
      // Construct unit square QUAD4
      std::vector<Point> pts = {Point(0, 0, 0), Point(1, 0, 0), Point(1, 1, 0), Point(0, 1, 0)};
      auto [elem, nodes] = this->construct_elem(pts, QUAD4);
      libmesh_ignore(nodes);

      // 1a) Reference Elem should have min angle == max angle == pi/2
      {
        Real min_angle = elem->quality(MIN_ANGLE);
        Real max_angle = elem->quality(MAX_ANGLE);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(/*expected=*/90, /*actual=*/min_angle, TOLERANCE);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(/*expected=*/90, /*actual=*/max_angle, TOLERANCE);
      }

      // 1b) Rotate all points about x-axis by 90 degrees
      Real cost = std::cos(.5*libMesh::pi);
      Real sint = std::sin(.5*libMesh::pi);
      RealTensorValue Rx(1, 0, 0,
                         0, cost, sint,
                         0, sint, cost);

      for (auto & pt : pts)
        pt = Rx * pt;

      {
        Real min_angle = elem->quality(MIN_ANGLE);
        Real max_angle = elem->quality(MAX_ANGLE);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(/*expected=*/90, /*actual=*/min_angle, TOLERANCE);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(/*expected=*/90, /*actual=*/max_angle, TOLERANCE);
      }

      // 1c) Rotate all points about z-axis by 90 degrees
      RealTensorValue Rz(cost, -sint, 0,
                         sint,  cost, 0,
                         0,        0, 1);

      for (auto & pt : pts)
        pt = Rz * pt;

      {
        Real min_angle = elem->quality(MIN_ANGLE);
        Real max_angle = elem->quality(MAX_ANGLE);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(/*expected=*/90, /*actual=*/min_angle, TOLERANCE);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(/*expected=*/90, /*actual=*/max_angle, TOLERANCE);
      }

      // 1d) Rotate all points about y-axis by 270 degrees
      RealTensorValue Ry(cost,  0, sint,
                         0,     1, 0,
                         -sint, 0, cost);

      for (int cnt=0; cnt<3; ++cnt)
        for (auto & pt : pts)
          pt = Ry * pt;

      {
        Real min_angle = elem->quality(MIN_ANGLE);
        Real max_angle = elem->quality(MAX_ANGLE);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(/*expected=*/90, /*actual=*/min_angle, TOLERANCE);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(/*expected=*/90, /*actual=*/max_angle, TOLERANCE);
      }
    }

    // Case 2: Rhombus QUAD4. This case should have an interior min
    // angle of theta and an interior max angle of pi-theta, where
    // "theta" is the specified amount that we "sheared" the element
    // by on creation
    {
      // Helper lambda function that constructs a rhombus quad with
      // interior acute angle theta.
      auto test_rhombus_quad = [this](Real theta)
      {
        Real ct = std::cos(theta);
        Real st = std::sin(theta);
        std::vector<Point> pts = {
          Point(0, 0, 0),
          Point(1, 0, 0),
          Point(1. + ct, st, 0),
          Point(     ct, st, 0)};
        auto [elem, nodes] = this->construct_elem(pts, QUAD4);
        libmesh_ignore(nodes);

        Real min_angle = elem->quality(MIN_ANGLE);
        Real max_angle = elem->quality(MAX_ANGLE);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(/*expected=*/Real(180) / libMesh::pi * theta,                 /*actual=*/min_angle, TOLERANCE);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(/*expected=*/Real(180) / libMesh::pi * (libMesh::pi - theta), /*actual=*/max_angle, TOLERANCE);
      };

      // 2a) Rhombus with interior angle theta=pi/6.
      test_rhombus_quad(libMesh::pi / 6);

      // 2b) Rhombus with interior angle theta=pi/3.
      test_rhombus_quad(libMesh::pi / 3);
    }
  }

  void testQuad4Jacobian()
  {
    LOG_UNIT_TEST;

    // Rhombus QUAD4. This case should have a JACOBIAN and
    // SCALED_JACOBIAN value of sin(theta), where "theta" is the
    // specified amount that we "sheared" the element by on creation.
    // The JACOBIAN and SCALED_JACOBIAN are the same for this element
    // because the edge lengths are all = 1.
    {
      // Helper lambda function that constructs a rhombus quad with
      // interior acute angle theta.
      auto test_rhombus_quad = [this](Real theta)
      {
        Real ct = std::cos(theta);
        Real st = std::sin(theta);
        std::vector<Point> pts = {
          Point(0, 0, 0),
          Point(1, 0, 0),
          Point(1. + ct, st, 0),
          Point(     ct, st, 0)};
        auto [elem, nodes] = this->construct_elem(pts, QUAD4);
        libmesh_ignore(nodes);

        Real jac = elem->quality(JACOBIAN);
        Real scaled_jac = elem->quality(SCALED_JACOBIAN);

        // Debugging
        // libMesh::out << "jac = " << jac << std::endl;
        // libMesh::out << "scaled_jac = " << scaled_jac << std::endl;

        CPPUNIT_ASSERT_DOUBLES_EQUAL(/*expected=*/std::abs(std::sin(theta)), /*actual=*/jac, TOLERANCE);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(/*expected=*/std::abs(std::sin(theta)), /*actual=*/scaled_jac, TOLERANCE);
      };

      // 2a) Rhombus with interior angle theta=pi/6.
      test_rhombus_quad(libMesh::pi / 6);

      // 2b) Rhombus with interior angle theta=pi/3.
      test_rhombus_quad(libMesh::pi / 3);
    }
  }

  void testTet4DihedralAngle()
  {
    LOG_UNIT_TEST;

    // Construct a tetrahedron whose projection into the x-y plane looks like the unit square,
    // and where the apex node is just slightly out of the x-y plane. This element should have
    // reasonable MIN,MAX_ANGLE values but MIN,MAX_DIHEDRAL angles of nearly 0. This is an
    // example that demonstrates one should not solely use edge angles to determine Elem quality
    // in 3D.
    std::vector<Point> pts = {Point(0, 0, 0), Point(1, 0, 0), Point(0, 1, 0), Point(1, 1, 0.01)};
    auto [elem, nodes] = this->construct_elem(pts, TET4);
    libmesh_ignore(nodes);

    Real min_angle = elem->quality(MIN_ANGLE);
    Real max_angle = elem->quality(MAX_ANGLE);
    Real min_dihedral_angle = elem->quality(MIN_DIHEDRAL_ANGLE);
    Real max_dihedral_angle = elem->quality(MAX_DIHEDRAL_ANGLE);

    // Debugging
    // libMesh::out << "Squashed Tet4 min_angle = " << min_angle << " degrees" << std::endl;
    // libMesh::out << "Squashed Tet4 max_angle = " << max_angle << " degrees"  << std::endl;
    // libMesh::out << "Squashed Tet4 min_dihedral_angle = " << min_dihedral_angle << " degrees" << std::endl;
    // libMesh::out << "Squashed Tet4 max_dihedral_angle = " << max_dihedral_angle << " degrees" << std::endl;

    CPPUNIT_ASSERT_DOUBLES_EQUAL(/*expected=*/44.9985676771277, /*actual=*/min_angle, TOLERANCE);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(/*expected=*/90, /*actual=*/max_angle, TOLERANCE);

    // Assert that both min and max dihedral angles are less than 1 degree (a very low quality element)
    CPPUNIT_ASSERT_LESS(1.0, min_dihedral_angle);
    CPPUNIT_ASSERT_LESS(1.0, max_dihedral_angle);
  }

  void testTet4Jacobian()
  {
    LOG_UNIT_TEST;

    {
      // Same element as in testTet4DihedralAngle(). Here we verify that
      // the element is low quality since the SCALED_JACOBIAN is < O(h)
      // as h -> 0, and h can be arbitrarily small in the squashed element.
      Real h = 0.01;
      std::vector<Point> pts = {Point(0, 0, 0), Point(1, 0, 0), Point(0, 1, 0), Point(1, 1, h)};
      auto [elem, nodes] = this->construct_elem(pts, TET4);
      libmesh_ignore(nodes);

      Real jac = elem->quality(JACOBIAN);
      Real scaled_jac = elem->quality(SCALED_JACOBIAN);

      // Debugging
      // libMesh::out << "Squashed Tet4 jac = " << jac << std::endl;
      // libMesh::out << "Squashed Tet4 scaled_jac = " << scaled_jac << std::endl;

      CPPUNIT_ASSERT_DOUBLES_EQUAL(/*expected=*/h, /*actual=*/jac, TOLERANCE);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(/*expected=*/h/std::sqrt(h*h + 2), /*actual=*/scaled_jac, TOLERANCE);
    }

    {
      // A representative Tet generated by calling build_cube(). This
      // element has minimum interior edge angle of approximately
      // 35.26 degrees, so it it is not a simple refinement of the
      // reference element. It is not located at the origin, but still
      // has a simple to compute value for the JACOBIAN and
      // SCALED_JACOBIAN quality metrics.
      std::vector<Point> pts = {
        Point(2.5, 2.5, 4.5),
        Point(3.5, 1.5, 5.5),
        Point(1.5, 1.5, 5.5),
        Point(2.5, 2.5, 5.5)
      };
      auto [elem, nodes] = this->construct_elem(pts, TET4);
      libmesh_ignore(nodes);

      Real jac = elem->quality(JACOBIAN);
      Real scaled_jac = elem->quality(SCALED_JACOBIAN);

      // Debugging
      // libMesh::out << "build_cube() Tet4 jac = " << jac << std::endl;
      // libMesh::out << "build_cube() Tet4 scaled_jac = " << scaled_jac << std::endl;

      CPPUNIT_ASSERT_DOUBLES_EQUAL(/*expected=*/2,                          /*actual=*/jac, TOLERANCE);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(/*expected=*/Real(1)/std::sqrt(Real(6)), /*actual=*/scaled_jac, TOLERANCE);
    }
  }


protected:

  // Helper function that is called by test_elem_invertible() to build an Elem
  // of the requested elem_type from the provided Points. Note: the
  // Nodes which are constructed in order to construct the Elem are
  // also returned since
  std::pair<std::unique_ptr<Elem>, std::vector<std::unique_ptr<Node>>>
  construct_elem(const std::vector<Point> & pts,
                 ElemType elem_type)
  {
    const unsigned int n_points = pts.size();

    // Create Nodes
    std::vector<std::unique_ptr<Node>> nodes(n_points);
    for (unsigned int i=0; i<n_points; i++)
      nodes[i] = Node::build(pts[i], /*id*/ i);

    // Create Elem, assign nodes
    std::unique_ptr<Elem> elem = Elem::build(elem_type, /*parent*/ nullptr);

    // Make sure we were passed consistent input to build this type of Elem
    libmesh_error_msg_if(elem->n_nodes() != n_points,
                         "Wrong number of points "
                         << n_points
                         << " provided to build a "
                         << Utility::enum_to_string(elem_type));

    for (unsigned int i=0; i<n_points; i++)
      elem->set_node(i) = nodes[i].get();

    // Return Elem and Nodes we created
    return std::make_pair(std::move(elem), std::move(nodes));
  }

  // Helper function that builds the specified type of Elem from a
  // vector of Points and returns the value of has_invertible_map()
  // for that Elem.
  bool test_elem_invertible(const std::vector<Point> & pts,
                            ElemType elem_type)
  {
    // Construct Elem of desired type
    auto [elem, nodes] = this->construct_elem(pts, elem_type);
    libmesh_ignore(nodes);

    // Return whether or not this Elem has an invertible map
    return elem->has_invertible_map();
  }

  // Helper function for testing true_centroid() and volume() implementations for
  // 3D elements
  void test_true_centroid_and_volume(ElemType elem_type)
  {
    // Check that derived class true_centroid() gives same result as
    // the base class Elem::true_centroid() on a 2x2x2 mesh of 10%
    // distorted elements.
    {
      ReplicatedMesh mesh(*TestCommWorld);

      MeshTools::Generation::build_cube(mesh,
                                        /*nelem=*/2, /*nelem=*/2, /*nelem=*/2,
                                        /*xmin=*/-1, /*xmax=*/1,
                                        /*ymin=*/-1, /*ymax=*/1,
                                        /*zmin=*/-1, /*zmax=*/1,
                                        elem_type);

      MeshTools::Modification::distort(mesh,
                                       /*factor=*/0.1,
                                       /*perturb_boundary=*/false);

      for (const auto & elem : mesh.element_ptr_range())
        {
          Point derived_centroid = elem->true_centroid();
          Point base_centroid = elem->Elem::true_centroid();

          // Debugging: check results in detail
          // auto flags = libMesh::out.flags();
          // libMesh::out << std::scientific << std::setprecision(16);
          // libMesh::out << "derived_centroid = " << derived_centroid << std::endl;
          // libMesh::out << "base_centroid = " << base_centroid << std::endl;
          // libMesh::out.flags(flags);

          CPPUNIT_ASSERT(derived_centroid.absolute_fuzzy_equals(base_centroid, TOLERANCE*TOLERANCE));

          // Make sure that base class and "optimized" routines for computing the cell volume agree
          Real derived_volume = elem->volume();
          Real base_volume = elem->Elem::volume();
          LIBMESH_ASSERT_FP_EQUAL(base_volume, derived_volume, TOLERANCE*TOLERANCE);
        }
    }
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION( VolumeTest );
