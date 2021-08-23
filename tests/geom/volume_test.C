// libmesh includes
#include <libmesh/elem.h>
#include <libmesh/enum_elem_type.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_modification.h>
#include <libmesh/mesh.h>
#include <libmesh/reference_elem.h>
#include <libmesh/node.h>
#include <libmesh/enum_to_string.h>
#include <libmesh/tensor_value.h>

// unit test includes
#include "test_comm.h"
#include "libmesh_cppunit.h"

// C++ includes
#include <iomanip>

using namespace libMesh;

class VolumeTest : public CppUnit::TestCase
{

public:
  CPPUNIT_TEST_SUITE( VolumeTest );
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

  void testHex8TrueCentroid() { test_true_centroid_and_volume(HEX8); }
  void testPrism6TrueCentroid() { test_true_centroid_and_volume(PRISM6); }

  void testHex20PLevelTrueCentroid()
  {
    // Test that Elem base class true_centroid() implementation works
    // for an elevated p_level HEX20
    {
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
    }
  }

  void testEdge3Volume()
  {
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
    // 1.) This is the original test which started the investigation
    // of determining invertibility.  In this test, the actual
    // midpoint of nodes 0 and 1 is 0.5*(1.100328e2 + 1.176528e2) =
    // 113.8428, so we can see that the middle node is closer to the
    // left endpoint. In this case, it is too close and the element is
    // not invertible.
    bool invertible = test_elem({
      Point(-3.566160e1, -6.690970e-1, 1.100328e2),
      Point(-3.566160e1, -6.690970e-1, 1.176528e2),
      Point(-3.566160e1, -6.690970e-1, 1.115568e2)}, EDGE3);
    CPPUNIT_ASSERT(!invertible);

    // 2.) Just like case 1, but now node 2 is at the midpoint, so
    // this case is invertible.
    invertible = test_elem({
      Point(-3.566160e1, -6.690970e-1, 1.100328e2),
      Point(-3.566160e1, -6.690970e-1, 1.176528e2),
      Point(-3.566160e1, -6.690970e-1, 113.8428)}, EDGE3);
    CPPUNIT_ASSERT(invertible);

    // 3.) Non-collinear case where the mid-edge node is "above" and "way
    // past" the right endpoint. This case is not invertible
    invertible = test_elem({Point(0, 0, 0), Point(1, 0, 0), Point(3.5, 1.5, 0)}, EDGE3);
    CPPUNIT_ASSERT(!invertible);
  }

  void testEdge4Invertible()
  {
    // Reference Elem should be invertible
    {
      const Elem & edge4 = ReferenceElem::get(EDGE4);
      CPPUNIT_ASSERT(edge4.has_invertible_map());
    }

    // If node 2 goes to the left past -5/9 = -.555, the element becomes non-invertible
    {
      // x2 > -5/9, the map is still invertible
      bool invertible =
        test_elem({Point(-1, 0, 0), Point(1, 0, 0), Point(-0.5, 0, 0), Point(Real(1)/3, 0, 0)},
                  EDGE4);
      CPPUNIT_ASSERT(invertible);

      // x2 < -5/9, it is too close to x0 now
      invertible =
        test_elem({Point(-1, 0, 0), Point(1, 0, 0), Point(-0.57, 0, 0), Point(Real(1)/3, 0, 0)},
                  EDGE4);
      CPPUNIT_ASSERT(!invertible);
    }

    // If node 2 goes to the right past 5/21 ~ 0.2381, the element becomes non-invertible
    {
      // x2 < 5/21, the map should still be invertible
      bool invertible =
        test_elem({Point(-1, 0, 0), Point(1, 0, 0), Point(Real(3)/21, 0, 0), Point(Real(1)/3, 0, 0)},
                  EDGE4);
      CPPUNIT_ASSERT(invertible);

      // x2 > 5/21, x2 is too close to x3 now
      invertible =
        test_elem({Point(-1, 0, 0), Point(1, 0, 0), Point(Real(6)/21, 0, 0), Point(Real(1)/3, 0, 0)},
                  EDGE4);
      CPPUNIT_ASSERT(!invertible);
    }
  }

  void testQuad4Invertible()
  {
    // Case 1: Test that rigid body rotations have no effect on the
    // invertibility of the reference element
    {
      // 1a) The reference element rotated into various different different planes.
      std::vector<Point> pts = {Point(0, 0, 0), Point(1, 0, 0), Point(1, 1, 0), Point(0, 1, 0)};
      bool invertible = test_elem(pts, QUAD4);
      CPPUNIT_ASSERT(invertible);

      // 1b) Rotate all points about x-axis by 90 degrees
      Real cost = std::cos(.5*libMesh::pi);
      Real sint = std::sin(.5*libMesh::pi);
      RealTensorValue Rx(1, 0, 0,
                         0, cost, sint,
                         0, sint, cost);

      for (auto & pt : pts)
        pt = Rx * pt;

      invertible = test_elem(pts, QUAD4);
      CPPUNIT_ASSERT(invertible);

      // 1c) Rotate all points about z-axis by 90 degrees
      RealTensorValue Rz(cost, -sint, 0,
                         sint,  cost, 0,
                         0,        0, 1);

      for (auto & pt : pts)
        pt = Rz * pt;

      invertible = test_elem(pts, QUAD4);
      CPPUNIT_ASSERT(invertible);

      // 1d) Rotate all points about y-axis by 270 degrees
      RealTensorValue Ry(cost,  0, sint,
                         0,     1, 0,
                         -sint, 0, cost);

      for (int cnt=0; cnt<3; ++cnt)
        for (auto & pt : pts)
          pt = Ry * pt;

      invertible = test_elem(pts, QUAD4);
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
        test_elem({Point(0, 0, 0), Point(1, 0, 0), Point(alpha, alpha, 0), Point(0, 1, 0)}, QUAD4);

      CPPUNIT_ASSERT(!invertible);
    }

    // Case 3) Top right corner is moved to (alpha, 1, 0). Element
    // becomes non-invertible when alpha < 0.
    {
      const Real alpha = -0.25;

      bool invertible =
        test_elem({Point(0, 0, 0), Point(1, 0, 0), Point(alpha, 1, 0), Point(0, 1, 0)}, QUAD4);

      CPPUNIT_ASSERT(!invertible);
    }

    // Case 4) Degenerate case - all 4 points at same location. This
    // zero-volume element does not have an invertible map.
    {
      const Real alpha = std::log(2);

      bool invertible =
        test_elem({Point(alpha, alpha, alpha),
                   Point(alpha, alpha, alpha),
                   Point(alpha, alpha, alpha),
                   Point(alpha, alpha, alpha)}, QUAD4);

      CPPUNIT_ASSERT(!invertible);
    }
  }

protected:

  // Helper function that builds the specified type of Elem from a
  // vector of Points and returns the value of has_invertible_map()
  // for that Elem.
  bool test_elem(const std::vector<Point> & pts,
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
