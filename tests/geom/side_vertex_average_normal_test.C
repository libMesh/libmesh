// libmesh includes
#include <libmesh/cell_c0polyhedron.h>
#include <libmesh/elem.h>
#include <libmesh/enum_elem_type.h>
#include <libmesh/face_c0polygon.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_modification.h>
#include <libmesh/mesh.h>
#include <libmesh/reference_elem.h>
#include <libmesh/replicated_mesh.h>
#include <libmesh/node.h>
#include <libmesh/enum_to_string.h>
#include <libmesh/tensor_value.h>
#include <libmesh/enum_elem_quality.h>
#include <libmesh/fe_base.h>
#include <libmesh/quadrature_gauss.h>

// unit test includes
#include "test_comm.h"
#include "libmesh_cppunit.h"

// C++ includes
#include <iomanip>

using namespace libMesh;

class SideVertexAverageNormalTest : public CppUnit::TestCase
{

public:
  LIBMESH_CPPUNIT_TEST_SUITE( SideVertexAverageNormalTest );
  CPPUNIT_TEST( testEdge2 );
  CPPUNIT_TEST( testTri3 );
  CPPUNIT_TEST( testQuad4 );
  CPPUNIT_TEST( testPyramid5 );
  CPPUNIT_TEST( testPrism6 );
  CPPUNIT_TEST( testHex8 );
  CPPUNIT_TEST( testC0Polygon );
  CPPUNIT_TEST( testC0Polyhedron );
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp()
  {
  }

  void tearDown()
  {
  }

  void testEdge2()
  {
    LOG_UNIT_TEST;

    {
      // Reference
      const Elem & edge2 = ReferenceElem::get(EDGE2);
      const Point n1 = edge2.side_vertex_average_normal(0);
      LIBMESH_ASSERT_FP_EQUAL(-1, n1(0), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(0, n1(1), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(0, n1(2), TOLERANCE*TOLERANCE);
      const Point n2 = edge2.side_vertex_average_normal(1);
      LIBMESH_ASSERT_FP_EQUAL(1, n2(0), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(0, n2(1), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(0, n2(2), TOLERANCE*TOLERANCE);
    }
    {
      // Oriented
      std::vector<Point> pts = {Point(1, 0, 0), Point(1, 3, 0)};
      auto [edge2, nodes] = this->construct_elem(pts, EDGE2);
      const Point n1 = edge2->side_vertex_average_normal(0);
      LIBMESH_ASSERT_FP_EQUAL(0, n1(0), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(-1, n1(1), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(0, n1(2), TOLERANCE*TOLERANCE);
      const Point n2 = edge2->side_vertex_average_normal(1);
      LIBMESH_ASSERT_FP_EQUAL(0, n2(0), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(1, n2(1), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(0, n2(2), TOLERANCE*TOLERANCE);
    }
  }

  void testTri3()
  {
    LOG_UNIT_TEST;

    {
      // Reference
      const Elem & tri3 = ReferenceElem::get(TRI3);
      const Point n1 = tri3.side_vertex_average_normal(0);
      LIBMESH_ASSERT_FP_EQUAL(0, n1(0), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(-1, n1(1), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(0, n1(2), TOLERANCE*TOLERANCE);
      const Point n2 = tri3.side_vertex_average_normal(1);
      LIBMESH_ASSERT_FP_EQUAL(sqrt(2) / 2, n2(0), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(sqrt(2) / 2, n2(1), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(0, n2(2), TOLERANCE*TOLERANCE);
      const Point n3 = tri3.side_vertex_average_normal(2);
      LIBMESH_ASSERT_FP_EQUAL(-1, n3(0), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(0, n3(1), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(0, n3(2), TOLERANCE*TOLERANCE);
    }
    {
      // General shape
      std::vector<Point> pts = {Point(1, 0, 0), Point(1, 1, 0), Point(0, 3, 1)};
      auto [tri3, nodes] = this->construct_elem(pts, TRI3);
      const Point n1 = tri3->side_vertex_average_normal(0);
      LIBMESH_ASSERT_FP_EQUAL(sqrt(2) / 2, n1(0), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(0, n1(1), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(-sqrt(2) / 2, n1(2), TOLERANCE*TOLERANCE);
      const Point n2 = tri3->side_vertex_average_normal(1);
      LIBMESH_ASSERT_FP_EQUAL(1. / sqrt(3), n2(0), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(1. / sqrt(3), n2(1), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(-1. / sqrt(3), n2(2), TOLERANCE*TOLERANCE);
      const Point n3 = tri3->side_vertex_average_normal(2);
      LIBMESH_ASSERT_FP_EQUAL(-3. / sqrt(22), n3(0), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(-2. / sqrt(22), n3(1), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(3. / sqrt(22), n3(2), TOLERANCE*TOLERANCE);
    }
  }

  void testQuad4()
  {
    LOG_UNIT_TEST;

    {
      // Reference
      const Elem & quad4 = ReferenceElem::get(QUAD4);
      const Point n1 = quad4.side_vertex_average_normal(0);
      LIBMESH_ASSERT_FP_EQUAL(0, n1(0), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(-1, n1(1), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(0, n1(2), TOLERANCE*TOLERANCE);
      const Point n2 = quad4.side_vertex_average_normal(1);
      LIBMESH_ASSERT_FP_EQUAL(1, n2(0), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(0, n2(1), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(0, n2(2), TOLERANCE*TOLERANCE);
      const Point n3 = quad4.side_vertex_average_normal(2);
      LIBMESH_ASSERT_FP_EQUAL(0, n3(0), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(1, n3(1), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(0, n3(2), TOLERANCE*TOLERANCE);
      const Point n4 = quad4.side_vertex_average_normal(3);
      LIBMESH_ASSERT_FP_EQUAL(-1, n4(0), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(0, n4(1), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(0, n4(2), TOLERANCE*TOLERANCE);
    }

    {
      // Planar, general shape
      std::vector<Point> pts = {Point(1, 0, 0), Point(1, 3, 0), Point(-1, -1, 0), Point(0, -1, 0)};
      auto [quad4, nodes] = this->construct_elem(pts, QUAD4);
      const Point n1 = quad4->side_vertex_average_normal(0);
      LIBMESH_ASSERT_FP_EQUAL(1, n1(0), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(0, n1(1), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(0, n1(2), TOLERANCE*TOLERANCE);
      const Point n2 = quad4->side_vertex_average_normal(1);
      LIBMESH_ASSERT_FP_EQUAL(-2 / sqrt(5), n2(0), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(1 / sqrt(5), n2(1), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(0, n2(2), TOLERANCE*TOLERANCE);
      const Point n3 = quad4->side_vertex_average_normal(2);
      LIBMESH_ASSERT_FP_EQUAL(0, n3(0), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(-1, n3(1), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(0, n3(2), TOLERANCE*TOLERANCE);
      const Point n4 = quad4->side_vertex_average_normal(3);
      LIBMESH_ASSERT_FP_EQUAL(1 / sqrt(2), n4(0), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(-1 / sqrt(2), n4(1), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(0, n4(2), TOLERANCE*TOLERANCE);
    }

    {
      // Non-planar
      std::vector<Point> pts = {Point(0, 0, 0), Point(1, 0, 1), Point(1, 1, 0), Point(0, 1, 1)};
      auto [quad4, nodes] = this->construct_elem(pts, QUAD4);
      const Point n1 = quad4->side_vertex_average_normal(0);
      LIBMESH_ASSERT_FP_EQUAL(0, n1(0), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(-1, n1(1), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(0, n1(2), TOLERANCE*TOLERANCE);
      const Point n2 = quad4->side_vertex_average_normal(1);
      LIBMESH_ASSERT_FP_EQUAL(1, n2(0), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(0, n2(1), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(0, n2(2), TOLERANCE*TOLERANCE);
      const Point n3 = quad4->side_vertex_average_normal(2);
      LIBMESH_ASSERT_FP_EQUAL(0, n3(0), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(1, n3(1), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(0, n3(2), TOLERANCE*TOLERANCE);
      const Point n4 = quad4->side_vertex_average_normal(3);
      LIBMESH_ASSERT_FP_EQUAL(-1, n4(0), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(0, n4(1), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(0, n4(2), TOLERANCE*TOLERANCE);
    }

    {
      // Non-planar, general
      std::vector<Point> pts = {Point(0, 0, 0), Point(1, 0, 3), Point(1, 1, 0), Point(0, 2, 1)};
      auto [quad4, nodes] = this->construct_elem(pts, QUAD4);
      const std::unique_ptr<const Elem> face = quad4->build_side_ptr(0);
      std::unique_ptr<libMesh::FEBase> fe(libMesh::FEBase::build(2, libMesh::FEType(1)));
      libMesh::QGauss qface(1, libMesh::CONSTANT);
      fe->attach_quadrature_rule(&qface);
      const std::vector<Point> & normals = fe->get_normals();
      for (const auto s : make_range(quad4->n_sides()))
      {
        const std::unique_ptr<const Elem> face = quad4->build_side_ptr(s);
        fe->attach_quadrature_rule(&qface);
        fe->reinit(quad4.get(), s, TOLERANCE);
        const Point n1 = quad4->side_vertex_average_normal(s);
        LIBMESH_ASSERT_FP_EQUAL(normals[0](0), n1(0), TOLERANCE*TOLERANCE);
        LIBMESH_ASSERT_FP_EQUAL(normals[0](1), n1(1), TOLERANCE*TOLERANCE);
        LIBMESH_ASSERT_FP_EQUAL(normals[0](2), n1(2), TOLERANCE*TOLERANCE);
      }
    }
  }

  void testPyramid5()
  {
    LOG_UNIT_TEST;

    {
      // Reference
      const Elem & pyr5 = ReferenceElem::get(PYRAMID5);
      const Point n1 = pyr5.side_vertex_average_normal(0);
      LIBMESH_ASSERT_FP_EQUAL(0, n1(0), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(-sqrt(2) / 2, n1(1), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(sqrt(2) / 2, n1(2), TOLERANCE*TOLERANCE);
      const Point n2 = pyr5.side_vertex_average_normal(1);
      LIBMESH_ASSERT_FP_EQUAL(sqrt(2) / 2, n2(0), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(0, n2(1), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(sqrt(2) / 2, n2(2), TOLERANCE*TOLERANCE);
      const Point n3 = pyr5.side_vertex_average_normal(2);
      LIBMESH_ASSERT_FP_EQUAL(0, n3(0), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(sqrt(2) / 2, n3(1), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(sqrt(2) / 2, n3(2), TOLERANCE*TOLERANCE);
      const Point n4 = pyr5.side_vertex_average_normal(3);
      LIBMESH_ASSERT_FP_EQUAL(-sqrt(2) / 2, n4(0), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(0, n4(1), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(sqrt(2) / 2, n4(2), TOLERANCE*TOLERANCE);
      const Point n5 = pyr5.side_vertex_average_normal(4);
      LIBMESH_ASSERT_FP_EQUAL(0, n5(0), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(0, n5(1), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(-1, n5(2), TOLERANCE*TOLERANCE);
    }
  }

  void testPrism6()
  {
    LOG_UNIT_TEST;

    {
      // Reference
      const Elem & pri6 = ReferenceElem::get(PRISM6);
      const Point n1 = pri6.side_vertex_average_normal(0);
      LIBMESH_ASSERT_FP_EQUAL(0, n1(0), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(0, n1(1), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(-1, n1(2), TOLERANCE*TOLERANCE);
      const Point n2 = pri6.side_vertex_average_normal(1);
      LIBMESH_ASSERT_FP_EQUAL(0, n2(0), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(-1, n2(1), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(0, n2(2), TOLERANCE*TOLERANCE);
      const Point n3 = pri6.side_vertex_average_normal(2);
      LIBMESH_ASSERT_FP_EQUAL(sqrt(2) / 2, n3(0), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(sqrt(2) / 2, n3(1), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(0, n3(2), TOLERANCE*TOLERANCE);
      const Point n4 = pri6.side_vertex_average_normal(3);
      LIBMESH_ASSERT_FP_EQUAL(-1, n4(0), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(0, n4(1), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(0, n4(2), TOLERANCE*TOLERANCE);
      const Point n5 = pri6.side_vertex_average_normal(4);
      LIBMESH_ASSERT_FP_EQUAL(0, n5(0), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(0, n5(1), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(1, n5(2), TOLERANCE*TOLERANCE);
    }
  }

  void testHex8()
  {
    LOG_UNIT_TEST;

    {
      // Reference
      const Elem & hex8 = ReferenceElem::get(HEX8);
      const Point n1 = hex8.side_vertex_average_normal(0);
      LIBMESH_ASSERT_FP_EQUAL(0, n1(0), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(0, n1(1), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(-1, n1(2), TOLERANCE*TOLERANCE);
      const Point n2 = hex8.side_vertex_average_normal(1);
      LIBMESH_ASSERT_FP_EQUAL(0, n2(0), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(-1, n2(1), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(0, n2(2), TOLERANCE*TOLERANCE);
      const Point n3 = hex8.side_vertex_average_normal(2);
      LIBMESH_ASSERT_FP_EQUAL(1, n3(0), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(0, n3(1), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(0, n3(2), TOLERANCE*TOLERANCE);
      const Point n4 = hex8.side_vertex_average_normal(3);
      LIBMESH_ASSERT_FP_EQUAL(0, n4(0), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(1, n4(1), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(0, n4(2), TOLERANCE*TOLERANCE);
      const Point n5 = hex8.side_vertex_average_normal(4);
      LIBMESH_ASSERT_FP_EQUAL(-1, n5(0), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(0, n5(1), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(0, n5(2), TOLERANCE*TOLERANCE);
      const Point n6 = hex8.side_vertex_average_normal(5);
      LIBMESH_ASSERT_FP_EQUAL(0, n6(0), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(0, n6(1), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(1, n6(2), TOLERANCE*TOLERANCE);
    }

    {
      // Non-planar quad faces
      std::vector<Point> pts = {Point(0, 0, 0), Point(1, 0, 3), Point(1, 1, 0), Point(0, 2, 1),
                                Point(-2.2, 0, 4), Point(-1, 0, 5), Point(-0.5, 1, 4), Point(-2.2, 2, 6)};
      auto [hex8, nodes] = this->construct_elem(pts, HEX8);
      std::unique_ptr<libMesh::FEBase> fe(libMesh::FEBase::build(3, libMesh::FEType(1)));
      libMesh::QGauss qface(2, libMesh::CONSTANT);
      const std::vector<Point> & normals = fe->get_normals();
      for (const auto s : make_range(hex8->n_sides()))
      {
        const std::unique_ptr<const Elem> face = hex8->build_side_ptr(s);
        fe->attach_quadrature_rule(&qface);
        fe->reinit(hex8.get(), s, TOLERANCE);
        const Point n1 = hex8->side_vertex_average_normal(s);
        LIBMESH_ASSERT_FP_EQUAL(normals[0](0), n1(0), TOLERANCE*TOLERANCE);
        LIBMESH_ASSERT_FP_EQUAL(normals[0](1), n1(1), TOLERANCE*TOLERANCE);
        LIBMESH_ASSERT_FP_EQUAL(normals[0](2), n1(2), TOLERANCE*TOLERANCE);
      }
    }
  }

  void testC0Polygon()
  {
    LOG_UNIT_TEST;
    {
      // Square
      std::vector<Point> pts = {Point(0, 0, 0), Point(1, 0, 0), Point(1, 1, 0), Point(0, 1, 0)};
      auto [square, nodes] = this->construct_elem(pts, C0POLYGON);
      const Point n1 = square->side_vertex_average_normal(0);
      LIBMESH_ASSERT_FP_EQUAL(0, n1(0), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(-1, n1(1), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(0, n1(2), TOLERANCE*TOLERANCE);
      const Point n2 = square->side_vertex_average_normal(1);
      LIBMESH_ASSERT_FP_EQUAL(1, n2(0), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(0, n2(1), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(0, n2(2), TOLERANCE*TOLERANCE);
      const Point n3 = square->side_vertex_average_normal(2);
      LIBMESH_ASSERT_FP_EQUAL(0, n3(0), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(1, n3(1), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(0, n3(2), TOLERANCE*TOLERANCE);
      const Point n4 = square->side_vertex_average_normal(3);
      LIBMESH_ASSERT_FP_EQUAL(-1, n4(0), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(0, n4(1), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(0, n4(2), TOLERANCE*TOLERANCE);
    }
    {
      // Hexagon (not the ref one but an easy one to draw)
      std::vector<Point> pts = {Point(0, 0, 0), Point(1, 0, 0), Point(1.5, 1, 0), Point(1, 2, 0), Point(0, 2, 0) , Point(-0.5, 1, 0)};
      auto [hexagon, nodes] = this->construct_elem(pts, C0POLYGON);
      const Point n1 = hexagon->side_vertex_average_normal(0);
      LIBMESH_ASSERT_FP_EQUAL(0, n1(0), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(-1, n1(1), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(0, n1(2), TOLERANCE*TOLERANCE);
      const Point n2 = hexagon->side_vertex_average_normal(1);
      LIBMESH_ASSERT_FP_EQUAL(2 / sqrt(5), n2(0), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(-1 / sqrt(5), n2(1), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(0, n2(2), TOLERANCE*TOLERANCE);
      const Point n3 = hexagon->side_vertex_average_normal(2);
      LIBMESH_ASSERT_FP_EQUAL(2 / sqrt(5), n3(0), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(1 / sqrt(5), n3(1), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(0, n3(2), TOLERANCE*TOLERANCE);
      const Point n4 = hexagon->side_vertex_average_normal(3);
      LIBMESH_ASSERT_FP_EQUAL(0, n4(0), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(1, n4(1), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(0, n4(2), TOLERANCE*TOLERANCE);
      const Point n5 = hexagon->side_vertex_average_normal(4);
      LIBMESH_ASSERT_FP_EQUAL(-2 / sqrt(5), n5(0), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(1 / sqrt(5), n5(1), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(0, n5(2), TOLERANCE*TOLERANCE);
      const Point n6 = hexagon->side_vertex_average_normal(5);
      LIBMESH_ASSERT_FP_EQUAL(-2 / sqrt(5), n6(0), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(-1 / sqrt(5), n6(1), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(0, n6(2), TOLERANCE*TOLERANCE);
    }
    {
      // Non-planar quad (see quad unit test)
      // Note that we don't necessarily allow non-planar polys at this time
      std::vector<Point> pts = {Point(0, 0, 0), Point(1, 0, 3), Point(1, 1, 0), Point(0, 2, 1)};
      auto [poly4, nodes] = this->construct_elem(pts, C0POLYGON);
      // Use a quad4 to get the normals
      auto [quad4, nodes2] = this->construct_elem(pts, QUAD4);
      const std::unique_ptr<const Elem> face = quad4->build_side_ptr(0);
      std::unique_ptr<libMesh::FEBase> fe(libMesh::FEBase::build(2, libMesh::FEType(1)));
      libMesh::QGauss qface(1, libMesh::CONSTANT);
      fe->attach_quadrature_rule(&qface);
      const std::vector<Point> & normals = fe->get_normals();
      for (const auto s : make_range(quad4->n_sides()))
      {
        const std::unique_ptr<const Elem> face = quad4->build_side_ptr(s);
        fe->attach_quadrature_rule(&qface);
        fe->reinit(quad4.get(), s, TOLERANCE);
        const Point n1 = quad4->side_vertex_average_normal(s);
        LIBMESH_ASSERT_FP_EQUAL(normals[0](0), n1(0), TOLERANCE*TOLERANCE);
        LIBMESH_ASSERT_FP_EQUAL(normals[0](1), n1(1), TOLERANCE*TOLERANCE);
        LIBMESH_ASSERT_FP_EQUAL(normals[0](2), n1(2), TOLERANCE*TOLERANCE);
      }
    }
  }

  void testC0Polyhedron()
  {
    LOG_UNIT_TEST;
    {
      // Cube
      std::vector<Point> points = {Point(0, 0, 0), Point(1, 0, 0), Point(1, 1, 0), Point(0, 1, 0),
                                   Point(0, 0, 1), Point(1, 0, 1), Point(1, 1, 1),  Point(0, 1, 1)};

      // See notes in elem_test.h
      const std::vector<std::vector<unsigned int>> nodes_on_side =
        { {0, 1, 2, 3},   // min z
          {0, 1, 5, 4},   // min y
          {2, 6, 5, 1},   // max x
          {2, 3, 7, 6},   // max y
          {0, 4, 7, 3},   // min x
          {5, 6, 7, 4} }; // max z

      // Create Nodes
      std::vector<std::unique_ptr<Node>> nodes(points.size());
      for (const auto i : index_range(points))
        nodes[i] = Node::build(points[i], /*id*/ i);

      // Build all the sides of the cube
      std::vector<std::shared_ptr<Polygon>> sides(nodes_on_side.size());

      for (auto s : index_range(nodes_on_side))
      {
        const auto & nodes_on_s = nodes_on_side[s];
        sides[s] = std::make_shared<C0Polygon>(nodes_on_s.size());
        for (auto i : index_range(nodes_on_s))
          sides[s]->set_node(i, nodes[nodes_on_s[i]].get());
      }

      std::unique_ptr<Elem> polyhedron = std::make_unique<C0Polyhedron>(sides);
      const Point n1 = polyhedron->side_vertex_average_normal(0);
      LIBMESH_ASSERT_FP_EQUAL(0, n1(0), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(0, n1(1), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(-1, n1(2), TOLERANCE*TOLERANCE);
      const Point n2 = polyhedron->side_vertex_average_normal(1);
      LIBMESH_ASSERT_FP_EQUAL(0, n2(0), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(-1, n2(1), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(0, n2(2), TOLERANCE*TOLERANCE);
      const Point n3 = polyhedron->side_vertex_average_normal(2);
      LIBMESH_ASSERT_FP_EQUAL(1, n3(0), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(0, n3(1), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(0, n3(2), TOLERANCE*TOLERANCE);
      const Point n4 = polyhedron->side_vertex_average_normal(3);
      LIBMESH_ASSERT_FP_EQUAL(0, n4(0), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(1, n4(1), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(0, n4(2), TOLERANCE*TOLERANCE);
      const Point n5 = polyhedron->side_vertex_average_normal(4);
      LIBMESH_ASSERT_FP_EQUAL(-1, n5(0), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(0, n5(1), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(0, n5(2), TOLERANCE*TOLERANCE);
      const Point n6 = polyhedron->side_vertex_average_normal(5);
      LIBMESH_ASSERT_FP_EQUAL(0, n6(0), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(0, n6(1), TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(1, n6(2), TOLERANCE*TOLERANCE);
    }

    {
      // Deformed cube, but still with planar faces
      // Note that we don't necessarily allow polyhedras with non-planar sides at this time
      std::vector<Point> pts = {Point(0, 0, 0), Point(1, 0, 3), Point(1, 1, 0), Point(0, 2, -1),
                                Point(0, 0, 4), Point(1, 0, 5), Point(1, 1, 4), Point(0, 2, 6)};

      // See notes in elem_test.h
      const std::vector<std::vector<unsigned int>> nodes_on_side =
        { {0, 1, 2, 3},   // min z
          {0, 1, 5, 4},   // min y
          {2, 6, 5, 1},   // max x
          {2, 3, 7, 6},   // max y
          {0, 4, 7, 3},   // min x
          {5, 6, 7, 4} }; // max z

      // Create Nodes
      std::vector<std::unique_ptr<Node>> nodes(pts.size());
      for (const auto i : index_range(pts))
        nodes[i] = Node::build(pts[i], /*id*/ i);

      // Build all the sides of the cube
      std::vector<std::shared_ptr<Polygon>> sides(nodes_on_side.size());

      for (auto s : index_range(nodes_on_side))
      {
        const auto & nodes_on_s = nodes_on_side[s];
        sides[s] = std::make_shared<C0Polygon>(nodes_on_s.size());
        for (auto i : index_range(nodes_on_s))
          sides[s]->set_node(i, nodes[nodes_on_s[i]].get());
      }

      std::unique_ptr<Elem> polyhedron = std::make_unique<C0Polyhedron>(sides);

      // Get a hex8 for normal comparisons
      auto [hex8, nodes2] = this->construct_elem(pts, HEX8);
      std::unique_ptr<libMesh::FEBase> fe(libMesh::FEBase::build(3, libMesh::FEType(1)));
      libMesh::QGauss qface(2, libMesh::CONSTANT);
      const std::vector<Point> & normals = fe->get_normals();
      for (const auto s : make_range(hex8->n_sides()))
      {
        const std::unique_ptr<const Elem> face = hex8->build_side_ptr(s);
        fe->attach_quadrature_rule(&qface);
        fe->reinit(hex8.get(), s, TOLERANCE);
        const Point n1 = polyhedron->side_vertex_average_normal(s);
        LIBMESH_ASSERT_FP_EQUAL(normals[0](0), n1(0), TOLERANCE*TOLERANCE);
        LIBMESH_ASSERT_FP_EQUAL(normals[0](1), n1(1), TOLERANCE*TOLERANCE);
        LIBMESH_ASSERT_FP_EQUAL(normals[0](2), n1(2), TOLERANCE*TOLERANCE);
      }
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
    std::unique_ptr<Elem> elem;
    if (elem_type != C0POLYGON)
      elem = Elem::build(elem_type, /*parent*/ nullptr);
    else
      elem = std::make_unique<C0Polygon>(n_points);

    // Make sure we were passed consistent input to build this type of Elem
    libmesh_error_msg_if(elem->n_nodes() != n_points,
                          "Wrong number of points "
                          << n_points
                          << " provided to build a "
                          << Utility::enum_to_string(elem_type));

    for (unsigned int i=0; i<n_points; i++)
      elem->set_node(i, nodes[i].get());

    // Return Elem and Nodes we created
    return std::make_pair(std::move(elem), std::move(nodes));
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION( SideVertexAverageNormalTest );
