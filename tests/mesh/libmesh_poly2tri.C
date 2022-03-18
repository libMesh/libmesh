
#include "libmesh_cppunit.h"

#ifdef LIBMESH_HAVE_POLY2TRI
#  include "poly2tri/poly2tri.h"
#endif

#include <numeric>

class LibMeshPoly2TriTest : public CppUnit::TestCase {
public:
  CPPUNIT_TEST_SUITE( LibMeshPoly2TriTest );

#ifdef LIBMESH_HAVE_POLY2TRI
  CPPUNIT_TEST( testLibMeshPoly2Tri );
  CPPUNIT_TEST( testLibMeshPoly2TriHole );
#endif

  CPPUNIT_TEST_SUITE_END();

private:

public:
  void setUp()
  {}

  void tearDown()
  {}

#ifdef LIBMESH_HAVE_POLY2TRI
  void testLibMeshPoly2Tri ()
  {
    std::vector<p2t::Point> pentagon {{0,0},{1,0},{1,1},{.5,1.5},{0,1}};

    std::vector<p2t::Point *> api_shim(pentagon.size());
    std::iota(api_shim.begin(), api_shim.end(), pentagon.data());

    p2t::CDT cdt(api_shim);

    cdt.Triangulate();

    auto tris = cdt.GetTriangles();

    // We mostly wanted to make sure this compiled, but might as well
    // make sure it gave us the expected triangle count while we're at
    // it.
    CPPUNIT_ASSERT_EQUAL(tris.size(), std::size_t(3));
  }

  void testLibMeshPoly2TriHole ()
  {
//    const double eps = 0;  // This (or negative eps) succeeds
//    const double eps = 1e-12;  // This (or larger eps) succeeds
    const double eps = 1e-15; // This gave EdgeEvent - null triangle with older Poly2Tri

    std::vector<p2t::Point> outer_boundary
      {{0,0},{0.5,eps},{1,0},{1-eps,0.836541},
       {1,2},{.46,1.46+eps},{0,1},{eps,0.5}};

    std::vector<p2t::Point *> api_shim(outer_boundary.size());
    std::iota(api_shim.begin(), api_shim.end(), outer_boundary.data());

    const double r2o4 = std::sqrt(2.)/4;
    std::vector<p2t::Point> hole_boundary
      {{0.5+r2o4,0.5},{0.5,0.5+r2o4},{0.5-r2o4,0.5},{0.5-eps,0.5-r2o4}};

    std::vector<p2t::Point *> hole_shim(hole_boundary.size());
    std::iota(hole_shim.begin(), hole_shim.end(), hole_boundary.data());

    p2t::CDT cdt(api_shim);
    cdt.AddHole(hole_shim);

    std::vector<p2t::Point> interior_points
      {{0.21,0.79},{0.21,0.21},{0.79,0.21}};
    for (auto & p : interior_points)
      cdt.AddPoint(&p);

    cdt.Triangulate();

    auto tris = cdt.GetTriangles();

    std::cout << "With perturbation " << eps << " we got " << tris.size() << " triangles!" << std::endl;

    // We mostly wanted to make sure this compiled, but might as well
    // make sure it gave us the expected triangle count while we're at
    // it.
//    CPPUNIT_ASSERT_EQUAL(tris.size(), std::size_t(3));
  }

#endif
};

CPPUNIT_TEST_SUITE_REGISTRATION( LibMeshPoly2TriTest );
