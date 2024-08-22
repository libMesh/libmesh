
#include "libmesh_cppunit.h"

#ifdef LIBMESH_HAVE_NETGEN
namespace nglib {
#include "netgen/nglib/nglib.h"
}
#endif

#include <numeric>

using namespace nglib;

class LibMeshNetgenTest : public CppUnit::TestCase {
public:
  LIBMESH_CPPUNIT_TEST_SUITE( LibMeshNetgenTest );

#ifdef LIBMESH_HAVE_NETGEN
  CPPUNIT_TEST( testLibMeshNetgen );
  // CPPUNIT_TEST( testBadHole );
#endif

  CPPUNIT_TEST_SUITE_END();

private:

public:

  void setUp()
  {}

  void tearDown()
  {}

#ifdef LIBMESH_HAVE_NETGEN
  void testLibMeshNetgen ()
  {
    LOG_UNIT_TEST;

    auto ngmesh = Ng_NewMesh();
    Ng_DeleteMesh(ngmesh);
  }

  void testBadHole ()
  {
    LOG_UNIT_TEST;

    auto ngmesh = Ng_NewMesh();

    Ng_Meshing_Parameters params;
    params.uselocalh = false;
    params.minh = 0;
    params.elementsperedge = 1;
    params.elementspercurve = 1;
    params.closeedgeenable = false;
    params.closeedgefact = 0;
    params.minedgelenenable = false;
    params.minedgelen = 0;
    params.maxh = std::pow(0.4, 1./3.);

    auto add_point = [ngmesh](double x, double y, double z)
    {
      std::array<double, 3> point_val {x,y,z};
      Ng_AddPoint(ngmesh, point_val.data());
    };

    auto add_elem = [ngmesh](int n1, int n2, int n3)
    {
      std::array<int, 3> elem_nodes {n1, n2, n3};
      Ng_AddSurfaceElement(ngmesh, NG_TRIG, elem_nodes.data());
    };

    add_point(-4, 4, -4);
    add_point(-4, 4, 4);
    add_point(4, 4, 4);
    add_elem(1, 2, 3);
    add_point(4, -4, 4);
    add_elem(4, 3, 2);
    add_point(4, 4, -4);
    add_elem(5, 1, 3);
    add_elem(5, 3, 4);
    add_point(4, -4, -4);
    add_elem(1, 5, 6);
    add_elem(6, 5, 4);
    add_point(-4, -4, 4);
    add_elem(6, 4, 7);
    add_elem(7, 4, 2);
    add_point(-4, -4, -4);
    add_elem(8, 1, 6);
    add_elem(8, 6, 7);
    add_elem(1, 8, 2);
    add_elem(8, 7, 2);

    add_point(3, -3, -3);
    add_point(-3, -3, 3);
    add_point(-3, -3, -3);
    add_elem(11, 10, 9);
    add_point(-3, 3, -3);
    add_elem(11, 9, 12);
    add_point(-3, 3, 3);
    add_elem(10, 11, 13);
    add_point(3, 3, 3);
    add_elem(13, 12, 14);
    add_elem(12, 13, 11);
    add_point(3, -3, 3);
    add_elem(9, 10, 15);
    add_elem(13, 14, 15);
    add_elem(13, 15, 10);
    add_point(3, 3, -3);
    add_elem(15, 14, 16);
    add_elem(16, 14, 12);
    add_elem(15, 16, 9);
    add_elem(9, 16, 12);

    Ng_GenerateVolumeMesh(ngmesh, &params);
    const int n_elem = Ng_GetNE(ngmesh);
    const int n_points = Ng_GetNP(ngmesh);
    CPPUNIT_ASSERT(n_points >= 16);

    std::set<int> nodes_seen;
    for (int i = 0; i != n_elem; ++i)
      {
        // Avoid segfault even if ngtype isn't a tet4
        int ngnodes[11];
        Ng_Volume_Element_Type ngtype =
          Ng_GetVolumeElement(ngmesh, i+1, ngnodes);
        CPPUNIT_ASSERT(ngtype == NG_TET);
        for (int n = 0; n != 4; ++n)
          if (ngnodes[n] < 17 && ngnodes[n] > 0)
            nodes_seen.insert(ngnodes[n]);
      }

    CPPUNIT_ASSERT_EQUAL(nodes_seen.size(), std::size_t(16));

    Ng_DeleteMesh(ngmesh);
  }
#endif
};

CPPUNIT_TEST_SUITE_REGISTRATION( LibMeshNetgenTest );
