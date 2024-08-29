
#include "libmesh_cppunit.h"

#ifdef LIBMESH_HAVE_NETGEN
namespace nglib {
#include "netgen/nglib/nglib.h"
}

using namespace nglib;
#endif

#include <numeric>

class LibMeshNetgenTest : public CppUnit::TestCase {
public:
  LIBMESH_CPPUNIT_TEST_SUITE( LibMeshNetgenTest );

#ifdef LIBMESH_HAVE_NETGEN
  CPPUNIT_TEST( testLibMeshNetgen );
  CPPUNIT_TEST( testBadHole );

  // Too expensive to do regularly...
  // CPPUNIT_TEST( testAllPermutations );
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

  std::size_t testPermutedHole (const std::vector<std::size_t> & permutation)
  {
    auto ngmesh = Ng_NewMesh();

    Ng_Meshing_Parameters params;
    params.elementsperedge = 1;
    params.elementspercurve = 1;

    std::vector<std::array<double, 3>> point_coords =
      {{-4, 4, -4},
      {-4, 4, 4},
      {4, 4, 4},
      {4, -4, 4},
      {4, 4, -4},
      {4, -4, -4},
      {-4, -4, 4},
      {-4, -4, -4},

      {3, -3, -3},
      {-3, -3, 3},
      {-3, -3, -3},
      {-3, 3, -3},
      {-3, 3, 3},
      {3, 3, 3},
      {3, -3, 3},
      {3, 3, -3}};

    auto n_input_points = point_coords.size();

    CPPUNIT_ASSERT_EQUAL(permutation.size(), n_input_points);

    std::vector<std::size_t> inverse_permutation (n_input_points);

    for (std::size_t i=0; i != n_input_points; ++i)
      {
        inverse_permutation[permutation[i]] = i;
        Ng_AddPoint(ngmesh, point_coords[permutation[i]].data());
      }

    std::vector<std::array<int, 3>> tri_nodes =
      {{1, 2, 3},
      {4, 3, 2},
      {5, 1, 3},
      {5, 3, 4},
      {1, 5, 6},
      {6, 5, 4},
      {6, 4, 7},
      {7, 4, 2},
      {8, 1, 6},
      {8, 6, 7},
      {1, 8, 2},
      {8, 7, 2},

      {11, 10, 9},
      {11, 9, 12},
      {10, 11, 13},
      {13, 12, 14},
      {12, 13, 11},
      {9, 10, 15},
      {13, 14, 15},
      {13, 15, 10},
      {15, 14, 16},
      {16, 14, 12},
      {15, 16, 9},
      {9, 16, 12}};

    for (auto & nodes : tri_nodes)
      {
        std::array<int, 3> permuted_nodes = nodes;
        for (auto & node_i : permuted_nodes)
          node_i = inverse_permutation[node_i-1]+1;
        Ng_AddSurfaceElement(ngmesh, NG_TRIG, permuted_nodes.data());
      }

    Ng_GenerateVolumeMesh(ngmesh, &params);
    const std::size_t n_elem = Ng_GetNE(ngmesh);
    const std::size_t n_points = Ng_GetNP(ngmesh);
    CPPUNIT_ASSERT_GREATEREQUAL(n_input_points, n_points);

    std::set<int> nodes_seen;
    for (std::size_t i = 0; i != n_elem; ++i)
      {
        // Avoid segfault even if ngtype isn't a tet4
        int ngnodes[11];
        Ng_Volume_Element_Type ngtype =
          Ng_GetVolumeElement(ngmesh, i+1, ngnodes);
        CPPUNIT_ASSERT(ngtype == NG_TET);
        for (int n = 0; n != 4; ++n)
          if (ngnodes[n] <= int(n_input_points) && ngnodes[n] > 0)
            nodes_seen.insert(ngnodes[n]);
      }

    Ng_DeleteMesh(ngmesh);

    return nodes_seen.size();
  }

  void testBadHole ()
  {
    LOG_UNIT_TEST;

    std::vector<std::size_t> permutation (16);
    std::iota(permutation.begin(), permutation.end(), 0);

    int n_original_nodes = testPermutedHole(permutation);

    CPPUNIT_ASSERT_EQUAL(n_original_nodes, 16);
  }

  void testAllPermutations ()
  {
    LOG_UNIT_TEST;

    const std::size_t n_input_points = 16;

    std::vector<std::size_t> permutation (n_input_points);
    std::iota(permutation.begin(), permutation.end(), 0);

    int fails=0, successes = 0;

    auto run_test = [this, &fails, &successes, &permutation]()
    {
      int n_original_nodes = testPermutedHole(permutation);
      fails += (n_original_nodes != n_input_points);
      successes += (n_original_nodes == n_input_points);
      std::cout << "Found " << n_original_nodes << "/" << n_input_points << " original nodes\n";
      std::cout << "Fails = " << fails << ", successes = " << successes << std::endl;
    };

    // Heap's Algorithm, non-recursive version
    run_test();

    std::vector<std::size_t> c(permutation.size());

    for (std::size_t i = 1; i != n_input_points ; ++i)
      {
        if (c[i] < i)
          {
            if (i%2)
              std::swap(permutation[c[i]], permutation[i]);
            else
              std::swap(permutation[0], permutation[i]);
            run_test();

            ++c[i];
            i = 0; // Yeah, this runs in factorial time
          }
        else
          c[i] = 0;
      }
  }
#endif
};

CPPUNIT_TEST_SUITE_REGISTRATION( LibMeshNetgenTest );
