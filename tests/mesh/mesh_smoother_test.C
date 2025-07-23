#include <libmesh/boundary_info.h>
#include <libmesh/elem.h>
#include <libmesh/enum_elem_type.h>
#include <libmesh/function_base.h>
#include <libmesh/libmesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_modification.h>
#include <libmesh/mesh_smoother_laplace.h>
#include <libmesh/mesh_smoother_vsmoother.h>
#include <libmesh/mesh_tools.h>
#include <libmesh/node.h>
#include <libmesh/reference_elem.h>
#include <libmesh/replicated_mesh.h>
#include <libmesh/system.h> // LIBMESH_HAVE_SOLVER define

#include "test_comm.h"
#include "libmesh_cppunit.h"

namespace
{
using namespace libMesh;

// Distortion function that doesn't distort boundary nodes
// 2D only, use for LaplaceMeshSmoother
class DistortSquare : public FunctionBase<Real>
{
  std::unique_ptr<FunctionBase<Real>> clone() const override
  {
    return std::make_unique<DistortSquare>();
  }

  Real operator()(const Point &, const Real = 0.) override
  {
    libmesh_not_implemented();
  } // scalar-only API

  // Skew inward based on a cubic function
  void operator()(const Point & p, const Real, DenseVector<Real> & output)
  {
    output.resize(3);
    const Real eta = 2 * p(0) - 1;
    const Real zeta = 2 * p(1) - 1;
    output(0) = p(0) + (std::pow(eta, 3) - eta) * p(1) * (1 - p(1));
    output(1) = p(1) + (std::pow(zeta, 3) - zeta) * p(0) * (1 - p(0));
    output(2) = 0;
  }
};

// Distortion function supporting 1D, 2D, and 3D
// Use for VariationalMeshSmoother
class DistortHyperCube : public FunctionBase<Real>
{
public:
  DistortHyperCube(const unsigned int dim) : _dim(dim) {}

private:
  std::unique_ptr<FunctionBase<Real>> clone() const override
  {
    return std::make_unique<DistortHyperCube>(_dim);
  }

  Real operator()(const Point &, const Real = 0.) override { libmesh_not_implemented(); }

  void operator()(const Point & p, const Real, DenseVector<Real> & output) override
  {
    output.resize(3);
    output.zero();

    // Count how many coordinates are exactly on the boundary (0 or 1)
    unsigned int boundary_dims = 0;
    std::array<bool, 3> is_on_boundary = {false, false, false};
    for (unsigned int i = 0; i < _dim; ++i)
      {
        if (std::abs(p(i)) < TOLERANCE || std::abs(p(i) - 1.) < TOLERANCE)
          {
            ++boundary_dims;
            is_on_boundary[i] = true;
          }
      }

    // If all coordinates are on the boundary, treat as vertex — leave unchanged
    if (boundary_dims == _dim)
      {
        for (unsigned int i = 0; i < _dim; ++i)
          output(i) = p(i);
        return;
      }

    // Distort only those directions not fixed on the boundary
    for (unsigned int i = 0; i < _dim; ++i)
      {
        if (!is_on_boundary[i]) // only distort free dimensions
          {
            Real xi = 2. * p(i) - 1.;
            Real modulation = 0.3; // This value constrols the strength of the distortion
            for (unsigned int j = 0; j < _dim; ++j)
              {
                if (j != i)
                  {
                    Real pj = std::clamp(p(j), 0., 1.);         // ensure numeric safety
                    modulation *= (pj - 0.5) * (pj - 0.5) * 4.; // quadratic bump centered at 0.5
                  }
              }
            output(i) = p(i) + (std::pow(xi, 3) - xi) * modulation;
          }
        else
          {
            output(i) = p(i); // dimension on boundary remains unchanged
          }
      }
  }

  const unsigned int _dim;
};

class SquareToParallelogram : public FunctionBase<Real>
{
  std::unique_ptr<FunctionBase<Real>> clone() const override
  {
    return std::make_unique<SquareToParallelogram>();
  }

  Real operator()(const Point &, const Real = 0.) override
  {
    libmesh_not_implemented();
  } // scalar-only API

  // Has the effect that a square, meshed into right triangles with diagonals
  // rising from lower-left to upper-right, is transformed into a left-leaning
  // parallelogram of eqilateral triangles.
  void operator()(const Point & p, const Real, DenseVector<Real> & output)
  {
    output.resize(3);
    output(0) = p(0) - 0.5 * p(1);
    output(1) = 0.5 * std::sqrt(3) * p(1);
    output(2) = 0;
  }
};

class ParallelogramToSquare : public FunctionBase<Real>
{
  std::unique_ptr<FunctionBase<Real>> clone() const override
  {
    return std::make_unique<ParallelogramToSquare>();
  }

  Real operator()(const Point &, const Real = 0.) override
  {
    libmesh_not_implemented();
  } // scalar-only API

  // Has the effect that a left-leaning parallelogram of equilateral triangles is transformed
  // into a square of right triangle with diagonals rising from lower-left to upper-right.
  // This is the inversion of the SquareToParallelogram mapping.
  void operator()(const Point & p, const Real, DenseVector<Real> & output)
  {
    output.resize(3);
    output(0) = p(0) + p(1) / std::sqrt(3.);
    output(1) = (2. / std::sqrt(3)) * p(1);
    output(2) = 0;
  }
};

}

using namespace libMesh;

class MeshSmootherTest : public CppUnit::TestCase
{
  /**
   * The goal of this test is to verify proper operation of the
   * MeshSmoother subclasses.
   */
public:
  LIBMESH_CPPUNIT_TEST_SUITE(MeshSmootherTest);

#if LIBMESH_DIM > 2
  CPPUNIT_TEST(testLaplaceQuad);
  CPPUNIT_TEST(testLaplaceTri);
#if defined(LIBMESH_ENABLE_VSMOOTHER) && defined(LIBMESH_HAVE_SOLVER)
  CPPUNIT_TEST(testVariationalEdge);
  CPPUNIT_TEST(testVariationalEdgeMultipleSubdomains);
  CPPUNIT_TEST(testVariationalQuad);
  CPPUNIT_TEST(testVariationalQuadMultipleSubdomains);
  CPPUNIT_TEST(testVariationalTri);
  CPPUNIT_TEST(testVariationalTriMultipleSubdomains);
  CPPUNIT_TEST(testVariationalHex);
  CPPUNIT_TEST(testVariationalHexMultipleSubdomains);
#endif // LIBMESH_ENABLE_VSMOOTHER
#endif

  CPPUNIT_TEST_SUITE_END();

public:
  void setUp() {}

  void tearDown() {}

  void testLaplaceSmoother(ReplicatedMesh & mesh, MeshSmoother & smoother, ElemType type)
  {
    LOG_UNIT_TEST;

    constexpr unsigned int n_elems_per_side = 4;

    MeshTools::Generation::build_square(
        mesh, n_elems_per_side, n_elems_per_side, 0., 1., 0., 1., type);

    // Move it around so we have something that needs smoothing
    DistortSquare ds;
    MeshTools::Modification::redistribute(mesh, ds);

    // Assert the distortion is as expected
    auto center_distortion_is =
        [](const Node & node, int d, bool distortion, Real distortion_tol = TOLERANCE) {
          const Real r = node(d);
          const Real R = r * n_elems_per_side;
          CPPUNIT_ASSERT_GREATER(-TOLERANCE * TOLERANCE, r);
          CPPUNIT_ASSERT_GREATER(-TOLERANCE * TOLERANCE, 1 - r);

          // If we're at the boundaries we should *never* be distorted
          if (std::abs(node(0)) < TOLERANCE * TOLERANCE ||
              std::abs(node(0) - 1) < TOLERANCE * TOLERANCE)
            {
              const Real R1 = node(1) * n_elems_per_side;
              CPPUNIT_ASSERT_LESS(TOLERANCE * TOLERANCE, std::abs(R1 - std::round(R1)));
              return true;
            }

          if (std::abs(node(1)) < TOLERANCE * TOLERANCE ||
              std::abs(node(1) - 1) < TOLERANCE * TOLERANCE)
            {
              const Real R0 = node(0) * n_elems_per_side;
              CPPUNIT_ASSERT_LESS(TOLERANCE * TOLERANCE, std::abs(R0 - std::round(R0)));

              return true;
            }

          // If we're at the center we're fine
          if (std::abs(r - 0.5) < TOLERANCE * TOLERANCE)
            return true;

          return ((std::abs(R - std::round(R)) > distortion_tol) == distortion);
        };

    for (auto node : mesh.node_ptr_range())
      {
        CPPUNIT_ASSERT(center_distortion_is(*node, 0, true));
        CPPUNIT_ASSERT(center_distortion_is(*node, 1, true));
      }

    // Enough iterations to mostly fix us up.  Laplace seems to be at 1e-3
    // tolerance by iteration 6, so hopefully everything is there on any
    // system by 8.
    for (unsigned int i = 0; i != 8; ++i)
      smoother.smooth();

    // Make sure we're not too distorted anymore.
    for (auto node : mesh.node_ptr_range())
      {
        CPPUNIT_ASSERT(center_distortion_is(*node, 0, false, 1e-3));
        CPPUNIT_ASSERT(center_distortion_is(*node, 1, false, 1e-3));
      }
  }

  void testVariationalSmoother(ReplicatedMesh & mesh,
                               MeshSmoother & smoother,
                               const ElemType type,
                               const bool multiple_subdomains = false)
  {
    LOG_UNIT_TEST;

    const auto dim = ReferenceElem::get(type).dim();

    unsigned int n_elems_per_side = 5;

    // If n_elems_per_side is even, then some sliding boundary nodes will have a
    // coordinante with value 0.5, which is already the optimal position,
    // causing distortion_is(node, true) to return false when evaluating the
    // distorted mesh. To avoid this, we require n_elems_per_side to be odd.
    libmesh_error_msg_if(n_elems_per_side % 2 != 1, "n_elems_per_side should be odd.");

    switch (dim)
      {
        case 1:
          MeshTools::Generation::build_line(mesh, n_elems_per_side, 0., 1., type);
          break;
        case 2:
          MeshTools::Generation::build_square(
              mesh, n_elems_per_side, n_elems_per_side, 0., 1., 0., 1., type);
          break;

        case 3:
          MeshTools::Generation::build_cube(mesh,
                                            n_elems_per_side,
                                            n_elems_per_side,
                                            n_elems_per_side,
                                            0.,
                                            1.,
                                            0.,
                                            1.,
                                            0.,
                                            1.,
                                            type);
          break;

        default:
          libmesh_error_msg("Unsupported dimension " << dim);
      }

    // Move it around so we have something that needs smoothing
    DistortHyperCube dh(dim);
    MeshTools::Modification::redistribute(mesh, dh);

    // Add multiple subdomains if requested
    std::unordered_map<dof_id_type, Point> subdomain_boundary_node_id_to_point;
    if (multiple_subdomains)
      {
        // Modify the subdomain ids in an interesting way
        for (auto * elem : mesh.active_element_ptr_range())
          {
            unsigned int subdomain_id = 0;
            for (const auto d : make_range(dim))
              if (elem->vertex_average()(d) > 0.5)
                ++subdomain_id;
            elem->subdomain_id() += subdomain_id;
          }

        // This loop should NOT be combined with the one above because we need to
        // finish checking and updating subdomain ids for all elements before
        // recording the final subdomain boundary.
        for (auto * elem : mesh.active_element_ptr_range())
          for (const auto & s : elem->side_index_range())
            {
              const auto * neighbor = elem->neighbor_ptr(s);
              if (neighbor == nullptr)
                continue;

              if (elem->subdomain_id() != neighbor->subdomain_id())
                // This side is part of a subdomain boundary, record the
                // corresponding node locations
                for (const auto & n : elem->nodes_on_side(s))
                  subdomain_boundary_node_id_to_point[elem->node_id(n)] =
                      Point(*(elem->get_nodes()[n]));
            }
      }

    // Function to assert the distortion is as expected
    const auto & boundary_info = mesh.get_boundary_info();
    auto distortion_is = [&n_elems_per_side, &dim, &boundary_info](
                             const Node & node, bool distortion, Real distortion_tol = TOLERANCE) {
      // Get boundary ids associated with the node
      std::vector<boundary_id_type> boundary_ids;
      boundary_info.boundary_ids(&node, boundary_ids);

      // This tells us what type of node we are: internal, sliding, or fixed
      const auto num_dofs = dim - boundary_ids.size();
      /*
       * The following cases of num_dofs are possible, ASSUMING all boundaries
       * are non-overlapping
       * 3D: 3-0,     3-1,     3-2,     3-3
       *   = 3        2        1        0
       *     internal sliding, sliding, fixed
       * 2D: 2-0,     2-1,     2-2
       *   = 2        1        0
       *     internal sliding, fixed
       * 1D: 1-0,     1-1
       *   = 1        0
       *     internal fixed
       *
       * We expect that R is an integer in [0, n_elems_per_side] for
       * num_dofs of the node's cooridinantes, while the remaining coordinates
       * are fixed to the boundary with value 0 or 1. In other words, at LEAST
       * dim - num_dofs coordinantes should be 0 or 1.
       */

      size_t num_zero_or_one = 0;

      bool distorted = false;
      for (const auto d : make_range(dim))
        {
          const Real r = node(d);
          const Real R = r * n_elems_per_side;
          CPPUNIT_ASSERT_GREATER(-distortion_tol * distortion_tol, r);
          CPPUNIT_ASSERT_GREATER(-distortion_tol * distortion_tol, 1 - r);

          // Due to the type of distortion used, nodes on the x, y, or z plane
          // of symmetry do not have their respective x, y, or z node adjusted.
          // Just continue to the next dimension.
          if (std::abs(r - 0.5) < distortion_tol * distortion_tol)
            continue;

          const bool d_distorted = std::abs(R - std::round(R)) > distortion_tol;
          distorted |= d_distorted;
          num_zero_or_one += (absolute_fuzzy_equals(r, 0.) || absolute_fuzzy_equals(r, 1.));
        }

      CPPUNIT_ASSERT_GREATEREQUAL(dim - num_dofs, num_zero_or_one);

      // We can never expect a fixed node to be distorted
      if (num_dofs == 0)
        // if (num_dofs < dim)
        return true;
      return distorted == distortion;
    };

    // Function to check if a given node has changed based on previous mapping
    auto is_subdomain_boundary_node_the_same = [&subdomain_boundary_node_id_to_point](
                                                   const Node & node) {
      auto it = subdomain_boundary_node_id_to_point.find(node.id());
      if (it != subdomain_boundary_node_id_to_point.end())
        return (relative_fuzzy_equals(Point(node), subdomain_boundary_node_id_to_point[node.id()]));
      else
        // node is not a subdomain boundary node, just return true
        return true;
    };

    // Make sure our DistortSquare transformation has distorted the mesh
    for (auto node : mesh.node_ptr_range())
      CPPUNIT_ASSERT(distortion_is(*node, true));

    // Transform the square mesh of triangles to a parallelogram mesh of
    // triangles. This will allow the Variational Smoother to smooth the mesh
    // to the optimal case of equilateral triangles
    if (type == TRI3)
      {
        SquareToParallelogram stp;
        MeshTools::Modification::redistribute(mesh, stp);
      }

    smoother.smooth();

    // Transform the parallelogram mesh back to a square mesh. In the case of
    // the Variational Smoother, equilateral triangular elements will be
    // transformed into right triangular elements that align with the original
    // undistorted mesh.
    if (type == TRI3)
      {
        ParallelogramToSquare pts;
        MeshTools::Modification::redistribute(mesh, pts);
      }

    // Make sure we're not too distorted anymore OR that interval subdomain boundary nodes did not
    // change.
    for (auto node : mesh.node_ptr_range())
      {
        if (multiple_subdomains)
          CPPUNIT_ASSERT(is_subdomain_boundary_node_the_same(*node));
        else
          CPPUNIT_ASSERT(distortion_is(*node, false, 1e-3));
      }
  }

  void testLaplaceQuad()
  {
    ReplicatedMesh mesh(*TestCommWorld);
    LaplaceMeshSmoother laplace(mesh);

    testLaplaceSmoother(mesh, laplace, QUAD4);
  }

  void testLaplaceTri()
  {
    ReplicatedMesh mesh(*TestCommWorld);
    LaplaceMeshSmoother laplace(mesh);

    testLaplaceSmoother(mesh, laplace, TRI3);
  }

#ifdef LIBMESH_ENABLE_VSMOOTHER
  void testVariationalEdge()
  {
    ReplicatedMesh mesh(*TestCommWorld);
    VariationalMeshSmoother variational(mesh);

    testVariationalSmoother(mesh, variational, EDGE2);
  }

  void testVariationalEdgeMultipleSubdomains()
  {
    ReplicatedMesh mesh(*TestCommWorld);
    VariationalMeshSmoother variational(mesh);

    testVariationalSmoother(mesh, variational, EDGE2, true);
  }

  void testVariationalQuad()
  {
    ReplicatedMesh mesh(*TestCommWorld);
    VariationalMeshSmoother variational(mesh);

    testVariationalSmoother(mesh, variational, QUAD4);
  }

  void testVariationalQuadMultipleSubdomains()
  {
    ReplicatedMesh mesh(*TestCommWorld);
    VariationalMeshSmoother variational(mesh);

    testVariationalSmoother(mesh, variational, QUAD4, true);
  }

  void testVariationalTri()
  {
    ReplicatedMesh mesh(*TestCommWorld);
    VariationalMeshSmoother variational(mesh);

    testVariationalSmoother(mesh, variational, TRI3);
  }

  void testVariationalTriMultipleSubdomains()
  {
    ReplicatedMesh mesh(*TestCommWorld);
    VariationalMeshSmoother variational(mesh);

    testVariationalSmoother(mesh, variational, TRI3, true);
  }

  void testVariationalHex()
  {
    ReplicatedMesh mesh(*TestCommWorld);
    VariationalMeshSmoother variational(mesh);

    testVariationalSmoother(mesh, variational, HEX8);
  }

  void testVariationalHexMultipleSubdomains()
  {
    ReplicatedMesh mesh(*TestCommWorld);
    VariationalMeshSmoother variational(mesh);

    testVariationalSmoother(mesh, variational, HEX8, true);
  }
#endif // LIBMESH_ENABLE_VSMOOTHER
};

CPPUNIT_TEST_SUITE_REGISTRATION(MeshSmootherTest);
