#include <libmesh/libmesh.h>
#include <libmesh/enum_elem_type.h>
#include <libmesh/function_base.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_modification.h>
#include <libmesh/mesh_smoother_laplace.h>
#include <libmesh/mesh_smoother_vsmoother.h>
#include <libmesh/mesh_tools.h>
#include <libmesh/node.h>
#include <libmesh/elem.h>
#include <libmesh/replicated_mesh.h>
#include <libmesh/boundary_info.h>

#include "test_comm.h"
#include "libmesh_cppunit.h"

namespace {
using namespace libMesh;

class DistortSquare : public FunctionBase<Real>
{
  std::unique_ptr<FunctionBase<Real>> clone () const override
  { return std::make_unique<DistortSquare>(); }

  Real operator() (const Point &,
                   const Real = 0.) override
  { libmesh_not_implemented(); } // scalar-only API

  // Skew inward based on a cubic function
  void operator() (const Point & p,
                   const Real,
                   DenseVector<Real> & output)
  {
    output.resize(3);
    const Real eta = 2*p(0)-1;
    const Real zeta = 2*p(1)-1;
    output(0) = p(0) + (std::pow(eta,3)-eta)*p(1)*(1-p(1));
    output(1) = p(1) + (std::pow(zeta,3)-zeta)*p(0)*(1-p(0));
    output(2) = 0;
  }
};

class SquareToParallelogram : public FunctionBase<Real>
{
  std::unique_ptr<FunctionBase<Real>> clone () const override
  { return std::make_unique<SquareToParallelogram>(); }

  Real operator() (const Point &,
                   const Real = 0.) override
  { libmesh_not_implemented(); } // scalar-only API

  // Has the effect that a square, meshed into right triangles with diagonals
  // rising from lower-left to upper-right, is transformed into a left-leaning
  // parallelogram of eqilateral triangles.
  void operator() (const Point & p,
                   const Real,
                   DenseVector<Real> & output)
  {
    output.resize(3);
    output(0) = p(0) - 0.5 * p(1);
    output(1) = 0.5 * std::sqrt(3) * p(1);
    output(2) = 0;
  }
};

class ParallelogramToSquare : public FunctionBase<Real>
{
  std::unique_ptr<FunctionBase<Real>> clone () const override
  { return std::make_unique<ParallelogramToSquare>(); }

  Real operator() (const Point &,
                   const Real = 0.) override
  { libmesh_not_implemented(); } // scalar-only API

  // Has the effect that a left-leaning parallelogram of equilateral triangles is transformed
  // into a square of right triangle with diagonals rising from lower-left to upper-right.
  // This is the inversion of the SquareToParallelogram mapping.
  void operator() (const Point & p,
                   const Real,
                   DenseVector<Real> & output)
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
  LIBMESH_CPPUNIT_TEST_SUITE( MeshSmootherTest );

#if LIBMESH_DIM > 2
  CPPUNIT_TEST( testLaplaceQuad );
  CPPUNIT_TEST( testLaplaceTri );
#  ifdef LIBMESH_ENABLE_VSMOOTHER
  CPPUNIT_TEST( testVariationalQuad );
  CPPUNIT_TEST( testVariationalTri );
  CPPUNIT_TEST( testVariationalQuadMultipleSubdomains );
#  endif // LIBMESH_ENABLE_VSMOOTHER
#endif

  CPPUNIT_TEST_SUITE_END();

public:
  void setUp() {}

  void tearDown() {}

  void testSmoother(ReplicatedMesh & mesh, MeshSmoother & smoother, const ElemType type, const bool multiple_subdomains=false)
  {
    LOG_UNIT_TEST;

    constexpr unsigned int n_elems_per_side = 4;

    MeshTools::Generation::build_square(mesh, n_elems_per_side, n_elems_per_side,
                                        0.,1.,0.,1., type);

    // Move it around so we have something that needs smoothing
    DistortSquare ds;
    MeshTools::Modification::redistribute(mesh, ds);

    std::unordered_map<dof_id_type, Point> subdomain_boundary_node_id_to_point;
    if (multiple_subdomains)
    {
      // Increment the subdomain id on the right half by 1
      for (auto * elem : mesh.active_element_ptr_range())
        if (elem->vertex_average()(0) > 0.5)
          ++elem->subdomain_id();

      // This loop should NOT be combined with the one above because we need to
      // finish checking and updating subdomain ids for all elements before
      // recording the final subdomain boundary.
      for (auto * elem : mesh.active_element_ptr_range())
        for (const auto & s : elem->side_index_range())
        {
          const auto* neighbor = elem->neighbor_ptr(s);
          if (neighbor == nullptr)
            continue;

          if (elem->subdomain_id() != neighbor->subdomain_id())
            // This side is part of a subdomain boundary, record the
            // corresponding node locations
            for (const auto & n : elem->nodes_on_side(s))
              subdomain_boundary_node_id_to_point[elem->node_id(n)] = Point(*(elem->get_nodes()[n]));
        }
      }

    const auto & boundary_info = mesh.get_boundary_info();

    // Function to assert the distortion is as expected
    auto center_distortion_is = [&boundary_info]
      (const Node & node, int d, bool distortion,
       Real distortion_tol=TOLERANCE) {
      const Real r = node(d);
      const Real R = r * n_elems_per_side;
      CPPUNIT_ASSERT_GREATER(-distortion_tol*distortion_tol, r);
      CPPUNIT_ASSERT_GREATER(-distortion_tol*distortion_tol, 1-r);

      // If we're at the center we're fine
      if (std::abs(r-0.5) < distortion_tol*distortion_tol)
        return true;

      // Boundary nodes are allowed to slide along the boundary.
      // However, nodes that are part of more than one boundary (i.e., corners) should remain fixed.

      // Get boundary ids associated with the node
      std::vector<boundary_id_type> boundary_ids;
      boundary_info.boundary_ids(&node, boundary_ids);

      switch (boundary_ids.size())
      {
        // Internal node
        case 0:
          return ((std::abs(R-std::round(R)) > distortion_tol) == distortion);
          break;
        // Sliding boundary node
        case 1:
          // Since sliding boundary nodes may or may not already be in the optimal
          // position, they may or may not be different from the originally distorted
          // mesh. Return true here to avoid issues.
          return true;
          break;
        // Fixed boundary node, should not have moved
        case 2:
          if (std::abs(node(0)) < distortion_tol*distortion_tol ||
              std::abs(node(0)-1) < distortion_tol*distortion_tol)
            {
              const Real R1 = node(1) * n_elems_per_side;
              CPPUNIT_ASSERT_LESS(distortion_tol*distortion_tol, std::abs(R1-std::round(R1)));
              return true;
            }

          if (std::abs(node(1)) < distortion_tol*distortion_tol ||
              std::abs(node(1)-1) < distortion_tol*distortion_tol)
            {
              const Real R0 = node(0) * n_elems_per_side;
              CPPUNIT_ASSERT_LESS(distortion_tol*distortion_tol, std::abs(R0-std::round(R0)));

              return true;
            }
            return false;
          break;
          default:
            libmesh_error_msg("Node has unsupported number of boundary ids = " << boundary_ids.size());
      }
    };

    // Function to check if a given node has changed based on previous mapping
    auto is_internal_subdomain_boundary_node_the_same = [&subdomain_boundary_node_id_to_point]
      (const Node & node) {
      auto it = subdomain_boundary_node_id_to_point.find(node.id());
      if (it != subdomain_boundary_node_id_to_point.end())
        return (Point(node) == subdomain_boundary_node_id_to_point[node.id()]);
      else
        // node is not an internal subdomain boundary node, just return true
        return true;
    };

    // Make sure our DistortSquare transformation has distorted the mesh
    for (auto node : mesh.node_ptr_range())
      {
        CPPUNIT_ASSERT(center_distortion_is(*node, 0, true));
        CPPUNIT_ASSERT(center_distortion_is(*node, 1, true));
      }

    // Transform the square mesh of triangles to a parallelogram mesh of triangles.
    // This will allow the Variational Smoother to smooth the mesh to the optimal case
    // of equilateral triangles
   const bool is_variational_smoother_type = (dynamic_cast<VariationalMeshSmoother*>(&smoother) != nullptr);
   if (type == TRI3 && is_variational_smoother_type)
    {
      SquareToParallelogram stp;
      MeshTools::Modification::redistribute(mesh, stp);
    }

    // Enough iterations to mostly fix us up.  Laplace seems to be at 1e-3
    // tolerance by iteration 6, so hopefully everything is there on any
    // system by 8.
    for (unsigned int i=0; i != 8; ++i)
      smoother.smooth();

    // Transform the parallelogram mesh back to a square mesh. In the case of the
    // Variational Smoother, equilateral triangular elements will be transformed
    // into right triangular elements that align with the original undistorted mesh.
    if (type == TRI3 && is_variational_smoother_type)
    {
      ParallelogramToSquare pts;
      MeshTools::Modification::redistribute(mesh, pts);
    }

    // Make sure we're not too distorted anymore OR that interval subdomain boundary nodes did not change.
    for (auto node : mesh.node_ptr_range())
      {
        if (multiple_subdomains)
        {
          CPPUNIT_ASSERT(is_internal_subdomain_boundary_node_the_same(*node));
        }
        else
        {
          CPPUNIT_ASSERT(center_distortion_is(*node, 0, false, 1e-3));
          CPPUNIT_ASSERT(center_distortion_is(*node, 1, false, 1e-3));
        }
      }
  }


  void testLaplaceQuad()
  {
    ReplicatedMesh mesh(*TestCommWorld);
    LaplaceMeshSmoother laplace(mesh);

    testSmoother(mesh, laplace, QUAD4);
  }


  void testLaplaceTri()
  {
    ReplicatedMesh mesh(*TestCommWorld);
    LaplaceMeshSmoother laplace(mesh);

    testSmoother(mesh, laplace, TRI3);
  }


#ifdef LIBMESH_ENABLE_VSMOOTHER
  void testVariationalQuad()
  {
    ReplicatedMesh mesh(*TestCommWorld);
    VariationalMeshSmoother variational(mesh);

    testSmoother(mesh, variational, QUAD4);
  }


  void testVariationalTri()
  {
    ReplicatedMesh mesh(*TestCommWorld);
    VariationalMeshSmoother variational(mesh);

    testSmoother(mesh, variational, TRI3);
  }

  void testVariationalQuadMultipleSubdomains()
  {
    ReplicatedMesh mesh(*TestCommWorld);
    VariationalMeshSmoother variational(mesh);

    testSmoother(mesh, variational, QUAD4, true);
  }
#endif // LIBMESH_ENABLE_VSMOOTHER
};


CPPUNIT_TEST_SUITE_REGISTRATION( MeshSmootherTest );
