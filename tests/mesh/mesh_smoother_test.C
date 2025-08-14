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
#include "libmesh/face_tri.h"
#include "libmesh/utility.h"
#include "libmesh/enum_to_string.h"
#include "libmesh/parallel_ghost_sync.h"

#include "test_comm.h"
#include "libmesh_cppunit.h"

#include <random>

namespace
{
using namespace libMesh;
using Utility::pow;

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
    output(0) = p(0) + (pow<3>(eta) - eta) * p(1) * (1 - p(1));
    output(1) = p(1) + (pow<3>(zeta) - zeta) * p(0) * (1 - p(0));
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
    std::array<bool, 3> is_on_boundary = {false, false, false};
    for (unsigned int i = 0; i < _dim; ++i)
      if (std::abs(p(i)) < TOLERANCE || std::abs(p(i) - 1.) < TOLERANCE)
        is_on_boundary[i] = true;

    // Distort only those directions not fixed on the boundary
    for (unsigned int i = 0; i < _dim; ++i)
      {
        if (!is_on_boundary[i]) // only distort free dimensions
          {
            // This value controls the strength of the distortion
            Real modulation = 0.3;
            Real xi = 2. * p(i) - 1.;
            for (unsigned int j = 0; j < _dim; ++j)
              {
                if (j != i)
                  {
                    Real pj = std::clamp(p(j), 0., 1.);         // ensure numeric safety
                    modulation *= (pj - 0.5) * (pj - 0.5) * 4.; // quadratic bump centered at 0.5
                  }
              }
            const auto delta = (pow<3>(xi) - xi) * modulation;
            // Check for delta = 0 to make sure we perturb every point
            output(i) = (delta > TOLERANCE) ? p(i) + delta : 1.05 * p(i);
          }
        else
          output(i) = p(i); // dimension on boundary remains unchanged
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
    output(1) = 0.5 * std::sqrt(Real(3)) * p(1);
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
    output(0) = p(0) + p(1) / std::sqrt(Real(3));
    output(1) = (2. / std::sqrt(Real(3))) * p(1);
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
  CPPUNIT_TEST(testVariationalEdge2);
  CPPUNIT_TEST(testVariationalEdge3);
  CPPUNIT_TEST(testVariationalEdge3MultipleSubdomains);

  CPPUNIT_TEST(testVariationalQuad);
  CPPUNIT_TEST(testVariationalQuadMultipleSubdomains);
  CPPUNIT_TEST(testVariationalQuadTangled);

  CPPUNIT_TEST(testVariationalTri3);
  CPPUNIT_TEST(testVariationalTri6);
  CPPUNIT_TEST(testVariationalTri6MultipleSubdomains);

  CPPUNIT_TEST(testVariationalHex8);
  CPPUNIT_TEST(testVariationalHex20);
  CPPUNIT_TEST(testVariationalHex27);
  CPPUNIT_TEST(testVariationalHex27Tangled);
  CPPUNIT_TEST(testVariationalHex27MultipleSubdomains);
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
                               VariationalMeshSmoother & smoother,
                               const ElemType type,
                               const bool multiple_subdomains = false,
                               const bool tangle_mesh=false)
  {
    LOG_UNIT_TEST;

    if (multiple_subdomains && tangle_mesh)
      libmesh_not_implemented_msg(
          "Arbitrary mesh tangling with multiple subdomains is not supported.");

    // Get mesh dimension, determine whether type is triangular
    const auto * ref_elem = &(ReferenceElem::get(type));
    const auto dim = ref_elem->dim();
    const bool type_is_tri = Utility::enum_to_string(type).compare(0, 3, "TRI") == 0;

    unsigned int n_elems_per_side = 4;

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

    const auto & boundary_info = mesh.get_boundary_info();

    if (tangle_mesh)
      {
        bool swapped_nodes = false;
        for (auto * elem : mesh.active_local_element_ptr_range())
          {
            // The solver has trouble if the mesh is too tangled, so only tangle
            // some elements
            const auto centroid = elem->vertex_average();
            if ((centroid - Point(0.5, (dim > 1) ? 0.5 : 0.0, (dim > 2) ? 0.5 : 0.)).norm() < 0.25)
              {
                // Identify the first two non-boundary nodes in elem
                dof_id_type node1_id = DofObject::invalid_id;
                dof_id_type node2_id = DofObject::invalid_id;
                for (const auto & node_id : make_range(elem->n_nodes()))
                  {
                    // Get boundary ids associated with the node
                    std::vector<boundary_id_type> boundary_ids;
                    const auto * node_ptr = elem->node_ptr(node_id);
                    boundary_info.boundary_ids(node_ptr, boundary_ids);

                    // Skip if boundary node
                    if (boundary_ids.size())
                      continue;

                    if (node1_id == DofObject::invalid_id)
                      node1_id = node_id;
                    else if (node2_id == DofObject::invalid_id)
                      {
                        node2_id = node_id;
                        break;
                      }
                  }

                // Swap the nodes
                if ((node1_id != DofObject::invalid_id) && (node2_id != DofObject::invalid_id))
                  {
                    // Make sure we haven't swapped any nodes yet on this or other processors
                    mesh.comm().max(swapped_nodes);
                    if (swapped_nodes)
                      break;

                    auto & node1 = elem->node_ref(node1_id);
                    auto & node2 = elem->node_ref(node2_id);

                    for (const auto & d : make_range(dim))
                      {
                        auto tmp = node1(d);
                        node1(d) = node2(d);
                        node2(d) = tmp;
                      }

                    // Once we have swapped two elements of the mesh, do not swap
                    // any more. Again, we still have issues with "too tangled"
                    // meshes.
                    swapped_nodes = true;
                  }
              }
          }

        SyncNodalPositions sync_object(mesh);
        Parallel::sync_dofobject_data_by_id (mesh.comm(), mesh.nodes_begin(), mesh.nodes_end(), sync_object);
      }

    // Add multiple subdomains if requested
    std::unordered_map<subdomain_id_type, Real> distorted_subdomain_volumes;
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
        // recording the final subdomain volumes.
        for (auto * elem : mesh.active_element_ptr_range())
          distorted_subdomain_volumes[elem->subdomain_id()] += elem->volume();
      }


      // Get the mesh order
      const auto &elem_orders = mesh.elem_default_orders();
      libmesh_error_msg_if(
          elem_orders.size() != 1,
          "The variational smoother cannot be used for mixed-order meshes!");
      const auto fe_order = *elem_orders.begin();

      // Function to assert the distortion is as expected
      auto distortion_is = [
        &n_elems_per_side, &dim, &boundary_info, &fe_order
      ](const Node &node, bool distortion, Real distortion_tol = TOLERANCE)
      {
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

        std::size_t num_zero_or_one = 0;

        bool distorted = false;
        for (const auto d : make_range(dim))
          {
            const Real r = node(d);
            const Real R = r * n_elems_per_side * fe_order;

            const bool d_distorted = std::abs(R - std::round(R)) > distortion_tol;
            distorted |= d_distorted;
            num_zero_or_one +=
                (absolute_fuzzy_equals(r, 0.) || absolute_fuzzy_equals(r, 1.));
          }

        CPPUNIT_ASSERT_GREATEREQUAL(dim - num_dofs, num_zero_or_one);

        // We can never expect a fixed node to be distorted
        if (num_dofs == 0)
          // if (num_dofs < dim)
          return true;
        return distorted == distortion;

    };

    // Make sure our DistortHyperCube transformation has distorted the mesh
    for (auto node : mesh.node_ptr_range())
      CPPUNIT_ASSERT(distortion_is(*node, true));

    // Make sure our DistortHyperCube transformation has tangled the mesh
    if (tangle_mesh)
    {
      smoother.setup(); // Need this to create the system we are about to access
      const auto & unsmoothed_info = smoother.get_mesh_info();
      CPPUNIT_ASSERT(unsmoothed_info.mesh_is_tangled);
    }

    // Transform the square mesh of triangles to a parallelogram mesh of
    // triangles. This will allow the Variational Smoother to smooth the mesh
    // to the optimal case of equilateral triangles
    if (type_is_tri)
      {
        SquareToParallelogram stp;
        MeshTools::Modification::redistribute(mesh, stp);
      }

    smoother.smooth();

    // Transform the parallelogram mesh back to a square mesh. In the case of
    // the Variational Smoother, equilateral triangular elements will be
    // transformed into right triangular elements that align with the original
    // undistorted mesh.
    if (type_is_tri)
      {
        ParallelogramToSquare pts;
        MeshTools::Modification::redistribute(mesh, pts);
      }

    if (multiple_subdomains)
      {
        // Make sure the subdomain volumes didn't change. Although this doesn't
        // guarantee that the subdomain didn't change, it is a good indicator
        // and is friedly to sliding subdomain boundary nodes.
        std::unordered_map<subdomain_id_type, Real> smoothed_subdomain_volumes;
        for (auto * elem : mesh.active_element_ptr_range())
          smoothed_subdomain_volumes[elem->subdomain_id()] += elem->volume();

        for (const auto & pair : distorted_subdomain_volumes)
          {
            const auto & subdomain_id = pair.first;
            CPPUNIT_ASSERT(
                relative_fuzzy_equals(libmesh_map_find(distorted_subdomain_volumes, subdomain_id),
                                      libmesh_map_find(smoothed_subdomain_volumes, subdomain_id)));
          }
      }
    else
      // Make sure we're not too distorted anymore
      for (auto node : mesh.node_ptr_range())
        CPPUNIT_ASSERT(distortion_is(*node, false, 1e-2));
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
  void testVariationalEdge2()
  {
    ReplicatedMesh mesh(*TestCommWorld);
    VariationalMeshSmoother variational(mesh);

    testVariationalSmoother(mesh, variational, EDGE2);
  }

  void testVariationalEdge3()
  {
    ReplicatedMesh mesh(*TestCommWorld);
    VariationalMeshSmoother variational(mesh);

    testVariationalSmoother(mesh, variational, EDGE3);
  }

  void testVariationalEdge3MultipleSubdomains()
  {
    ReplicatedMesh mesh(*TestCommWorld);
    VariationalMeshSmoother variational(mesh);

    testVariationalSmoother(mesh, variational, EDGE3, true);
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

  void testVariationalQuadTangled()
  {
    ReplicatedMesh mesh(*TestCommWorld);
    VariationalMeshSmoother variational(mesh);

    testVariationalSmoother(mesh, variational, QUAD4, false, true);
  }

  void testVariationalTri3()
  {
    ReplicatedMesh mesh(*TestCommWorld);
    VariationalMeshSmoother variational(mesh);

    testVariationalSmoother(mesh, variational, TRI3);
  }

  void testVariationalTri6()
  {
    ReplicatedMesh mesh(*TestCommWorld);
    VariationalMeshSmoother variational(mesh);

    testVariationalSmoother(mesh, variational, TRI6);
  }

  void testVariationalTri6MultipleSubdomains()
  {
    ReplicatedMesh mesh(*TestCommWorld);
    VariationalMeshSmoother variational(mesh);

    testVariationalSmoother(mesh, variational, TRI3, true);
  }

  void testVariationalHex8()
  {
    ReplicatedMesh mesh(*TestCommWorld);
    VariationalMeshSmoother variational(mesh);

    testVariationalSmoother(mesh, variational, HEX8);
  }

  void testVariationalHex20()
  {
    ReplicatedMesh mesh(*TestCommWorld);
    VariationalMeshSmoother variational(mesh);

    testVariationalSmoother(mesh, variational, HEX20);
  }

  void testVariationalHex27()
  {
    ReplicatedMesh mesh(*TestCommWorld);
    VariationalMeshSmoother variational(mesh);

    testVariationalSmoother(mesh, variational, HEX27);
  }

  void testVariationalHex27MultipleSubdomains()
  {
    ReplicatedMesh mesh(*TestCommWorld);
    VariationalMeshSmoother variational(mesh);

    testVariationalSmoother(mesh, variational, HEX27, true);
  }

  void testVariationalHex27Tangled()
  {
    ReplicatedMesh mesh(*TestCommWorld);
    VariationalMeshSmoother variational(mesh);

    testVariationalSmoother(mesh, variational, HEX27, false, true);
  }
#endif // LIBMESH_ENABLE_VSMOOTHER
};

CPPUNIT_TEST_SUITE_REGISTRATION(MeshSmootherTest);
