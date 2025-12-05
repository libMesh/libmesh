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
#include <libmesh/mesh.h>
#include <libmesh/system.h> // LIBMESH_HAVE_SOLVER define
#include "libmesh/face_tri.h"
#include "libmesh/cell_prism.h"
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
            output(i) = (delta > TOLERANCE) ? p(i) + delta : 1.03 * p(i);
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
    output(2) = p(2);
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
    output(2) = p(2);
  }
};

class CubeToParallelepiped : public FunctionBase<Real>
{
  std::unique_ptr<FunctionBase<Real>> clone() const override
  {
    return std::make_unique<CubeToParallelepiped>();
  }

  Real operator()(const Point &, const Real = 0.) override
  {
    libmesh_not_implemented();
  } // scalar-only API

  // Has the effect that a cube, meshed into right triangular prisms with diagonals
  // rising from lower-right to upper-left, is transformed into a right-leaning
  // parallelepiped of equilateral triangular prisms. Additionally, the z direction
  // is scaled to ensure element height to base aspect ratios match the target element.
  // Without the correct aspect ratio, the smoothed mesh nodes are off and the
  // distortion_is assertions fail.
  void operator()(const Point & p, const Real, DenseVector<Real> & output)
  {
    output.resize(3);
    output(0) = p(0) + 0.5 * p(1);
    output(1) = 0.5 * std::sqrt(3) * p(1);
    // Adjusting z to get the triangular area to element height aspect ratio To
    // match the target element
    output(2) = p(2) * 0.25 * std::sqrt(3.);
  }
};

class ParallelepipedToCube : public FunctionBase<Real>
{
  std::unique_ptr<FunctionBase<Real>> clone() const override
  {
    return std::make_unique<ParallelepipedToCube>();
  }

  Real operator()(const Point &, const Real = 0.) override
  {
    libmesh_not_implemented();
  } // scalar-only API

  // Has the effect that a right-leaning parallelepiped of equilaterial triangular
  // prisms with is transformed into a cube of right triangular prisms with
  // diagonals rising from lower-right to upper-left. Additionally, the z direction
  // is scaled to ensure element heights align with the values expected by the
  // distortion_is function.
  // This is the inversion of the CubeToParallelepiped mapping.
  void operator()(const Point & p, const Real, DenseVector<Real> & output)
  {
    output.resize(3);
    output(0) = p(0) - p(1) / std::sqrt(3.);
    output(1) = (2. / std::sqrt(3)) * p(1);
    // Undoing the aspect ratio adjustment to get the z divisions to work with
    // distortion_is
    output(2) = p(2) * 4. / std::sqrt(3.);
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

  CPPUNIT_TEST(testVariationalQuad4);
  CPPUNIT_TEST(testVariationalQuad4MultipleSubdomains);
  CPPUNIT_TEST(testVariationalQuad8);
  CPPUNIT_TEST(testVariationalQuad9);
  CPPUNIT_TEST(testVariationalQuad9MultipleSubdomains);
  CPPUNIT_TEST(testVariationalQuad4Tangled);

  CPPUNIT_TEST(testVariationalTri3);
  CPPUNIT_TEST(testVariationalTri6);
  CPPUNIT_TEST(testVariationalTri6MultipleSubdomains);

  CPPUNIT_TEST(testVariationalHex8);
  CPPUNIT_TEST(testVariationalHex20);
  CPPUNIT_TEST(testVariationalHex27);
  CPPUNIT_TEST(testVariationalHex27Tangled);
  CPPUNIT_TEST(testVariationalHex27MultipleSubdomains);

  CPPUNIT_TEST(testVariationalPyramid5);
  CPPUNIT_TEST(testVariationalPyramid13);
  CPPUNIT_TEST(testVariationalPyramid14);
  CPPUNIT_TEST(testVariationalPyramid18);
  CPPUNIT_TEST(testVariationalPyramid18MultipleSubdomains);

  CPPUNIT_TEST(testVariationalPrism6);
  CPPUNIT_TEST(testVariationalPrism15);
  CPPUNIT_TEST(testVariationalPrism18);
  CPPUNIT_TEST(testVariationalPrism20);
  CPPUNIT_TEST(testVariationalPrism21);
  CPPUNIT_TEST(testVariationalPrism21MultipleSubdomains);

  CPPUNIT_TEST(testVariationalTet4);
  CPPUNIT_TEST(testVariationalTet10);
  CPPUNIT_TEST(testVariationalTet14);
  CPPUNIT_TEST(testVariationalTet14MultipleSubdomains);

#if defined(LIBMESH_HAVE_GZSTREAM)
  CPPUNIT_TEST(testVariationalMixed2D);
  CPPUNIT_TEST(testVariationalMixed3D);
#endif

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

  // Helper function to determine how many dimensions of a given point lie at
  // the center and faces of a sub-cube
  std::pair<unsigned int, unsigned int>
  numCenteredAndFacedDimensions(const Point & point, const Real & side_length, const Real & tol)
  {
    const auto half_length = 0.5 * side_length;
    unsigned int num_centered = 0;
    unsigned int num_face = 0;
    for (const auto d : make_range(3))
      {
        const auto num_half_lengths = point(d) / half_length;
        if (std::abs(num_half_lengths - std::round(num_half_lengths)) > tol)
          // This coordinante does not occur at a sub-cube face or center
          continue;

        if (std::lround(num_half_lengths) % 2)
          // An odd number of half lengths means point(d) is at a sub-cube center
          ++num_centered;
        else
          // An odd number of half lengths means point(d) is on a sub-cube face
          ++num_face;
      }

    return std::make_pair(num_centered, num_face);
  }

  // Helper function to determine whether a given point lies at the center of
  // a sub-cube
  bool pointIsCubeCenter(const Point & point, const Real & side_length, const Real & tol = TOLERANCE)
  {
    return numCenteredAndFacedDimensions(point, side_length, tol).first == 3;
  }

  // Helper function to determine whether a given point lies at the vertex of
  // a sub-cube
  bool pointIsCubeVertex(const Point & point, const Real & side_length, const Real & tol = TOLERANCE)
  {
    return numCenteredAndFacedDimensions(point, side_length, tol).second == 3;
  }

  // Helper function to determine whether a given point lies at the center of
  // a sub-cube face
  bool pointIsCubeFaceCenter(const Point & point, const Real & side_length, const Real & tol = TOLERANCE)
  {
    const auto result = numCenteredAndFacedDimensions(point, side_length, tol);
    return (result.first == 2) && (result.second == 1);
  }

  void testVariationalSmoother(Mesh &mesh,
                               VariationalMeshSmoother &smoother,
                               const ElemType type,
                               const bool multiple_subdomains = false,
                               const bool tangle_mesh = false,
                               const Real tangle_damping_factor = 1.0)
  {
    LOG_UNIT_TEST;

    if (multiple_subdomains && tangle_mesh)
      libmesh_not_implemented_msg(
          "Arbitrary mesh tangling with multiple subdomains is not supported.");

    // Get mesh dimension, determine whether type is triangular
    const auto * ref_elem = &(ReferenceElem::get(type));
    const auto dim = ref_elem->dim();
    const bool type_is_tri = Utility::enum_to_string(type).compare(0, 3, "TRI") == 0;
    const bool type_is_pyramid = Utility::enum_to_string(type).compare(0, 7, "PYRAMID") == 0;
    const bool type_is_prism = Utility::enum_to_string(type).compare(0, 5, "PRISM") == 0;
    const bool type_is_tet = Utility::enum_to_string(type).compare(0, 3, "TET") == 0;

    // Used fewer elems for higher order types, as extra midpoint nodes will add
    // enough complexity
    unsigned int n_elems_per_side = 4 / Elem::type_to_default_order_map[type];
    const auto side_length = 1.0 / n_elems_per_side;

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
    const auto & proc_id = mesh.processor_id();

    if (tangle_mesh)
      {
        // We tangle the mesh by (partially) swapping the locations of 2 nodes.
        // If the n-dimentional hypercube is represented as a grid with
        // (undistorted) divisions occuring at integer numbers, we will
        // (partially) swap the nodes nearest p1 = (1,1,1) and p2 = (2,1,2).

        // Define p1 and p2 given the integer indices
        // Undistorted elem side length
        const Real dr = 1. / (n_elems_per_side * Elem::type_to_default_order_map[type]);
        const Point p1 = Point(dr, dim > 1 ? dr : 0, dim > 2 ? dr : 0);
        const Point p2 = Point(2. * dr, dim > 1 ? dr : 0, dim > 2 ? 2. * dr : 0);

        // We need to find the nodes of the mesh that are closest to p1 and p2.
        // This will require some parallel communication for distributed meshes.

        // Distances between local mesh nodes and p1
        std::map<dof_id_type, Real> dist1_map;
        // Distances between local mesh nodes and p2
        std::map<dof_id_type, Real> dist2_map;

        for (const auto * node : mesh.local_node_ptr_range())
          {
          dist1_map[node->id()] = (*node - p1).norm();
          dist2_map[node->id()] = (*node - p2).norm();
          }

          // Helper function to analyze a dist_map accross all processors and
          // return the node id and spatial coordinates for the node with the
          // smallest distance
          auto get_closet_point_accross_all_procs =
              [&mesh, &proc_id](const std::map<dof_id_type, Real> &dist_map) {

                processor_id_type broadcasting_proc = 0;
                Real d_min_local = std::numeric_limits<Real>::max();
                dof_id_type node_id = DofObject::invalid_id;

                // Only execute this if this proc owns any nodes. Otherwise,
                // use default values defined above
                if (dist_map.size())
                {
                  // Iterator to the map entry with the smallest distance
                  auto min_it =
                      std::min_element(dist_map.begin(), dist_map.end(),
                                       [](const auto &a, const auto &b) {
                                         return a.second < b.second;
                                       });

                  node_id = min_it->first;
                  d_min_local = min_it->second;
                }

                // Get the smallest distance accross all procs
                auto d_min_global = d_min_local;
                mesh.comm().min(d_min_global);

                // Find the proc id, node_id, and spatial coordinantes for the
                // node with the smallest distance
                Point node;
                if (d_min_local == d_min_global) {
                  broadcasting_proc = proc_id;
                  node = mesh.node_ref(node_id);
                }

                mesh.comm().max(broadcasting_proc);
                // Broadcast information about this node to all procs
                mesh.comm().broadcast(node_id, broadcasting_proc);
                mesh.comm().broadcast(node, broadcasting_proc);

                return std::pair(node_id, node);
              };

          const auto [node_id1, node1] =
              get_closet_point_accross_all_procs(dist1_map);
          const auto [node_id2, node2] =
              get_closet_point_accross_all_procs(dist2_map);

          // modify the nodes if they are on this proc
          const auto displacement = tangle_damping_factor * (node1 - node2);
          if (mesh.query_node_ptr(node_id1))
            mesh.node_ref(node_id1) -= displacement;
          if (mesh.query_node_ptr(node_id2))
            mesh.node_ref(node_id2) += displacement;

          SyncNodalPositions sync_object(mesh);
          Parallel::sync_dofobject_data_by_id(mesh.comm(), mesh.nodes_begin(),
                                              mesh.nodes_end(), sync_object);

          // Make sure the mesh is tangled
          // Need this to create the system we are about to access
          smoother.setup();
          const auto &unsmoothed_info = smoother.get_mesh_info();
          CPPUNIT_ASSERT(unsmoothed_info.mesh_is_tangled);
      }

    // Add multiple subdomains if requested
    std::unordered_map<subdomain_id_type, Real> distorted_subdomain_volumes;
    subdomain_id_type highest_subdomain_id = 0;
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
          for (auto *elem : mesh.active_element_ptr_range()) {
            const auto sub_id = elem->subdomain_id();
            // Don't double count an element's volume on multiple procs
            if (elem->processor_id() != proc_id)
              continue;
            distorted_subdomain_volumes[sub_id] += elem->volume();
            highest_subdomain_id = std::max(sub_id, highest_subdomain_id);
          }

          // Parallel communication
          mesh.comm().max(highest_subdomain_id);
          for (const auto sub_id : make_range(highest_subdomain_id + 1)) {
            // Make sure sub_id exists in the map on this proc
            if (distorted_subdomain_volumes.find(sub_id) ==
                distorted_subdomain_volumes.end())
              distorted_subdomain_volumes[sub_id] = 0.;

            mesh.comm().sum(distorted_subdomain_volumes[sub_id]);
          }

        // We've just invalidated the get_mesh_subdomains() cache by
        // adding a new one; fix it.
        mesh.cache_elem_data();
      }

    // Get the mesh order
    const auto & elem_orders = mesh.elem_default_orders();
    libmesh_error_msg_if(elem_orders.size() != 1,
                         "The variational smoother cannot be used for mixed-order meshes!");
    // For pyramids, the factor 2 accounts for the account that cubes of side
    // length s are divided into pyramids of base side length s and height s/2.
    // The factor 4 lets us catch triangular face midpoint nodes for PYRAMID18 elements.
    // Similar reasoning for tets.
    const auto scale_factor = *elem_orders.begin() * ((type_is_pyramid || type_is_tet) ? 2 * 4 : 1);

    // Function to assert the node distortion is as expected
    auto node_distortion_is = [&n_elems_per_side, &dim, &boundary_info, &scale_factor, &type_is_prism](
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

      std::size_t num_zero_or_one = 0;

      bool distorted = false;
      for (const auto d : make_range(dim))
        {
          const Real r = node(d);
          const Real R = r * n_elems_per_side * scale_factor;
          CPPUNIT_ASSERT_GREATER(-distortion_tol * distortion_tol, r);
          CPPUNIT_ASSERT_GREATER(-distortion_tol * distortion_tol, 1 - r);

          bool d_distorted = std::abs(R - std::round(R)) > distortion_tol;
          if (type_is_prism && (scale_factor == 3))
            {
              // Adjustment required for triangular face nodes of PRISM20/21 elements.
              // These elements have fe_order 3. This takes care of nodes occuring at
              // thirds, but not halves.
              const Real R_prism = R / scale_factor * 2;
              const bool d_distorted_prism =
                  std::abs(R_prism - std::round(R_prism)) > distortion_tol;
              d_distorted &= d_distorted_prism;
            }
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

    // Make sure our DistortHyperCube transformation has distorted the mesh
    for (auto node : mesh.node_ptr_range())
      CPPUNIT_ASSERT(node_distortion_is(*node, true));

    // Transform the square mesh of triangles to a parallelogram mesh of
    // triangles. This will allow the Variational Smoother to smooth the mesh
    // to the optimal case of equilateral triangles
    if (type_is_tri)
      {
        // Transform the square mesh of triangles to a parallelogram mesh of
        // triangles. This will allow the Variational Smoother to smooth the mesh
        // to the optimal case of equilateral triangles
        SquareToParallelogram stp;
        MeshTools::Modification::redistribute(mesh, stp);
      }
    else if (type_is_prism)
      {
        // Transform the cube mesh of prisms to a parallelepiped mesh of
        // prisms. This will allow the Variational Smoother to smooth the mesh
        // to the optimal case of equilateral triangular prisms
        CubeToParallelepiped ctp;
        MeshTools::Modification::redistribute(mesh, ctp);
      }

    smoother.smooth();

    if (type_is_tri)
      {
        // Transform the parallelogram mesh back to a square mesh. In the case of
        // the Variational Smoother, equilateral triangules will be transformed
        // into right triangules that align with the original undistorted mesh.
        ParallelogramToSquare pts;
        MeshTools::Modification::redistribute(mesh, pts);
      }
    else if (type_is_prism)
      {
        // Transform the parallelepiped mesh back to a cube mesh. In the case of
        // the Variational Smoother, equilateral triangular prisms will be transformed
        // into right triangular prisms that align with the original undistorted mesh.
        ParallelepipedToCube ptc;
        MeshTools::Modification::redistribute(mesh, ptc);
      }

    if (multiple_subdomains)
      {
        // Make sure the subdomain volumes didn't change. Although this doesn't
        // guarantee that the subdomain didn't change, it is a good indicator
        // and is friedly to sliding subdomain boundary nodes.
        std::unordered_map<subdomain_id_type, Real> smoothed_subdomain_volumes;
        for (auto *elem : mesh.active_element_ptr_range()) {
          // Don't double count an element's volume on multiple procs
          if (elem->processor_id() != proc_id)
            continue;
          smoothed_subdomain_volumes[elem->subdomain_id()] += elem->volume();
        }

        // Parallel communication
        for (const auto sub_id : make_range(highest_subdomain_id + 1)) {
          // Make sure sub_id exists in the map on this proc
          if (smoothed_subdomain_volumes.find(sub_id) ==
              smoothed_subdomain_volumes.end())
            smoothed_subdomain_volumes[sub_id] = 0.;

          mesh.comm().sum(smoothed_subdomain_volumes[sub_id]);
        }

        for (const auto sub_id : make_range(highest_subdomain_id + 1))
          CPPUNIT_ASSERT(relative_fuzzy_equals(
              libmesh_map_find(distorted_subdomain_volumes, sub_id),
              libmesh_map_find(smoothed_subdomain_volumes, sub_id), TOLERANCE));
      }
    else
      {
        std::unordered_map<dof_id_type, std::vector<const Elem *>> nodes_to_elem_map;
        MeshTools::build_nodes_to_elem_map(mesh, nodes_to_elem_map);

        // Make sure we're not too distorted anymore, using the given
        // tolerance.
        std::set<dof_id_type> nodes_checked;
        const Real tol = TOLERANCE;

        for (const auto * elem : mesh.active_element_ptr_range())
          {
            for (const auto local_node_id : make_range(elem->n_nodes()))
              {
                const auto & node = elem->node_ref(local_node_id);
                if (nodes_checked.find(node.id()) != nodes_checked.end())
                  continue;

                nodes_checked.insert(node.id());

                // Check for special case where pyramidal base-to-apex midpoint
                // nodes are not smoothed to the actual midpoints.
                if (type_is_pyramid)
                  {
                    if (local_node_id > 8 && local_node_id < 13)
                      {
                        // Base-Apex midpoint nodes
                        //
                        // Due to the nature of the pyramidal target element and the
                        // cubic nature of the test mesh, for a dilation weight o
                        // 0.5, smoothed base-apex midpoint nodes are smoothed to
                        // the point base + x * (apex - base) instead of
                        // base + 0.5 (apex - base), where x is in (0,1) and
                        // depends on the number of nodes in the element.
                        // We hard-code this check below.

                        const auto & base = elem->node_ref(local_node_id - 9);
                        const auto & apex = elem->node_ref(4);
                        const Real x = (type == PYRAMID18) ? 0.56646084 : 0.54985875;

                        CPPUNIT_ASSERT(node.absolute_fuzzy_equals(base + x * (apex - base), tol));
                        continue;
                      }
                    else if (local_node_id > 13)
                      {
                        // Triangular face midpoint nodes
                        //
                        // Due to the nature of the pyramidal target element and the
                        // cubic nature of the test mesh, for a dilation weight o
                        // 0.5, smoothed triangular face midpoint nodes are
                        // smoothed a weighted average of the constituent
                        // vertices instead of an unweighted average.
                        // We hard-code this check below.
                        const auto & base1 = elem->node_ref(local_node_id - 14);
                        const auto & base2 = elem->node_ref((local_node_id - 13) % 4);
                        const auto & apex = elem->node_ref(4);

                        const auto node_approx =
                            (0.31401599 * base1 + 0.31401599 * base2 + 0.37196802 * apex);
                        CPPUNIT_ASSERT(node.absolute_fuzzy_equals(node_approx, tol));
                        continue;
                      }
                  }

                // Check for special case where some tet midpoint nodes are not
                // smoothed to the actual midpoints.
                else if (type_is_tet && !elem->is_vertex(local_node_id))
                  {
                    // We have a non-vertex node. Determine what "type" of
                    // midpoint node with respect to the mesh geometry.
                    // First, get the nodes that neighbor this node
                    std::vector<const Node *> neighbors;
                    MeshTools::find_nodal_or_face_neighbors(
                        mesh, node, nodes_to_elem_map, neighbors);

                    switch (neighbors.size())
                      {
                          case 2: {
                            // Midpoint node
                            // Determine what "type" of midpoint node

                            // Is one of the neighbors at the center of a sub-cube?
                            const auto is_0_cube_center =
                                pointIsCubeCenter(*neighbors[0], side_length);
                            const auto is_1_cube_center =
                                pointIsCubeCenter(*neighbors[1], side_length);
                            libmesh_assert(!(is_0_cube_center && is_1_cube_center));

                            if (is_0_cube_center || is_1_cube_center)
                              {
                                const auto & cube_center =
                                    is_0_cube_center ? *neighbors[0] : *neighbors[1];
                                const auto & other =
                                    is_0_cube_center ? *neighbors[1] : *neighbors[0];

                                // Since one neighbor is a sub-cube center, the
                                // other will either be a sub-cube face center or
                                // a sub-cube vertex. Determine which one.
                                if (pointIsCubeFaceCenter(other, side_length))
                                  {
                                    const Real x = (type == TET10) ? 0.42895041 : 0.41486385;
                                    CPPUNIT_ASSERT(node.absolute_fuzzy_equals(
                                        other + x * (cube_center - other), tol));
                                  }

                                else if (pointIsCubeVertex(other, side_length))
                                  {
                                    const Real x = (type == TET10) ? 0.55388920 : 0.58093516;
                                    CPPUNIT_ASSERT(node.absolute_fuzzy_equals(
                                        other + x * (cube_center - other), tol));
                                  }
                              }

                            else
                              {
                                // The remaining possibilities are
                                // cube-vertex-to-cube-vertex midpoints and
                                // cube-vertex-to-cube-face-center midpoints.
                                const auto is_0_cube_vertex =
                                    pointIsCubeVertex(*neighbors[0], side_length);
                                const auto is_1_cube_vertex =
                                    pointIsCubeVertex(*neighbors[1], side_length);
                                const auto is_0_cube_face_center =
                                    pointIsCubeFaceCenter(*neighbors[0], side_length);
                                const auto is_1_cube_face_center =
                                    pointIsCubeFaceCenter(*neighbors[1], side_length);

                                if (is_0_cube_vertex && is_1_cube_vertex)
                                  // Due to symmetry, this type of midpoint is
                                  // the geometric midpoint. Let the
                                  // node_distortion_is function check it
                                  CPPUNIT_ASSERT(node_distortion_is(node, false));

                                else
                                  {
                                    libmesh_error_msg_if(
                                        (is_0_cube_center || is_0_cube_face_center) &&
                                            (is_1_cube_center || is_1_cube_face_center),
                                        "We should never get here!");
                                    const auto & cube_vertex =
                                        is_0_cube_vertex ? *neighbors[0] : *neighbors[1];
                                    const auto & cube_face_center =
                                        is_0_cube_face_center ? *neighbors[0] : *neighbors[1];
                                    const Real x = (type == TET10) ? 0.61299101 : 0.65125580;
                                    CPPUNIT_ASSERT(node.absolute_fuzzy_equals(
                                        cube_vertex + x * (cube_face_center - cube_vertex), tol));
                                  }
                              }

                            continue;
                            break;
                          }

                          case 6: {
                            // Face node
                            // Only keep the face vertex nodes
                            neighbors.erase(std::remove_if(neighbors.begin(),
                                                           neighbors.end(),
                                                           [&elem](const Node * n) {
                                                             return !elem->is_vertex(
                                                                 elem->local_node(n->id()));
                                                           }),
                                            neighbors.end());

                            // There are three types of faces:
                            //
                            // 1) Isoscelese triangle with vertices at two
                            //    sub-cube vertices and the sub-cube center,
                            // 2) Isosceles triangle with vertices at two
                            //    sub-cube vertices and a sub-cube face center,
                            // 3) Scalene triangle with vertices at a sub-cube
                            //    center, a sub-cube vertex, and a sub-cube face
                            //    center.
                            //
                            // Determine which kind of face this node belongs to based on the vertex
                            // neighbors

                            unsigned int num_vertices_at_cube_center = 0;
                            unsigned int num_vertices_at_cube_vertices = 0;
                            unsigned int num_vertices_at_cube_face_centers = 0;
                            for (const auto * neighbor : neighbors)
                              {
                                if (pointIsCubeCenter(*neighbor, side_length))
                                  num_vertices_at_cube_center += 1;
                                else if (pointIsCubeVertex(*neighbor, side_length))
                                  num_vertices_at_cube_vertices += 1;
                                else if (pointIsCubeFaceCenter(*neighbor, side_length))
                                  num_vertices_at_cube_face_centers += 1;
                              }

                            // We will express the face node at a linear
                            // combination of the face vertices
                            Node node_approx = Node(0, 0, 0);

                            // Isosceles triangular face
                            if (num_vertices_at_cube_vertices == 2)
                              {
                                // Isosceles triangular face (type 1)
                                if (num_vertices_at_cube_center)
                                  for (const auto * neighbor : neighbors)
                                    {
                                      Real weight;
                                      if (pointIsCubeVertex(*neighbor, side_length))
                                        weight = 0.30600747;
                                      else if (pointIsCubeCenter(*neighbor, side_length))
                                        weight = 0.38798506;
                                      else
                                        libmesh_error_msg("We should never get here!");

                                      node_approx += weight * (*neighbor);
                                    }

                                // Isosceles triangular face (type 2)
                                else if (num_vertices_at_cube_face_centers)
                                  for (const auto * neighbor : neighbors)
                                    {
                                      Real weight;
                                      if (pointIsCubeVertex(*neighbor, side_length))
                                        weight = 0.28078090;
                                      else if (pointIsCubeFaceCenter(*neighbor, side_length))
                                        weight = 0.43843820;
                                      else
                                        libmesh_error_msg("We should never get here!");

                                      node_approx += weight * (*neighbor);
                                    }

                                else
                                  libmesh_error_msg("We should never get here!");
                              }

                            // Scalene triangular face (type 3)
                            else if (num_vertices_at_cube_center && num_vertices_at_cube_vertices &&
                                     num_vertices_at_cube_face_centers)
                              for (const auto * neighbor : neighbors)
                                {
                                  Real weight;
                                  if (pointIsCubeCenter(*neighbor, side_length))
                                    weight = 0.33102438;
                                  else if (pointIsCubeVertex(*neighbor, side_length))
                                    weight = 0.27147230;
                                  else if (pointIsCubeFaceCenter(*neighbor, side_length))
                                    weight = 0.39750332;
                                  else
                                    libmesh_error_msg("We should never get here!");

                                  node_approx += weight * (*neighbor);
                                }

                            else
                              libmesh_error_msg("We should never get here!");

                            CPPUNIT_ASSERT(node.absolute_fuzzy_equals(node_approx, tol));

                            continue;
                            break;
                          }
                          default: {
                            libmesh_error_msg(node << " has unexpected number of neighbors ("
                                                   << neighbors.size() << ")");
                            break;
                          }
                      } // switch (neighbors.size())
                  } // if (type_is_tet)

                CPPUNIT_ASSERT(node_distortion_is(node, false));

              }
          }
      }
  }

  // Function that distorts and smooths a mesh and checks that the result is
  // equal to the original "gold" mesh
  void testVariationalSmootherRegression(const ReplicatedMesh & gold_mesh)
  {
    // Make a copy of the gold mesh to distort and then smooth
    ReplicatedMesh mesh(gold_mesh);

    // Move it around so we have something that needs smoothing
    DistortHyperCube dh(mesh.mesh_dimension());
    MeshTools::Modification::redistribute(mesh, dh);

    // Function to assert the distortion is as expected
    auto mesh_distortion_is = [](const ReplicatedMesh & gold_mesh,
                                 const ReplicatedMesh & mesh,
                                 const bool distortion,
                                 Real distortion_tol = TOLERANCE) {

      CPPUNIT_ASSERT(gold_mesh.n_nodes() ==  mesh.n_nodes());
      for (const auto node_id : make_range(gold_mesh.n_nodes()))
        {
          const auto & gold_node = gold_mesh.node_ref(node_id);
          const auto & node = mesh.node_ref(node_id);
          for (const auto d : make_range(gold_mesh.mesh_dimension()))
            {
              const bool d_distorted = std::abs(gold_node(d) - node(d)) > distortion_tol;
              if (d_distorted)
                // mesh is distorted from gold_mesh
                return distortion;
            }
        }

      // mesh is not distorted from gold_mesh
      return !distortion;
    };

    // Make sure the mesh has been distorted
    CPPUNIT_ASSERT(mesh_distortion_is(gold_mesh, mesh, true));

    // Turn off subdomain boundary preservation because DistortHyperCube
    // does not preserve subdomain boundaries
    // Also set dilation coefficient to 0 because the reference volume of the
    // distorted mesh won't be equal to the reference volume of the smoothed
    // mesh and will smooth to a slightly different solution.
    VariationalMeshSmoother smoother(mesh, 0.0, false);
    smoother.smooth();

    // Make sure the mesh has been smoothed to the gold mesh
    CPPUNIT_ASSERT(mesh_distortion_is(gold_mesh, mesh, false));
  }

  void testLaplaceQuad()
  {
    ReplicatedMesh mesh(*TestCommWorld);
    LaplaceMeshSmoother laplace(mesh, 8);

    testLaplaceSmoother(mesh, laplace, QUAD4);
  }

  void testLaplaceTri()
  {
    ReplicatedMesh mesh(*TestCommWorld);
    LaplaceMeshSmoother laplace(mesh, 8);

    testLaplaceSmoother(mesh, laplace, TRI3);
  }

#ifdef LIBMESH_ENABLE_VSMOOTHER
  void testVariationalEdge2()
  {
    Mesh mesh(*TestCommWorld);
    VariationalMeshSmoother variational(mesh);

    testVariationalSmoother(mesh, variational, EDGE2);
  }

  void testVariationalEdge3()
  {
    Mesh mesh(*TestCommWorld);
    VariationalMeshSmoother variational(mesh);

    testVariationalSmoother(mesh, variational, EDGE3);
  }

  void testVariationalEdge3MultipleSubdomains()
  {
    Mesh mesh(*TestCommWorld);
    VariationalMeshSmoother variational(mesh);

    testVariationalSmoother(mesh, variational, EDGE3, true);
  }

  void testVariationalQuad4()
  {
    Mesh mesh(*TestCommWorld);
    VariationalMeshSmoother variational(mesh, 0.5, true, TOLERANCE * TOLERANCE, TOLERANCE * TOLERANCE, 100);

    // High verbosity for a 3D mesh to increase code coverage, but silence the output
    // If we want to save the output for processing later, send it somewhere
    // besides nullptr
    std::streambuf * out_buf = libMesh::out.rdbuf(nullptr);

    testVariationalSmoother(mesh, variational, QUAD4);

    // Reset stdout
    libMesh::out.rdbuf(out_buf);
  }

  void testVariationalQuad4MultipleSubdomains()
  {
    Mesh mesh(*TestCommWorld);
    VariationalMeshSmoother variational(mesh);

    testVariationalSmoother(mesh, variational, QUAD4, true);
  }

  void testVariationalQuad8()
  {
    Mesh mesh(*TestCommWorld);
    VariationalMeshSmoother variational(mesh);

    testVariationalSmoother(mesh, variational, QUAD8);
  }

  void testVariationalQuad9()
  {
    Mesh mesh(*TestCommWorld);
    VariationalMeshSmoother variational(mesh);

    testVariationalSmoother(mesh, variational, QUAD9);
  }

  void testVariationalQuad9MultipleSubdomains()
  {
    Mesh mesh(*TestCommWorld);
    VariationalMeshSmoother variational(mesh);

    testVariationalSmoother(mesh, variational, QUAD9, true);
  }

  void testVariationalQuad4Tangled()
  {
    Mesh mesh(*TestCommWorld);
    VariationalMeshSmoother variational(mesh);

    testVariationalSmoother(mesh, variational, QUAD4, false, true, 0.65);
  }

  void testVariationalTri3()
  {
    Mesh mesh(*TestCommWorld);
    VariationalMeshSmoother variational(mesh);

    testVariationalSmoother(mesh, variational, TRI3);
  }

  void testVariationalTri6()
  {
    Mesh mesh(*TestCommWorld);
    VariationalMeshSmoother variational(mesh);

    testVariationalSmoother(mesh, variational, TRI6);
  }

  void testVariationalTri6MultipleSubdomains()
  {
    Mesh mesh(*TestCommWorld);
    VariationalMeshSmoother variational(mesh);

    testVariationalSmoother(mesh, variational, TRI3, true);
  }

  void testVariationalHex8()
  {
    Mesh mesh(*TestCommWorld);
    VariationalMeshSmoother variational(mesh);

    testVariationalSmoother(mesh, variational, HEX8);
  }

  void testVariationalHex20()
  {
    Mesh mesh(*TestCommWorld);
    VariationalMeshSmoother variational(mesh);

    testVariationalSmoother(mesh, variational, HEX20);
  }

  void testVariationalHex27()
  {
    Mesh mesh(*TestCommWorld);
    VariationalMeshSmoother variational(mesh);

    testVariationalSmoother(mesh, variational, HEX27);
  }

  void testVariationalHex27MultipleSubdomains()
  {
    Mesh mesh(*TestCommWorld);
    VariationalMeshSmoother variational(mesh);

    testVariationalSmoother(mesh, variational, HEX27, true);
  }

  void testVariationalHex27Tangled()
  {
    Mesh mesh(*TestCommWorld);
    VariationalMeshSmoother variational(mesh);

    testVariationalSmoother(mesh, variational, HEX27, false, true, 0.55);
  }

  void testVariationalPyramid5()
  {
    Mesh mesh(*TestCommWorld);
    VariationalMeshSmoother variational(mesh);

    testVariationalSmoother(mesh, variational, PYRAMID5);
  }

  void testVariationalPyramid13()
  {
    Mesh mesh(*TestCommWorld);
    VariationalMeshSmoother variational(mesh);

    testVariationalSmoother(mesh, variational, PYRAMID13);
  }

  void testVariationalPyramid14()
  {
    Mesh mesh(*TestCommWorld);
    VariationalMeshSmoother variational(mesh);

    testVariationalSmoother(mesh, variational, PYRAMID14);
  }

  void testVariationalPyramid18()
  {
    Mesh mesh(*TestCommWorld);
    VariationalMeshSmoother variational(mesh);

    testVariationalSmoother(mesh, variational, PYRAMID18);
  }

  void testVariationalPyramid18MultipleSubdomains()
  {
    Mesh mesh(*TestCommWorld);
    VariationalMeshSmoother variational(mesh, 0.0);

    testVariationalSmoother(mesh, variational, PYRAMID18, true);
  }

  void testVariationalPrism6()
  {
    Mesh mesh(*TestCommWorld);
    VariationalMeshSmoother variational(mesh);

    testVariationalSmoother(mesh, variational, PRISM6);
  }

  void testVariationalPrism15()
  {
    Mesh mesh(*TestCommWorld);
    VariationalMeshSmoother variational(mesh);

    testVariationalSmoother(mesh, variational, PRISM15);
  }

  void testVariationalPrism18()
  {
    Mesh mesh(*TestCommWorld);
    VariationalMeshSmoother variational(
        mesh, 0.5, true, TOLERANCE * TOLERANCE, 10 * TOLERANCE * TOLERANCE);

    testVariationalSmoother(mesh, variational, PRISM18);
  }

  void testVariationalPrism20()
  {
    Mesh mesh(*TestCommWorld);
    VariationalMeshSmoother variational(
        mesh, 0.5, true, TOLERANCE * TOLERANCE, 10 * TOLERANCE * TOLERANCE);

    testVariationalSmoother(mesh, variational, PRISM20);
  }

  void testVariationalPrism21()
  {
    Mesh mesh(*TestCommWorld);
    VariationalMeshSmoother variational(
        mesh, 0.5, true, TOLERANCE * TOLERANCE, 10 * TOLERANCE * TOLERANCE);

    testVariationalSmoother(mesh, variational, PRISM21);
  }

  void testVariationalPrism21MultipleSubdomains()
  {
    Mesh mesh(*TestCommWorld);
    mesh.allow_renumbering(false);
    VariationalMeshSmoother variational(mesh);

    testVariationalSmoother(mesh, variational, PRISM21, true);
  }

  void testVariationalTet4()
  {
    Mesh mesh(*TestCommWorld);
    VariationalMeshSmoother variational(mesh, 0.5, true, TOLERANCE * TOLERANCE, TOLERANCE * TOLERANCE, 100);

    // High verbosity for a 3D mesh to increase code coverage, but silence the output
    // If we want to save the output for processing later, send it somewhere
    // besides nullptr
    std::streambuf * out_buf = libMesh::out.rdbuf(nullptr);

    testVariationalSmoother(mesh, variational, TET4);

    // Reset stdout
    libMesh::out.rdbuf(out_buf);
  }

  void testVariationalTet10()
  {
    Mesh mesh(*TestCommWorld);
    VariationalMeshSmoother variational(mesh);

    testVariationalSmoother(mesh, variational, TET10);
  }

  void testVariationalTet14()
  {
    Mesh mesh(*TestCommWorld);
    VariationalMeshSmoother variational(mesh);

    testVariationalSmoother(mesh, variational, TET14);
  }

  void testVariationalTet14MultipleSubdomains()
  {
    Mesh mesh(*TestCommWorld);
    VariationalMeshSmoother variational(mesh);

    testVariationalSmoother(mesh, variational, TET14, true);
  }

  void testVariationalMixed2D()
  {
    ReplicatedMesh mesh(*TestCommWorld);
    mesh.read("meshes/quad4_tri3_smoothed.xda.gz");

    testVariationalSmootherRegression(mesh);
  }

  void testVariationalMixed3D()
  {
    ReplicatedMesh mesh(*TestCommWorld);
    mesh.read("meshes/hex8_prism6_smoothed.xda.gz");

    testVariationalSmootherRegression(mesh);
  }

#endif // LIBMESH_ENABLE_VSMOOTHER
};

CPPUNIT_TEST_SUITE_REGISTRATION(MeshSmootherTest);
