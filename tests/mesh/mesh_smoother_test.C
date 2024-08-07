#include <libmesh/libmesh.h>
#include <libmesh/enum_elem_type.h>
#include <libmesh/function_base.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_modification.h>
#include <libmesh/mesh_smoother_laplace.h>
#include <libmesh/mesh_smoother_vsmoother.h>
#include <libmesh/mesh_tools.h>
#include <libmesh/node.h>
#include <libmesh/replicated_mesh.h>

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
// Yeah we'd like to but this segfaults
//     CPPUNIT_TEST( testVariationalQuad );
//     CPPUNIT_TEST( testVariationalTri );
#  endif // LIBMESH_ENABLE_VSMOOTHER
#endif

  CPPUNIT_TEST_SUITE_END();

public:
  void setUp() {}

  void tearDown() {}

  void testSmoother(ReplicatedMesh & mesh, MeshSmoother & smoother, ElemType type)
  {
    LOG_UNIT_TEST;

    constexpr unsigned int n_elems_per_side = 4;

    MeshTools::Generation::build_square(mesh, n_elems_per_side, n_elems_per_side,
                                        0.,1.,0.,1., type);

    // Move it around so we have something that needs smoothing
    DistortSquare ds;
    MeshTools::Modification::redistribute(mesh, ds);

    // Assert the distortion is as expected
    auto center_distortion_is = []
      (const Node & node, int d, bool distortion,
       Real distortion_tol=TOLERANCE) {
      const Real r = node(d);
      const Real R = r * n_elems_per_side;
      CPPUNIT_ASSERT_GREATER(-TOLERANCE*TOLERANCE, r);
      CPPUNIT_ASSERT_GREATER(-TOLERANCE*TOLERANCE, 1-r);

      // If we're at the boundaries we should *never* be distorted
      if (std::abs(node(0)) < TOLERANCE*TOLERANCE ||
          std::abs(node(0)-1) < TOLERANCE*TOLERANCE)
        {
          const Real R1 = node(1) * n_elems_per_side;
          CPPUNIT_ASSERT_LESS(TOLERANCE*TOLERANCE, std::abs(R1-std::round(R1)));
          return true;
        }

      if (std::abs(node(1)) < TOLERANCE*TOLERANCE ||
          std::abs(node(1)-1) < TOLERANCE*TOLERANCE)
        {
          const Real R0 = node(0) * n_elems_per_side;
          CPPUNIT_ASSERT_LESS(TOLERANCE*TOLERANCE, std::abs(R0-std::round(R0)));

          return true;
        }

      // If we're at the center we're fine
      if (std::abs(r-0.5) < TOLERANCE*TOLERANCE)
        return true;

      return ((std::abs(R-std::round(R)) > distortion_tol) == distortion);
    };

    for (auto node : mesh.node_ptr_range())
      {
        CPPUNIT_ASSERT(center_distortion_is(*node, 0, true));
        CPPUNIT_ASSERT(center_distortion_is(*node, 1, true));
      }

    // Enough iterations to mostly fix us up.  Laplace seems to be at 1e-3
    // tolerance by iteration 6, so hopefully everything is there on any
    // system by 8.
    for (unsigned int i=0; i != 8; ++i)
      smoother.smooth();

    // Make sure we're not too distorted anymore.
    for (auto node : mesh.node_ptr_range())
      {
        CPPUNIT_ASSERT(center_distortion_is(*node, 0, false, 1e-3));
        CPPUNIT_ASSERT(center_distortion_is(*node, 1, false, 1e-3));
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
#endif // LIBMESH_ENABLE_VSMOOTHER
};


CPPUNIT_TEST_SUITE_REGISTRATION( MeshSmootherTest );
