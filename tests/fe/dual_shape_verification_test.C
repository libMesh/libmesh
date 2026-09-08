// libmesh includes
#include "libmesh/libmesh.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/edge_edge2.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/mesh.h"

// unit test includes
#include "test_comm.h"

#include "libmesh_cppunit.h"

using namespace libMesh;

/**
 * This class is for unit testing dual coefficient and shape function values
 */
class DualShapeTest : public CppUnit::TestCase
{
public:
  LIBMESH_CPPUNIT_TEST_SUITE( DualShapeTest );
  CPPUNIT_TEST( testEdge2Lagrange );
  CPPUNIT_TEST( testQuad8Lagrange );
  CPPUNIT_TEST( testTri6Lagrange );
  CPPUNIT_TEST_SUITE_END();

private:
  std::unique_ptr<Mesh> _mesh;
  std::unique_ptr<QGauss> _qrule;
  std::unique_ptr<FEBase> _fe;
  Elem * _elem;

public:

  void testEdge2Lagrange ()
  {
    LOG_UNIT_TEST;
    if (!_elem)
      return;

    _fe->reinit(_elem);

    const DenseMatrix<Real> & dual_coeff = _fe->get_dual_coeff();

    CPPUNIT_ASSERT_EQUAL(dual_coeff.m(), unsigned(2));
    CPPUNIT_ASSERT_EQUAL(dual_coeff.n(), unsigned(2));

    // TOLERANCE*TOLERANCE works with double but not float128
    Real my_tol = TOLERANCE*std::sqrt(TOLERANCE);

    LIBMESH_ASSERT_FP_EQUAL(2, dual_coeff(0,0), my_tol);
    LIBMESH_ASSERT_FP_EQUAL(-1, dual_coeff(0,1), my_tol);
    LIBMESH_ASSERT_FP_EQUAL(-1, dual_coeff(1,0), my_tol);
    LIBMESH_ASSERT_FP_EQUAL(2, dual_coeff(1,1), my_tol);

    const auto & dual_phi = _fe->get_dual_phi();

    CPPUNIT_ASSERT_EQUAL(std::size_t(2), dual_phi.size());

    const auto & qpoints = _qrule->get_points();

    CPPUNIT_ASSERT_EQUAL(qpoints.size(), dual_phi[0].size());

    for (auto qp : index_range(dual_phi[0]))
      LIBMESH_ASSERT_FP_EQUAL(1./2. * (1. - 3.*qpoints[qp](0)), dual_phi[0][qp],
        my_tol);

    CPPUNIT_ASSERT_EQUAL(qpoints.size(), dual_phi[1].size());

    for (auto qp : index_range(dual_phi[1]))
      LIBMESH_ASSERT_FP_EQUAL(1./2. * (1. + 3.*qpoints[qp](0)), dual_phi[1][qp],
        my_tol);
  }

  /**
   * Check what the transformed dual basis on QUAD8/TRI6 guarantees, none of which the untransformed
   * construction gives: strictly positive weights, reproduction of constants, and biorthogonality
   * against Ntilde = T N rather than N. \p vertex_over_mid pins down alpha.
   */
  void testTransformedDual (const ElemType elem_type, const Real vertex_over_mid)
  {
    Mesh mesh(*TestCommWorld);
    MeshTools::Generation::build_square(mesh, 1, 1, -1., 1., -1., 1., elem_type);

    auto rng = mesh.active_local_element_ptr_range();
    if (rng.begin() == rng.end())
      return;
    const Elem * elem = *(rng.begin());

    FEType fe_type(SECOND, LAGRANGE);
    std::unique_ptr<FEBase> fe = FEBase::build(2, fe_type);
    const auto & JxW = fe->get_JxW();
    const auto & phi = fe->get_phi();
    const auto & dual_phi = fe->get_dual_phi();
    QGauss qrule(2, fe_type.default_quadrature_order());
    fe->attach_quadrature_rule(&qrule);
    fe->reinit(elem);

    const unsigned int n = elem->n_nodes();
    CPPUNIT_ASSERT_EQUAL(std::size_t(n), phi.size());
    CPPUNIT_ASSERT_EQUAL(std::size_t(n), dual_phi.size());

    // TOLERANCE*TOLERANCE works with double but not float128
    const Real my_tol = TOLERANCE*std::sqrt(TOLERANCE);

    // dtilde_j = \int dual_phi_j; these are the weights that must be positive.
    std::vector<Real> dtilde(n, 0);
    for (const auto j : make_range(n))
      for (const auto qp : index_range(JxW))
        dtilde[j] += JxW[qp]*dual_phi[j][qp];

    Real sum_dtilde = 0, volume = 0;
    for (const auto j : make_range(n))
      {
        CPPUNIT_ASSERT(dtilde[j] > 0);
        sum_dtilde += dtilde[j];
      }
    for (const auto qp : index_range(JxW))
      volume += JxW[qp];

    // Partition of unity is preserved, so constants are reproduced and the weights sum to the volume.
    LIBMESH_ASSERT_FP_EQUAL(volume, sum_dtilde, my_tol);
    for (const auto qp : index_range(JxW))
      {
        Real sum = 0;
        for (const auto j : make_range(n))
          sum += dual_phi[j][qp];
        LIBMESH_ASSERT_FP_EQUAL(1., sum, my_tol);
      }

    // Scale-free signature of alpha.
    LIBMESH_ASSERT_FP_EQUAL(vertex_over_mid,
                            dtilde[0]/dtilde[elem->n_vertices()], my_tol);

    // T from the element's own adjacency, not an assumed index convention.
    const Real alpha = Real(1)/5;
    DenseMatrix<Real> T(n, n);
    for (const auto i : make_range(n))
      T(i,i) = 1.;
    for (const auto m : make_range(elem->n_vertices(), n))
      {
        T(m,m) = 1. - 2.*alpha;
        for (const auto v : make_range(elem->n_second_order_adjacent_vertices(m)))
          T(elem->second_order_adjacent_vertex(m,v), m) += alpha;
      }

    for (const auto k : make_range(n))
      for (const auto j : make_range(n))
        {
          Real integral = 0;
          for (const auto qp : index_range(JxW))
            {
              Real ntilde_k = 0;
              for (const auto l : make_range(n))
                ntilde_k += T(k,l)*phi[l][qp];
              integral += JxW[qp]*dual_phi[j][qp]*ntilde_k;
            }
          LIBMESH_ASSERT_FP_EQUAL(k == j ? dtilde[j] : 0., integral, my_tol);
        }
  }

  void testQuad8Lagrange ()
  {
    LOG_UNIT_TEST;
    // QUAD8 weights 1/5 (corner) and 4/5 (mid-edge) on the reference element
    testTransformedDual(QUAD8, Real(1)/4);
  }

  void testTri6Lagrange ()
  {
    LOG_UNIT_TEST;
    // TRI6 weights 1/15 (vertex) and 1/10 (mid-edge) on the reference element
    testTransformedDual(TRI6, Real(2)/3);
  }

  void setUp()
    {
      FEType fe_type(FIRST, LAGRANGE);
      _fe = FEBase::build(1, fe_type);
      _fe->get_phi();
      _fe->get_dual_phi();

      _mesh = std::make_unique<Mesh>(*TestCommWorld);

      MeshTools::Generation::build_line(*_mesh, 1, -1, 1, EDGE2);

      auto rng = _mesh->active_local_element_ptr_range();
      _elem = rng.begin() == rng.end() ? nullptr : *(rng.begin());

      _qrule = std::make_unique<QGauss>(1, fe_type.default_quadrature_order());
      _fe->attach_quadrature_rule(_qrule.get());
    }

  void tearDown() {}

};

CPPUNIT_TEST_SUITE_REGISTRATION( DualShapeTest );
