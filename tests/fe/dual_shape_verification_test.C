// libmesh includes
#include "libmesh/libmesh.h"
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
  CPPUNIT_TEST_SUITE_END();

private:
  Mesh * _mesh;
  QGauss * _qrule;
  FEBase * _fe;
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

    LIBMESH_ASSERT_FP_EQUAL(2, dual_coeff(0,0), TOLERANCE*TOLERANCE);
    LIBMESH_ASSERT_FP_EQUAL(-1, dual_coeff(0,1), TOLERANCE*TOLERANCE);
    LIBMESH_ASSERT_FP_EQUAL(-1, dual_coeff(1,0), TOLERANCE*TOLERANCE);
    LIBMESH_ASSERT_FP_EQUAL(2, dual_coeff(1,1), TOLERANCE*TOLERANCE);

    const auto & dual_phi = _fe->get_dual_phi();

    CPPUNIT_ASSERT_EQUAL(std::size_t(2), dual_phi.size());

    const auto & qpoints = _qrule->get_points();

    CPPUNIT_ASSERT_EQUAL(qpoints.size(), dual_phi[0].size());

    for (auto qp : index_range(dual_phi[0]))
      LIBMESH_ASSERT_FP_EQUAL(1./2. * (1. - 3.*qpoints[qp](0)), dual_phi[0][qp],
        TOLERANCE*TOLERANCE);

    CPPUNIT_ASSERT_EQUAL(qpoints.size(), dual_phi[1].size());

    for (auto qp : index_range(dual_phi[1]))
      LIBMESH_ASSERT_FP_EQUAL(1./2. * (1. + 3.*qpoints[qp](0)), dual_phi[1][qp],
        TOLERANCE*TOLERANCE);
  }

  void setUp()
    {
      FEType fe_type(FIRST, LAGRANGE);
      _fe = FEBase::build(1, fe_type).release();
      _fe->get_phi();
      _fe->get_dual_phi();

      _mesh = new Mesh(*TestCommWorld);

      MeshTools::Generation::build_line(*_mesh, 1, -1, 1, EDGE2);

      auto rng = _mesh->active_local_element_ptr_range();
      _elem = rng.begin() == rng.end() ? nullptr : *(rng.begin());

      _qrule = new QGauss(1, fe_type.default_quadrature_order());
      _fe->attach_quadrature_rule(_qrule);
    }

  void tearDown()
    {
      delete _mesh;
      delete _fe;
      delete _qrule;
    }

};

CPPUNIT_TEST_SUITE_REGISTRATION( DualShapeTest );
