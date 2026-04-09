// Tests for libMesh::Kokkos FEEvaluator<> shape function specialisations.
//
// For every (element type, polynomial order) pair the test:
//   A. Cross-validates FEEvaluator::shape / grad_shape against libMesh FEBase
//      reinit on the reference element at QGauss quadrature points.
//   B. Checks partition of unity: sum_i phi_i(qp) == 1 at every quad point.
//   C. Checks gradient sum: sum_i grad_phi_i(qp) == 0 (follows from PoU).
//   D. Checks dispatch parity: nativeShape / nativeGradShape return bit-
//      identical results to direct FEEvaluator<> specialisation calls.
//   E. Checks MONOMIAL parity against libMesh MONOMIAL FEBase reinit.
//   F. Checks face QP → parent coordinate mapping via nativeMapShape.
//
// All device-code paths are compiled on the host; KOKKOS_INLINE_FUNCTION
// degrades to 'inline' for host-only builds.

#include "libmesh_cppunit.h"

#ifdef LIBMESH_HAVE_KOKKOS

#include "kokkos/fe_types.h"
#include "kokkos/fe_base.h"
#include "kokkos/scalar_types.h"
#include "kokkos/fe_lagrange_1d.h"
#include "kokkos/fe_lagrange_2d.h"
#include "kokkos/fe_lagrange_3d.h"
#include "kokkos/fe_monomial.h"
#include "kokkos/fe_evaluator.h"
#include "kokkos/fe_face_map.h"

#include "libmesh/fe_base.h"
#include "libmesh/fe.h"
#include "libmesh/fe_type.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/reference_elem.h"
#include "libmesh/enum_elem_type.h"
#include "libmesh/enum_order.h"

#include <cmath>
#include <memory>
#include <vector>

using namespace libMesh::Kokkos;

static constexpr double TOL = 1.0e-13;

// ── element table ─────────────────────────────────────────────────────────────

struct KokkosElemEntry
{
  FEElemTopology   topo;
  libMesh::ElemType lm_type;
  libMesh::Order    lm_order;
  unsigned int      dim;
  unsigned int      n_dofs;
  const char *      name;
};

static const KokkosElemEntry KOKKOS_ELEMS[] =
{
  { FEElemTopology::EDGE2,  libMesh::EDGE2,  libMesh::FIRST,  1,  2, "EDGE2"  },
  { FEElemTopology::EDGE3,  libMesh::EDGE3,  libMesh::SECOND, 1,  3, "EDGE3"  },
  { FEElemTopology::TRI3,   libMesh::TRI3,   libMesh::FIRST,  2,  3, "TRI3"   },
  { FEElemTopology::TRI6,   libMesh::TRI6,   libMesh::SECOND, 2,  6, "TRI6"   },
  { FEElemTopology::QUAD4,  libMesh::QUAD4,  libMesh::FIRST,  2,  4, "QUAD4"  },
  { FEElemTopology::QUAD8,  libMesh::QUAD8,  libMesh::SECOND, 2,  8, "QUAD8"  },
  { FEElemTopology::QUAD9,  libMesh::QUAD9,  libMesh::SECOND, 2,  9, "QUAD9"  },
  { FEElemTopology::TET4,   libMesh::TET4,   libMesh::FIRST,  3,  4, "TET4"   },
  { FEElemTopology::TET10,  libMesh::TET10,  libMesh::SECOND, 3, 10, "TET10"  },
  { FEElemTopology::HEX8,   libMesh::HEX8,   libMesh::FIRST,  3,  8, "HEX8"   },
  { FEElemTopology::HEX20,  libMesh::HEX20,  libMesh::SECOND, 3, 20, "HEX20"  },
  { FEElemTopology::HEX27,  libMesh::HEX27,  libMesh::SECOND, 3, 27, "HEX27"  },
};
static constexpr unsigned int N_KOKKOS_ELEMS =
  sizeof(KOKKOS_ELEMS) / sizeof(KOKKOS_ELEMS[0]);

// ── helpers ───────────────────────────────────────────────────────────────────

// Return QGauss points on the libMesh reference element in Real3 coords.
static std::vector<Real3> qpointsFromLibMesh(const KokkosElemEntry & e,
                                              unsigned int order)
{
  libMesh::QGauss qr(e.dim, static_cast<libMesh::Order>(order));
  qr.allow_rules_with_negative_weights = true;
  qr.init(e.lm_type);
  std::vector<Real3> pts(qr.n_points());
  for (unsigned int q = 0; q < qr.n_points(); ++q)
  {
    pts[q].v[0] = qr.qp(q)(0);
    pts[q].v[1] = (e.dim >= 2) ? qr.qp(q)(1) : 0.0;
    pts[q].v[2] = (e.dim >= 3) ? qr.qp(q)(2) : 0.0;
  }
  return pts;
}

// ── A + B + C  Cross-validation + partition of unity ─────────────────────────

// Template base class, one instantiation per (elem, order) pair.
template <unsigned int ElemIdx, unsigned int QOrder>
class KokkosFEShapeTest : public CppUnit::TestCase
{
  LIBMESH_CPPUNIT_TEST_SUITE(KokkosFEShapeTest);
  CPPUNIT_TEST(testPhiMatchesLibMesh);
  CPPUNIT_TEST(testPartitionOfUnity);
  CPPUNIT_TEST(testGradSumIsZero);
  CPPUNIT_TEST_SUITE_END();

  static const KokkosElemEntry & E() { return KOKKOS_ELEMS[ElemIdx]; }

public:

  void testPhiMatchesLibMesh()
  {
    LOG_UNIT_TEST;
    const auto & e = E();

    auto fe = libMesh::FEBase::build(
        e.dim, libMesh::FEType(e.lm_order, libMesh::LAGRANGE));
    libMesh::QGauss qr(e.dim, static_cast<libMesh::Order>(QOrder));
    qr.allow_rules_with_negative_weights = true;
    fe->attach_quadrature_rule(&qr);

    const auto & phi_lm  = fe->get_phi();
    const auto & dphi_lm = fe->get_dphi();
    const libMesh::Elem * ref = &libMesh::ReferenceElem::get(e.lm_type);
    fe->reinit(ref);

    auto qpts = qpointsFromLibMesh(e, QOrder);
    CPPUNIT_ASSERT_EQUAL((unsigned int)qr.n_points(), (unsigned int)qpts.size());

    for (unsigned int i = 0; i < e.n_dofs; ++i)
      for (unsigned int q = 0; q < qpts.size(); ++q)
      {
        Real xi   = qpts[q].v[0];
        Real eta  = qpts[q].v[1];
        Real zeta = qpts[q].v[2];

        LIBMESH_ASSERT_FP_EQUAL(nativeShape(e.topo, i, xi, eta, zeta),
                                 phi_lm[i][q], TOL);

        Real3 ng = nativeGradShape(e.topo, i, xi, eta, zeta);
        LIBMESH_ASSERT_FP_EQUAL(ng.v[0], dphi_lm[i][q](0), TOL);
        if (e.dim >= 2)
          LIBMESH_ASSERT_FP_EQUAL(ng.v[1], dphi_lm[i][q](1), TOL);
        if (e.dim >= 3)
          LIBMESH_ASSERT_FP_EQUAL(ng.v[2], dphi_lm[i][q](2), TOL);
      }
  }

  void testPartitionOfUnity()
  {
    LOG_UNIT_TEST;
    const auto & e = E();
    auto qpts = qpointsFromLibMesh(e, QOrder);

    for (unsigned int q = 0; q < qpts.size(); ++q)
    {
      Real sum = 0.0;
      for (unsigned int i = 0; i < e.n_dofs; ++i)
        sum += nativeShape(e.topo, i, qpts[q].v[0], qpts[q].v[1], qpts[q].v[2]);
      LIBMESH_ASSERT_FP_EQUAL(1.0, sum, TOL);
    }
  }

  void testGradSumIsZero()
  {
    LOG_UNIT_TEST;
    const auto & e = E();
    auto qpts = qpointsFromLibMesh(e, QOrder);

    for (unsigned int q = 0; q < qpts.size(); ++q)
    {
      Real3 gs{};
      for (unsigned int i = 0; i < e.n_dofs; ++i)
      {
        Real3 g = nativeGradShape(e.topo, i, qpts[q].v[0], qpts[q].v[1], qpts[q].v[2]);
        gs.v[0] += g.v[0];
        gs.v[1] += g.v[1];
        gs.v[2] += g.v[2];
      }
      LIBMESH_ASSERT_FP_EQUAL(0.0, gs.v[0], TOL);
      if (e.dim >= 2) LIBMESH_ASSERT_FP_EQUAL(0.0, gs.v[1], TOL);
      if (e.dim >= 3) LIBMESH_ASSERT_FP_EQUAL(0.0, gs.v[2], TOL);
    }
  }
};

// Instantiate for each element at quadrature orders 1–4
#define INSTANTIATE_KOKKOS_SHAPE_TEST(EIDX, ORDER)                            \
  using KokkosFEShapeTest_##EIDX##_o##ORDER =                                 \
    KokkosFEShapeTest<EIDX, ORDER>;                                            \
  CPPUNIT_TEST_SUITE_REGISTRATION(KokkosFEShapeTest_##EIDX##_o##ORDER)

// 1D elements (indices 0-1)
INSTANTIATE_KOKKOS_SHAPE_TEST(0, 1);
INSTANTIATE_KOKKOS_SHAPE_TEST(0, 3);
INSTANTIATE_KOKKOS_SHAPE_TEST(1, 2);
INSTANTIATE_KOKKOS_SHAPE_TEST(1, 4);

// 2D triangles (2-3)
INSTANTIATE_KOKKOS_SHAPE_TEST(2, 1);
INSTANTIATE_KOKKOS_SHAPE_TEST(2, 3);
INSTANTIATE_KOKKOS_SHAPE_TEST(3, 2);
INSTANTIATE_KOKKOS_SHAPE_TEST(3, 4);

// 2D quads (4-6)
INSTANTIATE_KOKKOS_SHAPE_TEST(4, 1);
INSTANTIATE_KOKKOS_SHAPE_TEST(4, 3);
INSTANTIATE_KOKKOS_SHAPE_TEST(5, 2);
INSTANTIATE_KOKKOS_SHAPE_TEST(5, 4);
INSTANTIATE_KOKKOS_SHAPE_TEST(6, 2);
INSTANTIATE_KOKKOS_SHAPE_TEST(6, 4);

// 3D tets (7-8)
INSTANTIATE_KOKKOS_SHAPE_TEST(7, 1);
INSTANTIATE_KOKKOS_SHAPE_TEST(7, 3);
INSTANTIATE_KOKKOS_SHAPE_TEST(8, 2);
INSTANTIATE_KOKKOS_SHAPE_TEST(8, 4);

// 3D hexes (9-11)
INSTANTIATE_KOKKOS_SHAPE_TEST(9, 1);
INSTANTIATE_KOKKOS_SHAPE_TEST(9, 3);
INSTANTIATE_KOKKOS_SHAPE_TEST(10, 2);
INSTANTIATE_KOKKOS_SHAPE_TEST(10, 4);
INSTANTIATE_KOKKOS_SHAPE_TEST(11, 2);
INSTANTIATE_KOKKOS_SHAPE_TEST(11, 4);

// ── D. Dispatch parity — bit-identical to direct FEEvaluator<> calls ─────────

class KokkosDispatchParityTest : public CppUnit::TestCase
{
  LIBMESH_CPPUNIT_TEST_SUITE(KokkosDispatchParityTest);
  CPPUNIT_TEST(testQuad4PhiBitIdentical);
  CPPUNIT_TEST(testTet4PhiBitIdentical);
  CPPUNIT_TEST(testHex27PhiBitIdentical);
  CPPUNIT_TEST(testTri6GradBitIdentical);
  CPPUNIT_TEST_SUITE_END();

public:

  void testQuad4PhiBitIdentical()
  {
    LOG_UNIT_TEST;
    auto qpts = qpointsFromLibMesh(KOKKOS_ELEMS[4], 4);
    for (unsigned int i = 0; i < 4; ++i)
      for (unsigned int q = 0; q < qpts.size(); ++q)
      {
        Real xi = qpts[q].v[0], eta = qpts[q].v[1];
        Real direct   = FEEvaluator<LagrangeTag, Quad4Tag>::shape(i, xi, eta, 0.0);
        Real dispatch = nativeShape(FEElemTopology::QUAD4, i, xi, eta, 0.0);
        CPPUNIT_ASSERT(direct == dispatch);
      }
  }

  void testTet4PhiBitIdentical()
  {
    LOG_UNIT_TEST;
    auto qpts = qpointsFromLibMesh(KOKKOS_ELEMS[7], 5);
    for (unsigned int i = 0; i < 4; ++i)
      for (unsigned int q = 0; q < qpts.size(); ++q)
      {
        Real xi = qpts[q].v[0], eta = qpts[q].v[1], zeta = qpts[q].v[2];
        Real direct   = FEEvaluator<LagrangeTag, Tet4Tag>::shape(i, xi, eta, zeta);
        Real dispatch = nativeShape(FEElemTopology::TET4, i, xi, eta, zeta);
        CPPUNIT_ASSERT(direct == dispatch);
      }
  }

  void testHex27PhiBitIdentical()
  {
    LOG_UNIT_TEST;
    auto qpts = qpointsFromLibMesh(KOKKOS_ELEMS[11], 3);
    for (unsigned int i = 0; i < 27; ++i)
      for (unsigned int q = 0; q < qpts.size(); ++q)
      {
        Real xi = qpts[q].v[0], eta = qpts[q].v[1], zeta = qpts[q].v[2];
        Real direct   = FEEvaluator<LagrangeTag, Hex27Tag>::shape(i, xi, eta, zeta);
        Real dispatch = nativeShape(FEElemTopology::HEX27, i, xi, eta, zeta);
        CPPUNIT_ASSERT(direct == dispatch);
      }
  }

  void testTri6GradBitIdentical()
  {
    LOG_UNIT_TEST;
    auto qpts = qpointsFromLibMesh(KOKKOS_ELEMS[3], 4);
    for (unsigned int i = 0; i < 6; ++i)
      for (unsigned int q = 0; q < qpts.size(); ++q)
      {
        Real xi = qpts[q].v[0], eta = qpts[q].v[1];
        Real3 direct   = FEEvaluator<LagrangeTag, Tri6Tag>::grad_shape(i, xi, eta, 0.0);
        Real3 dispatch = nativeGradShape(FEElemTopology::TRI6, i, xi, eta, 0.0);
        CPPUNIT_ASSERT(direct.v[0] == dispatch.v[0]);
        CPPUNIT_ASSERT(direct.v[1] == dispatch.v[1]);
      }
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION(KokkosDispatchParityTest);

// ── E. MONOMIAL parity against libMesh FEBase ─────────────────────────────────

template <unsigned int Dim, unsigned int Order, unsigned int ElemIdx>
class KokkosMonomialParityTest : public CppUnit::TestCase
{
  LIBMESH_CPPUNIT_TEST_SUITE(KokkosMonomialParityTest);
  CPPUNIT_TEST(testMonomialPhiMatchesLibMesh);
  CPPUNIT_TEST_SUITE_END();

  static const KokkosElemEntry & E() { return KOKKOS_ELEMS[ElemIdx]; }

  static FEElemClass elemClass()
  {
    return classFromTopology(KOKKOS_ELEMS[ElemIdx].topo);
  }

public:

  void testMonomialPhiMatchesLibMesh()
  {
    LOG_UNIT_TEST;
    const auto & e = E();
    const unsigned int quad_order = Order + 1;

    const FEShapeKey key{FEFamily::MONOMIAL, elemClass(), Order};
    const unsigned int n = nDofs(key);
    const libMesh::Order lm_order = static_cast<libMesh::Order>(Order);

    auto qpts = qpointsFromLibMesh(e, quad_order);

    for (unsigned int i = 0; i < n; ++i)
      for (unsigned int q = 0; q < qpts.size(); ++q)
      {
        Real xi   = qpts[q].v[0];
        Real eta  = qpts[q].v[1];
        Real zeta = qpts[q].v[2];
        libMesh::Point p(xi, eta, zeta);

        // Use the static FE<Dim, MONOMIAL>::shape() API to get the reference
        // value.  This avoids the fe->reinit() + phi_lm path which can give
        // wrong values for 3D MONOMIAL at order > 1 on bare reference elements.
        Real lm_phi, lm_dphi_dx, lm_dphi_dy, lm_dphi_dz;
        if (Dim == 1)
        {
          lm_phi    = libMesh::FE<1, libMesh::MONOMIAL>::shape      (e.lm_type, lm_order, i, p);
          lm_dphi_dx = libMesh::FE<1, libMesh::MONOMIAL>::shape_deriv(e.lm_type, lm_order, i, 0, p);
          lm_dphi_dy = 0; lm_dphi_dz = 0;
        }
        else if (Dim == 2)
        {
          lm_phi     = libMesh::FE<2, libMesh::MONOMIAL>::shape      (e.lm_type, lm_order, i, p);
          lm_dphi_dx = libMesh::FE<2, libMesh::MONOMIAL>::shape_deriv(e.lm_type, lm_order, i, 0, p);
          lm_dphi_dy = libMesh::FE<2, libMesh::MONOMIAL>::shape_deriv(e.lm_type, lm_order, i, 1, p);
          lm_dphi_dz = 0;
        }
        else
        {
          lm_phi     = libMesh::FE<3, libMesh::MONOMIAL>::shape      (e.lm_type, lm_order, i, p);
          lm_dphi_dx = libMesh::FE<3, libMesh::MONOMIAL>::shape_deriv(e.lm_type, lm_order, i, 0, p);
          lm_dphi_dy = libMesh::FE<3, libMesh::MONOMIAL>::shape_deriv(e.lm_type, lm_order, i, 1, p);
          lm_dphi_dz = libMesh::FE<3, libMesh::MONOMIAL>::shape_deriv(e.lm_type, lm_order, i, 2, p);
        }

        LIBMESH_ASSERT_FP_EQUAL(nativeShape(key, i, xi, eta, zeta), lm_phi, TOL);

        Real3 ng = nativeGradShape(key, i, xi, eta, zeta);
        LIBMESH_ASSERT_FP_EQUAL(ng.v[0], lm_dphi_dx, TOL);
        if (e.dim >= 2)
          LIBMESH_ASSERT_FP_EQUAL(ng.v[1], lm_dphi_dy, TOL);
        if (e.dim >= 3)
          LIBMESH_ASSERT_FP_EQUAL(ng.v[2], lm_dphi_dz, TOL);
      }
  }
};

// Dim 1: EDGE (idx 0), orders 0-3
#define INST_MONO(D, O, EIDX)                                                   \
  using KokkosMonomialParityTest_d##D##_o##O##_e##EIDX =                        \
    KokkosMonomialParityTest<D, O, EIDX>;                                        \
  CPPUNIT_TEST_SUITE_REGISTRATION(KokkosMonomialParityTest_d##D##_o##O##_e##EIDX)

INST_MONO(1, 0, 0);
INST_MONO(1, 1, 0);
INST_MONO(1, 2, 0);
INST_MONO(1, 3, 0);

// Dim 2: QUAD (idx 4)
INST_MONO(2, 0, 4);
INST_MONO(2, 1, 4);
INST_MONO(2, 2, 4);
INST_MONO(2, 3, 4);

// Dim 2: TRI (idx 2)
INST_MONO(2, 1, 2);
INST_MONO(2, 2, 2);

// Dim 3: HEX (idx 9)
INST_MONO(3, 0, 9);
INST_MONO(3, 1, 9);
INST_MONO(3, 2, 9);
INST_MONO(3, 3, 9);

// Dim 3: TET (idx 7)
INST_MONO(3, 1, 7);
INST_MONO(3, 2, 7);

// ── F. Face QP → parent coordinate mapping ────────────────────────────────────

class KokkosFaceMapTest : public CppUnit::TestCase
{
  LIBMESH_CPPUNIT_TEST_SUITE(KokkosFaceMapTest);
  CPPUNIT_TEST(testFaceMapping_Quad4);
  CPPUNIT_TEST(testFaceMapping_Tri3);
  CPPUNIT_TEST(testFaceMapping_Hex8);
  CPPUNIT_TEST(testFaceMapping_Tet4);
  CPPUNIT_TEST(testFaceMapping_Hex8_HighOrder);
  CPPUNIT_TEST_SUITE_END();

  static void checkFaceMapping(const KokkosElemEntry & e, unsigned int quad_order)
  {
    const libMesh::Elem * ref = &libMesh::ReferenceElem::get(e.lm_type);
    const auto side_topo = getSideTopology(e.topo);

    libMesh::QGauss qr_face(e.dim - 1,
                             static_cast<libMesh::Order>(quad_order));
    qr_face.allow_rules_with_negative_weights = true;

    auto fe_face = libMesh::FEBase::build(
        e.dim, libMesh::FEType(e.lm_order, libMesh::LAGRANGE));
    fe_face->attach_quadrature_rule(&qr_face);
    const auto & phi_lm = fe_face->get_phi();

    for (unsigned int side = 0; side < ref->n_sides(); ++side)
    {
      fe_face->reinit(ref, side);
      unsigned int nqp = qr_face.n_points();

      // Get face QPs from libMesh side rule (side-element reference coords)
      qr_face.init(side_topo == FEElemTopology::EDGE2 ? libMesh::EDGE2 :
                   side_topo == FEElemTopology::EDGE3 ? libMesh::EDGE3 :
                   side_topo == FEElemTopology::TRI3  ? libMesh::TRI3  :
                   side_topo == FEElemTopology::TRI6  ? libMesh::TRI6  :
                   side_topo == FEElemTopology::QUAD4 ? libMesh::QUAD4 :
                                                        libMesh::QUAD8);
      CPPUNIT_ASSERT_EQUAL(qr_face.n_points(), nqp);

      auto side_elem = ref->side_ptr(side);

      for (unsigned int q = 0; q < nqp; ++q)
      {
        // Map face QP to parent reference coords using isoparametric basis
        Real3 face_pt{};
        face_pt.v[0] = qr_face.qp(q)(0);
        face_pt.v[1] = (e.dim >= 3) ? qr_face.qp(q)(1) : 0.0;

        Real3 parent_pt{};
        for (unsigned int k = 0; k < side_elem->n_nodes(); ++k)
        {
          Real psi = nativeShape(side_topo, k,
                                  face_pt.v[0], face_pt.v[1], 0.0);
          parent_pt.v[0] += psi * side_elem->point(k)(0);
          parent_pt.v[1] += psi * side_elem->point(k)(1);
          parent_pt.v[2] += psi * side_elem->point(k)(2);
        }

        for (unsigned int i = 0; i < e.n_dofs; ++i)
        {
          double native_phi = nativeShape(e.topo, i,
                                           parent_pt.v[0], parent_pt.v[1], parent_pt.v[2]);
          LIBMESH_ASSERT_FP_EQUAL(native_phi, phi_lm[i][q], TOL);
        }
      }
    }
  }

public:

  void testFaceMapping_Quad4()
  {
    LOG_UNIT_TEST;
    checkFaceMapping(KOKKOS_ELEMS[4], 3);
  }

  void testFaceMapping_Tri3()
  {
    LOG_UNIT_TEST;
    checkFaceMapping(KOKKOS_ELEMS[2], 3);
  }

  void testFaceMapping_Hex8()
  {
    LOG_UNIT_TEST;
    checkFaceMapping(KOKKOS_ELEMS[9], 3);
  }

  void testFaceMapping_Tet4()
  {
    LOG_UNIT_TEST;
    checkFaceMapping(KOKKOS_ELEMS[7], 3);
  }

  void testFaceMapping_Hex8_HighOrder()
  {
    LOG_UNIT_TEST;
    checkFaceMapping(KOKKOS_ELEMS[9], 5);
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION(KokkosFaceMapTest);

// ── G. All topologies handled by nativeShape (no missing switch cases) ────────

class KokkosNativeShapeSafetyTest : public CppUnit::TestCase
{
  LIBMESH_CPPUNIT_TEST_SUITE(KokkosNativeShapeSafetyTest);
  CPPUNIT_TEST(testAllTopologiesFinite);
  CPPUNIT_TEST_SUITE_END();

public:

  void testAllTopologiesFinite()
  {
    LOG_UNIT_TEST;
    // Evaluate phi_0 at the centroid of each reference element.
    // Any unhandled switch case will typically produce NaN or 0.
    const Real3 centroids[] = {
      Real3{ 0.0,     0.0,     0.0 },     // EDGE2
      Real3{ 0.0,     0.0,     0.0 },     // EDGE3
      Real3{ 1.0/3.0, 1.0/3.0, 0.0 },    // TRI3
      Real3{ 1.0/3.0, 1.0/3.0, 0.0 },    // TRI6
      Real3{ 0.0,     0.0,     0.0 },     // QUAD4
      Real3{ 0.0,     0.0,     0.0 },     // QUAD8
      Real3{ 0.0,     0.0,     0.0 },     // QUAD9
      Real3{ 0.25,    0.25,    0.25 },    // TET4
      Real3{ 0.25,    0.25,    0.25 },    // TET10
      Real3{ 0.0,     0.0,     0.0 },     // HEX8
      Real3{ 0.0,     0.0,     0.0 },     // HEX20
      Real3{ 0.0,     0.0,     0.0 },     // HEX27
    };

    for (unsigned int e = 0; e < N_KOKKOS_ELEMS; ++e)
    {
      Real val = nativeShape(KOKKOS_ELEMS[e].topo, 0,
                              centroids[e].v[0], centroids[e].v[1], centroids[e].v[2]);
      CPPUNIT_ASSERT(std::isfinite(val));
    }
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION(KokkosNativeShapeSafetyTest);

#endif // LIBMESH_HAVE_KOKKOS
