// Tests for libMesh::Kokkos FE type conversion helpers (CPU-side, no device
// compiler required).  Mirrors the MOOSE KokkosFETypesTest unit tests.
//
// Run:
//   make check   (or)   ./unit_tests-opt --test KokkosFETypesTest

#include "libmesh_cppunit.h"

#ifdef LIBMESH_HAVE_KOKKOS

#include "kokkos/fe_types.h"

#include "libmesh/enum_elem_type.h"
#include "libmesh/enum_fe_family.h"

using namespace libMesh::Kokkos;

class KokkosFETypesTest : public CppUnit::TestCase
{
  LIBMESH_CPPUNIT_TEST_SUITE(KokkosFETypesTest);

  // toKokkosTopology
  CPPUNIT_TEST(testToTopology_Edge);
  CPPUNIT_TEST(testToTopology_2D);
  CPPUNIT_TEST(testToTopology_3D);

  // toKokkosFamily
  CPPUNIT_TEST(testToFamily_Standard);
  CPPUNIT_TEST(testToFamily_Monomial);
  CPPUNIT_TEST(testToFamily_UnknownSentinel);

  // getSideTopology
  CPPUNIT_TEST(testSideTopology_1D);
  CPPUNIT_TEST(testSideTopology_2D);
  CPPUNIT_TEST(testSideTopology_3D);

  // nDofs(FEFamily, FEElemTopology)
  CPPUNIT_TEST(testNDofs_Lagrange_1D);
  CPPUNIT_TEST(testNDofs_Lagrange_2D);
  CPPUNIT_TEST(testNDofs_Lagrange_3D);

  // classFromTopology
  CPPUNIT_TEST(testClassFromTopology_1D);
  CPPUNIT_TEST(testClassFromTopology_2D);
  CPPUNIT_TEST(testClassFromTopology_3D);

  // nDofs(FEShapeKey) — LAGRANGE
  CPPUNIT_TEST(testNDofs_Key_Lagrange_1D);
  CPPUNIT_TEST(testNDofs_Key_Lagrange_2D);
  CPPUNIT_TEST(testNDofs_Key_Lagrange_3D);

  // nDofs(FEShapeKey) — MONOMIAL
  CPPUNIT_TEST(testNDofs_Key_Monomial_1D);
  CPPUNIT_TEST(testNDofs_Key_Monomial_2D);
  CPPUNIT_TEST(testNDofs_Key_Monomial_3D);

  // Enum contiguity
  CPPUNIT_TEST(testEnumContiguity);

  CPPUNIT_TEST_SUITE_END();

public:

  // ── toKokkosTopology ───────────────────────────────────────────────────────

  void testToTopology_Edge()
  {
    LOG_UNIT_TEST;
    CPPUNIT_ASSERT_EQUAL(toKokkosTopology(libMesh::EDGE2), FEElemTopology::EDGE2);
    CPPUNIT_ASSERT_EQUAL(toKokkosTopology(libMesh::EDGE3), FEElemTopology::EDGE3);
  }

  void testToTopology_2D()
  {
    LOG_UNIT_TEST;
    CPPUNIT_ASSERT_EQUAL(toKokkosTopology(libMesh::TRI3),  FEElemTopology::TRI3);
    CPPUNIT_ASSERT_EQUAL(toKokkosTopology(libMesh::TRI6),  FEElemTopology::TRI6);
    CPPUNIT_ASSERT_EQUAL(toKokkosTopology(libMesh::QUAD4), FEElemTopology::QUAD4);
    CPPUNIT_ASSERT_EQUAL(toKokkosTopology(libMesh::QUAD8), FEElemTopology::QUAD8);
    CPPUNIT_ASSERT_EQUAL(toKokkosTopology(libMesh::QUAD9), FEElemTopology::QUAD9);
  }

  void testToTopology_3D()
  {
    LOG_UNIT_TEST;
    CPPUNIT_ASSERT_EQUAL(toKokkosTopology(libMesh::TET4),  FEElemTopology::TET4);
    CPPUNIT_ASSERT_EQUAL(toKokkosTopology(libMesh::TET10), FEElemTopology::TET10);
    CPPUNIT_ASSERT_EQUAL(toKokkosTopology(libMesh::HEX8),  FEElemTopology::HEX8);
    CPPUNIT_ASSERT_EQUAL(toKokkosTopology(libMesh::HEX20), FEElemTopology::HEX20);
    CPPUNIT_ASSERT_EQUAL(toKokkosTopology(libMesh::HEX27), FEElemTopology::HEX27);
  }

  // ── toKokkosFamily ─────────────────────────────────────────────────────────

  void testToFamily_Standard()
  {
    LOG_UNIT_TEST;
    CPPUNIT_ASSERT_EQUAL(toKokkosFamily(libMesh::LAGRANGE),     FEFamily::LAGRANGE);
    CPPUNIT_ASSERT_EQUAL(toKokkosFamily(libMesh::LAGRANGE_VEC), FEFamily::LAGRANGE_VEC);
    CPPUNIT_ASSERT_EQUAL(toKokkosFamily(libMesh::HERMITE),      FEFamily::HERMITE);
  }

  void testToFamily_Monomial()
  {
    LOG_UNIT_TEST;
    CPPUNIT_ASSERT_EQUAL(toKokkosFamily(libMesh::MONOMIAL),     FEFamily::MONOMIAL);
    CPPUNIT_ASSERT_EQUAL(toKokkosFamily(libMesh::MONOMIAL_VEC), FEFamily::MONOMIAL_VEC);
  }

  void testToFamily_UnknownSentinel()
  {
    LOG_UNIT_TEST;
    // Unrecognised families must return UNKNOWN rather than abort.
    CPPUNIT_ASSERT_EQUAL(toKokkosFamily(libMesh::HIERARCHIC), FEFamily::UNKNOWN);
    CPPUNIT_ASSERT(toKokkosFamily(libMesh::HIERARCHIC) != FEFamily::LAGRANGE);
    CPPUNIT_ASSERT(toKokkosFamily(libMesh::HIERARCHIC) != FEFamily::MONOMIAL);
  }

  // ── getSideTopology ────────────────────────────────────────────────────────

  void testSideTopology_1D()
  {
    LOG_UNIT_TEST;
    CPPUNIT_ASSERT_EQUAL(getSideTopology(FEElemTopology::EDGE2), FEElemTopology::EDGE2);
    CPPUNIT_ASSERT_EQUAL(getSideTopology(FEElemTopology::EDGE3), FEElemTopology::EDGE2);
  }

  void testSideTopology_2D()
  {
    LOG_UNIT_TEST;
    // First-order 2D → linear edge sides
    CPPUNIT_ASSERT_EQUAL(getSideTopology(FEElemTopology::TRI3),  FEElemTopology::EDGE2);
    CPPUNIT_ASSERT_EQUAL(getSideTopology(FEElemTopology::QUAD4), FEElemTopology::EDGE2);
    // Second-order 2D → quadratic edge sides
    CPPUNIT_ASSERT_EQUAL(getSideTopology(FEElemTopology::TRI6),  FEElemTopology::EDGE3);
    CPPUNIT_ASSERT_EQUAL(getSideTopology(FEElemTopology::QUAD8), FEElemTopology::EDGE3);
    CPPUNIT_ASSERT_EQUAL(getSideTopology(FEElemTopology::QUAD9), FEElemTopology::EDGE3);
  }

  void testSideTopology_3D()
  {
    LOG_UNIT_TEST;
    CPPUNIT_ASSERT_EQUAL(getSideTopology(FEElemTopology::TET4),  FEElemTopology::TRI3);
    CPPUNIT_ASSERT_EQUAL(getSideTopology(FEElemTopology::HEX8),  FEElemTopology::QUAD4);
    CPPUNIT_ASSERT_EQUAL(getSideTopology(FEElemTopology::TET10), FEElemTopology::TRI6);
    CPPUNIT_ASSERT_EQUAL(getSideTopology(FEElemTopology::HEX20), FEElemTopology::QUAD8);
    CPPUNIT_ASSERT_EQUAL(getSideTopology(FEElemTopology::HEX27), FEElemTopology::QUAD9);
  }

  // ── nDofs(FEFamily, FEElemTopology) ───────────────────────────────────────

  void testNDofs_Lagrange_1D()
  {
    LOG_UNIT_TEST;
    CPPUNIT_ASSERT_EQUAL(nDofs(FEFamily::LAGRANGE, FEElemTopology::EDGE2), 2u);
    CPPUNIT_ASSERT_EQUAL(nDofs(FEFamily::LAGRANGE, FEElemTopology::EDGE3), 3u);
  }

  void testNDofs_Lagrange_2D()
  {
    LOG_UNIT_TEST;
    CPPUNIT_ASSERT_EQUAL(nDofs(FEFamily::LAGRANGE, FEElemTopology::TRI3),  3u);
    CPPUNIT_ASSERT_EQUAL(nDofs(FEFamily::LAGRANGE, FEElemTopology::TRI6),  6u);
    CPPUNIT_ASSERT_EQUAL(nDofs(FEFamily::LAGRANGE, FEElemTopology::QUAD4), 4u);
    CPPUNIT_ASSERT_EQUAL(nDofs(FEFamily::LAGRANGE, FEElemTopology::QUAD8), 8u);
    CPPUNIT_ASSERT_EQUAL(nDofs(FEFamily::LAGRANGE, FEElemTopology::QUAD9), 9u);
  }

  void testNDofs_Lagrange_3D()
  {
    LOG_UNIT_TEST;
    CPPUNIT_ASSERT_EQUAL(nDofs(FEFamily::LAGRANGE, FEElemTopology::TET4),  4u);
    CPPUNIT_ASSERT_EQUAL(nDofs(FEFamily::LAGRANGE, FEElemTopology::TET10), 10u);
    CPPUNIT_ASSERT_EQUAL(nDofs(FEFamily::LAGRANGE, FEElemTopology::HEX8),  8u);
    CPPUNIT_ASSERT_EQUAL(nDofs(FEFamily::LAGRANGE, FEElemTopology::HEX20), 20u);
    CPPUNIT_ASSERT_EQUAL(nDofs(FEFamily::LAGRANGE, FEElemTopology::HEX27), 27u);
  }

  // ── classFromTopology ──────────────────────────────────────────────────────

  void testClassFromTopology_1D()
  {
    LOG_UNIT_TEST;
    CPPUNIT_ASSERT_EQUAL(classFromTopology(FEElemTopology::EDGE2), FEElemClass::EDGE);
    CPPUNIT_ASSERT_EQUAL(classFromTopology(FEElemTopology::EDGE3), FEElemClass::EDGE);
  }

  void testClassFromTopology_2D()
  {
    LOG_UNIT_TEST;
    CPPUNIT_ASSERT_EQUAL(classFromTopology(FEElemTopology::TRI3),  FEElemClass::TRI);
    CPPUNIT_ASSERT_EQUAL(classFromTopology(FEElemTopology::TRI6),  FEElemClass::TRI);
    CPPUNIT_ASSERT_EQUAL(classFromTopology(FEElemTopology::QUAD4), FEElemClass::QUAD);
    CPPUNIT_ASSERT_EQUAL(classFromTopology(FEElemTopology::QUAD8), FEElemClass::QUAD);
    CPPUNIT_ASSERT_EQUAL(classFromTopology(FEElemTopology::QUAD9), FEElemClass::QUAD);
  }

  void testClassFromTopology_3D()
  {
    LOG_UNIT_TEST;
    CPPUNIT_ASSERT_EQUAL(classFromTopology(FEElemTopology::TET4),  FEElemClass::TET);
    CPPUNIT_ASSERT_EQUAL(classFromTopology(FEElemTopology::TET10), FEElemClass::TET);
    CPPUNIT_ASSERT_EQUAL(classFromTopology(FEElemTopology::HEX8),  FEElemClass::HEX);
    CPPUNIT_ASSERT_EQUAL(classFromTopology(FEElemTopology::HEX20), FEElemClass::HEX);
    CPPUNIT_ASSERT_EQUAL(classFromTopology(FEElemTopology::HEX27), FEElemClass::HEX);
  }

  // ── nDofs(FEShapeKey) — LAGRANGE ──────────────────────────────────────────

  void testNDofs_Key_Lagrange_1D()
  {
    LOG_UNIT_TEST;
    CPPUNIT_ASSERT_EQUAL(nDofs(FEShapeKey{FEFamily::LAGRANGE, FEElemClass::EDGE, 1}), 2u);
    CPPUNIT_ASSERT_EQUAL(nDofs(FEShapeKey{FEFamily::LAGRANGE, FEElemClass::EDGE, 2}), 3u);
  }

  void testNDofs_Key_Lagrange_2D()
  {
    LOG_UNIT_TEST;
    CPPUNIT_ASSERT_EQUAL(nDofs(FEShapeKey{FEFamily::LAGRANGE, FEElemClass::TRI,  1}),  3u);
    CPPUNIT_ASSERT_EQUAL(nDofs(FEShapeKey{FEFamily::LAGRANGE, FEElemClass::TRI,  2}),  6u);
    CPPUNIT_ASSERT_EQUAL(nDofs(FEShapeKey{FEFamily::LAGRANGE, FEElemClass::QUAD, 1}),  4u);
    CPPUNIT_ASSERT_EQUAL(nDofs(FEShapeKey{FEFamily::LAGRANGE, FEElemClass::QUAD, 2}),  9u);
  }

  void testNDofs_Key_Lagrange_3D()
  {
    LOG_UNIT_TEST;
    CPPUNIT_ASSERT_EQUAL(nDofs(FEShapeKey{FEFamily::LAGRANGE, FEElemClass::TET, 1}),  4u);
    CPPUNIT_ASSERT_EQUAL(nDofs(FEShapeKey{FEFamily::LAGRANGE, FEElemClass::TET, 2}), 10u);
    CPPUNIT_ASSERT_EQUAL(nDofs(FEShapeKey{FEFamily::LAGRANGE, FEElemClass::HEX, 1}),  8u);
    CPPUNIT_ASSERT_EQUAL(nDofs(FEShapeKey{FEFamily::LAGRANGE, FEElemClass::HEX, 2}), 27u);
  }

  // ── nDofs(FEShapeKey) — MONOMIAL ──────────────────────────────────────────
  // 1D: p+1  |  2D: (p+1)(p+2)/2  |  3D: (p+1)(p+2)(p+3)/6

  void testNDofs_Key_Monomial_1D()
  {
    LOG_UNIT_TEST;
    for (unsigned int p = 0; p <= 5; ++p)
      CPPUNIT_ASSERT_EQUAL(nDofs(FEShapeKey{FEFamily::MONOMIAL, FEElemClass::EDGE, p}), p + 1);
  }

  void testNDofs_Key_Monomial_2D()
  {
    LOG_UNIT_TEST;
    for (unsigned int p = 0; p <= 5; ++p)
    {
      const unsigned int expected = (p + 1) * (p + 2) / 2;
      CPPUNIT_ASSERT_EQUAL(nDofs(FEShapeKey{FEFamily::MONOMIAL, FEElemClass::TRI,  p}), expected);
      CPPUNIT_ASSERT_EQUAL(nDofs(FEShapeKey{FEFamily::MONOMIAL, FEElemClass::QUAD, p}), expected);
    }
  }

  void testNDofs_Key_Monomial_3D()
  {
    LOG_UNIT_TEST;
    for (unsigned int p = 0; p <= 5; ++p)
    {
      const unsigned int expected = (p + 1) * (p + 2) * (p + 3) / 6;
      CPPUNIT_ASSERT_EQUAL(nDofs(FEShapeKey{FEFamily::MONOMIAL, FEElemClass::TET, p}), expected);
      CPPUNIT_ASSERT_EQUAL(nDofs(FEShapeKey{FEFamily::MONOMIAL, FEElemClass::HEX, p}), expected);
    }
  }

  // ── Enum contiguity (values usable as array indices) ──────────────────────

  void testEnumContiguity()
  {
    LOG_UNIT_TEST;
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(FEElemTopology::EDGE2),   0u);
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(FEElemTopology::EDGE3),   1u);
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(FEElemTopology::TRI3),    2u);
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(FEElemTopology::TRI6),    3u);
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(FEElemTopology::QUAD4),   4u);
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(FEElemTopology::QUAD8),   5u);
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(FEElemTopology::QUAD9),   6u);
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(FEElemTopology::TET4),    7u);
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(FEElemTopology::TET10),   8u);
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(FEElemTopology::HEX8),    9u);
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(FEElemTopology::HEX20),  10u);
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(FEElemTopology::HEX27),  11u);
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(FEElemTopology::N_TYPES), 12u);
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION(KokkosFETypesTest);

#endif // LIBMESH_HAVE_KOKKOS
