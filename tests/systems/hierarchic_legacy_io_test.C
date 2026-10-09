#include <libmesh/dof_map.h>
#include <libmesh/elem.h>
#include <libmesh/enum_xdr_mode.h>
#include <libmesh/equation_systems.h>
#include <libmesh/explicit_system.h>
#include <libmesh/fe.h>
#include <libmesh/fe_interface.h>
#include <libmesh/libmesh_version.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/numeric_vector.h>
#include <libmesh/replicated_mesh.h>

#include "test_comm.h"
#include "libmesh_cppunit.h"

#include <cmath>
#include <fstream>
#include <sstream>

using namespace libMesh;

namespace
{

// A smooth field that excites every mode of a projection onto it
Number field(const Point & p, const Parameters &, const std::string &, const std::string & var)
{
  const Real w = (var == "u") ? 1 : (var == "v") ? 2 : 3;
  return std::sin(w*p(0) + 0.5) * std::cos(0.7*p(1) - 0.2*w) * std::exp(0.3*p(2)) +
    0.1*w*p(0)*p(1)*p(1);
}

// The one-dimensional hierarchic function of order i as libMesh defined it before its bubbles
// were rescaled: the vertex functions, then (xi^i - 1)/i! for even i and (xi^i - xi)/i! for odd
Real legacy_1D_shape(const unsigned int i, const Real xi)
{
  if (i == 0)
    return (1 - xi)/2;
  if (i == 1)
    return (1 + xi)/2;

  Real factorial = 1;
  for (const auto n : make_range(2u, i + 1))
    factorial *= n;

  return (std::pow(xi, i) - ((i % 2) ? xi : 1)) / factorial;
}

}


class HierarchicLegacyIOTest : public CppUnit::TestCase {
public:
  LIBMESH_CPPUNIT_TEST_SUITE( HierarchicLegacyIOTest );

#if LIBMESH_DIM > 1
  CPPUNIT_TEST( testEdgeRatio );
  CPPUNIT_TEST( testQuadRatio );
  CPPUNIT_TEST( testReadLegacyQuad9 );
  CPPUNIT_TEST( testReadLegacyTri6 );
  CPPUNIT_TEST( testReadLegacyQuad9L2 );
  CPPUNIT_TEST( testReadLegacyQuad9Side );
#endif
#if LIBMESH_DIM > 2
  CPPUNIT_TEST( testReadLegacyHex27 );
  CPPUNIT_TEST( testReadLegacyTet14 );
  CPPUNIT_TEST( testReadLegacyPrism21 );
  CPPUNIT_TEST( testReadLegacyHex27Side );
#endif
  CPPUNIT_TEST( testRejectNewerFormat );
#if LIBMESH_DIM > 1
  CPPUNIT_TEST( testCurrentRoundTrip );
#endif

  CPPUNIT_TEST_SUITE_END();

public:
  void setUp() {}

  void tearDown() {}

  // On an edge each shape function is a single one-dimensional function, so the ratio is the old
  // shape function divided by the new one wherever the new one is nonzero
  void testEdgeRatio()
  {
    LOG_UNIT_TEST;

    ReplicatedMesh mesh(*TestCommWorld);
    MeshTools::Generation::build_line(mesh, 1, -1, 1, EDGE3);
    const Elem & elem = mesh.elem_ref(0);

    const unsigned int order = 9;
    const FEType fe_type(order, HIERARCHIC);
    const Point p(0.37);

    for (const auto i : make_range(order + 1))
      {
        const Real new_shape = FEInterface::shape(fe_type, &elem, i, p);
        const Real old_shape = legacy_1D_shape(i, p(0));
        LIBMESH_ASSERT_FP_EQUAL
          (old_shape / new_shape,
           fe_hierarchic_legacy_coefficient_ratio(HIERARCHIC, elem, Order(order), i),
           TOLERANCE*TOLERANCE);
      }
  }

  // Likewise for the tensor-product quadrilateral, whose shape functions are products of two
  // one-dimensional functions up to an orientation sign that both bases share
  void testQuadRatio()
  {
    LOG_UNIT_TEST;

    ReplicatedMesh mesh(*TestCommWorld);
    MeshTools::Generation::build_square(mesh, 1, 1, -1, 1, -1, 1, QUAD9);
    const Elem & elem = mesh.elem_ref(0);

    const unsigned int order = 6;
    const FEType fe_type(order, HIERARCHIC);
    const Point p(0.37, -0.61);

    for (const auto i : make_range((order + 1)*(order + 1)))
      {
        const auto [i0, i1, f] = fe_hierarchic_quad_tensor_indices(&elem, order, i);
        const Real new_shape = FEInterface::shape(fe_type, &elem, i, p);
        const Real old_shape = f * legacy_1D_shape(i0, p(0)) * legacy_1D_shape(i1, p(1));
        LIBMESH_ASSERT_FP_EQUAL
          (old_shape / new_shape,
           fe_hierarchic_legacy_coefficient_ratio(HIERARCHIC, elem, Order(order), i),
           TOLERANCE*TOLERANCE);
      }
  }

  // A file written with I/O compatibility version 1.7.0, before the bubbles were rescaled, must
  // read back as the same finite element function. Alongside each file is a table of the
  // function's values at sample points, which the library that wrote it evaluated in its own basis.
  void checkLegacyFile(const std::string & name)
  {
    const std::string base = "solutions/legacy_hierarchic/" + name;

    ReplicatedMesh mesh(*TestCommWorld);
    mesh.read(base + "_mesh.xda");

    EquationSystems es(mesh);
    es.read(base + ".xda", READ,
            EquationSystems::READ_HEADER | EquationSystems::READ_DATA);

    const System & sys = es.get_system("sys");

    std::ifstream values(base + "_values.txt");
    CPPUNIT_ASSERT(values.good());

    // Each processor checks the rows on the elements whose coefficients it holds, which between
    // them cover every row
    dof_id_type elem_id;
    Point p;
    unsigned int n_rows = 0, n_checked = 0;
    while (values >> elem_id >> p(0) >> p(1) >> p(2))
      {
        const Elem & elem = mesh.elem_ref(elem_id);
        const bool evaluable = sys.get_dof_map().is_evaluable(elem);
        for (const auto var : make_range(sys.n_vars()))
          {
            Real expected;
            values >> expected;

            if (evaluable)
              LIBMESH_ASSERT_FP_EQUAL(expected, libmesh_real(sys.point_value(var, p, elem)),
                                      TOLERANCE*TOLERANCE);
          }
        ++n_rows;
        n_checked += (elem.processor_id() == mesh.processor_id());
      }

    // Every element in the mesh has rows, and every row was checked on its element's owner
    CPPUNIT_ASSERT(n_rows >= mesh.n_elem());
    mesh.comm().sum(n_checked);
    CPPUNIT_ASSERT_EQUAL(n_rows, n_checked);
  }

  void testReadLegacyQuad9()      { LOG_UNIT_TEST; checkLegacyFile("quad9_hier4"); }
  void testReadLegacyTri6()       { LOG_UNIT_TEST; checkLegacyFile("tri6_hier5"); }
  void testReadLegacyQuad9L2()    { LOG_UNIT_TEST; checkLegacyFile("quad9_l2hier5"); }
  void testReadLegacyQuad9Side()  { LOG_UNIT_TEST; checkLegacyFile("quad9_side4"); }
  void testReadLegacyHex27()      { LOG_UNIT_TEST; checkLegacyFile("hex27_hier4"); }
  void testReadLegacyTet14()      { LOG_UNIT_TEST; checkLegacyFile("tet14_hier3"); }
  void testReadLegacyPrism21()    { LOG_UNIT_TEST; checkLegacyFile("prism21_hier4"); }
  void testReadLegacyHex27Side()  { LOG_UNIT_TEST; checkLegacyFile("hex27_side3"); }

  // A file written in the current format holds current-basis coefficients and must read back
  // without conversion
  void testCurrentRoundTrip()
  {
    LOG_UNIT_TEST;

    ReplicatedMesh mesh(*TestCommWorld);
    MeshTools::Generation::build_square(mesh, 2, 2, -1, 1, -1, 1, QUAD9);

    EquationSystems es(mesh);
    ExplicitSystem & sys = es.add_system<ExplicitSystem>("sys");
    sys.add_variable("u", FIFTH, HIERARCHIC);
    es.init();
    sys.project_solution(field, nullptr, es.parameters);
    es.write("hierarchic_round_trip.xda", WRITE, EquationSystems::WRITE_DATA);

    EquationSystems es2(mesh);
    es2.read("hierarchic_round_trip.xda", READ,
             EquationSystems::READ_HEADER | EquationSystems::READ_DATA);
    const System & sys2 = es2.get_system("sys");

    std::unique_ptr<NumericVector<Number>> diff = sys2.solution->clone();
    *diff -= *sys.solution;
    LIBMESH_ASSERT_FP_EQUAL(0, diff->linfty_norm(), TOLERANCE*TOLERANCE);
  }

  // A file claiming a format version newer than this library's must be refused
  void testRejectNewerFormat()
  {
    LOG_UNIT_TEST;

    std::istringstream iss(get_io_compatibility_version());
    int major = 0, minor = 0, patch = 0;
    char dot;
    iss >> major >> dot >> minor >> dot >> patch;

    // This library's own version is accepted, with or without the qualifiers that follow it
    const std::string current = "libMesh-" + get_io_compatibility_version();
    CPPUNIT_ASSERT_EQUAL(LIBMESH_VERSION_ID(major, minor, patch),
                         parse_io_compatibility_version(current + " parallel"));

#ifdef LIBMESH_ENABLE_EXCEPTIONS
    const std::string newer =
      "libMesh-" + std::to_string(major) + "." + std::to_string(minor + 1) + ".0";
    CPPUNIT_ASSERT_THROW_MESSAGE
      ("A newer file format was not refused",
       parse_io_compatibility_version(newer),
       libMesh::LogicError);
#endif
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION( HierarchicLegacyIOTest );
