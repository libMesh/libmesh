#include <libmesh/boundary_info.h>
#include <libmesh/elem.h>
#include <libmesh/fe_map.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_modification.h>
#include <libmesh/mesh_tools.h>

#include "test_comm.h"
#include "libmesh_cppunit.h"


using namespace libMesh;

class AllRBBTest : public CppUnit::TestCase
{
  /**
   * The goal of this test is to verify proper operation of the
   * all_rbb() mesh modification, by converting meshes of classic
   * "isogeometric" types from Lagrange interpolants to Rational
   * Bezier-Bernstein splines, then verifying that the converted
   * geometry is exact to within floating-point error (for point
   * radii) or quadrature error (for element volume).
   */
public:
  LIBMESH_CPPUNIT_TEST_SUITE( AllRBBTest );

  CPPUNIT_TEST( testAllRBBNodeElem );
  CPPUNIT_TEST( testAllRBBEdge );
  CPPUNIT_TEST( testAllRBBEdge3 );

  // 2D tests
#if LIBMESH_DIM > 1
  CPPUNIT_TEST( testAllRBBTri );
  CPPUNIT_TEST( testAllRBBTri6 );
  CPPUNIT_TEST( testAllRBBQuad );
  CPPUNIT_TEST( testAllRBBQuad8 );
  CPPUNIT_TEST( testAllRBBQuad9 );

  CPPUNIT_TEST( testAllRBBCircle4 );
  CPPUNIT_TEST( testAllRBBCircle8 );
  CPPUNIT_TEST( testAllRBBCircle16 );

  CPPUNIT_TEST( testAllRBBDisk5 );
  CPPUNIT_TEST( testAllRBBDisk20 );
  CPPUNIT_TEST( testAllRBBDisk80 );

  CPPUNIT_TEST( testAllRBBTri6Disk10 );
  CPPUNIT_TEST( testAllRBBTri6Disk40 );
  CPPUNIT_TEST( testAllRBBTri6Disk160 );
#endif

  // 3D tests
#if LIBMESH_DIM > 2
  CPPUNIT_TEST( testAllRBBTet );
  CPPUNIT_TEST( testAllRBBTet10 );
  CPPUNIT_TEST( testAllRBBHex );
  CPPUNIT_TEST( testAllRBBHex20 );
  CPPUNIT_TEST( testAllRBBHex27 );
  CPPUNIT_TEST( testAllRBBPrism6 );
  CPPUNIT_TEST( testAllRBBPrism18 );
#endif

  CPPUNIT_TEST_SUITE_END();

protected:
  // We don't do anything interesting to affine elements in all_rbb(),
  // but we can verify that we're not screwing them up.
  void test_box(ElemType elem_type)
  {
    Mesh mesh(*TestCommWorld);

    auto dim = Elem::type_to_dim_map[elem_type];

    MeshTools::Generation::build_cube(mesh,
                                        dim > 0 ? 2 : 0, dim > 1 ? 1 : 0, dim > 2 ? 1 : 0,
                                        0., 1.,
                                        0., 1.,
                                        0., 1.,
                                        elem_type);

    const auto n_orig_elem = mesh.n_elem();

    MeshTools::Modification::all_rbb(mesh);

    CPPUNIT_ASSERT_EQUAL(n_orig_elem, mesh.n_elem());

    unsigned char weight_index = mesh.default_mapping_data();

    for (auto & elem : mesh.element_ptr_range())
    {
      CPPUNIT_ASSERT_EQUAL(elem->mapping_type(), RATIONAL_BERNSTEIN_MAP);
      CPPUNIT_ASSERT(elem->has_affine_map());
    }

    for (auto & node : mesh.node_ptr_range())
    {
      const Real w = node->get_extra_datum<Real>(weight_index);

      CPPUNIT_ASSERT_EQUAL(Real(1), w);
    }
  }

  void test_circle(unsigned int n_refinements)
  {
    Mesh interior_mesh(*TestCommWorld),
         boundary_mesh(*TestCommWorld);

    const Real radius = 1;
    const Real circumference = 2 * pi * radius;
    const Real tol = TOLERANCE*TOLERANCE;

    // Build a filled circle
    MeshTools::Generation::build_sphere (interior_mesh, radius,
                                         n_refinements, QUAD9);

    // Get just the outer EDGE3 circle mesh
    interior_mesh.get_boundary_info().sync(boundary_mesh);

    const dof_id_type n_edges = 4 << n_refinements;

    CPPUNIT_ASSERT_EQUAL(boundary_mesh.n_elem(), n_edges);
    CPPUNIT_ASSERT_EQUAL(boundary_mesh.n_elem(), n_edges);

    for (auto & node : boundary_mesh.node_ptr_range())
      {
        const Point p = *node;
        LIBMESH_ASSERT_FP_EQUAL(p.norm(), radius, tol);
      }

    // We just did Lagrange interpolation, so our mesh measure
    // shouldn't be *quite* right.  Empirically, we converge from
    // beneath, and our error looks like Ch^4.
    const Real max_lagrange_error =
      radius * 5e-2 / (1 << (4*n_refinements));
    LIBMESH_ASSERT_FP_EQUAL(MeshTools::volume(boundary_mesh),
                            circumference, max_lagrange_error);

    MeshTools::Modification::all_rbb(boundary_mesh);

    for (auto & elem : boundary_mesh.element_ptr_range())
      {
        CPPUNIT_ASSERT_EQUAL(RATIONAL_BERNSTEIN_MAP, elem->mapping_type());

        // We can no longer assert that each Node is at a specified
        // radius from the circle center, because these are now spline
        // control nodes, but we can assert that physical points
        // within the element are at the desired radius.
        constexpr int n_intervals = 4;
        Point master_pt;
        for (master_pt(0) = -1; master_pt(0) <= 1 + TOLERANCE;
             master_pt(0) += Real(2)/n_intervals)
          {
            const Point p = FEMap::map(elem->dim(), elem, master_pt);
            LIBMESH_ASSERT_FP_EQUAL(radius, p.norm(), tol);
          }
      }

    // We're using quadrature for volume approximation, so we still
    // have error, but our quadrature error looks something like Ch^6
    // with a much smaller C.
    const Real max_rbb_error =
      radius * 1e-3 / (1 << (6*n_refinements));
    LIBMESH_ASSERT_FP_EQUAL(MeshTools::volume(boundary_mesh),
                            circumference, max_rbb_error);
  }

  void test_disk(unsigned int n_refinements, const ElemType type = QUAD9)
  {
    Mesh mesh(*TestCommWorld);

    const Real radius = 1;
    const Real area = pi * radius * radius;
    const Real tol = TOLERANCE*TOLERANCE;

    // Build a filled circle
    MeshTools::Generation::build_sphere (mesh, radius,
                                         n_refinements, type);

    const dof_id_type n_elem =
      (5 << (n_refinements*2)) * (type == QUAD9 ? 1 : 2);

    CPPUNIT_ASSERT_EQUAL(mesh.n_elem(), n_elem);

    // We just did Lagrange interpolation, so our mesh measure
    // shouldn't be *quite* right.  Empirically, we converge from
    // beneath, and our error looks like Ch^4.
    const Real max_lagrange_error =
      radius * 5e-2 / (1 << (4*n_refinements));

    LIBMESH_ASSERT_FP_EQUAL(MeshTools::volume(mesh),
                            area, max_lagrange_error);

    MeshTools::Modification::all_rbb(mesh);

    for (const Elem * elem : mesh.element_ptr_range())
      {
        CPPUNIT_ASSERT_EQUAL(RATIONAL_BERNSTEIN_MAP, elem->mapping_type());

        // We can no longer assert that each Node is at a specified
        // radius from the circle center, because these are now spline
        // control nodes, but we can assert that physical points
        // within the element are at the desired radius.
        for (auto s : make_range(elem->n_sides()))
          {
            if (elem->neighbor_ptr(s))
              continue;

            constexpr int n_intervals = 4;
            Point master_pt = elem->master_point(s);
            const Point step =
              (elem->master_point((s+1)%elem->n_sides()) - master_pt)
              / n_intervals;
            for (auto i : make_range(n_intervals+1))
              {
                libmesh_ignore(i);
                const Point p = FEMap::map(elem->dim(), elem, master_pt);
                LIBMESH_ASSERT_FP_EQUAL(radius, p.norm(), tol);
                master_pt += step;
              }
          }
      }

    // We're using quadrature for volume approximation, so we still
    // have error, but our quadrature error looks like Ch^6 with a
    // much smaller C.
    const Real max_rbb_error =
      radius * 2e-3 / (1 << (6*n_refinements));
    LIBMESH_ASSERT_FP_EQUAL(MeshTools::volume(mesh),
                            area, max_rbb_error);
  }

public:
  void setUp() {}

  void tearDown() {}

  void testAllRBBNodeElem() { LOG_UNIT_TEST; test_box(NODEELEM); }
  void testAllRBBEdge() { LOG_UNIT_TEST; test_box(EDGE2); }
  void testAllRBBEdge3() { LOG_UNIT_TEST; test_box(EDGE3); }
  void testAllRBBTri() { LOG_UNIT_TEST; test_box(TRI3); }
  void testAllRBBTri6() { LOG_UNIT_TEST; test_box(TRI6); }
  void testAllRBBQuad() { LOG_UNIT_TEST; test_box(QUAD4); }
  void testAllRBBQuad8() { LOG_UNIT_TEST; test_box(QUAD8); }
  void testAllRBBQuad9() { LOG_UNIT_TEST; test_box(QUAD9); }
  void testAllRBBTet() { LOG_UNIT_TEST; test_box(TET4); }
  void testAllRBBTet10() { LOG_UNIT_TEST; test_box(TET10); }
  void testAllRBBHex() { LOG_UNIT_TEST; test_box(HEX8); }
  void testAllRBBHex20() { LOG_UNIT_TEST; test_box(HEX20); }
  void testAllRBBHex27() { LOG_UNIT_TEST; test_box(HEX27); }
  void testAllRBBPrism6() { LOG_UNIT_TEST; test_box(PRISM6); }
  void testAllRBBPrism18() { LOG_UNIT_TEST; test_box(PRISM18); }

  // We still don't support general Polys, Tri7s or anything with them
  // as faces, infinite elements, or anything above quadratic.

  // 0 refinements of our default circle gives us 4 RBB edges
  void testAllRBBCircle4() { LOG_UNIT_TEST; test_circle(0); }
  void testAllRBBCircle8() { LOG_UNIT_TEST; test_circle(1); }
  void testAllRBBCircle16() { LOG_UNIT_TEST; test_circle(2); }

  // 0 refinements of our default disk gives us 5 RBB quads or 10 RBB
  // triangles
  void testAllRBBDisk5() { LOG_UNIT_TEST; test_disk(0); }
  void testAllRBBDisk20() { LOG_UNIT_TEST; test_disk(1); }
  void testAllRBBDisk80() { LOG_UNIT_TEST; test_disk(2); }

  void testAllRBBTri6Disk10() { LOG_UNIT_TEST; test_disk(0, TRI6); }
  void testAllRBBTri6Disk40() { LOG_UNIT_TEST; test_disk(1, TRI6); }
  void testAllRBBTri6Disk160() { LOG_UNIT_TEST; test_disk(2, TRI6); }
};


CPPUNIT_TEST_SUITE_REGISTRATION( AllRBBTest );
