#include <libmesh/boundary_info.h>
#include <libmesh/elem.h>
#include <libmesh/enum_quadrature_type.h>
#include <libmesh/fe_map.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_modification.h>
#include <libmesh/mesh_tools.h>
#include <libmesh/quadrature.h>

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

  // We use AMR when generating circles
#  ifdef LIBMESH_ENABLE_AMR
  CPPUNIT_TEST( testAllRBBCircle4 );
  CPPUNIT_TEST( testAllRBBCircle8 );
  CPPUNIT_TEST( testAllRBBCircle16 );

  CPPUNIT_TEST( testAllRBBDisk5 );
  CPPUNIT_TEST( testAllRBBDisk20 );
  CPPUNIT_TEST( testAllRBBDisk80 );

  CPPUNIT_TEST( testAllRBBTri6Disk10 );
  CPPUNIT_TEST( testAllRBBTri6Disk40 );
  CPPUNIT_TEST( testAllRBBTri6Disk160 );
#  endif
#endif

  // 3D tests
#if LIBMESH_DIM > 2
  CPPUNIT_TEST( testAllRBBTet );
  CPPUNIT_TEST( testAllRBBTet10 );
  CPPUNIT_TEST( testAllRBBHex );
  CPPUNIT_TEST( testAllRBBHex20 );
  CPPUNIT_TEST( testAllRBBHex27 );

// We need to add BERNSTEIN support for Prisms!
//  CPPUNIT_TEST( testAllRBBPrism6 );
//  CPPUNIT_TEST( testAllRBBPrism18 );

  // We use AMR when generating circles or spheres
#  ifdef LIBMESH_ENABLE_AMR
  CPPUNIT_TEST( testAllRBBCylinder10 );
  CPPUNIT_TEST( testAllRBBCylinder80 );

  CPPUNIT_TEST( testAllRBBSphere6 );
  CPPUNIT_TEST( testAllRBBSphere24 );
  CPPUNIT_TEST( testAllRBBSphere96 );

  CPPUNIT_TEST( testAllRBBTri6Sphere12 );
  CPPUNIT_TEST( testAllRBBTri6Sphere48 );
  CPPUNIT_TEST( testAllRBBTri6Sphere192 );
#  endif

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
    CPPUNIT_ASSERT_EQUAL(n_orig_elem, mesh.max_elem_id());

    const Real orig_volume = MeshTools::volume(mesh) / n_orig_elem;

    std::vector<Real> orig_hmin(n_orig_elem), orig_hmax(n_orig_elem);

    // In the tet case our elements all have the same volume but
    // they're stretched differently and have different hmin/hmax
    for (auto & elem : mesh.element_ptr_range())
    {
      orig_hmin[elem->id()] = elem->hmin();
      orig_hmax[elem->id()] = elem->hmax();
    }

    MeshTools::Modification::all_rbb(mesh);

    CPPUNIT_ASSERT_EQUAL(n_orig_elem, mesh.n_elem());
    CPPUNIT_ASSERT_EQUAL(n_orig_elem, mesh.max_elem_id());

    unsigned char weight_index = mesh.default_mapping_data();

    for (auto & elem : mesh.element_ptr_range())
    {
      CPPUNIT_ASSERT_EQUAL(elem->mapping_type(),
                           RATIONAL_BERNSTEIN_MAP);
      CPPUNIT_ASSERT(elem->has_affine_map());

      // Tri6 has this much FP error??
      LIBMESH_ASSERT_FP_EQUAL(elem->volume(), orig_volume,
                              20*TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(elem->hmax(), orig_hmax[elem->id()],
                              TOLERANCE*TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(elem->hmin(), orig_hmin[elem->id()],
                              TOLERANCE*TOLERANCE);
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
    CPPUNIT_ASSERT_EQUAL(boundary_mesh.n_nodes(), n_edges*2);

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

  void test_cylinder (unsigned int n_refinements, const ElemType type = HEX27)
  {
    Mesh disk_mesh(*TestCommWorld), mesh(*TestCommWorld);

    const Real radius = 1;
    const Real height = 3;
    const Real volume = pi * radius * radius * height;
    const Real tol = TOLERANCE*TOLERANCE;

    // We're extruding a circle from side 0 up
    const ElemType side_type = Elem::build(type)->side_type(0);

    // Build a filled circle
    MeshTools::Generation::build_sphere (disk_mesh, radius,
                                         n_refinements, side_type);

    // Then extrude it into a cylinder
    const unsigned int nz = 2;
    const RealVectorValue extrusion_vector{0,0,height};
    MeshTools::Generation::build_extrusion (mesh, disk_mesh, nz, extrusion_vector);

    const dof_id_type n_elem =
      (5 << (n_refinements*2)) * (side_type == QUAD9 ? 1 : 2) * nz;

    CPPUNIT_ASSERT_EQUAL(mesh.n_elem(), n_elem);

    // We just did Lagrange interpolation, so our mesh measure
    // shouldn't be *quite* right, but we should converge similarly to
    // how we did with the filled disk.
    const Real max_lagrange_error =
      radius * 5e-2 / (1 << (4*n_refinements)) *
      radius * radius * height;

    LIBMESH_ASSERT_FP_EQUAL(MeshTools::volume(mesh),
                            volume, max_lagrange_error);

    MeshTools::Modification::all_rbb(mesh);

    std::unique_ptr<const Elem> elem_side;
    constexpr int n_intervals = 4;
    auto qrule = QBase::build(QGRID, /*dim=*/2, Order(n_intervals));

    for (const Elem * elem : mesh.element_ptr_range())
      {
        CPPUNIT_ASSERT_EQUAL(RATIONAL_BERNSTEIN_MAP, elem->mapping_type());

        // We cannot assert that each Node is at a specified radius
        // from the cylinder axis, because these are now spline
        // control nodes, but we can assert that physical points
        // within the rounded surface are at the desired radius.
        for (auto s : make_range(elem->n_sides()))
          {
            if (elem->neighbor_ptr(s))
              continue;

            const Point side_normal =
              elem->side_vertex_average_normal(s);

            // We ought to either be on the rounded surface or the end
            // caps
            if (std::abs(side_normal(2)) > TOLERANCE*TOLERANCE)
              {
                LIBMESH_ASSERT_FP_EQUAL(side_normal(0), 0, TOLERANCE*TOLERANCE);
                LIBMESH_ASSERT_FP_EQUAL(side_normal(1), 0, TOLERANCE*TOLERANCE);
                continue;
              }

            elem->build_side_ptr(elem_side, s);
            qrule->init(*elem_side);

            for (auto i : make_range(qrule->n_points()))
              {
                Point p = FEMap::map(elem_side->dim(), elem_side.get(), qrule->qp(i));
                p(2) = 0; // Just look at r in cylindrical coordinates
                LIBMESH_ASSERT_FP_EQUAL(radius, p.norm(), tol);
              }
          }
      }

    // We're using quadrature for volume approximation, so we still
    // have error, but our quadrature error looks like Ch^6 with a
    // much smaller C.
    const Real max_rbb_error =
      radius * 2e-3 / (1 << (6*n_refinements)) *
      radius * radius * height;
    LIBMESH_ASSERT_FP_EQUAL(MeshTools::volume(mesh),
                            volume, max_rbb_error);
  }

  void test_sphere(unsigned int n_refinements, const ElemType type = QUAD9)
  {
    Mesh interior_mesh(*TestCommWorld),
         boundary_mesh(*TestCommWorld);

    const Real radius = 1;
    const Real surface_area = 4 * pi * radius * radius;
    const Real tol = TOLERANCE*TOLERANCE;

    // Build a filled sphere.  We're going to avoid using the `flat`
    // direct 2D option here because that only currently supports
    // TRI3.
    MeshTools::Generation::build_sphere (interior_mesh, radius,
                                         n_refinements, HEX27);

    // Get just the outer QUAD9 sphere mesh
    interior_mesh.get_boundary_info().sync(boundary_mesh);

    const dof_id_type n_faces = 6 << (2*n_refinements);

    CPPUNIT_ASSERT_EQUAL(boundary_mesh.n_elem(), n_faces);

    auto check_radii = [&boundary_mesh, radius, type](Real radius_tol) {
      constexpr int n_intervals = 40;

      Real max_radius_error = 0;
      for (auto & elem : boundary_mesh.element_ptr_range())
        {
          // We can't necessarily assert that each Node is at a
          // specified radius from the circle center, because these
          // may be spline control nodes.  We want physical points
          // within an element to all be at the desired radius, but on
          // rational quadratics that's only possible for "latitude /
          // longitude" quad edges, so we need a non-trivial tolerance
          // here.
          Point master_pt;
          if (type == TRI6)
            {
              for (master_pt(0) = 0; master_pt(0) <= 1 + TOLERANCE;
                   master_pt(0) += Real(1)/n_intervals)
                {
                  for (master_pt(1) = 0; master_pt(1) <= 1 - master_pt(0) + TOLERANCE;
                       master_pt(1) += Real(1)/n_intervals)
                    {
                      const Point p = FEMap::map(elem->dim(), elem, master_pt);
                      max_radius_error = std::max(max_radius_error, std::abs(radius-p.norm()));
                    }
                }
            }
          else
            {
              libmesh_assert_equal_to(type, QUAD9);
              for (master_pt(0) = -1; master_pt(0) <= 1 + TOLERANCE;
                   master_pt(0) += Real(2)/n_intervals)
                {
                  for (master_pt(1) = -1; master_pt(1) <= 1 + TOLERANCE;
                       master_pt(1) += Real(2)/n_intervals)
                    {
                      const Point p = FEMap::map(elem->dim(), elem, master_pt);
                      max_radius_error = std::max(max_radius_error, std::abs(radius-p.norm()));
                    }
                }
            }
          LIBMESH_ASSERT_FP_EQUAL(max_radius_error, 0, radius_tol);
        }
    };

    auto verify_equispaced_midnodes = [&boundary_mesh, type]()
    {
      std::unique_ptr<Elem> side_ptr;
      for (auto & elem : boundary_mesh.element_ptr_range())
        {
          // We expect midnodes to be and to remain equispaced on
          // quad edges.  But our QUAD9 elements can be curved
          // trapezoids, so we don't expect "diagonal" edges to have
          // an equispaced mid-face node.
          if (type != TRI6)
            for (auto s : elem->side_index_range())
              {
                elem->build_side_ptr(side_ptr, s);

                auto c02 = side_ptr->point(2) - side_ptr->point(0);
                auto c12 = side_ptr->point(2) - side_ptr->point(1);
                LIBMESH_ASSERT_FP_EQUAL
                  (c02.norm_sq(), c12.norm_sq(), TOLERANCE*TOLERANCE);
              }
        }
    };

    verify_equispaced_midnodes();

    if (type == TRI6)
      {
        MeshTools::Modification::all_tri(boundary_mesh);
        verify_equispaced_midnodes();
      }
    else
      libmesh_assert_equal_to(type, QUAD9);

    // For our sphere construction, all of our Lagrange nodes should
    // be at exactly the right radius.
    for (auto & node : boundary_mesh.node_ptr_range())
      {
        const Point p = *node;
        LIBMESH_ASSERT_FP_EQUAL(p.norm(), radius, tol);
      }

    // But Lagrange isn't isogeometric, so non-nodes should have
    // radius error.  Empirically, our error looks like Ch^3
    const Real max_lagrange_rad_error =
      radius * 0.2 / (1 << (3*n_refinements));
    check_radii(max_lagrange_rad_error);

    // We just did Lagrange interpolation, so our mesh measure
    // shouldn't be *quite* right.  Empirically, we converge from
    // beneath, and our error looks like Ch^4.
    const Real max_lagrange_vol_error =
      radius * radius * 1.5 / (1 << (4*n_refinements));
    LIBMESH_ASSERT_FP_EQUAL(MeshTools::volume(boundary_mesh),
                            surface_area, max_lagrange_vol_error);

    MeshTools::Modification::all_rbb(boundary_mesh);

    // All of our vertices should still be at exactly the right
    // radius; for these the control point is the point.
    for (auto & elem : boundary_mesh.element_ptr_range())
      {
        CPPUNIT_ASSERT_EQUAL(RATIONAL_BERNSTEIN_MAP, elem->mapping_type());

        for (auto v : make_range(elem->n_vertices()))
          {
            const Point & p = elem->point(v);
            LIBMESH_ASSERT_FP_EQUAL(p.norm(), radius, tol);
          }
      }

    // We're not building a stereographic mesh on the sphere, so we
    // still have error ... of only slightly better magnitude?
    const Real max_rbb_rad_error =
      radius * 0.2 / (1 << (3*n_refinements));
    check_radii(max_rbb_rad_error);

    // And we also get about the same order on volume():
    const Real max_rbb_vol_error =
      radius * radius * 1.5 / (1 << (4*n_refinements));
    LIBMESH_ASSERT_FP_EQUAL(MeshTools::volume(boundary_mesh),
                            surface_area, max_rbb_vol_error);
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

  // We need to add BERNSTEIN support for Prisms!
  // void testAllRBBPrism6() { LOG_UNIT_TEST; test_box(PRISM6); }
  // void testAllRBBPrism18() { LOG_UNIT_TEST; test_box(PRISM18); }

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

  // Extruding a disk into 2 layers gives us twice as many hexes as we
  // had quads
  void testAllRBBCylinder10() { LOG_UNIT_TEST; test_cylinder(0); }
  void testAllRBBCylinder80() { LOG_UNIT_TEST; test_cylinder(1); }

  // 0 refinements of our default sphere gives us 6 RBB quad faces or
  // 12 RBB triangles
  void testAllRBBSphere6() { LOG_UNIT_TEST; test_sphere(0); }
  void testAllRBBSphere24() { LOG_UNIT_TEST; test_sphere(1); }
  void testAllRBBSphere96() { LOG_UNIT_TEST; test_sphere(2); }

  void testAllRBBTri6Sphere12() { LOG_UNIT_TEST; test_sphere(0, TRI6); }
  void testAllRBBTri6Sphere48() { LOG_UNIT_TEST; test_sphere(1, TRI6); }
  void testAllRBBTri6Sphere192() { LOG_UNIT_TEST; test_sphere(2, TRI6); }
};


CPPUNIT_TEST_SUITE_REGISTRATION( AllRBBTest );
