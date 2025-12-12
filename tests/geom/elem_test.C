#include "elem_test.h"

#include <libmesh/boundary_info.h>
#include <libmesh/enum_elem_quality.h>
#include <libmesh/elem_side_builder.h>
#include <libmesh/mesh_modification.h>
#include <libmesh/mesh_refinement.h>
#include <libmesh/parallel_implementation.h>
#include <libmesh/enum_to_string.h>

using namespace libMesh;

template <ElemType elem_type>
class ElemTest : public PerElemTest<elem_type> {
public:
  void test_bounding_box()
  {
    LOG_UNIT_TEST;

    for (const auto & elem :
         this->_mesh->active_local_element_ptr_range())
      {
        const BoundingBox bbox = elem->loose_bounding_box();

        // The "loose" bounding box should actually be pretty tight
        // in most of these cases, but for weirdly aligned triangles
        // (such as occur in pyramid elements) it won't be, so we'll
        // just test against a widened bounding box.
        BoundingBox wide_bbox(elem->point(0), elem->point(0));

        for (unsigned int n = 0; n != elem->n_nodes(); ++n)
          {
            const Point & p = elem->point(n);

            if (!elem->infinite())
              CPPUNIT_ASSERT(bbox.contains_point(p));

            wide_bbox.union_with
              (BoundingBox(elem->point(n), elem->point(n)));
          }

        wide_bbox.scale(1. / 3.);

        if (!elem->infinite() && elem->dim())
          {
            CPPUNIT_ASSERT(!bbox.contains_point(wide_bbox.min()));
            CPPUNIT_ASSERT(!bbox.contains_point(wide_bbox.max()));
          }
      }
  }

  void test_quality()
  {
    LOG_UNIT_TEST;

    for (const auto & elem : this->_mesh->active_local_element_ptr_range())
      {
        // EDGE_LENGTH_RATIO is one metric that is defined on all elements
        const Real edge_length_ratio = elem->quality(EDGE_LENGTH_RATIO);

        // We use "0" to mean infinity rather than inf or NaN, and
        // every quality other than that should be 1 or larger (worse)
        CPPUNIT_ASSERT_LESSEQUAL(edge_length_ratio, Real(1)); // 1 <= edge_length_ratio

        // We're building isotropic meshes, where even elements
        // dissected from cubes ought to have tolerable quality.
        //
        // Worst I see is 2 on tets, but let's add a little tolerance
        // in case we decide to play with rotated meshes here later
        CPPUNIT_ASSERT_LESSEQUAL(Real(2+TOLERANCE), edge_length_ratio); // edge_length_ratio <= 2

        // The MIN_ANGLE and MAX_ANGLE quality metrics are also defined on all elements
        const Real min_angle = elem->quality(MIN_ANGLE);
        const Real max_angle = elem->quality(MAX_ANGLE);

        // Reference Quads/Hexes have maximum internal angles of 90
        // degrees, but pyramids actually have an obtuse interior
        // angle of acos(-1/3) ~ 109.47 deg formed by edge pairs
        // {(0, 4), (2, 4)} and {(1,4), (3,4)}.
        if (elem->type() != C0POLYGON)
          CPPUNIT_ASSERT_LESSEQUAL((std::acos(Real(-1)/3) * 180 / libMesh::pi) + TOLERANCE, max_angle);

        // We're doing some polygon testing with a squashed pentagon
        // that has a 135 degree angle
        CPPUNIT_ASSERT_LESSEQUAL(135 + TOLERANCE, max_angle);

        // Notes on minimum angle we expect to see:
        // 1.) 1D Elements don't have interior angles, so the base
        //     class implementation currently returns 0 for those
        //     elements.
        // 2.) Reference triangles/tetrahedra have min interior angle
        //     of 45 degrees, however, here we are checking Tets in a
        //     build_cube() mesh which have angles as small as
        //     acos(2/sqrt(6)) so we use that as our lower bound here.
        if (elem->dim() > 1)
          CPPUNIT_ASSERT_GREATEREQUAL((std::acos(Real(2)/std::sqrt(Real(6))) * 180 / libMesh::pi) - TOLERANCE, min_angle);

        // MIN,MAX_DIHEDRAL_ANGLE are implemented for all 3D elements
        if (elem->dim() > 2)
          {
            const Real min_dihedral_angle = elem->quality(MIN_DIHEDRAL_ANGLE);
            const Real max_dihedral_angle = elem->quality(MAX_DIHEDRAL_ANGLE);

            // Debugging
            // libMesh::out << "Elem type: " << Utility::enum_to_string(elem->type())
            //              << ", min_dihedral_angle = " << min_dihedral_angle
            //              << ", max_dihedral_angle = " << max_dihedral_angle
            //              << std::endl;

            // Assert that we match expected values for build_cube() meshes.
            // * build_cube() meshes of hexes, tetrahedra, and prisms have max_dihedral_angle == 90
            // * build_cube() meshes of pyramids have max_dihedral_angle == 60 between adjacent triangular faces
            CPPUNIT_ASSERT_LESSEQUAL   (90 + TOLERANCE, max_dihedral_angle);

            // * build_cube() meshes of hexes have min_dihedral_angle == 90
            // * build_cube() meshes of prisms, tets, and pyramids have min_dihedral_angle == 45
            // * For the InfPrism tests, we construct a single
            //   InfPrism by hand with interior angle ~53.13 deg. so that
            //   is the minimum dihedral angle we expect in that case.
            CPPUNIT_ASSERT_GREATEREQUAL(45 - TOLERANCE, min_dihedral_angle);
          }

        // The JACOBIAN and SCALED_JACOBIAN quality metrics are also defined on all elements

        // The largest non-infinite elements in this test (Hexes and
        // Prisms) have a nodal Jacobian (volume) of 8, e.g. the 2x2x2
        // reference Hex. The infinite elements that we construct for
        // this testing have a max nodal area of 64 since they are
        // created (see setUp()) with max side lengths of 4 (4*4*4=64).
        const Real jac = elem->quality(JACOBIAN);
        if (elem->infinite())
          CPPUNIT_ASSERT_LESSEQUAL   (64 + TOLERANCE, jac);
        else
          CPPUNIT_ASSERT_LESSEQUAL   (8 + TOLERANCE, jac);

        // The smallest 2D/3D nodal areas for regular elements in this
        // test are found in tetrahedra, which have a minimum value of
        // 2.0. However, we return a default value of 1.0 for 0D and
        // 1D elements here, and we have a custom distorted-pentagon
        // C0Polygon with a 0.5 at 2 nodes and a custom C0Polyhedron
        // with a 1, so we handle those cases too.
        if (elem->dim() < 2 || elem->type() == C0POLYHEDRON)
          CPPUNIT_ASSERT_GREATEREQUAL(1 - TOLERANCE, jac);
        else if (elem->type() == C0POLYGON)
          CPPUNIT_ASSERT_GREATEREQUAL(0.5 - TOLERANCE, jac);
        else
          CPPUNIT_ASSERT_GREATEREQUAL(2 - TOLERANCE, jac);

        // The scaled Jacobian should always be <= 1. The minimum
        // scaled Jacobian value I observed was ~0.408248 for the
        // tetrahedral meshes. This is consistent with the way that
        // we generate Tet meshes with build_cube(), as we don't
        // simply refine the reference tetrahedron in that case.
        const Real scaled_jac = elem->quality(SCALED_JACOBIAN);
        CPPUNIT_ASSERT_LESSEQUAL   (1 + TOLERANCE, scaled_jac);
        CPPUNIT_ASSERT_GREATEREQUAL(    Real(0.4), scaled_jac);

        // Debugging
        // libMesh::out << "elem->type() = " << Utility::enum_to_string(elem->type())
        //              << ", jac = " << jac
        //              << ", scaled_jac = " << scaled_jac
        //              << std::endl;
      }
  }

  void test_maps()
  {
    LOG_UNIT_TEST;

    for (const auto & elem :
         this->_mesh->active_local_element_ptr_range())
      {
        for (const auto edge : elem->edge_index_range())
          for (const auto side_on_edge : elem->sides_on_edge(edge))
            for (const auto node_on_edge : elem->nodes_on_edge(edge))
              CPPUNIT_ASSERT(elem->is_node_on_side(node_on_edge, side_on_edge));

        for (const auto side : elem->side_index_range())
          for (const auto node_on_side : elem->nodes_on_side(side))
            CPPUNIT_ASSERT(elem->is_node_on_side(node_on_side, side));

        for (const auto edge : elem->edge_index_range())
          for (const auto node_on_edge : elem->nodes_on_edge(edge))
            CPPUNIT_ASSERT(elem->is_node_on_edge(node_on_edge, edge));

        for (const auto edge : elem->edge_index_range())
          for (const auto side_on_edge : elem->sides_on_edge(edge))
            CPPUNIT_ASSERT(elem->is_edge_on_side(edge, side_on_edge));
      }
  }

  void test_static_data()
  {
    LOG_UNIT_TEST;

    for (const auto & elem :
         this->_mesh->active_local_element_ptr_range())
      {
        CPPUNIT_ASSERT(elem->n_nodes() <= Elem::max_n_nodes);

        const ElemType etype = elem->type();

        CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(elem->dim()),
                             Elem::type_to_dim_map[etype]);
        CPPUNIT_ASSERT_EQUAL(elem->default_order(),
                             Elem::type_to_default_order_map[etype]);

        // If we have an element type with topology defined solely by
        // the type, then that should match the runtime topology.  If
        // not, then we should be aware of it.
        if (!elem->runtime_topology())
          {
            CPPUNIT_ASSERT_EQUAL(elem->n_nodes(), Elem::type_to_n_nodes_map[etype]);
            CPPUNIT_ASSERT_EQUAL(elem->n_sides(), Elem::type_to_n_sides_map[etype]);
            CPPUNIT_ASSERT_EQUAL(elem->n_edges(), Elem::type_to_n_edges_map[etype]);
          }
        else
          {
            CPPUNIT_ASSERT_EQUAL(invalid_uint, Elem::type_to_n_nodes_map[etype]);
            CPPUNIT_ASSERT_EQUAL(invalid_uint, Elem::type_to_n_sides_map[etype]);
            CPPUNIT_ASSERT_EQUAL(invalid_uint, Elem::type_to_n_edges_map[etype]);
          }
      }
  }

  void test_contains_point_node()
  {
    LOG_UNIT_TEST;

    for (const auto & elem :
         this->_mesh->active_local_element_ptr_range())
      {
        if (elem->infinite())
          continue;

        for (const auto n : elem->node_index_range())
#ifndef LIBMESH_ENABLE_EXCEPTIONS
          // If this node has a singular Jacobian, we need exceptions in order
          // to catch the failed inverse_map solve and return the singular
          // master point. Therefore, if we don't have exceptions and we're
          // at a singular node, we can't test this. As of the writing of
          // this comment, this issue exists for only Pyramid elements at
          // the apex.
          if (elem->local_singular_node(elem->point(n), TOLERANCE*TOLERANCE) == invalid_uint)
#endif
            CPPUNIT_ASSERT(elem->contains_point(elem->point(n)));
      }
  }

  void test_permute()
  {
    LOG_UNIT_TEST;

    for (const auto & elem :
         this->_mesh->active_local_element_ptr_range())
      {
        if (elem->infinite())
          continue;

        const Point centroid = elem->true_centroid();
        const Point vertex_avg = elem->vertex_average();
        Point quasicc;
        if (elem->dim() < 3)
          quasicc = elem->quasicircumcenter();

        for (const auto p : IntRange<unsigned int>(0, elem->n_permutations()))
          {
            elem->permute(p);
            CPPUNIT_ASSERT(elem->has_invertible_map());
            const Point new_centroid = elem->true_centroid();
            const Point new_vertex_avg = elem->vertex_average();
            Point new_quasicc;
            if (elem->dim() < 3)
              new_quasicc = elem->quasicircumcenter();
            for (const auto d : make_range(LIBMESH_DIM))
              {
                // Getting a little FP error from Pyramid18
                LIBMESH_ASSERT_FP_EQUAL(centroid(d), new_centroid(d),
                                        TOLERANCE*std::sqrt(TOLERANCE));
                LIBMESH_ASSERT_FP_EQUAL(vertex_avg(d), new_vertex_avg(d),
                                        TOLERANCE*TOLERANCE);
                LIBMESH_ASSERT_FP_EQUAL(quasicc(d), new_quasicc(d),
                                        TOLERANCE*TOLERANCE);
              }
          }
      }
  }

  void test_flip()
  {
    LOG_UNIT_TEST;

    BoundaryInfo & boundary_info = this->_mesh->get_boundary_info();

    for (const auto & elem :
         this->_mesh->active_local_element_ptr_range())
      {
        if (elem->infinite())
          continue;

        if (elem->type() == C0POLYHEDRON)
          continue;

        const Point vertex_avg = elem->vertex_average();

        const unsigned int n_sides = elem->n_sides();
        std::vector<std::set<Point*>> side_nodes(n_sides);
        std::vector<Elem*> neighbors(n_sides);
        std::vector<std::vector<boundary_id_type>> bcids(n_sides);
        for (auto s : make_range(n_sides))
          {
            for (auto n : elem->nodes_on_side(s))
              side_nodes[s].insert(elem->node_ptr(n));
            neighbors[s] = elem->neighbor_ptr(s);
            boundary_info.boundary_ids(elem, s, bcids[s]);
          }

        elem->flip(&boundary_info);

        // We should just be flipped, not twisted, so our map should
        // still be affine.
        // ... except for stupid singular pyramid maps
        // ... or the polygons we're deliberately testing non-affine
        if ((elem->dim() < 3 ||
             elem->n_vertices() != 5) &&
            elem_type != C0POLYGON)
          CPPUNIT_ASSERT(elem->has_affine_map());
        else if (elem_type == C0POLYGON)
          CPPUNIT_ASSERT(!elem->has_affine_map());

        // The neighbors and bcids should have flipped in a way
        // consistently with the nodes (unless this is a 0D NodeElem)
        bool something_changed = false;
        for (auto s : make_range(n_sides))
          {
            std::set<Point*> new_side_nodes;
            for (auto n : elem->nodes_on_side(s))
              new_side_nodes.insert(elem->node_ptr(n));

            std::vector<boundary_id_type> new_bcids;
            boundary_info.boundary_ids(elem, s, new_bcids);

            unsigned int old_side = libMesh::invalid_uint;
            for (auto os : make_range(n_sides))
              if (new_side_nodes == side_nodes[os])
                old_side = os;

            if (old_side != s)
              something_changed = true;

            CPPUNIT_ASSERT(old_side != libMesh::invalid_uint);

            CPPUNIT_ASSERT(neighbors[old_side] ==
                           elem->neighbor_ptr(s));

            CPPUNIT_ASSERT(bcids[old_side] == new_bcids);
          }
        CPPUNIT_ASSERT(!elem->dim() || something_changed);

        const Point new_vertex_avg = elem->vertex_average();
        for (const auto d : make_range(LIBMESH_DIM))
          LIBMESH_ASSERT_FP_EQUAL(vertex_avg(d), new_vertex_avg(d),
                                  TOLERANCE*TOLERANCE);
      }
  }

  void test_orient()
  {
    LOG_UNIT_TEST;

    BoundaryInfo & boundary_info = this->_mesh->get_boundary_info();

    for (const auto & elem :
         this->_mesh->active_local_element_ptr_range())
      {
        if (elem->infinite())
          continue;

        const Point vertex_avg = elem->vertex_average();

        const unsigned int n_sides = elem->n_sides();
        std::vector<std::set<Point*>> side_nodes(n_sides);
        std::vector<Elem*> neighbors(n_sides);
        std::vector<std::vector<boundary_id_type>> bcids(n_sides);
        for (auto s : make_range(n_sides))
          {
            for (auto n : elem->nodes_on_side(s))
              side_nodes[s].insert(elem->node_ptr(n));
            neighbors[s] = elem->neighbor_ptr(s);
            boundary_info.boundary_ids(elem, s, bcids[s]);
          }

        CPPUNIT_ASSERT(!elem->is_flipped());

        if (elem->id()%2)
          {
            elem->flip(&boundary_info);
            CPPUNIT_ASSERT(elem->is_flipped());
          }

        elem->orient(&boundary_info);
        CPPUNIT_ASSERT(!elem->is_flipped());

        // Our map should still be affine.
        // ... except for stupid singular pyramid maps
        // ... or the polygons we're deliberately testing non-affine
        // ... or the polyhedra we deliberately don't define "affine
        // map" for.
        if ((elem->dim() < 3 ||
             elem->n_vertices() != 5) &&
            elem_type != C0POLYGON &&
            elem_type != C0POLYHEDRON)
          CPPUNIT_ASSERT(elem->has_affine_map());
        else if (elem_type == C0POLYGON &&
                 elem_type != C0POLYHEDRON)
          CPPUNIT_ASSERT(!elem->has_affine_map());

        // The neighbors and bcids should have flipped back to where
        // they were.
        for (auto s : make_range(n_sides))
          {
            std::set<Point*> new_side_nodes;
            for (auto n : elem->nodes_on_side(s))
              new_side_nodes.insert(elem->node_ptr(n));

            std::vector<boundary_id_type> new_bcids;
            boundary_info.boundary_ids(elem, s, new_bcids);

            CPPUNIT_ASSERT(side_nodes[s] ==
                           new_side_nodes);

            CPPUNIT_ASSERT(neighbors[s] ==
                           elem->neighbor_ptr(s));

            CPPUNIT_ASSERT(bcids[s] == new_bcids);
          }

        const Point new_vertex_avg = elem->vertex_average();
        for (const auto d : make_range(LIBMESH_DIM))
          LIBMESH_ASSERT_FP_EQUAL(vertex_avg(d), new_vertex_avg(d),
                                  TOLERANCE*TOLERANCE);
      }
  }

  void test_orient_elements()
  {
    LOG_UNIT_TEST;

    const Mesh old_mesh {*this->_mesh};

    BoundaryInfo & boundary_info = this->_mesh->get_boundary_info();
    const BoundaryInfo & old_boundary_info = old_mesh.get_boundary_info();
    CPPUNIT_ASSERT(&boundary_info != &old_boundary_info);

    for (const auto & elem :
         this->_mesh->active_local_element_ptr_range())
      {
        if (elem->infinite())
          continue;

        if (elem->id()%2)
          {
            elem->flip(&boundary_info);
            CPPUNIT_ASSERT(elem->is_flipped());
          }
      }

    MeshTools::Modification::orient_elements(*this->_mesh);

    // I should really create a MeshBase::operator==()...
    for (const auto & elem :
         this->_mesh->active_local_element_ptr_range())
      {
        const Elem & old_elem = old_mesh.elem_ref(elem->id());

        CPPUNIT_ASSERT(!elem->is_flipped());

        // Elem::operator==() uses node ids to compare
        CPPUNIT_ASSERT(*elem == old_elem);

        const unsigned int n_sides = elem->n_sides();
        for (auto s : make_range(n_sides))
          {
            std::vector<boundary_id_type> bcids, old_bcids;
            boundary_info.boundary_ids(elem, s, bcids);
            old_boundary_info.boundary_ids(&old_elem, s, old_bcids);
            CPPUNIT_ASSERT(bcids == old_bcids);

            if (elem->neighbor_ptr(s))
              {
                CPPUNIT_ASSERT(old_elem.neighbor_ptr(s));
                CPPUNIT_ASSERT_EQUAL(elem->neighbor_ptr(s)->id(),
                                     old_elem.neighbor_ptr(s)->id());
              }
            else
              CPPUNIT_ASSERT(!old_elem.neighbor_ptr(s));
          }
      }
  }

  void test_center_node_on_side()
  {
    LOG_UNIT_TEST;

    for (const auto & elem :
         this->_mesh->active_local_element_ptr_range())
      for (const auto s : elem->side_index_range())
        {
          if (elem->type() == EDGE2 || elem->type() == EDGE3 || elem->type() == EDGE4)
            CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(s), elem->center_node_on_side(s));
          else if (elem->type() == TRI6 || elem->type() == TRI7)
            CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(s + 3), elem->center_node_on_side(s));
          else if (elem->type() == QUAD8 || elem->type() == QUAD9 ||
                   elem->type() == QUADSHELL8 || elem->type() == QUADSHELL9)
            CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(s + 4), elem->center_node_on_side(s));
          else if (elem->type() == HEX27)
            CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(s + 20), elem->center_node_on_side(s));
          else if (elem->type() == PRISM18 && s >= 1 && s <= 3)
            CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(s + 14), elem->center_node_on_side(s));
          else if ((elem->type() == PRISM20 ||
                    elem->type() == PRISM21) && s >= 1 && s <= 3)
            CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(s + 14), elem->center_node_on_side(s));
          else if (elem->type() == PRISM20 ||
                   elem->type() == PRISM21)
            CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(18 + (s == 4)), elem->center_node_on_side(s));
          else if (elem->type() == PYRAMID14 && s == 4)
            CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(13), elem->center_node_on_side(s));
          else if (elem->type() == PYRAMID18)
            {
              if (s < 4)
                CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(s + 14), elem->center_node_on_side(s));
              else
                CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(13), elem->center_node_on_side(s));
            }
          else
            CPPUNIT_ASSERT_EQUAL(invalid_uint, elem->center_node_on_side(s));
        }
  }

  void test_side_type()
  {
    LOG_UNIT_TEST;

    for (const auto & elem :
         this->_mesh->active_local_element_ptr_range())
      for (const auto s : elem->side_index_range())
        CPPUNIT_ASSERT_EQUAL(elem->build_side_ptr(s)->type(), elem->side_type(s));
  }

  void test_side_subdomain()
  {
    LOG_UNIT_TEST;

    for (const auto & elem :
         this->_mesh->active_local_element_ptr_range())
      {
        std::unique_ptr<const Elem> side;

        for (const auto s : elem->side_index_range())
          {
            CPPUNIT_ASSERT_EQUAL(elem->build_side_ptr(s)->subdomain_id(), elem->subdomain_id());

            elem->build_side_ptr(side, s);
            CPPUNIT_ASSERT_EQUAL(side->subdomain_id(), elem->subdomain_id());

            // We don't require that the "lightweight" side_ptr work
            // the same way
            // CPPUNIT_ASSERT_EQUAL(elem->side_ptr(s)->subdomain_id(), elem->subdomain_id());
          }
      }
  }

  void test_elem_side_builder()
  {
    LOG_UNIT_TEST;

    ElemSideBuilder cache;
    for (auto & elem : this->_mesh->active_local_element_ptr_range())
      for (const auto s : elem->side_index_range())
      {
        const auto side = elem->build_side_ptr(s);

        auto & cached_side = cache(*elem, s);
        CPPUNIT_ASSERT_EQUAL(side->type(), cached_side.type());
        for (const auto n : side->node_index_range())
          CPPUNIT_ASSERT_EQUAL(side->node_ref(n), cached_side.node_ref(n));

        const auto & const_cached_side = cache(const_cast<const Elem &>(*elem), s);
        CPPUNIT_ASSERT_EQUAL(side->type(), const_cached_side.type());
        for (const auto n : side->node_index_range())
          CPPUNIT_ASSERT_EQUAL(side->node_ref(n), const_cached_side.node_ref(n));
      }
  }

  void test_n_refinements(unsigned int n)
  {
#ifdef LIBMESH_ENABLE_AMR
    // We don't support refinement of all element types
    if (elem_type == EDGE4 ||
        elem_type == PRISM20 ||
        elem_type == PYRAMID5 ||
        elem_type == PYRAMID13 ||
        elem_type == PYRAMID14 ||
        elem_type == PYRAMID18 ||
        elem_type == C0POLYGON ||
        elem_type == C0POLYHEDRON)
      return;

    auto refining_mesh = this->_mesh->clone();

    MeshRefinement mr(*refining_mesh);
    mr.uniformly_refine(n);

    std::set<std::pair<dof_id_type, unsigned int>> parent_node_was_touched;
    std::set<std::pair<dof_id_type, unsigned int>> parent_child_was_touched;

    for (const Elem * elem : refining_mesh->active_element_ptr_range())
      {
        CPPUNIT_ASSERT_EQUAL(elem->level(), n);
        CPPUNIT_ASSERT(!elem->ancestor());
        CPPUNIT_ASSERT(elem->active());
        CPPUNIT_ASSERT(!elem->subactive());
        CPPUNIT_ASSERT(!elem->has_children());
        CPPUNIT_ASSERT(!elem->has_ancestor_children());
        CPPUNIT_ASSERT(!elem->interior_parent());

        const Elem * parent = elem->parent();
        CPPUNIT_ASSERT(parent);
        CPPUNIT_ASSERT(parent->ancestor());
        CPPUNIT_ASSERT(!parent->active());
        CPPUNIT_ASSERT(!parent->subactive());
        CPPUNIT_ASSERT(parent->has_children());
        CPPUNIT_ASSERT(!parent->has_ancestor_children());
        CPPUNIT_ASSERT(!parent->interior_parent());
        if (n == 1)
          {
            CPPUNIT_ASSERT_EQUAL(parent, elem->top_parent());
            CPPUNIT_ASSERT_EQUAL(parent, parent->top_parent());
          }
        else
          {
            CPPUNIT_ASSERT(parent != elem->top_parent());
            CPPUNIT_ASSERT(parent != parent->top_parent());
            CPPUNIT_ASSERT_EQUAL(elem->top_parent(), parent->top_parent());
          }

        CPPUNIT_ASSERT(parent->is_ancestor_of(elem));
        const unsigned int c = parent->which_child_am_i(elem);
        CPPUNIT_ASSERT(c < parent->n_children());
        CPPUNIT_ASSERT_EQUAL(elem, parent->child_ptr(c));
        parent_child_was_touched.emplace(parent->id(), c);

        CPPUNIT_ASSERT_EQUAL(parent->n_nodes_in_child(c), elem->n_nodes());
        for (auto n : make_range(elem->n_nodes()))
          {
            CPPUNIT_ASSERT_EQUAL(parent->is_vertex_on_child(c, n), elem->is_vertex(n));

            auto pn = parent->as_parent_node(c, n);
            CPPUNIT_ASSERT_EQUAL(pn, parent->get_node_index(elem->node_ptr(n)));
            if (pn == libMesh::invalid_uint)
              continue;
            CPPUNIT_ASSERT_EQUAL(parent->is_vertex_on_parent(c, n), parent->is_vertex(pn));
            parent_node_was_touched.emplace(parent->id(), pn);
          }

        for (auto s : make_range(parent->n_sides()))
          {
            if (parent->is_child_on_side(c,s))
              {
                auto parent_side = parent->build_side_ptr(s);

                // Implicitly assuming here that s is the child side
                // too - we support that now and hopefully won't have
                // to change it later
                auto child_side  = elem->build_side_ptr(s);

                // 2D Inf FE inverse_map not yet implemented?
                if (!parent_side->infinite())
                  for (const Node & node : child_side->node_ref_range())
                    CPPUNIT_ASSERT(parent_side->contains_point(node));

                // We can at least check for sharing some node with
                // a side, on both finite and infinite elements, for
                // those children that aren't just the middle elements
                // of triangular sides.
                bool shares_a_side_node = false;
                bool shares_a_vertex = false;
                for (const Node & node : child_side->node_ref_range())
                  if (parent_side->local_node(node.id()) != invalid_uint)
                    shares_a_side_node = true;
                for (const Node & node : elem->node_ref_range())
                  if (parent->local_node(node.id()) < parent->n_vertices())
                    shares_a_vertex = true;

                CPPUNIT_ASSERT(shares_a_side_node || !shares_a_vertex);

                if (elem->neighbor_ptr(s) && !elem->neighbor_ptr(s)->is_remote())
                  CPPUNIT_ASSERT_EQUAL(parent->child_neighbor(elem->neighbor_ptr(s)), elem);
              }
          }

        for (auto e : make_range(parent->n_edges()))
          {
            if (parent->is_child_on_edge(c,e))
              {
                auto parent_edge = parent->build_edge_ptr(e);

                // Implicitly assuming here that e is the child edge
                // too - we support that now and hopefully won't have
                // to change it later
                auto child_edge  = elem->build_edge_ptr(e);

                // 1D Inf FE inverse_map not yet implemented?
                if (!parent_edge->infinite())
                  for (const Node & node : child_edge->node_ref_range())
                    CPPUNIT_ASSERT(parent_edge->contains_point(node));

                // We can at least check for sharing some node with
                // an edge, on both finite and infinite elements
                bool shares_an_edge_node = false;
                for (const Node & node : child_edge->node_ref_range())
                  if (parent_edge->local_node(node.id()) != invalid_uint)
                    shares_an_edge_node = true;

                CPPUNIT_ASSERT(shares_an_edge_node);
              }
          }

        if (parent->has_affine_map())
          CPPUNIT_ASSERT(elem->has_affine_map());
      }

    // It's possible for a parent element on a distributed mesh to not
    // have all its children available on any one processor
    TestCommWorld->set_union(parent_child_was_touched);
    TestCommWorld->set_union(parent_node_was_touched);

    for (const Elem * elem : refining_mesh->local_element_ptr_range())
      {
        if (elem->active())
          continue;

        // With only one layer of refinement the family tree methods
        // should have the full number of elements, even if some are
        // remote.
        if (n == 1)
          {
            std::vector<const Elem *> family;
            elem->family_tree(family);
            CPPUNIT_ASSERT_EQUAL(family.size(),
                                 std::size_t(elem->n_children() + 1));

            family.clear();
            elem->total_family_tree(family);
            CPPUNIT_ASSERT_EQUAL(family.size(),
                                 std::size_t(elem->n_children() + 1));

            family.clear();
            elem->active_family_tree(family);
            CPPUNIT_ASSERT_EQUAL(family.size(),
                                 std::size_t(elem->n_children()));

            for (auto s : make_range(elem->n_sides()))
              {
                family.clear();
                elem->active_family_tree_by_side(family,s);
                if (!elem->build_side_ptr(s)->infinite())
                  CPPUNIT_ASSERT_EQUAL(double(family.size()),
                                       std::pow(2.0, int(elem->dim()-1)));
                else
                  CPPUNIT_ASSERT_EQUAL(double(family.size()),
                                       std::pow(2.0, int(elem->dim()-2)));
                for (const Elem * child : family)
                  {
                    if (child->is_remote())
                      continue;

                    unsigned int c = elem->which_child_am_i(child);
                    CPPUNIT_ASSERT(elem->is_child_on_side(c, s));
                  }
              }
          }

        if (elem->level() + 1 == n)
          {
            for (auto c : make_range(elem->n_children()))
              {
                auto it = parent_child_was_touched.find(std::make_pair(elem->id(), c));
                CPPUNIT_ASSERT(it != parent_child_was_touched.end());
              }

            for (auto n : make_range(elem->n_nodes()))
              {
                auto it = parent_node_was_touched.find(std::make_pair(elem->id(), n));
                CPPUNIT_ASSERT(it != parent_node_was_touched.end());
              }
          }
      }
#endif
  }

  void test_refinement()
  {
    LOG_UNIT_TEST;

    test_n_refinements(1);
  }

  void test_double_refinement()
  {
    LOG_UNIT_TEST;

    test_n_refinements(2);
  }

  void test_is_internal()
  {
    LOG_UNIT_TEST;

    for (const auto & elem :
         this->_mesh->active_local_element_ptr_range())
      for (const auto nd : elem->node_index_range())
        {
          if ((elem->type() == EDGE3 || elem->type() == EDGE4) && nd >= 2)
            CPPUNIT_ASSERT(elem->is_internal(nd));
          else if (elem->type() == HEX27 && nd == 26)
            CPPUNIT_ASSERT(elem->is_internal(nd));
          else if (elem->type() == PRISM21 && nd == 20)
            CPPUNIT_ASSERT(elem->is_internal(nd));
          else if ((elem->type() == QUAD9 || elem->type() == QUADSHELL9) && nd == 8)
            CPPUNIT_ASSERT(elem->is_internal(nd));
          else if (elem->type() == TRI7 && nd == 6)
            CPPUNIT_ASSERT(elem->is_internal(nd));
          else if (elem->type() == INFHEX18 && nd == 17)
            CPPUNIT_ASSERT(elem->is_internal(nd));
          else if (elem->type() == INFQUAD6 && nd == 5)
            CPPUNIT_ASSERT(elem->is_internal(nd));
          else
            CPPUNIT_ASSERT(!elem->is_internal(nd));
        }
  }

  void test_node_edge_map_consistency()
  {
    LOG_UNIT_TEST;

    for (const auto & elem : this->_mesh->active_local_element_ptr_range())
      {
        for (const auto nd : elem->node_index_range())
          {
            auto adjacent_edge_ids = elem->edges_adjacent_to_node(nd);

            if (elem->dim() < 2)
              {
                // 0D elements don't have edges.
                // 1D elements *are* "edges", but we don't consider
                // them to *have* edges, so assert that here.
                CPPUNIT_ASSERT(adjacent_edge_ids.empty());
              }
            else
              {
                // For 2D and 3D elements, on each edge which is
                // claimed to be adjacent, check that "nd" is indeed
                // on it
                for (const auto & edge_id : adjacent_edge_ids)
                  {
                    auto node_ids_on_edge = elem->nodes_on_edge(edge_id);
                    CPPUNIT_ASSERT(std::find(node_ids_on_edge.begin(), node_ids_on_edge.end(), nd) != node_ids_on_edge.end());
                  }
              }
          }
      }
  }

};

#define ELEMTEST                                \
  CPPUNIT_TEST( test_bounding_box );            \
  CPPUNIT_TEST( test_quality );                 \
  CPPUNIT_TEST( test_node_edge_map_consistency ); \
  CPPUNIT_TEST( test_maps );                    \
  CPPUNIT_TEST( test_static_data );             \
  CPPUNIT_TEST( test_permute );                 \
  CPPUNIT_TEST( test_flip );                    \
  CPPUNIT_TEST( test_orient );                  \
  CPPUNIT_TEST( test_orient_elements );         \
  CPPUNIT_TEST( test_contains_point_node );     \
  CPPUNIT_TEST( test_center_node_on_side );     \
  CPPUNIT_TEST( test_side_type );               \
  CPPUNIT_TEST( test_side_subdomain );          \
  CPPUNIT_TEST( test_elem_side_builder );       \
  CPPUNIT_TEST( test_refinement );              \
  CPPUNIT_TEST( test_double_refinement );       \
  CPPUNIT_TEST( test_is_internal )

#define INSTANTIATE_ELEMTEST(elemtype)                          \
  class ElemTest_##elemtype : public ElemTest<elemtype> {       \
  public:                                                       \
  ElemTest_##elemtype() :                                       \
    ElemTest<elemtype>() {                                      \
    if (unitlog->summarized_logs_enabled())                     \
      this->libmesh_suite_name = "ElemTest";                    \
    else                                                        \
      this->libmesh_suite_name = "ElemTest_" #elemtype;         \
  }                                                             \
  CPPUNIT_TEST_SUITE( ElemTest_##elemtype );                    \
  ELEMTEST;                                                     \
  CPPUNIT_TEST_SUITE_END();                                     \
  };                                                            \
                                                                \
  CPPUNIT_TEST_SUITE_REGISTRATION( ElemTest_##elemtype )

INSTANTIATE_ELEMTEST(NODEELEM);

INSTANTIATE_ELEMTEST(EDGE2);
INSTANTIATE_ELEMTEST(EDGE3);
INSTANTIATE_ELEMTEST(EDGE4);

#if LIBMESH_DIM > 1
INSTANTIATE_ELEMTEST(TRI3);
INSTANTIATE_ELEMTEST(TRISHELL3);
INSTANTIATE_ELEMTEST(TRI6);
INSTANTIATE_ELEMTEST(TRI7);

INSTANTIATE_ELEMTEST(QUAD4);
INSTANTIATE_ELEMTEST(QUADSHELL4);
INSTANTIATE_ELEMTEST(QUAD8);
INSTANTIATE_ELEMTEST(QUADSHELL8);
INSTANTIATE_ELEMTEST(QUAD9);
INSTANTIATE_ELEMTEST(QUADSHELL9);

// This just tests with one pentagon, but better than nothing
INSTANTIATE_ELEMTEST(C0POLYGON);

// And this is just one truncated cube; also better than nothing
INSTANTIATE_ELEMTEST(C0POLYHEDRON);

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
INSTANTIATE_ELEMTEST(INFQUAD4);
INSTANTIATE_ELEMTEST(INFQUAD6);
#endif
#endif // LIBMESH_DIM > 1

#if LIBMESH_DIM > 2
INSTANTIATE_ELEMTEST(TET4);
INSTANTIATE_ELEMTEST(TET10);
INSTANTIATE_ELEMTEST(TET14);

INSTANTIATE_ELEMTEST(HEX8);
INSTANTIATE_ELEMTEST(HEX20);
INSTANTIATE_ELEMTEST(HEX27);

INSTANTIATE_ELEMTEST(PRISM6);
INSTANTIATE_ELEMTEST(PRISM15);
INSTANTIATE_ELEMTEST(PRISM18);
INSTANTIATE_ELEMTEST(PRISM20);
INSTANTIATE_ELEMTEST(PRISM21);

// These tests use PointLocator, which uses contains_point(), which
// uses inverse_map(), which doesn't play nicely on Pyramids unless we
// have exceptions support
#ifdef LIBMESH_ENABLE_EXCEPTIONS
INSTANTIATE_ELEMTEST(PYRAMID5);
INSTANTIATE_ELEMTEST(PYRAMID13);
INSTANTIATE_ELEMTEST(PYRAMID14);
INSTANTIATE_ELEMTEST(PYRAMID18);
#endif

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
INSTANTIATE_ELEMTEST(INFHEX8);
INSTANTIATE_ELEMTEST(INFHEX16);
INSTANTIATE_ELEMTEST(INFHEX18);

INSTANTIATE_ELEMTEST(INFPRISM6);
INSTANTIATE_ELEMTEST(INFPRISM12);
#endif
#endif // LIBMESH_DIM > 2
