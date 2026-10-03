// Tests that the shape functions of a degree of freedom on a vertex, edge, or face follow only
// that entity's orientation index (Elem::edge_orientation(), Elem::face_orientation()), and that
// elements sharing a side agree on the shape functions of the degrees of freedom they share.

#include "fe_test.h"

#include <libmesh/enum_elem_type.h>
#include <libmesh/enum_fe_family.h>
#include <libmesh/enum_order.h>
#include <libmesh/fe_map.h>
#include <libmesh/fe_type.h>
#include <libmesh/parallel_implementation.h>
#include <libmesh/remote_elem.h>

#include <algorithm>
#include <cmath>
#include <map>
#include <memory>
#include <utility>
#include <vector>

using namespace libMesh;

/// The entity of an element that owns a degree of freedom, as a kind and an index within that kind
enum EntityKind : unsigned int { VERTEX = 0, EDGE = 1, FACE = 2, INTERIOR = 3 };
typedef std::pair<unsigned int, unsigned int> Entity;

template <ElemType elem_type, Order order, FEFamily family>
class FEOrientationTest : public CppUnit::TestCase
{
protected:
  std::unique_ptr<Mesh> _mesh;
  std::string libmesh_suite_name;

public:
  void setUp()
  {
    _mesh = std::make_unique<Mesh>(*TestCommWorld);
    build_mesh();

    // An affine skew, so that no element symmetry masks a shape function following the wrong
    // entity, and so that inverse_map() is exact
    SkewFunc skew_func;
    MeshTools::Modification::redistribute(*_mesh, skew_func);
    _mesh->complete_preparation();
  }

  void tearDown() { _mesh.reset(); }

  /// Builds the mesh the tests run on, before it is skewed
  virtual void build_mesh()
  {
    const unsigned int dim = Elem::type_to_dim_map[elem_type];
    MeshTools::Generation::build_cube(*_mesh, 2, 2, 2*(dim > 2),
                                      0., 1., 0., 1., 0., 1., elem_type);
  }

  /// The entity of \p elem that owns the degrees of freedom sitting on node \p n
  Entity node_entity(const Elem & elem, const unsigned int n)
  {
    if (elem.is_vertex(n))
      return {VERTEX, n};

    if (elem.is_edge(n))
      for (const auto e : make_range(elem.n_edges()))
        if (elem.is_node_on_edge(n, e))
          return {EDGE, e};

    if (elem.is_face(n))
      for (const auto f : make_range(elem.n_faces()))
        if (elem.is_node_on_side(n, f))
          return {FACE, f};

    return {INTERIOR, 0};
  }

  /// The orientation index that the shape functions of \p entity follow
  unsigned int entity_orientation(const Elem & elem, const Entity & entity)
  {
    switch (entity.first)
      {
      case EDGE: return elem.edge_orientation(entity.second);
      case FACE: return elem.face_orientation(entity.second);
      default:   return 0;
      }
  }

  void test_orientation_determines_shapes()
  {
    LOG_UNIT_TEST;

    const FEType fe_type(order, family);

    // Number of comparisons between repeat visits to an entity orientation
    unsigned int n_compared = 0;

    for (auto * elem : _mesh->active_local_element_ptr_range())
      {
        QGauss qrule(elem->dim(), order);
        qrule.init(*elem);
        const std::vector<Point> & points = qrule.get_points();

        // Shape function values seen, per entity and orientation index
        std::map<Entity, std::map<unsigned int, std::vector<Real>>> seen;

        // Permuting renumbers the vertices without moving them, cycling the entity orientations
        for (const auto p : make_range(elem->n_permutations()))
          {
            elem->permute(p);

            // Element degrees of freedom are numbered by node, then interior
            unsigned int dof = 0;

            for (const auto n : elem->node_index_range())
              {
                const unsigned int n_dofs = FEInterface::n_dofs_at_node(fe_type, elem, n);
                if (!n_dofs)
                  continue;

                record(seen[node_entity(*elem, n)],
                       entity_orientation(*elem, node_entity(*elem, n)),
                       shapes(fe_type, *elem, dof, n_dofs, points),
                       n_compared);

                dof += n_dofs;
              }

            const unsigned int n_interior = FEInterface::n_dofs_per_elem(fe_type, elem);
            if (n_interior)
              record(seen[{INTERIOR, 0}], 0,
                     shapes(fe_type, *elem, dof, n_interior, points), n_compared);
          }
      }

    // Guard against a vacuous pass
    _mesh->comm().sum(n_compared);
    CPPUNIT_ASSERT(n_compared > 0);
  }

  void test_conformity_across_sides()
  {
    LOG_UNIT_TEST;

    const FEType fe_type(order, family);

    EquationSystems es(*_mesh);
    System & sys = es.add_system<System>("orientation");
    sys.add_variable("u", fe_type);
    es.init();
    const DofMap & dof_map = sys.get_dof_map();

    unsigned int n_compared = 0;

    for (const auto round : make_range(n_rounds()))
      {
        // Permute by element id so that every processor agrees on the numbering of a ghosted
        // element
        for (auto * elem : _mesh->element_ptr_range())
          elem->permute((round + elem->id()) % elem->n_permutations());

        for (const auto * elem : _mesh->active_local_element_ptr_range())
          {
            std::vector<dof_id_type> dofs;
            dof_map.dof_indices(elem, dofs);

            for (const auto s : make_range(elem->n_sides()))
              {
                const Elem * neighbor = elem->neighbor_ptr(s);
                if (!neighbor || neighbor == remote_elem)
                  continue;

                std::vector<dof_id_type> neighbor_dofs;
                dof_map.dof_indices(neighbor, neighbor_dofs);

                std::unique_ptr<const Elem> side = elem->build_side_ptr(s);
                QGauss side_rule(side->dim(), order);
                side_rule.init(*side);

                for (const auto & side_point : side_rule.get_points())
                  {
                    const Point physical = FEMap::map(side->dim(), side.get(), side_point);

                    // Map far tighter than the shape value tolerance, so mapping error cannot
                    // mask or cause a mismatch
                    const Point point =
                      FEMap::inverse_map(elem->dim(), elem, physical, map_tolerance);
                    const Point neighbor_point =
                      FEMap::inverse_map(neighbor->dim(), neighbor, physical, map_tolerance);

                    for (const auto i : index_range(dofs))
                      {
                        const auto it =
                          std::find(neighbor_dofs.begin(), neighbor_dofs.end(), dofs[i]);
                        if (it == neighbor_dofs.end())
                          continue;

                        const auto j = std::distance(neighbor_dofs.begin(), it);

                        LIBMESH_ASSERT_FP_EQUAL(
                          FEInterface::shape(fe_type, elem, i, point),
                          FEInterface::shape(fe_type, neighbor, j, neighbor_point),
                          TOLERANCE);
                        ++n_compared;
                      }
                  }
              }
          }
      }

    // Guard against a vacuous pass
    _mesh->comm().sum(n_compared);
    CPPUNIT_ASSERT(n_compared > 0);
  }

private:
  /// The tolerance for inverse_map()
  static constexpr Real map_tolerance = TOLERANCE * TOLERANCE;

  /// Rounds of permutation, which vary the relative numbering of shared sides
  unsigned int n_rounds() const
  {
    unsigned int rounds = 0;
    for (const auto * elem : _mesh->element_ptr_range())
      rounds = std::max(rounds, elem->n_permutations());
    _mesh->comm().max(rounds);
    return rounds;
  }

  /// The values of shape functions \p first through \p first + \p n_dofs - 1 at \p points
  std::vector<Real> shapes(const FEType & fe_type,
                           const Elem & elem,
                           const unsigned int first,
                           const unsigned int n_dofs,
                           const std::vector<Point> & points)
  {
    std::vector<Real> values;

    for (const auto & point : points)
      for (const auto i : make_range(n_dofs))
        values.push_back(FEInterface::shape(fe_type, &elem, first + i, point));

    return values;
  }

  /// Record \p values for \p orientation, or check them against those already recorded
  void record(std::map<unsigned int, std::vector<Real>> & seen,
              const unsigned int orientation,
              const std::vector<Real> & values,
              unsigned int & n_compared)
  {
    const auto [entry, inserted] = seen.emplace(orientation, values);
    if (inserted)
      return;

    CPPUNIT_ASSERT_EQUAL(entry->second.size(), values.size());
    for (const auto k : index_range(values))
      LIBMESH_ASSERT_FP_EQUAL(entry->second[k], values[k], TOLERANCE * TOLERANCE);

    if (seen.size() > 1)
      ++n_compared;
  }
};

/**
 * The orientation tests on a mesh mixing element types, where each pair of neighbors shares a
 * face that the two elements parametrize differently: a prism and a hexahedron share a
 * quadrilateral, two prisms extruded in different directions share a quadrilateral whose axes
 * they swap, and a prism and a tetrahedron share a triangle.
 */
template <Order order, FEFamily family>
class FEHybridOrientationTest : public FEOrientationTest<INVALID_ELEM, order, family>
{
public:
  void build_mesh() override
  {
    Mesh & mesh = *this->_mesh;

    // A prism extruded along z from the triangle (0,0), (1,0), (0,1), with a hexahedron below
    // y = 0, a prism extruded along y beyond x = 0, and a tetrahedron above z = 1
    const std::vector<Point> points =
      {{0,0,0}, {1,0,0}, {0,1,0}, {0,0,1}, {1,0,1}, {0,1,1},
       {0,-1,0}, {1,-1,0}, {0,-1,1}, {1,-1,1},
       {-1,0,0}, {-1,1,0},
       {0.25,0.25,2}};
    for (const auto i : index_range(points))
      mesh.add_point(points[i], i);

    const std::vector<std::pair<ElemType, std::vector<dof_id_type>>> elems =
      {{PRISM6, {0, 1, 2, 3, 4, 5}},
       {HEX8,   {6, 7, 1, 0, 8, 9, 4, 3}},
       {PRISM6, {0, 10, 3, 2, 11, 5}},
       {TET4,   {3, 4, 5, 12}}};
    for (const auto & [type, nodes] : elems)
      {
        Elem * elem = mesh.add_elem(Elem::build(type));
        for (const auto n : index_range(nodes))
          elem->set_node(n, mesh.node_ptr(nodes[n]));
      }

    mesh.prepare_for_use();
    mesh.all_complete_order();
  }
};

#define ORIENTATIONTEST                                 \
  CPPUNIT_TEST( test_orientation_determines_shapes );    \
  CPPUNIT_TEST( test_conformity_across_sides )

#define INSTANTIATE_FEORIENTATIONTEST(elemtype, order, family)           \
  class FEOrientationTest_##family##_##order##_##elemtype :              \
    public FEOrientationTest<elemtype, order, family> {                  \
  public:                                                               \
  FEOrientationTest_##family##_##order##_##elemtype() :                  \
    FEOrientationTest<elemtype, order, family>() {                       \
    if (unitlog->summarized_logs_enabled())                                     \
      this->libmesh_suite_name = "FEOrientationTest";                   \
    else                                                                \
      this->libmesh_suite_name =                                        \
        "FEOrientationTest_" #family "_" #order "_" #elemtype;          \
  }                                                                     \
  CPPUNIT_TEST_SUITE( FEOrientationTest_##family##_##order##_##elemtype ); \
  ORIENTATIONTEST;                                                      \
  CPPUNIT_TEST_SUITE_END();                                             \
  };                                                                    \
                                                                        \
  CPPUNIT_TEST_SUITE_REGISTRATION( FEOrientationTest_##family##_##order##_##elemtype )

// Side-only families have no interior shape values, so they get only the conformity test
#define INSTANTIATE_FECONFORMITYTEST(elemtype, order, family)             \
  class FEConformityTest_##family##_##order##_##elemtype :                \
    public FEOrientationTest<elemtype, order, family> {                   \
  public:                                                                 \
  FEConformityTest_##family##_##order##_##elemtype() :                    \
    FEOrientationTest<elemtype, order, family>() {                        \
    if (unitlog->summarized_logs_enabled())                               \
      this->libmesh_suite_name = "FEConformityTest";                      \
    else                                                                  \
      this->libmesh_suite_name =                                          \
        "FEConformityTest_" #family "_" #order "_" #elemtype;             \
  }                                                                       \
  CPPUNIT_TEST_SUITE( FEConformityTest_##family##_##order##_##elemtype ); \
  CPPUNIT_TEST( test_conformity_across_sides );                           \
  CPPUNIT_TEST_SUITE_END();                                               \
  };                                                                      \
                                                                          \
  CPPUNIT_TEST_SUITE_REGISTRATION( FEConformityTest_##family##_##order##_##elemtype )

INSTANTIATE_FEORIENTATIONTEST(TRI6,    THIRD,  HIERARCHIC);
INSTANTIATE_FEORIENTATIONTEST(QUAD9,   THIRD,  HIERARCHIC);
INSTANTIATE_FEORIENTATIONTEST(TET14,   THIRD,  HIERARCHIC);
INSTANTIATE_FEORIENTATIONTEST(TET14,   FOURTH, HIERARCHIC);
INSTANTIATE_FEORIENTATIONTEST(HEX27,   THIRD,  HIERARCHIC);
INSTANTIATE_FEORIENTATIONTEST(HEX27,   FOURTH, HIERARCHIC);
INSTANTIATE_FEORIENTATIONTEST(PRISM21, THIRD,  HIERARCHIC);
INSTANTIATE_FEORIENTATIONTEST(PRISM21, FOURTH, HIERARCHIC);
INSTANTIATE_FECONFORMITYTEST(HEX27,   THIRD,  SIDE_HIERARCHIC);
INSTANTIATE_FECONFORMITYTEST(HEX27,   FOURTH, SIDE_HIERARCHIC);
INSTANTIATE_FEORIENTATIONTEST(HEX27,   FIFTH,  HIERARCHIC);

#define INSTANTIATE_FEHYBRIDORIENTATIONTEST(order, family)                \
  class FEHybridOrientationTest_##family##_##order :                      \
    public FEHybridOrientationTest<order, family> {                       \
  public:                                                                 \
  FEHybridOrientationTest_##family##_##order() :                          \
    FEHybridOrientationTest<order, family>() {                            \
    if (unitlog->summarized_logs_enabled())                               \
      this->libmesh_suite_name = "FEHybridOrientationTest";               \
    else                                                                  \
      this->libmesh_suite_name =                                          \
        "FEHybridOrientationTest_" #family "_" #order;                    \
  }                                                                       \
  CPPUNIT_TEST_SUITE( FEHybridOrientationTest_##family##_##order );       \
  ORIENTATIONTEST;                                                        \
  CPPUNIT_TEST_SUITE_END();                                               \
  };                                                                      \
                                                                          \
  CPPUNIT_TEST_SUITE_REGISTRATION( FEHybridOrientationTest_##family##_##order )

INSTANTIATE_FEHYBRIDORIENTATIONTEST(THIRD,  HIERARCHIC);
INSTANTIATE_FEHYBRIDORIENTATIONTEST(FOURTH, HIERARCHIC);
