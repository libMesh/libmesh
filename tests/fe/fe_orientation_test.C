// A basis whose shape functions follow the order of an element's vertices, such as HIERARCHIC of
// order two and above, reads that order back from Elem::edge_orientation() and
// Elem::face_orientation(). Every degree of freedom of such a basis belongs to one vertex, edge, or
// face of the element, or to the element's interior, and its shape function follows the orientation
// of that one entity. These tests hold the families to that, which is what lets a cache of
// reference shape functions carry one entry per entity and orientation.
//
// The orientation indices of two elements sharing a side differ, since each element numbers the
// side's vertices its own way, and the indices exist so that both elements arrive at the same shape
// function for a degree of freedom they share. The conformity test holds the families to that.

#include <libmesh/dof_map.h>
#include <libmesh/elem.h>
#include <libmesh/enum_elem_type.h>
#include <libmesh/enum_fe_family.h>
#include <libmesh/enum_order.h>
#include <libmesh/equation_systems.h>
#include <libmesh/fe_interface.h>
#include <libmesh/fe_map.h>
#include <libmesh/fe_type.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/parallel_implementation.h>
#include <libmesh/quadrature_gauss.h>
#include <libmesh/remote_elem.h>
#include <libmesh/system.h>

#include "test_comm.h"
#include "libmesh_cppunit.h"

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

    // The geometric element the family needs, reached from the first-order type that the mesh
    // generator builds
    ElemType base_type = elem_type;
    bool second_order = false, complete_order = false;

    switch (elem_type)
      {
      case TRI6:     base_type = TRI3;    second_order = true;   break;
      case QUAD9:    base_type = QUAD4;   second_order = true;   break;
      case HEX27:    base_type = HEX8;    second_order = true;   break;
      case TET14:    base_type = TET4;    complete_order = true; break;
      case PRISM21:  base_type = PRISM6;  complete_order = true; break;
      default:       libmesh_error();
      }

    const unsigned int dim = Elem::type_to_dim_map[elem_type];

    if (dim == 2)
      MeshTools::Generation::build_square(*_mesh, 2, 2, 0., 1., 0., 1., base_type);
    else
      MeshTools::Generation::build_cube(*_mesh, 2, 2, 2, 0., 1., 0., 1., 0., 1., base_type);

    if (second_order)
      _mesh->all_second_order();
    if (complete_order)
      _mesh->all_complete_order();

    // Skew the mesh so that no symmetry of an element masks a shape function that follows the
    // wrong entity, and so that the vertices of an entity order differently from one element to
    // the next. The skew is linear and orientation preserving, which leaves every element affine
    // and valid, so that the map from a physical point back to a reference point is exact.
    for (auto * node : _mesh->node_ptr_range())
      {
        Node & p = *node;
        const Real y = p(1), z = (LIBMESH_DIM > 2) ? p(2) : 0.;

        p(0) += 0.3 * y + 0.17 * z;
        if (LIBMESH_DIM > 1)
          p(1) += 0.23 * z;
      }
  }

  void tearDown() { _mesh.reset(); }

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

    // How many entities were seen carrying more than one orientation, which is what the test
    // compares shape functions across
    unsigned int n_compared = 0;

    for (auto * elem : _mesh->active_local_element_ptr_range())
      {
        QGauss qrule(elem->dim(), order);
        qrule.init(*elem);
        const std::vector<Point> & points = qrule.get_points();

        // Per entity, the shape functions of its degrees of freedom, per orientation index
        std::map<Entity, std::map<unsigned int, std::vector<Real>>> seen;

        // Permuting an element leaves its vertices where they are and renumbers them, which walks
        // its entities through a range of orientations
        for (const auto p : make_range(elem->n_permutations()))
          {
            elem->permute(p);

            // libMesh numbers the degrees of freedom of an element by node, in node order, and
            // then those of the element's interior
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

    // A test that never saw an entity in two orientations would hold nothing
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

    // How many degrees of freedom shared across a side were compared
    unsigned int n_compared = 0;

    for (const auto round : make_range(n_rounds()))
      {
        // Permute by a function of the element id, so that every processor arrives at the same
        // numbering for an element it shares with another. A permutation moves an element's node
        // numbering and leaves its geometry and its degrees of freedom where they are.
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

                    // Solve the mapping to a tolerance well below the one the shape values are
                    // then held to, so that the comparison reads the shape functions rather than
                    // the accuracy of the mapping
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

    // A test that found no shared degree of freedom would hold nothing
    _mesh->comm().sum(n_compared);
    CPPUNIT_ASSERT(n_compared > 0);
  }

private:
  /// The tolerance the mapping from a physical point back to a reference point is solved to
  static constexpr Real map_tolerance = TOLERANCE * TOLERANCE;

  /**
   * How many rounds of permutation to walk the mesh through. Composing a permutation per round
   * carries the elements through a range of relative numberings of the sides they share.
   */
  unsigned int n_rounds() const { return Elem::build(elem_type)->n_permutations(); }

  /// The shape functions of \p n_dofs degrees of freedom starting at \p first, at every point
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

  /**
   * Hold the shape functions seen for an entity in one orientation against those seen for it in
   * that same orientation before.
   */
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

// A family whose degrees of freedom sit only on the sides of an element has no shape function to
// report at an interior point, so such a family is held to the conformity test alone.
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
