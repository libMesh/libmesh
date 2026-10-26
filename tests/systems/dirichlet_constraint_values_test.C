#include <libmesh/dirichlet_boundaries.h>
#include <libmesh/dof_map.h>
#include <libmesh/elem.h>
#include <libmesh/equation_systems.h>
#include <libmesh/function_base.h>
#include <libmesh/libmesh.h>
#include <libmesh/linear_implicit_system.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_refinement.h>
#include <libmesh/parallel.h>
#include <libmesh/periodic_boundary.h>
#include <libmesh/sparse_matrix.h>

#include "test_comm.h"
#include "libmesh_cppunit.h"

using namespace libMesh;

namespace
{
// build_square's boundary ids
const boundary_id_type bottom_id = 0, right_id = 1, top_id = 2, left_id = 3;

/**
 * A time-dependent boundary value, g = scale * (1 + t) * (1 + x + 2 y)
 */
class RampFunction : public FunctionBase<Number>
{
public:
  RampFunction(const Real scale) : _scale(scale) {}

  virtual std::unique_ptr<FunctionBase<Number>> clone() const override
  {
    return std::make_unique<RampFunction>(_scale);
  }

  virtual Number operator()(const Point & p, const Real time = 0.) override
  {
    return _scale * (1 + time) * (1 + p(0) + 2 * p(1));
  }

  virtual void operator()(const Point & p, const Real time, DenseVector<Number> & output) override
  {
    for (auto i : make_range(output.size()))
      output(i) = (*this)(p, time);
  }

private:
  const Real _scale;
};

/**
 * The processed constraint right-hand side of each constrained degree
 * of freedom this processor owns
 */
std::map<dof_id_type, Number> gather_local_rhs(DofMap & dof_map)
{
  std::map<dof_id_type, Number> rhs;
  const auto & values = dof_map.get_primal_constraint_values();
  for (const auto & [dof, row] : dof_map.get_dof_constraints())
    if (dof_map.local_index(dof))
      {
        libmesh_ignore(row);
        const auto it = values.find(dof);
        rhs[dof] = it == values.end() ? Number(0) : it->second;
      }
  return rhs;
}

/**
 * The processed constraint row of each constrained degree of freedom
 * this processor owns
 */
std::map<dof_id_type, DofConstraintRow> gather_local_rows(const DofMap & dof_map)
{
  std::map<dof_id_type, DofConstraintRow> rows;
  for (const auto & [dof, row] : dof_map.get_dof_constraints())
    if (dof_map.local_index(dof))
      rows[dof] = row;
  return rows;
}

/**
 * Checks that the processed constraints of \p dof_map match \p rows
 * and \p rhs exactly
 */
void assert_constraints_equal(const std::map<dof_id_type, DofConstraintRow> & rows,
                              const std::map<dof_id_type, Number> & rhs,
                              DofMap & dof_map)
{
  const auto other_rows = gather_local_rows(dof_map);
  const auto other_rhs = gather_local_rhs(dof_map);

  CPPUNIT_ASSERT_EQUAL(rows.size(), other_rows.size());
  for (const auto & [dof, row] : rows)
    {
      const auto it = other_rows.find(dof);
      CPPUNIT_ASSERT(it != other_rows.end());
      CPPUNIT_ASSERT_EQUAL(row.size(), it->second.size());
      for (const auto & [constraining_dof, coef] : row)
        {
          const auto jt = it->second.find(constraining_dof);
          CPPUNIT_ASSERT(jt != it->second.end());
          CPPUNIT_ASSERT_EQUAL(coef, jt->second);
        }
    }

  CPPUNIT_ASSERT_EQUAL(rhs.size(), other_rhs.size());
  for (const auto & [dof, value] : rhs)
    CPPUNIT_ASSERT_EQUAL(value, libmesh_map_find(other_rhs, dof));
}
}

class DirichletConstraintValuesTest : public CppUnit::TestCase
{
public:
  LIBMESH_CPPUNIT_TEST_SUITE(DirichletConstraintValuesTest);

#if LIBMESH_DIM > 1
#if defined(LIBMESH_ENABLE_DIRICHLET) && defined(LIBMESH_ENABLE_PERIODIC) && \
    defined(LIBMESH_ENABLE_AMR)
  CPPUNIT_TEST(testUpdateMatchesRebuild);
#endif
#endif
#if defined(LIBMESH_ENABLE_DIRICHLET) && defined(LIBMESH_HAVE_SOLVER)
  CPPUNIT_TEST(testSplineConstraintRows);
#endif

  CPPUNIT_TEST_SUITE_END();

private:
  void testUpdateMatchesRebuild()
  {
    LOG_UNIT_TEST;
    checkUpdate(4);
  }

  /**
   * Checks that System::reinit_dirichlet_constraint_values() gives
   * exactly the result of System::reinit_constraints() on a mesh with
   * spline constraint rows. The mesh is a line of 10 EDGE2 elements.
   * The constraint rows tie its nodes to 6 spline nodes. The Dirichlet
   * boundary is on a constrained node. Thus the processing of the
   * constraint rows projects the Dirichlet value onto the spline
   * nodes.
   */
  void testSplineConstraintRows()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
    EquationSystems es(mesh);
    auto & sys = es.add_system<LinearImplicitSystem>("sys");
    const auto u = sys.add_variable("u", FIRST);

    // The constraint operator below uses the node numbering of
    // build_line
    mesh.allow_renumbering(false);
    MeshTools::Generation::build_line(mesh, 10, 0., 1., EDGE2);

    // Each interior even node is a spline node. Each odd node takes
    // the mean of its two neighbors. The left node takes a combination
    // of a new spline node and the next spline node, and the right node
    // is a spline node. Processor 0 owns all rows of the operator.
    auto matrix = SparseMatrix<Number>::build(mesh.comm());
    const bool rank0 = !mesh.comm().rank();
    matrix->init(11, 6, rank0 ? 11 : 0, rank0 ? 6 : 0, rank0 ? 6 : 0, 0);
    if (rank0)
      {
        matrix->set(0, 0, 0.75);
        matrix->set(0, 1, 0.25);
        for (const auto i : make_range(5))
          {
            if (i)
              matrix->set(i*2, i, 1);
            matrix->set(i*2+1, i, 0.5);
            matrix->set(i*2+1, i+1, 0.5);
          }
        matrix->set(10, 5, 1);
      }
    matrix->close();
    mesh.copy_constraint_rows(*matrix);
    mesh.prepare_for_use();

    // build_line's left boundary id
    sys.get_dof_map().add_dirichlet_boundary(DirichletBoundary({0}, {u}, RampFunction(1.)));

    es.init();

    auto & dof_map = sys.get_dof_map();
    for (const Real time : {0., 0.5, 2.})
      {
        sys.time = time;
        sys.reinit_constraints();
        const auto rebuilt_rows = gather_local_rows(dof_map);
        const auto rebuilt_rhs = gather_local_rhs(dof_map);

        sys.reinit_dirichlet_constraint_values();
        assert_constraints_equal(rebuilt_rows, rebuilt_rhs, dof_map);
      }
  }

  /**
   * Checks that System::reinit_dirichlet_constraint_values() gives
   * exactly the result of System::reinit_constraints(). The mesh is an
   * n x n grid of QUAD9 elements with hanging nodes and a periodic
   * boundary. The periodic partners meet overlapping time-dependent
   * Dirichlet boundaries.
   */
  void checkUpdate(const unsigned int n)
  {
    Mesh mesh(*TestCommWorld);

    // Keep all neighbors until the periodic boundary is attached. The
    // periodic boundary needs them.
    mesh.allow_remote_element_removal(false);
    MeshTools::Generation::build_square(mesh, n, n, 0., 1., 0., 1., QUAD9);

    EquationSystems es(mesh);
    auto & sys = es.add_system<LinearImplicitSystem>("sys");
    const auto u = sys.add_variable("u", SECOND, HIERARCHIC);
    const auto v = sys.add_variable("v", SECOND, LAGRANGE);

    // The mesh is periodic in y. Thus periodic partners meet the
    // Dirichlet boundaries at the corners.
    PeriodicBoundary periodic(RealVectorValue(0., 1.));
    periodic.myboundary = bottom_id;
    periodic.pairedboundary = top_id;
    sys.get_dof_map().add_periodic_boundary(periodic);

    mesh.allow_remote_element_removal(true);
    mesh.prepare_for_use();

    // Add two overlapping time-dependent boundaries. The first claims u
    // and v on the left. The second claims u on the left and right.
    // The first boundary keeps the left side.
    sys.get_dof_map().add_dirichlet_boundary(
        DirichletBoundary({left_id}, {u, v}, RampFunction(1.)));
    sys.get_dof_map().add_dirichlet_boundary(
        DirichletBoundary({left_id, right_id}, {u}, RampFunction(-3.)));

    // Refine the elements near the left boundary. Then the hanging
    // nodes on the refined/unrefined interface reach the Dirichlet and
    // periodic boundaries.
    MeshRefinement refinement(mesh);
    for (auto & elem : mesh.active_element_ptr_range())
      {
        if (elem->vertex_average()(0) < 0.3)
          elem->set_refinement_flag(Elem::REFINE);
      }
    refinement.refine_elements();

    es.init();

    auto & dof_map = sys.get_dof_map();

    for (const Real time : {0., 0.25, 1., 2.5})
      {
        sys.time = time;

        sys.reinit_constraints();
        const auto rebuilt_rows = gather_local_rows(dof_map);
        const auto rebuilt_rhs = gather_local_rhs(dof_map);

        // The update must give the same result from any state that the
        // last processing left
        sys.reinit_dirichlet_constraint_values();
        assert_constraints_equal(rebuilt_rows, rebuilt_rhs, dof_map);
        const auto updated_rhs = gather_local_rhs(dof_map);

        // The values must change with time. Otherwise the comparison
        // above could pass on stale values.
        bool any_nonzero = false;
        for (const auto & [dof, value] : updated_rhs)
          {
            libmesh_ignore(dof);
            any_nonzero = any_nonzero || value != Number(0);
          }
        TestCommWorld->max(any_nonzero);
        CPPUNIT_ASSERT(any_nonzero);
      }

    // Two consecutive updates, with no full rebuild between them, must
    // also agree with a new full rebuild
    sys.time = 4.;
    sys.reinit_dirichlet_constraint_values();
    sys.reinit_dirichlet_constraint_values();
    const auto twice_updated_rows = gather_local_rows(dof_map);
    const auto twice_updated_rhs = gather_local_rhs(dof_map);
    sys.reinit_constraints();
    assert_constraints_equal(twice_updated_rows, twice_updated_rhs, dof_map);
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION(DirichletConstraintValuesTest);
