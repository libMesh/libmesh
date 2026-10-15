#include <libmesh/dirichlet_boundaries.h>
#include <libmesh/dof_map.h>
#include <libmesh/elem.h>
#include <libmesh/equation_systems.h>
#include <libmesh/function_base.h>
#include <libmesh/libmesh.h>
#include <libmesh/linear_implicit_system.h>
#include <libmesh/distributed_mesh.h>
#include <libmesh/mesh.h>
#include <libmesh/replicated_mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_refinement.h>
#include <libmesh/parallel.h>
#include <libmesh/periodic_boundary.h>

#include "test_comm.h"
#include "libmesh_cppunit.h"

#include <array>
#include <chrono>

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
}

class DirichletConstraintValuesTest : public CppUnit::TestCase
{
public:
  LIBMESH_CPPUNIT_TEST_SUITE(DirichletConstraintValuesTest);

#if LIBMESH_DIM > 1
#if defined(LIBMESH_ENABLE_DIRICHLET) && defined(LIBMESH_ENABLE_PERIODIC) && \
    defined(LIBMESH_ENABLE_AMR)
  CPPUNIT_TEST(testUpdateMatchesRebuild);
  CPPUNIT_TEST(testUpdateCost);
#endif
#endif

  CPPUNIT_TEST_SUITE_END();

private:
  void testUpdateMatchesRebuild()
  {
    LOG_UNIT_TEST;
    checkUpdate(4, /*report=*/false);
  }

  /**
   * Compares the cost of the update with the cost of a full rebuild,
   * phase by phase. The default grid is 64 x 64. That size gives
   * meaningful timings and keeps the suite fast. Options:
   * - --dirichlet-cost-n <n> sets the grid size, for a scaling study.
   * - --dirichlet-cost-interior-refinement refines a checkerboard over
   *   the whole domain. Then the number of hanging nodes grows with
   *   the volume, not with the boundary.
   * - --dirichlet-cost-distributed uses a DistributedMesh.
   */
  void testUpdateCost()
  {
    LOG_UNIT_TEST;
    checkUpdate(libMesh::command_line_next("--dirichlet-cost-n", 64u),
                /*report=*/true,
                libMesh::on_command_line("--dirichlet-cost-interior-refinement"),
                libMesh::on_command_line("--dirichlet-cost-distributed"));
  }

  /**
   * Checks that System::reinit_dirichlet_constraint_values() gives
   * exactly the result of System::reinit_constraints(). The mesh is an
   * n x n grid of QUAD9 elements with hanging nodes and a periodic
   * boundary. The periodic partners meet overlapping time-dependent
   * Dirichlet boundaries. If \p interior_refinement is true, the
   * hanging nodes cover the whole domain, not only the region near the
   * boundary. If \p distributed is true, the mesh is a DistributedMesh.
   */
  void checkUpdate(const unsigned int n,
                   const bool report,
                   const bool interior_refinement = false,
                   const bool distributed = false)
  {
    std::unique_ptr<UnstructuredMesh> mesh_ptr;
    if (distributed)
      mesh_ptr = std::make_unique<DistributedMesh>(*TestCommWorld);
    else
      mesh_ptr = std::make_unique<ReplicatedMesh>(*TestCommWorld);
    UnstructuredMesh & mesh = *mesh_ptr;

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

    // Refine the elements near the left boundary, or a checkerboard
    // over the whole domain. Then the hanging nodes on each
    // refined/unrefined interface reach the Dirichlet and periodic
    // boundaries.
    MeshRefinement refinement(mesh);
    for (auto & elem : mesh.active_element_ptr_range())
      {
        const Point c = elem->vertex_average();
        const bool checker = (unsigned(c(0) * n) + unsigned(c(1) * n)) % 2;
        if (interior_refinement ? checker : c(0) < 0.3)
          elem->set_refinement_flag(Elem::REFINE);
      }
    refinement.refine_elements();

    es.init();

    auto & dof_map = sys.get_dof_map();

    using clock = std::chrono::steady_clock;
    double rebuild_seconds = 0, update_seconds = 0;

    for (const Real time : {0., 0.25, 1., 2.5})
      {
        sys.time = time;

        auto start = clock::now();
        sys.reinit_constraints();
        rebuild_seconds += std::chrono::duration<double>(clock::now() - start).count();
        const auto rebuilt_rows = gather_local_rows(dof_map);
        const auto rebuilt_rhs = gather_local_rhs(dof_map);

        // The update must give the same result from any state that the
        // last processing left
        start = clock::now();
        sys.reinit_dirichlet_constraint_values();
        update_seconds += std::chrono::duration<double>(clock::now() - start).count();
        const auto updated_rows = gather_local_rows(dof_map);
        const auto updated_rhs = gather_local_rhs(dof_map);

        CPPUNIT_ASSERT_EQUAL(rebuilt_rows.size(), updated_rows.size());
        for (const auto & [dof, row] : rebuilt_rows)
          {
            const auto it = updated_rows.find(dof);
            CPPUNIT_ASSERT(it != updated_rows.end());
            CPPUNIT_ASSERT_EQUAL(row.size(), it->second.size());
            for (const auto & [constraining_dof, coef] : row)
              {
                const auto jt = it->second.find(constraining_dof);
                CPPUNIT_ASSERT(jt != it->second.end());
                CPPUNIT_ASSERT_EQUAL(coef, jt->second);
              }
          }

        CPPUNIT_ASSERT_EQUAL(rebuilt_rhs.size(), updated_rhs.size());
        for (const auto & [dof, value] : rebuilt_rhs)
          CPPUNIT_ASSERT_EQUAL(value, libmesh_map_find(updated_rhs, dof));

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
    const auto twice_updated_rhs = gather_local_rhs(dof_map);
    sys.reinit_constraints();
    const auto rebuilt_rhs = gather_local_rhs(dof_map);
    CPPUNIT_ASSERT_EQUAL(rebuilt_rhs.size(), twice_updated_rhs.size());
    for (const auto & [dof, value] : rebuilt_rhs)
      CPPUNIT_ASSERT_EQUAL(value, libmesh_map_find(twice_updated_rhs, dof));

    if (!report)
      return;

    // Time each phase of both paths separately. The update replaces
    // only the first phase. Thus the fraction of the update time in the
    // shared phases is the time that the update could also save if it
    // skipped process_constraints().
    const auto time_phases = [&](const bool update)
    {
      std::array<double, 3> seconds = {0, 0, 0};
      for (const Real time : {5., 6., 7., 8.})
        {
          sys.time = time;
          auto start = clock::now();
          if (update)
            dof_map.reinit_dirichlet_constraint_values(mesh, time);
          else
            dof_map.create_dof_constraints(mesh, time);
          sys.user_constrain();
          seconds[0] += std::chrono::duration<double>(clock::now() - start).count();

          start = clock::now();
          dof_map.process_constraints(mesh);
          seconds[1] += std::chrono::duration<double>(clock::now() - start).count();

          start = clock::now();
          dof_map.prepare_send_list();
          seconds[2] += std::chrono::duration<double>(clock::now() - start).count();
        }
      for (auto & s : seconds)
        TestCommWorld->max(s);
      return seconds;
    };

    const auto rebuild_phases = time_phases(false);
    const auto update_phases = time_phases(true);

    auto [cache_bytes, constraint_bytes] = dof_map.dirichlet_constraint_cache_memory();
    TestCommWorld->sum(cache_bytes);
    TestCommWorld->sum(constraint_bytes);
    TestCommWorld->max(rebuild_seconds);
    TestCommWorld->max(update_seconds);
    libMesh::out << "\nDirichletConstraintValuesTest: " << mesh.n_active_elem()
                 << " active elements, " << mesh.n_nodes() << " nodes, " << dof_map.n_dofs()
                 << " dofs, " << dof_map.n_constrained_dofs() << " constrained, "
                 << TestCommWorld->size() << " processors"
                 << "\n  time over 4 calls: reinit_constraints " << rebuild_seconds
                 << " s, reinit_dirichlet_constraint_values " << update_seconds << " s"
                 << "\n  phases over 4 calls (build raw | process_constraints | prepare_send_list):"
                 << "\n    rebuild " << rebuild_phases[0] << " | " << rebuild_phases[1] << " | "
                 << rebuild_phases[2] << " s"
                 << "\n    update  " << update_phases[0] << " | " << update_phases[1] << " | "
                 << update_phases[2] << " s"
                 << "\n  memory summed over processors: cache " << cache_bytes
                 << " B, processed dof constraints " << constraint_bytes << " B" << std::endl;
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION(DirichletConstraintValuesTest);
