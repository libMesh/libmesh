#include <libmesh/default_coupling.h>
#include <libmesh/dirichlet_boundaries.h>
#include <libmesh/distributed_mesh.h>
#include <libmesh/dof_map.h>
#include <libmesh/elem.h>
#include <libmesh/equation_systems.h>
#include <libmesh/exodusII_io.h>
#include <libmesh/function_base.h>
#include <libmesh/ghost_point_neighbors.h>
#include <libmesh/linear_implicit_system.h>
#include <libmesh/partitioner.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_refinement.h>
#include <libmesh/numeric_vector.h>
#include <libmesh/periodic_boundary.h>
#include <libmesh/quadrature_gauss.h>
#include <libmesh/sparse_matrix.h>
#include <libmesh/wrapped_function.h>

#include "test_comm.h"
#include "libmesh_cppunit.h"


using namespace libMesh;



Number quadratic_solution (const Point& p,
                           const Parameters&,
                           const std::string&,
                           const std::string&)
{
  const Real & y = LIBMESH_DIM > 1 ? p(1) : 0;

  // Discontinuity at y=1 from the forcing function there.
  // Periodic on -3 < y < 2
  return (y > 1) ? (y-3)*(y-3) : (y+1)*(y+1);
}



struct PeriodicQuadFunction : public FunctionBase<Number>
{
  PeriodicQuadFunction() = default;

  virtual std::unique_ptr<FunctionBase<Number>> clone () const
  { return std::make_unique<PeriodicQuadFunction>(); }

  // We only really need the vector-valued output for projections
  virtual Number operator() (const Point &,
                             const Real /*time*/ = 0.) override
  { libmesh_error(); }

  virtual void operator() (const Point & p,
                           const Real,
                           DenseVector<Number> & output) override
  {
    libmesh_assert_equal_to(output.size(), 1);
    Parameters params;
    output(0) = quadratic_solution(p, params, "", "");
  }

  Number component (unsigned int i,
                    const Point & p,
                    Real /* time */) override
  {
    Parameters params;
    switch (i) {
    case 0:
      return quadratic_solution(p, params, "", "");
    default:
      libmesh_error();
    }
    return 0;
  }
};


void periodic_bc_test_poisson(EquationSystems& es,
                              const std::string&)
{
  const MeshBase& mesh = es.get_mesh();
  LinearImplicitSystem& system = es.get_system<LinearImplicitSystem>("PBCSys");
  const DofMap& dof_map = system.get_dof_map();

  // Interior FE to integrate smooth part of poisson problem
  FEType fe_type = dof_map.variable_type(0);
  std::unique_ptr<FEBase> fe (FEBase::build(2, fe_type));
  QGauss qrule (2, fe_type.default_quadrature_order());
  fe->attach_quadrature_rule(&qrule);

  const std::vector<Real> & JxW = fe->get_JxW();
  const std::vector<std::vector<Real>> & phi = fe->get_phi();
  const std::vector<std::vector<RealGradient>> & dphi = fe->get_dphi();

  // Face FE to integrate delta function part of forcing
  std::unique_ptr<FEBase> fe_face (FEBase::build(2, fe_type));
  QGauss qface (1, fe_type.default_quadrature_order());
  fe_face->attach_quadrature_rule(&qface);

  const std::vector<Real> & JxW_face = fe_face->get_JxW();
  const std::vector<std::vector<Real>> & phi_face = fe_face->get_phi();
  std::unique_ptr<const Elem> elem_side;

  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;

  std::vector<dof_id_type> dof_indices;

  for (const Elem * elem : mesh.active_local_element_ptr_range())
    {
      dof_map.dof_indices (elem, dof_indices);
      const unsigned int n_dofs = dof_indices.size();

      Ke.resize (n_dofs, n_dofs);
      Fe.resize (n_dofs);

      fe->reinit (elem);

      // Integrate the poisson problem over the element interior, where we have
      // a nice f=2 forcing function aside from the discontinuity
      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {
          for (unsigned int i=0; i != n_dofs; i++)
            {
              for (unsigned int j=0; j != n_dofs; j++)
                Ke(i,j) += -JxW[qp]*(dphi[i][qp]*dphi[j][qp]);

              Fe(i) += JxW[qp]*phi[i][qp]*2;
            }
        }

      for(unsigned int s=0; s != elem->n_sides(); ++s)
      {
        elem->build_side_ptr(elem_side, s);
        Point centroid = elem_side->vertex_average();
        if (std::abs(centroid(1) - 1) < TOLERANCE)
          {
            fe_face->reinit(elem, s);

            for (unsigned int qp=0; qp<qface.n_points(); qp++)
              {
                // delta function divided by 2 since we hit it from
                // both sides
                for (unsigned int i=0; i != n_dofs; i++)
                  Fe(i) += JxW_face[qp]*phi_face[i][qp]*(-4);
              }
          }
      }

      dof_map.heterogenously_constrain_element_matrix_and_vector(Ke, Fe, dof_indices);
      system.matrix->add_matrix (Ke, dof_indices);
      system.rhs->add_vector    (Fe, dof_indices);
    }

  system.rhs->close();
  system.matrix->close();
}


// Assigns the element rows of an 8-row unit square to processors 0
// and 1.  At first processor 1 owns only row 2.  With two layers of
// ghosting it then holds the bottom row 0 but not the top row 7, which
// is the periodic image of row 0.  The repartitioned assignment also
// gives row 0 to processor 1.
class RowPartitioner : public Partitioner
{
public:
  RowPartitioner (bool repartitioned) : _repartitioned(repartitioned) {}

  std::unique_ptr<Partitioner> clone () const override
  { return std::make_unique<RowPartitioner>(*this); }

protected:
  void _do_partition (MeshBase & mesh, const unsigned int) override
  {
    for (auto & elem : mesh.active_element_ptr_range())
      {
        const auto row = static_cast<unsigned int>(elem->vertex_average()(1) * 8);
        const bool on_proc_1 = _repartitioned ? (row == 0 || row == 2) : (row == 2);
        elem->processor_id() = on_proc_1;
      }
  }

private:
  bool _repartitioned;
};


const boundary_id_type top_id    = 50;
const boundary_id_type bottom_id = 60;
const boundary_id_type side_id   = 70;


class PeriodicBCTest : public CppUnit::TestCase {
public:
  LIBMESH_CPPUNIT_TEST_SUITE( PeriodicBCTest );

#if LIBMESH_DIM > 1
#if defined(LIBMESH_HAVE_SOLVER) && defined(LIBMESH_HAVE_EXODUS_API) && defined(LIBMESH_HAVE_GZSTREAM)
  CPPUNIT_TEST( testPeriodicLagrange2 );
#endif
  CPPUNIT_TEST( testRedistributeToGhosts );
#ifdef LIBMESH_ENABLE_AMR
  CPPUNIT_TEST( testDistributedRefinedConstraints );
  CPPUNIT_TEST( testSecondSystemConstraints );
  CPPUNIT_TEST( testReplacedDefaultGhosting );
  CPPUNIT_TEST( testRefinementAcrossSets );
#endif
#endif // LIBMESH_DIM > 1

  CPPUNIT_TEST_SUITE_END();

private:

  void testPeriodicBC (FEType fe_type)
  {
    Mesh mesh(*TestCommWorld);
    BoundaryInfo & boundary = mesh.get_boundary_info();

    EquationSystems es(mesh);
    LinearImplicitSystem &sys =
      es.add_system<LinearImplicitSystem> ("PBCSys");

    // We need to change around BCIDs and set up the periodic BC first
    mesh.allow_remote_element_removal(false);

    mesh.read("meshes/shark_tooth_tri6.xda.gz");

    // Add some BCIDs on top and bottom for our periodic boundary
    for (const Elem * elem : mesh.active_element_ptr_range())
      {
        for (unsigned int s=0; s != 3; ++s)
          if (!elem->neighbor_ptr(s))
            {
              unsigned v2 = (s+1)%3;
              // Left-running edges are on bottom, because this xda
              // has a flipped mesh
              if (elem->point(s)(0) - elem->point(v2)(0) < -.1)
                boundary.add_side(elem, s, top_id);

              // Right-running edges are on top
              else if (elem->point(s)(0) - elem->point(v2)(0) > .1)
                boundary.add_side(elem, s, bottom_id);

              // Vertical edges are Dirichlet
              else
                boundary.add_side(elem, s, side_id);
            }
      }

    PeriodicBoundary vert(RealVectorValue(0., 4.0));
    vert.myboundary = bottom_id;
    vert.pairedboundary = top_id;
    sys.get_dof_map().add_periodic_boundary(vert);

    // Okay, *now* the mesh knows to save periodic neighbors
    mesh.allow_remote_element_removal(true);
    mesh.prepare_for_use();

    const unsigned int u_var = sys.add_variable("u", fe_type);
    sys.attach_assemble_function (periodic_bc_test_poisson);

    std::set<boundary_id_type> side_bdy { side_id };
    std::vector<unsigned int> all_vars (1, u_var);
    WrappedFunction<Number> exact_val(sys, quadratic_solution);
    DirichletBoundary exact_bc(side_bdy, all_vars, exact_val);
    sys.get_dof_map().add_dirichlet_boundary(exact_bc);

    es.init();

    sys.solve();

    ExodusII_IO(mesh).write_equation_systems("periodic.e", es);

    Parameters params;

    for (Real x = 0.1; x < 10; x += 0.2)
      for (Real ya = -2.9; ya < 1; ya += 0.2)
        {
          // Integrate the sharktooth shape
          const Real offset = 1-std::abs(std::fmod(x,2)-1);
          const Real y = ya + offset;
          Point p{x,y};
          const Number exact = quadratic_solution(p,params,"","");
          const Number approx = sys.point_value(0,p);
          LIBMESH_ASSERT_NUMBERS_EQUAL(exact, approx,
                                      TOLERANCE*TOLERANCE*20);
        }
  }


  // Repartition a distributed periodic mesh so that a processor takes
  // ownership of an element it only ghosted, before that element's
  // periodic neighbor arrives.
  void testRedistributeToGhosts()
  {
    LOG_UNIT_TEST;

    // Repartitioning moves elements only between processors 0 and 1
    if (TestCommWorld->size() < 2)
      return;

    DistributedMesh mesh(*TestCommWorld);
    mesh.partitioner() = std::make_unique<RowPartitioner>(false);

    EquationSystems es(mesh);
    System & sys = es.add_system<System>("PBCSys");
    sys.add_variable("u", FEType(FIRST, LAGRANGE));

    // build_square() numbers the bottom and top sides 0 and 2
    PeriodicBoundary vert(RealVectorValue(0., 1.));
    vert.myboundary = 0;
    vert.pairedboundary = 2;
    sys.get_dof_map().add_periodic_boundary(vert);

    // When processor 1 finds what processor 0 needs, the first layer
    // reaches from row 1 into row 0.  The second layer then asks for
    // the periodic neighbor of row 0, which processor 1 does not have.
    auto two_layers = std::make_shared<DefaultCoupling>();
    two_layers->set_mesh(&mesh);
    two_layers->set_n_levels(2);
    two_layers->set_periodic_boundaries(sys.get_dof_map().get_periodic_boundaries());
    mesh.add_ghosting_functor(two_layers);

    MeshTools::Generation::build_square(mesh, 4, 8, 0., 1., 0., 1., QUAD4);
    es.init();

    mesh.partitioner() = std::make_unique<RowPartitioner>(true);
    mesh.partition();

    CPPUNIT_ASSERT(!mesh.is_redistributing());
    for (const auto & elem : mesh.active_local_element_ptr_range())
      {
        const auto row = static_cast<unsigned int>(elem->vertex_average()(1) * 8);
        CPPUNIT_ASSERT_EQUAL(processor_id_type(row == 0 || row == 2), elem->processor_id());
      }
  }

  void testPeriodicLagrange2() { LOG_UNIT_TEST; testPeriodicBC(FEType(SECOND, LAGRANGE)); }

#ifdef LIBMESH_ENABLE_AMR
  /**
   * On a refined distributed mesh, each processor chooses which side of a periodic vertex to
   * constrain from the elements at the vertex and at its image. Missing ghost elements give
   * processors different choices, and both vertices are then constrained in terms of each other.
   *
   * \p n_systems_before systems with no periodic boundaries precede the periodic system. If
   * \p replace_default_ghosting is true, a GhostPointNeighbors that uses every periodic boundary
   * set on the mesh replaces the default ghosting.
   */
  void checkDistributedRefinedConstraints(const unsigned int n_systems_before,
                                          const bool replace_default_ghosting)
  {
    DistributedMesh mesh(*TestCommWorld);
    mesh.allow_remote_element_removal(false);
    MeshTools::Generation::build_square(mesh, 4, 4, 0., 1., 0., 1., QUAD9);

    EquationSystems es(mesh);
    for (const auto i : make_range(n_systems_before))
      es.add_system<LinearImplicitSystem>("Sys" + std::to_string(i)).add_variable("v", FIRST, LAGRANGE);

    auto & sys = es.add_system<LinearImplicitSystem>("PBCSys");
    sys.add_variable("u", SECOND, LAGRANGE);

    // build_square() puts boundary 0 at the bottom and boundary 2 at the top
    PeriodicBoundary vert(RealVectorValue(0., 1.));
    vert.myboundary = 0;
    vert.pairedboundary = 2;
    sys.get_dof_map().add_periodic_boundary(vert);

    if (replace_default_ghosting)
      {
        es.enable_default_ghosting(false);
        mesh.add_ghosting_functor(std::make_shared<GhostPointNeighbors>(mesh),
                                  /*use_all_pbcs=*/true);
      }

    mesh.allow_remote_element_removal(true);
    mesh.delete_remote_elements();

    // Refine one element, selected by position so that the mesh is independent of the
    // partitioning
    for (auto & elem : mesh.active_element_ptr_range())
      if (elem->vertex_average().absolute_fuzzy_equals(Point(0.375, 0.375)))
        elem->set_refinement_flag(Elem::REFINE);
    MeshRefinement(mesh).refine_elements();

    es.init();

    // Constraint rows are fully resolved, so no constrained dof appears in any row
    for (const auto s : make_range(es.n_systems()))
      {
        const DofMap & dof_map = es.get_system(s).get_dof_map();
        for (const auto & [dof, row] : dof_map.get_dof_constraints())
          for (const auto & [row_dof, coef] : row)
            {
              libmesh_ignore(dof, coef);
              CPPUNIT_ASSERT(!dof_map.is_constrained_dof(row_dof));
            }
      }
  }

  void testDistributedRefinedConstraints()
  {
    LOG_UNIT_TEST;
    checkDistributedRefinedConstraints(0, false);
  }

  // The default ghosting must use the periodic boundaries of a system other than the first
  void testSecondSystemConstraints()
  {
    LOG_UNIT_TEST;
    checkDistributedRefinedConstraints(1, false);
  }

  // A replacement for the default ghosting sees periodic boundaries added before it
  void testReplacedDefaultGhosting()
  {
    LOG_UNIT_TEST;
    checkDistributedRefinedConstraints(0, true);
  }

  // MeshRefinement keeps level one across the periodic boundaries of every system
  void testRefinementAcrossSets()
  {
    LOG_UNIT_TEST;

    DistributedMesh mesh(*TestCommWorld);
    mesh.allow_remote_element_removal(false);
    MeshTools::Generation::build_square(mesh, 4, 4, 0., 1., 0., 1., QUAD4);

    // build_square() numbers the sides 0 (bottom), 1 (right), 2 (top) and 3 (left)
    EquationSystems es(mesh);
    auto & vert_sys = es.add_system<System>("VertSys");
    vert_sys.add_variable("u", FIRST, LAGRANGE);
    PeriodicBoundary vert(RealVectorValue(0., 1.));
    vert.myboundary = 0;
    vert.pairedboundary = 2;
    vert_sys.get_dof_map().add_periodic_boundary(vert);

    auto & horiz_sys = es.add_system<System>("HorizSys");
    horiz_sys.add_variable("v", FIRST, LAGRANGE);
    PeriodicBoundary horiz(RealVectorValue(1., 0.));
    horiz.myboundary = 3;
    horiz.pairedboundary = 1;
    horiz_sys.get_dof_map().add_periodic_boundary(horiz);

    mesh.allow_remote_element_removal(true);
    mesh.delete_remote_elements();

    // Refine an element on the left boundary, then refine its child on that boundary. The
    // periodic image of the child's side is the unrefined element on the right boundary, which
    // the level-one rule must then refine too.
    MeshRefinement mesh_refinement(mesh);
    for (const Point & p : {Point(0.125, 0.375), Point(0.0625, 0.4375)})
      {
        for (auto & elem : mesh.active_element_ptr_range())
          if (elem->vertex_average().absolute_fuzzy_equals(p))
            elem->set_refinement_flag(Elem::REFINE);
        mesh_refinement.refine_elements();
      }

    CPPUNIT_ASSERT(mesh_refinement.test_level_one());

    bool image_refined = false;
    for (const auto & elem : mesh.local_level_element_ptr_range(0))
      if (elem->vertex_average().absolute_fuzzy_equals(Point(0.875, 0.375)))
        image_refined = elem->has_children();
    mesh.comm().max(image_refined);
    CPPUNIT_ASSERT(image_refined);
  }
#endif
};

CPPUNIT_TEST_SUITE_REGISTRATION( PeriodicBCTest );
