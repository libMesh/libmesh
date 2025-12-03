#include <libmesh/dirichlet_boundaries.h>
#include <libmesh/dof_map.h>
#include <libmesh/elem.h>
#include <libmesh/equation_systems.h>
#include <libmesh/exodusII_io.h>
#include <libmesh/function_base.h>
#include <libmesh/linear_implicit_system.h>
#include <libmesh/mesh.h>
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


  void testPeriodicLagrange2() { LOG_UNIT_TEST; testPeriodicBC(FEType(SECOND, LAGRANGE)); }
};

CPPUNIT_TEST_SUITE_REGISTRATION( PeriodicBCTest );
