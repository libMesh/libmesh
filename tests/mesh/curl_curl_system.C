#include "curl_curl_system.h"

#include "libmesh/dof_map.h"
#include "libmesh/fe_base.h"
#include "libmesh/fe_interface.h"
#include "libmesh/fem_context.h"
#include "libmesh/mesh.h"
#include "libmesh/quadrature.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/zero_function.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

CurlCurlSystem::CurlCurlSystem(EquationSystems & es,
                               const std::string & name_in,
                               const unsigned int number_in) :
  FEMSystem(es, name_in, number_in)
{}

void CurlCurlSystem::init_data ()
{
  // Add the solution variable
  u_var = this->add_variable ("u", FIRST, NEDELEC_ONE);
  v_var = this->add_variable ("v", FIRST, LAGRANGE);

  this->time_evolving(u_var);
  this->time_evolving(u_var);

  this->init_dirichlet_bcs();

  // Do the parent's initialization after variables and boundary constraints are defined
  FEMSystem::init_data();
}

void CurlCurlSystem::init_dirichlet_bcs()
{
  const boundary_id_type all_ids[4] = {0, 1, 2, 3};
  std::set<boundary_id_type> boundary_ids(all_ids, all_ids+4);

  std::vector<unsigned int> vars;
  vars.push_back(v_var);

  ZeroFunction<Number> func;

  this->get_dof_map().add_dirichlet_boundary(libMesh::DirichletBoundary(boundary_ids, vars, &func));
}

void CurlCurlSystem::init_context(DiffContext & context)
{
  FEMContext & c = cast_ref<FEMContext &>(context);

  // Get finite element object
  FEGenericBase<RealGradient> * fe;
  c.get_element_fe<RealGradient>(u_var, fe);

  // We should prerequest all the data
  // we will need to build the linear system.
  fe->get_JxW();
  fe->get_phi();
  fe->get_curl_phi();
  fe->get_xyz();

  // Get finite element object
  FEGenericBase<RealGradient> * side_fe;
  c.get_side_fe<RealGradient>(u_var, side_fe);

  side_fe->get_JxW();
  side_fe->get_phi();
  side_fe->get_xyz();

  // Now for scalar
  FEGenericBase<Real> * fes;
  c.get_element_fe<Real>(v_var, fes);

  // We should prerequest all the data
  // we will need to build the linear system.
  fes->get_JxW();
  fes->get_phi();
  fes->get_dphi();
  fes->get_xyz();
}


bool CurlCurlSystem::element_time_derivative (bool request_jacobian,
                                              DiffContext & context)
{
  FEMContext & c = cast_ref<FEMContext &>(context);

  // Get finite element objects
  FEGenericBase<RealGradient> * feu = libmesh_nullptr;
  c.get_element_fe<RealGradient>(u_var, feu);
  FEGenericBase<Real> * fev = libmesh_nullptr;
  c.get_element_fe<Real>(v_var, fev);

  // First we get some references to cell-specific data that
  // will be used to assemble the linear system.

  // Element Jacobian * quadrature weights for interior integration
  const std::vector<Real> & JxW = feu->get_JxW();

  const std::vector<std::vector<RealGradient> > & phiu = feu->get_phi();
  const std::vector<std::vector<RealGradient> > & curl_phiu = feu->get_curl_phi();
  const std::vector<std::vector<Real> > & phiv = fev->get_phi();
  const std::vector<std::vector<RealGradient> > & grad_phiv = fev->get_dphi();

  const std::vector<Point> & qpoint = feu->get_xyz();

  // The number of local degrees of freedom in each variable
  const unsigned int n_u_dofs = c.get_dof_indices(u_var).size();
  const unsigned int n_v_dofs = c.get_dof_indices(v_var).size();

  DenseSubMatrix<Number> & Kuu = c.get_elem_jacobian(u_var, u_var);
  DenseSubMatrix<Number> & Kvv = c.get_elem_jacobian(v_var, v_var);

  DenseSubVector<Number> & Fu = c.get_elem_residual(u_var);
  DenseSubVector<Number> & Fv = c.get_elem_residual(v_var);

  // Now we will build the element Jacobian and residual.
  // Constructing the residual requires the solution and its
  // gradient from the previous timestep.  This must be
  // calculated at each quadrature point by summing the
  // solution degree-of-freedom values by the appropriate
  // weight functions.
  const unsigned int n_qpoints = c.get_element_qrule().n_points();

  // Loop over quadrature points
  for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
      Gradient u;
      Gradient curl_u;
      Number v;
      Gradient grad_v;

      c.interior_value(u_var, qp, u);
      c.interior_curl(u_var, qp, curl_u);
      c.interior_value(v_var, qp, v);
      c.interior_gradient(v_var, qp, grad_v);

      // Value of the forcing function at this quadrature point
      RealGradient f = this->forcing(qpoint[qp]);

      // First, an i-loop over the degrees of freedom.
      for (unsigned int i=0; i != n_u_dofs; i++)
        {
          Fu(i) += (curl_u*curl_phiu[i][qp] + u*phiu[i][qp] - f*phiu[i][qp])*JxW[qp];

          // Matrix contributions for the uu and vv couplings.
          if (request_jacobian)
            for (unsigned int j=0; j != n_u_dofs; j++)
              Kuu(i,j) += (curl_phiu[j][qp]*curl_phiu[i][qp] +
                           phiu[j][qp]*phiu[i][qp])*JxW[qp];
        }
      for (unsigned int i=0; i != n_v_dofs; i++)
      {
        Fv(i) += (grad_v*grad_phiv[i][qp] - phiv[i][qp]) * JxW[qp];

        if (request_jacobian)
          for (unsigned j=0; j != n_v_dofs; j++)
            Kvv(i, j) += (grad_phiv[j][qp]*grad_phiv[i][qp]) * JxW[qp];
      }
    } // end of the quadrature point qp-loop

  return request_jacobian;
}


bool CurlCurlSystem::side_time_derivative (bool request_jacobian,
                                           DiffContext & context)
{
  FEMContext & c = cast_ref<FEMContext &>(context);

  // Get finite element object
  FEGenericBase<RealGradient> * side_fe = libmesh_nullptr;
  c.get_side_fe<RealGradient>(u_var, side_fe);

  // First we get some references to cell-specific data that
  // will be used to assemble the linear system.

  // Element Jacobian * quadrature weights for interior integration
  const std::vector<Real> & JxW = side_fe->get_JxW();

  // The velocity shape functions at interior quadrature points.
  const std::vector<std::vector<RealGradient> > & phi = side_fe->get_phi();

  // The number of local degrees of freedom in each variable
  const unsigned int n_u_dofs = c.get_dof_indices(u_var).size();

  const std::vector<Point> & normals = side_fe->get_normals();

  const std::vector<Point> & qpoint = side_fe->get_xyz();

  // The penalty value.  \frac{1}{\epsilon}
  // in the discussion above.
  const Real penalty = 1.e10;

  DenseSubMatrix<Number> & Kuu = c.get_elem_jacobian(u_var, u_var);
  DenseSubVector<Number> & Fu = c.get_elem_residual(u_var);

  const unsigned int n_qpoints = c.get_side_qrule().n_points();

  for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
      Gradient u;
      c.side_value(u_var, qp, u);

      RealGradient N(normals[qp](0), normals[qp](1));

      Gradient u_exact = this->exact_solution(qpoint[qp]);

      Gradient Ncu = (u - u_exact).cross(N);

      for (unsigned int i=0; i != n_u_dofs; i++)
        {
          Fu(i) += penalty*Ncu*(phi[i][qp].cross(N))*JxW[qp];

          if (request_jacobian)
            for (unsigned int j=0; j != n_u_dofs; j++)
              Kuu(i,j) += penalty*(phi[j][qp].cross(N))*(phi[i][qp].cross(N))*JxW[qp];
        }
    }

  return request_jacobian;
}

RealGradient CurlCurlSystem::forcing(const Point & p)
{
  Real x = p(0); Real y = p(1);

  return soln.forcing(x, y);
}

RealGradient CurlCurlSystem::exact_solution(const Point & p)
{
  Real x = p(0); Real y = p(1);

  return soln(x, y);
}
