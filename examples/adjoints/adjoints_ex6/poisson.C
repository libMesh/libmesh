/* This is where we define the assembly of the Laplace system */

// General libMesh includes
#include "libmesh/getpot.h"
#include "libmesh/boundary_info.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/dof_map.h"
#include "libmesh/fe_base.h"
#include "libmesh/fe_interface.h"
#include "libmesh/fem_context.h"
#include "libmesh/fem_system.h"
#include "libmesh/quadrature.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/mesh.h"
#include "libmesh/parallel.h"
#include "libmesh/quadrature.h"
#include "libmesh/fem_context.h"
#include "libmesh/zero_function.h"

// Local includes
#include "poisson.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

// Function to set the Dirichlet boundary function for the adjoint dirichlet boundary
class BdyFunction : public FunctionBase<Number>
{
public:
  BdyFunction (unsigned int T_var)
    : _T_var(T_var)
  { this->_initialized = true; }

  virtual Number operator() (const Point&, const Real = 0)
  { libmesh_not_implemented(); }

  virtual void operator() (const Point& p,
                           const Real,
                           DenseVector<Number>& output)
  {
    output.resize(2);
    output.zero();
    const Real x=p(0);
    // Set the parabolic weighting
    output(_T_var) = -x * (1 - x);
  }

  virtual UniquePtr<FunctionBase<Number> > clone() const
  { return UniquePtr<FunctionBase<Number> > (new BdyFunction(_T_var)); }

private:
  const unsigned int _T_var;
};

void PoissonSystem::init_data ()
{
  this->add_variable ("T", static_cast<Order>(_fe_order),
                      Utility::string_to_enum<FEFamily>(_fe_family));

  GetPot infile("poisson.in");
  exact_QoI[0] = infile("QoI_0", 0.0);
  alpha = infile("alpha", 100.0);

  // Now we will set the Dirichlet boundary conditions

  // Get the variable number of the variable for which we will set
  // the Dirichlet boundary conditions
  T_var = this->variable_number ("T");
  this->time_evolving(T_var);

  // Get boundary ids for all the boundaries
  const boundary_id_type all_bdry_id[4] = {0, 1, 2, 3};
  std::set<boundary_id_type> all_bdy(all_bdry_id, all_bdry_id+4);

  // For the adjoint problem, we will only set the bottom boundary to a non-zero value
  const boundary_id_type bottom_bdry_id = 0;
  std::set<boundary_id_type> bottom_bdry;
  bottom_bdry.insert(bottom_bdry_id);

  // The T_only identifier for setting the boundary conditions for T
  std::vector<unsigned int> T_only(1, T_var);

  // The zero function pointer for the primal all bdry bcs
  ZeroFunction<Number> zero;
  // Boundary function for bottom bdry adjoint condition
  BdyFunction bottom_adjoint(T_var);

  this->get_dof_map().add_dirichlet_boundary(DirichletBoundary (all_bdy, T_only, &zero));

  this->get_dof_map().add_adjoint_dirichlet_boundary(DirichletBoundary (bottom_bdry, T_only, &bottom_adjoint), 0);

  // Do the parent's initialization after variables are defined
  FEMSystem::init_data();
}

void PoissonSystem::init_context(DiffContext &context)
{
  FEMContext &c = cast_ref<FEMContext&>(context);

  // Now make sure we have requested all the data
  // we need to build the linear system.
  FEBase* elem_fe = NULL;
  c.get_element_fe( 0, elem_fe );
  elem_fe->get_JxW();
  elem_fe->get_phi();
  elem_fe->get_dphi();
  elem_fe->get_xyz();

  FEBase* side_fe = NULL;
  c.get_side_fe( 0, side_fe );

  side_fe->get_JxW();
  side_fe->get_phi();
  side_fe->get_dphi();
  side_fe->get_xyz();

}

#define optassert(X) {if (!(X)) libmesh_error_msg("Assertion " #X " failed.");}

// Assemble the element contributions to the stiffness matrix
bool PoissonSystem::element_time_derivative (bool request_jacobian,
                                             DiffContext &context)
{
  // Are the jacobians specified analytically ?
  bool compute_jacobian = request_jacobian && _analytic_jacobians;

  FEMContext &c = cast_ref<FEMContext&>(context);

  // First we get some references to cell-specific data that
  // will be used to assemble the linear system.
  FEBase* elem_fe = NULL;
  c.get_element_fe( 0, elem_fe );

  // Element Jacobian * quadrature weights for interior integration
  const std::vector<Real> &JxW = elem_fe->get_JxW();

  // Element basis functions
  const std::vector<std::vector<Real> >          &phi = elem_fe->get_phi();
  const std::vector<std::vector<RealGradient> > &dphi = elem_fe->get_dphi();

  // Quadrature point locations
  const std::vector<Point > &q_point = elem_fe->get_xyz();

  // The number of local degrees of freedom in each variable
  const unsigned int n_T_dofs = c.get_dof_indices(0).size();

  // The subvectors and submatrices we need to fill:
  DenseSubMatrix<Number> &K = c.get_elem_jacobian(0,0);
  DenseSubVector<Number> &F = c.get_elem_residual(0);

  // Now we will build the element Jacobian and residual.
  // Constructing the residual requires the solution and its
  // gradient from the previous timestep.  This must be
  // calculated at each quadrature point by summing the
  // solution degree-of-freedom values by the appropriate
  // weight functions.
  unsigned int n_qpoints = c.get_element_qrule().n_points();

  for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
      // Compute the solution gradient at the Newton iterate
      Gradient grad_T = c.interior_gradient(0, qp);

      // Location of the current qp
      const Real x = q_point[qp](0);
      const Real y = q_point[qp](1);

      // Forcing function
      Real f = -alpha * ( ( (- 4 * alpha * alpha) * exp(-alpha*x) * y * (1 - y) ) + ( -8 + ( 8 * exp(-alpha*x) ) + ( 8 * ( 1 - exp(-alpha) )* x) ) );

      // The residual contribution from this element
      for (unsigned int i=0; i != n_T_dofs; i++)
        F(i) += JxW[qp] * ( (f * phi[i][qp]) - alpha*(grad_T * dphi[i][qp]) ) ;
      if (compute_jacobian)
        for (unsigned int i=0; i != n_T_dofs; i++)
          for (unsigned int j=0; j != n_T_dofs; ++j)
            // The analytic jacobian
            K(i,j) += JxW[qp] * ( -alpha*(dphi[i][qp] * dphi[j][qp]) );
    } // end of the quadrature point qp-loop

  return compute_jacobian;
}

// Override the default DiffSystem postprocess function to compute the
// approximations to the QoIs
void PoissonSystem::postprocess()
{
  // Reset the array holding the computed QoIs
  computed_QoI[0] = 0.0;

  FEMSystem::postprocess();

  this->comm().sum(computed_QoI[0]);
}
