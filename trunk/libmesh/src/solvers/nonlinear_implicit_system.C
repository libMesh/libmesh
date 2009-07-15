// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
  
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
  
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA



// C++ includes

// Local includes
#include "nonlinear_implicit_system.h"
#include "equation_systems.h"
#include "libmesh_logging.h"
#include "linear_solver.h"
#include "nonlinear_solver.h"
#include "numeric_vector.h"
#include "sparse_matrix.h"

// ------------------------------------------------------------
// NonlinearImplicitSystem implementation
NonlinearImplicitSystem::NonlinearImplicitSystem (EquationSystems& es,
						  const std::string& name,
						  const unsigned int number) :
  
  Parent                    (es, name, number),
  nonlinear_solver          (NonlinearSolver<Number>::build(*this)),
  _n_nonlinear_iterations   (0),
  _final_nonlinear_residual (1.e20)
{
  // Set default parameters
  // These were chosen to match the Petsc defaults
  es.parameters.set<Real>        ("linear solver tolerance") = 1e-5;
  es.parameters.set<Real>        ("linear solver minimum tolerance") = 1e-5;
  es.parameters.set<unsigned int>("linear solver maximum iterations") = 10000;
  
  es.parameters.set<unsigned int>("nonlinear solver maximum iterations") = 50;
  es.parameters.set<unsigned int>("nonlinear solver maximum function evaluations") = 10000;

  es.parameters.set<Real>("nonlinear solver absolute residual tolerance") = 1e-50;
  es.parameters.set<Real>("nonlinear solver relative residual tolerance") = 1e-8;
  es.parameters.set<Real>("nonlinear solver absolute step tolerance") = 1e-8;
  es.parameters.set<Real>("nonlinear solver relative step tolerance") = 1e-8;
}



NonlinearImplicitSystem::~NonlinearImplicitSystem ()
{
  // Clear data
  this->clear();
}



void NonlinearImplicitSystem::clear ()
{
  // clear the nonlinear solver
  nonlinear_solver->clear();
  
  // clear the parent data
  Parent::clear();
}



void NonlinearImplicitSystem::reinit ()
{
  // re-initialize the nonlinear solver interface
  nonlinear_solver->clear();
  
  // initialize parent data
  Parent::reinit();  
}



void NonlinearImplicitSystem::set_solver_parameters ()
{
  // Get a reference to the EquationSystems
  const EquationSystems& es =
    this->get_equation_systems();  
  
  // Get the user-specifiied nonlinear solver tolerances
  const unsigned int maxits =
    es.parameters.get<unsigned int>("nonlinear solver maximum iterations");

  const unsigned int maxfuncs =
    es.parameters.get<unsigned int>("nonlinear solver maximum function evaluations");

  const Real abs_resid_tol =
    es.parameters.get<Real>("nonlinear solver absolute residual tolerance");

  const Real rel_resid_tol =
    es.parameters.get<Real>("nonlinear solver relative residual tolerance");

  const Real abs_step_tol =
    es.parameters.get<Real>("nonlinear solver absolute step tolerance");

  const Real rel_step_tol =
    es.parameters.get<Real>("nonlinear solver relative step tolerance");

  // Get the user-specified linear solver toleranaces
  const unsigned int maxlinearits =
    es.parameters.get<unsigned int>("linear solver maximum iterations");

  const Real linear_tol =
    es.parameters.get<Real>("linear solver tolerance");

  const Real linear_min_tol =
    es.parameters.get<Real>("linear solver minimum tolerance");

  // Set all the parameters on the NonlinearSolver
  nonlinear_solver->max_nonlinear_iterations = maxits;
  nonlinear_solver->max_function_evaluations = maxfuncs;
  nonlinear_solver->absolute_residual_tolerance = abs_resid_tol;
  nonlinear_solver->relative_residual_tolerance = rel_resid_tol;
  nonlinear_solver->absolute_step_tolerance = abs_step_tol;
  nonlinear_solver->relative_step_tolerance = rel_step_tol;
  nonlinear_solver->max_linear_iterations = maxlinearits;
  nonlinear_solver->initial_linear_tolerance = linear_tol;
  nonlinear_solver->minimum_linear_tolerance = linear_min_tol;
}



void NonlinearImplicitSystem::solve ()
{
  // Log how long the nonlinear solve takes.
  START_LOG("solve()", "System");
  
  this->set_solver_parameters();

  // Solve the nonlinear system.
  const std::pair<unsigned int, Real> rval =
    nonlinear_solver->solve (*matrix, *solution, *rhs, 
			     nonlinear_solver->relative_residual_tolerance,
                             nonlinear_solver->max_linear_iterations);

  // Store the number of nonlinear iterations required to
  // solve and the final residual.
  _n_nonlinear_iterations   = rval.first;
  _final_nonlinear_residual = rval.second;
    
  // Stop logging the nonlinear solve
  STOP_LOG("solve()", "System");

  // Update the system after the solve
  this->update();  
}



void NonlinearImplicitSystem::adjoint_solve ()
{
  // Log how long the nonlinear solve takes.
  START_LOG("adjoint_solve()", "System");

  // Adding an adjoint_solution vector, allocate an adjoint_solution if it doesn't already exist
  
  NumericVector<Number> & adjoint_solution = System::add_vector("adjoint_solution");
  
  // The nonlinear system should now already be solved.
  // Now assemble it's adjoint.

  if (this->assemble_before_solve)
    {
      // Build the Jacobian
      if (nonlinear_solver->jacobian != NULL)
        nonlinear_solver->jacobian (*current_local_solution.get(), *matrix);
      else if (nonlinear_solver->matvec != NULL)
        nonlinear_solver->matvec   (*current_local_solution.get(), NULL, matrix);
      else libmesh_error();

      // Take the discretization adjoint
      matrix->get_transpose(*matrix);

      // Reset and build the RHS from the QOI derivative
      rhs->zero();
      this->assemble_qoi_derivative();
    }

  // The adjoint problem is linear, so we only use the nonlinear
  // solver to get tolerances and the current jacobian from it.

  this->set_solver_parameters();

  AutoPtr<LinearSolver<Number> > linear_solver = LinearSolver<Number>::build();

  const std::pair<unsigned int, Real> rval =
    linear_solver->solve (*matrix, adjoint_solution, *rhs,
			   nonlinear_solver->relative_residual_tolerance,
                           nonlinear_solver->max_linear_iterations);

  // Store the number of nonlinear iterations required to
  // solve and the final residual.
  _n_nonlinear_iterations   = rval.first;
  _final_nonlinear_residual = rval.second;
    
  // Stop logging the nonlinear solve
  STOP_LOG("adjoint_solve()", "System");

  // Update the system after the solve
  this->update();  
}



void NonlinearImplicitSystem::qoi_parameter_sensitivity
  (std::vector<Number *>& parameters,
   std::vector<Number>&   sensitivities)
{
  // Get ready to fill in senstivities:
  sensitivities.clear();
  sensitivities.resize(parameters.size(), 0);

  // An introduction to the problem:
  //
  // Residual R(u(p),p) = 0
  // partial R / partial u = J = system matrix
  //
  // This implies that:
  // d/dp(R) = 0
  // (partial R / partial p) + 
  // (partial R / partial u) * (partial u / partial p) = 0

  // We first do an adjoint solve:
  // J^T * z = (partial q / partial u)

  this->adjoint_solve();

  // We use the identities:
  // dq/dp = (partial q / partial p) + (partial q / partial u) *
  //         (partial u / partial p)
  // dq/dp = (partial q / partial p) + (J^T * z) *
  //         (partial u / partial p)
  // dq/dp = (partial q / partial p) + z * J *
  //         (partial u / partial p)
 
  // Leading to our final formula:
  // dq/dp = (partial q / partial p) - z * (partial R / partial p)

  for (unsigned int i=0; i != parameters.size(); ++i)
    {
      // We currently get partial derivatives via central differencing
      Number delta_p = 1e-6;

      // (partial q / partial p) ~= (q(p+dp)-q(p-dp))/(2*dp)
      // (partial R / partial p) ~= (rhs(p+dp) - rhs(p-dp))/(2*dp)

      Number old_parameter = *parameters[i];
      Number old_qoi = this->qoi;

      *parameters[i] = old_parameter - delta_p;
      this->assemble_qoi();
      Number qoi_minus = this->qoi;

      this->assemble();
      this->rhs->close();
      AutoPtr<NumericVector<Number> > partialR_partialp = this->rhs->clone();
      *partialR_partialp *= -1;

      *parameters[i] = old_parameter + delta_p;
      this->assemble_qoi();
      Number qoi_plus = this->qoi;
      Number partialq_partialp = (qoi_plus - qoi_minus) / (2.*delta_p);

      this->assemble();
      this->rhs->close();
      *partialR_partialp -= *this->rhs;
      *partialR_partialp /= (2.*delta_p);

      // Don't leave the parameter changed
      *parameters[i] = old_parameter;

      sensitivities[i] = partialq_partialp +
			 partialR_partialp->dot(this->get_adjoint_solution());
    }

  // All parameters have been reset.
  // Don't leave the qoi or system changed - principle of least
  // surprise.
  this->assemble();
  this->assemble_qoi();
}
