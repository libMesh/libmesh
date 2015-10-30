// $Id: continuation_system.C,v 1.1 2007-08-02 17:27:19 jwpeterson Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2005  Benjamin S. Kirk, John W. Peterson
  
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

// LibMesh includes
#include "continuation_system.h"
#include "linear_solver.h"
#include "time_solver.h"
#include "newton_solver.h"


ContinuationSystem::ContinuationSystem (EquationSystems& es,
					const std::string& name,
					const unsigned int number)
  : Parent(es, name, number),
    continuation_parameter(NULL),
    quiet(true),
    ds(0.1),
    continuation_parameter_tolerance(1.e-6),
    solution_tolerance(1.e-6),
    rhs_mode(Residual),
    linear_solver(LinearSolver<Number>::build()),
    tangent_initialized(false),
    newton_solver(NULL),
    old_continuation_parameter(0.),
    dlambda_ds(1.0)
{
  // Warn about using untested code
  untested();
}




ContinuationSystem::~ContinuationSystem ()
{
  this->clear();
}




void ContinuationSystem::clear()
{
  Parent::clear();
}



void ContinuationSystem::init_data ()
{
  // Add a vector which stores the tangent "du/ds" to the system and save its pointer.
  du_ds = &(add_vector("du_ds"));

  // Add a vector to keep track of the previous nonlinear solution
  // at the old value of lambda.
  previous_u = &(add_vector("previous_u"));

  // Add a vector to keep track of the temporary solution "y" of Ay=G_{\lambda}.
  y = &(add_vector("y"));

  // Add a vector to keep track of the temporary solution "z" of Az=-G.
  z = &(add_vector("z"));

  // Add a vector to keep track of the Newton update during the constrained PDE solves.
  delta_u = &(add_vector("delta_u"));

  // Call the Parent's initialization routine.
  Parent::init_data();
}




void ContinuationSystem::solve()
{
  // Set the Residual RHS mode, and call the normal solve routine.
  rhs_mode      = Residual;
  DifferentiableSystem::solve();
}



// This is most of the "guts" of this class.  This is where we implement
// our custom Newton iterations and perform most of the solves.
void ContinuationSystem::continuation_solve()
{
  // Be sure the user has set the continuation parameter pointer
  if (!continuation_parameter)
    {
      std::cerr << "You must set the continuation_parameter pointer "
		<< "to a member variable of the derived class, preferably in the "
		<< "Derived class's init_data function.  This is how the ContinuationSystem "
		<< "updates the continuation parameter."
		<< std::endl;
      
      error();
    }

  
  // We can't start solving the augmented PDE system unless the tangent
  // vectors have been initialized.  This only needs to occur once.
  if (!tangent_initialized)
    initialize_tangent();


  // Set pointer to underlying Newton solver
  if (!newton_solver)
    newton_solver =
      dynamic_cast<NewtonSolver*> (this->time_solver->diff_solver().get());
  assert (newton_solver != NULL);
  
  // A pair for catching return values from linear system solves.
  std::pair<unsigned int, Real> rval;

  // The nonlinear loop
  unsigned int l=0;
  for (; l<newton_solver->max_nonlinear_iterations; ++l)
    {
      std::cout << "Starting Newton step " << l << std::endl;
      
      // Assemble the system matrix AND rhs, with rhs = G_{\lambda}
      rhs_mode = G_Lambda;

      // Assemble both rhs and Jacobian
      assembly(true, true);

      // Not sure if this is really necessary
      rhs->close();
      
      //std::cout << "||y||=" << y->l2_norm() << std::endl;
      
      // Solve G_u*y = G_{\lambda}
      // Initial guess?
      y->zero();
      y->close();
      rval =
	linear_solver->solve(*matrix,
			     *y,
			     *rhs,
			     1.e-12, // linear tolerance
			     newton_solver->max_linear_iterations);   // max linear iterations
      //here();

      if (!quiet)
	std::cout << "G_u*y = G_{lambda} solver converged at step "
		  << rval.first
		  << ", linear tolerance = "
		  << rval.second
		  << "."
		  << std::endl;

      //std::cout << "||y||=" << y->l2_norm() << std::endl;
      
      
      // Now, reassemble the system (skipping the Jacobian) with rhs = the normal residual.
      rhs_mode = Residual;
      assembly(true, false);
      rhs->close();
      
      // Solve the linear system G_u*z = G
      // Initial guess?
      z->zero(); // It seems to be extremely important to zero z here, otherwise the solver quits early.
      z->close();
      rval =
	linear_solver->solve(*matrix,
			     *z,
			     *rhs,
			     1.e-12, // linear tolerance
			     newton_solver->max_linear_iterations);   // max linear iterations
      //here();

      if (!quiet)
	std::cout << "G_u*z = G solver converged at step "
		  << rval.first
		  << " linear tolerance = "
		  << rval.second
		  << "."
		  << std::endl;
      
      // Note: need to scale z by -1 since our code always solves Jx=R
      // instead of Jx=-R.
      z->scale(-1.);
      z->close();
      //std::cout << "||z||=" << z->l2_norm() << std::endl;
      
      
      // Compute N, the residual of the arclength constraint eqn.
      // Note 1: N(u,lambda,s) := (u-u_{old}, du_ds) + (lambda-lambda_{old}, dlambda_ds) - _ds 
      // We temporarily use the delta_u vector as a temporary vector for this calculation.
      *delta_u = *solution;
      delta_u->add(-1., *previous_u);

      Real N =
	delta_u->dot(*du_ds) +
	((*continuation_parameter) - old_continuation_parameter)*dlambda_ds -
	ds;

      //std::cout << "N=" << N << std::endl;

      const Real duds_dot_z = du_ds->dot(*z);
      const Real duds_dot_y = du_ds->dot(*y);

      //std::cout << "duds_dot_z=" << duds_dot_z << std::endl;
      //std::cout << "duds_dot_y=" << duds_dot_y << std::endl;
      //std::cout << "dlambda_ds=" << dlambda_ds << std::endl;

      const Real delta_lambda_numerator   = -N -duds_dot_z;
      const Real delta_lambda_denominator = dlambda_ds  - duds_dot_y;

      assert (delta_lambda_denominator != 0.0);
      //std::cout << "delta_lambda_numerator=" << delta_lambda_numerator << std::endl;
      //std::cout << "delta_lambda_denominator=" << delta_lambda_denominator << std::endl;

      // Now, we are ready to compute the step delta_lambda
      // delta_lambda = (-N - dot(N_u,z)) / (N_{\lambda} - dot(N_u,y))
      Real delta_lambda = delta_lambda_numerator / delta_lambda_denominator;
      //	(-N          - duds_dot_z)  /  (dlambda_ds  - duds_dot_y);

      // Knowing delta_lambda, we are ready to update delta_u
      // delta_u = z - delta_lambda*y
      delta_u->zero();
      delta_u->add(1., *z);
      delta_u->add(-delta_lambda, *y);
      delta_u->close();


      // Update the system solution and the system parameter.
      //std::cout << "||solution||=" << solution->l2_norm() << std::endl;
      solution->add(1., *delta_u);
      solution->close();
      //std::cout << "||solution||=" << solution->l2_norm() << std::endl;
      *continuation_parameter += delta_lambda;
      
      // Check the convergence of the parameter and the solution.  If they are small
      // enough, we can break out of the Newton iteration loop.  
      const Real norm_delta_u = delta_u->l2_norm();
      std::cout << "  delta_lambda = " << delta_lambda << std::endl;
      std::cout << "  ||delta_u|| = " << norm_delta_u << std::endl;

      if ((fabs(delta_lambda) < continuation_parameter_tolerance) &&
	  (norm_delta_u       < solution_tolerance))
	{
	  if (!quiet)
	    std::cout << "Newton iterations converged!" << std::endl;
	  
	  break;
	}
    }

  //std::cout << "l=" << l << std::endl;
  if (l == newton_solver->max_nonlinear_iterations)
    {
      std::cout << "Newton iterations of augmented system did not converge!" << std::endl;
      error();
    }

  // Print converged solution control parameter and max value.
  std::cout << "lambda_current=" << *continuation_parameter << std::endl;
  std::cout << "u_max=" << solution->max() << std::endl;
  
  
  // Solve for the updated tangent du1/ds, d(lambda1)/ds
  solve_tangent();

  // Advance the solution and the parameter to the next value.
  update_solution();
}



// This function solves the tangent system:
// [ G_u          G_{lambda}        ][(du/ds)_new      ] = [  0 ]
// [ (du/ds)_old  (dlambda/ds)_old  ][(dlambda/ds)_new ]   [-N_s]
// The solution is given by:
// 1. Let G_u y = G_lambda
// 2. (du_ds)_new = -(dlambda/ds)_new * y
// 3. (dlambda/ds)_new = 1.0 / ( (dlambda/ds)_old - (du/ds)_old*y )
void ContinuationSystem::solve_tangent()
{
  // We shouldn't call this unless the current tangent already makes sense.
  assert (tangent_initialized);
  
  // Set pointer to underlying Newton solver
  if (!newton_solver)
    newton_solver =
      dynamic_cast<NewtonSolver*> (this->time_solver->diff_solver().get());
  assert (newton_solver != NULL);
  
  // Assemble the system matrix AND rhs, with rhs = G_{\lambda}
  this->rhs_mode = G_Lambda;
  
  // Assemble both rhs and Jacobian
  this->assembly(true, true);

  // Not sure if this is really necessary
  rhs->close();
      
  // Solve G_u*y =  G_{\lambda}
  std::pair<unsigned int, Real> rval = 
    linear_solver->solve(*matrix,
			 *y,
			 *rhs,
			 1.e-12, // linear tolerance
			 newton_solver->max_linear_iterations);   // max linear iterations
  
  if (!quiet)
    std::cout << "G_u*y = G_{lambda} solver converged at step "
	      << rval.first
	      << " linear tolerance = "
	      << rval.second
	      << "."
	      << std::endl;

  
  // Solve for the updated d(lambda)/ds
  // denom = N_{lambda}   - (du_ds)^t y
  //       = d(lambda)/ds - (du_ds)^t y
  Real denom = dlambda_ds - du_ds->dot(*y);

  //std::cout << "denom=" << denom << std::endl;
  assert (denom != 0.0);
  
  dlambda_ds = 1.0 / denom;
 
  // Compute the updated value of du/ds = -_dlambda_ds * y
  du_ds->zero();
  du_ds->add(-dlambda_ds, *y);
  du_ds->close();

  //std::cout << "||du_ds||=" << this->du_ds->l2_norm() << std::endl;
  
  //error();
}



void ContinuationSystem::initialize_tangent()
{
  // Be sure the tangent was not already initialized.
  assert (!tangent_initialized);
  
  // Compute delta_s_zero, the initial arclength travelled during the
  // first step.  Here we assume that previous_u and lambda_old store
  // the previous solution and control parameter.  You may need to
  // read in an old solution (or solve the non-continuation system)
  // first and call save_current_solution() before getting here.
  
  // Compute norms of the current and previous solutions
  Real norm_u          = solution->l2_norm();
  Real norm_previous_u = previous_u->l2_norm();

  if (norm_u == norm_previous_u)
    {
      std::cerr << "Warning, it appears u and previous_u are the "
		<< "same, are you sure this is correct?"
		<< "It's possible you forgot to set one or the other..."
		<< std::endl;
    }
  //std::cout << "||u_old||=" << this->solution->l2_norm() << std::endl;
  
  Real delta_s_zero = std::sqrt(
				(norm_u - norm_previous_u)*(norm_u - norm_previous_u) +
				(*continuation_parameter-old_continuation_parameter)*
				(*continuation_parameter-old_continuation_parameter)
				);

  if (!quiet)
    std::cout << "delta_s_zero=" << delta_s_zero << std::endl;

  // Now approximate the initial tangent d(lambda)/ds
  this->dlambda_ds = (*continuation_parameter-old_continuation_parameter) / delta_s_zero;

  if (!quiet)
    std::cout << "d(lambda_0)/ds=" << this->dlambda_ds << std::endl;

  // We can also approximate the deriv. wrt s by finite differences:
  // du/ds = (u1 - u0) / delta_s_zero.
  *du_ds = *solution;
  du_ds->add(-1., *previous_u);
  du_ds->scale(1./delta_s_zero);
  du_ds->close();
  
  // Set the flag which tells us the method has been initialized.
  tangent_initialized = true;
  
  // Advance the solution and the parameter to the next value.
  update_solution();
}




void ContinuationSystem::update_solution()
{
  // We must have a tangent that makes sense before we can update the solution.
  assert (tangent_initialized);

  // Save the current value of u and lambda before updating.
  save_current_solution();

  // Since we solved for the tangent vector, now we can compute
  // an initial guess for the new solution, and an initial guess for the new value of lambda.
  solution->add(ds, *du_ds);
  solution->close();

  //std::cout << "||u_new||=" << this->solution->l2_norm() << std::endl;
  //here();

  *continuation_parameter += ds*dlambda_ds;

  if (!quiet)
    std::cout << "lambda_new=" << *continuation_parameter << std::endl;
}




void ContinuationSystem::save_current_solution()
{
  // Save the old solution vector
  *previous_u = *solution;

  // Save the old value of lambda
  old_continuation_parameter = *continuation_parameter;
}
