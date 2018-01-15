// The libMesh Finite Element Library.
// Copyright (C) 2002-2018 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/equation_systems.h"
#include "libmesh/continuation_system.h"
#include "libmesh/linear_solver.h"
#include "libmesh/time_solver.h"
#include "libmesh/newton_solver.h"
#include "libmesh/sparse_matrix.h"

namespace libMesh
{

ContinuationSystem::ContinuationSystem (EquationSystems & es,
                                        const std::string & name_in,
                                        const unsigned int number_in) :
  Parent(es, name_in, number_in),
  continuation_parameter(libmesh_nullptr),
  quiet(true),
  continuation_parameter_tolerance(1.e-6),
  solution_tolerance(1.e-6),
  initial_newton_tolerance(0.01),
  old_continuation_parameter(0.),
  min_continuation_parameter(0.),
  max_continuation_parameter(0.),
  Theta(1.),
  Theta_LOCA(1.),
  n_backtrack_steps(5),
  n_arclength_reductions(5),
  ds_min(1.e-8),
  predictor(Euler),
  newton_stepgrowth_aggressiveness(1.),
  newton_progress_check(true),
  rhs_mode(Residual),
  linear_solver(LinearSolver<Number>::build(es.comm())),
  tangent_initialized(false),
  newton_solver(libmesh_nullptr),
  dlambda_ds(0.707),
  ds(0.1),
  ds_current(0.1),
  previous_dlambda_ds(0.),
  previous_ds(0.),
  newton_step(0)
{
  // Warn about using untested code
  libmesh_experimental();

  if (libMesh::on_command_line("--solver_system_names"))
    linear_solver->init((this->name()+"_").c_str());
  else
    linear_solver->init();
}




ContinuationSystem::~ContinuationSystem ()
{
  this->clear();
}




void ContinuationSystem::clear()
{
  // FIXME: Do anything here, e.g. zero vectors, etc?

  // Call the Parent's clear function
  Parent::clear();
}



void ContinuationSystem::init_data ()
{
  // Add a vector which stores the tangent "du/ds" to the system and save its pointer.
  du_ds = &(add_vector("du_ds"));

  // Add a vector which stores the tangent "du/ds" to the system and save its pointer.
  previous_du_ds = &(add_vector("previous_du_ds"));

  // Add a vector to keep track of the previous nonlinear solution
  // at the old value of lambda.
  previous_u = &(add_vector("previous_u"));

  // Add a vector to keep track of the temporary solution "y" of Ay=G_{\lambda}.
  y = &(add_vector("y"));

  // Add a vector to keep track of the "old value" of "y" which is the solution of Ay=G_{\lambda}.
  y_old = &(add_vector("y_old"));

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




void ContinuationSystem::initialize_tangent()
{
  // Be sure the tangent was not already initialized.
  libmesh_assert (!tangent_initialized);

  // Compute delta_s_zero, the initial arclength traveled during the
  // first step.  Here we assume that previous_u and lambda_old store
  // the previous solution and control parameter.  You may need to
  // read in an old solution (or solve the non-continuation system)
  // first and call save_current_solution() before getting here.

  // 1.) Compute delta_s_zero as ||u|| - ||u_old|| + ...
  // Compute norms of the current and previous solutions
  //   Real norm_u          = solution->l2_norm();
  //   Real norm_previous_u = previous_u->l2_norm();

  //   if (!quiet)
  //     {
  //       libMesh::out << "norm_u=" << norm_u << std::endl;
  //       libMesh::out << "norm_previous_u=" << norm_previous_u << std::endl;
  //     }

  //   if (norm_u == norm_previous_u)
  //     {
  //       libMesh::err << "Warning, it appears u and previous_u are the "
  //   << "same, are you sure this is correct?"
  //   << "It's possible you forgot to set one or the other..."
  //   << std::endl;
  //     }

  //   Real delta_s_zero = std::sqrt(
  //   (norm_u - norm_previous_u)*(norm_u - norm_previous_u) +
  //   (*continuation_parameter-old_continuation_parameter)*
  //   (*continuation_parameter-old_continuation_parameter)
  //   );

  //   // 2.) Compute delta_s_zero as ||u -u_old|| + ...
  //   *delta_u = *solution;
  //   delta_u->add(-1., *previous_u);
  //   delta_u->close();
  //   Real norm_delta_u = delta_u->l2_norm();
  //   Real norm_u          = solution->l2_norm();
  //   Real norm_previous_u = previous_u->l2_norm();

  //   // Scale norm_delta_u by the bigger of either norm_u or norm_previous_u
  //   norm_delta_u /= std::max(norm_u, norm_previous_u);

  //   if (!quiet)
  //     {
  //       libMesh::out << "norm_u=" << norm_u << std::endl;
  //       libMesh::out << "norm_previous_u=" << norm_previous_u << std::endl;
  //       //libMesh::out << "norm_delta_u=" << norm_delta_u << std::endl;
  //       libMesh::out << "norm_delta_u/max(|u|,|u_old|)=" << norm_delta_u << std::endl;
  //       libMesh::out << "|norm_u-norm_previous_u|=" << std::abs(norm_u - norm_previous_u) << std::endl;
  //     }

  //   const Real dlambda = *continuation_parameter-old_continuation_parameter;

  //   if (!quiet)
  //     libMesh::out << "dlambda=" << dlambda << std::endl;

  //   Real delta_s_zero = std::sqrt(
  //   (norm_delta_u*norm_delta_u) +
  //   (dlambda*dlambda)
  //   );

  //   if (!quiet)
  //     libMesh::out << "delta_s_zero=" << delta_s_zero << std::endl;

  // 1.) + 2.)
  //   // Now approximate the initial tangent d(lambda)/ds
  //   this->dlambda_ds = (*continuation_parameter-old_continuation_parameter) / delta_s_zero;


  //   // We can also approximate the deriv. wrt s by finite differences:
  //   // du/ds = (u1 - u0) / delta_s_zero.
  //   // FIXME: Use delta_u from above if we decide to keep that method.
  //   *du_ds = *solution;
  //   du_ds->add(-1., *previous_u);
  //   du_ds->scale(1./delta_s_zero);
  //   du_ds->close();


  // 3.) Treating (u-previous_u)/(lambda - lambda_old) as an approximation to du/d(lambda),
  // we follow the same technique as Carnes and Shadid.
  //   const Real dlambda = *continuation_parameter-old_continuation_parameter;
  //   libmesh_assert_greater (dlambda, 0.);

  //   // Use delta_u for temporary calculation of du/d(lambda)
  //   *delta_u = *solution;
  //   delta_u->add(-1., *previous_u);
  //   delta_u->scale(1. / dlambda);
  //   delta_u->close();

  //   // Determine initial normalization parameter
  //   const Real solution_size = std::max(solution->l2_norm(), previous_u->l2_norm());
  //   if (solution_size > 1.)
  //     {
  //       Theta = 1./solution_size;

  //       if (!quiet)
  // libMesh::out << "Setting Normalization Parameter Theta=" << Theta << std::endl;
  //     }

  //   // Compute d(lambda)/ds
  //   // The correct sign of d(lambda)/ds should be positive, since we assume that (lambda > lambda_old)
  //   // but we could always double-check that as well.
  //   Real norm_delta_u = delta_u->l2_norm();
  //   this->dlambda_ds = 1. / std::sqrt(1. + Theta*Theta*norm_delta_u*norm_delta_u);

  //   // Finally, compute du/ds = d(lambda)/ds * du/d(lambda)
  //   *du_ds = *delta_u;
  //   du_ds->scale(dlambda_ds);
  //   du_ds->close();


  // 4.) Use normalized arclength formula to estimate delta_s_zero
  //   // Determine initial normalization parameter
  //   set_Theta();

  //   // Compute (normalized) delta_s_zero
  //   *delta_u = *solution;
  //   delta_u->add(-1., *previous_u);
  //   delta_u->close();
  //   Real norm_delta_u = delta_u->l2_norm();

  //   const Real dlambda = *continuation_parameter-old_continuation_parameter;

  //   if (!quiet)
  //     libMesh::out << "dlambda=" << dlambda << std::endl;

  //   Real delta_s_zero = std::sqrt(
  //   (Theta_LOCA*Theta_LOCA*Theta*norm_delta_u*norm_delta_u) +
  //   (dlambda*dlambda)
  //   );
  //   *du_ds = *delta_u;
  //   du_ds->scale(1./delta_s_zero);
  //   dlambda_ds = dlambda / delta_s_zero;

  //   if (!quiet)
  //     {
  //       libMesh::out << "delta_s_zero=" << delta_s_zero << std::endl;
  //       libMesh::out << "initial d(lambda)/ds|_0 = " << dlambda_ds << std::endl;
  //       libMesh::out << "initial ||du_ds||_0 = " << du_ds->l2_norm() << std::endl;
  //     }

  //   // FIXME: Also store the initial finite-differenced approximation to -du/dlambda as y.
  //   // We stick to the convention of storing negative y, since that is what we typically
  //   // solve for anyway.
  //   *y = *delta_u;
  //   y->scale(-1./dlambda);
  //   y->close();



  // 5.) Assume dlambda/ds_0 ~ 1/sqrt(2) and determine the value of Theta_LOCA which
  // will satisfy this criterion

  // Initial change in parameter
  const Real dlambda = *continuation_parameter-old_continuation_parameter;
  libmesh_assert_not_equal_to (dlambda, 0.0);

  // Ideal initial value of dlambda_ds
  dlambda_ds = 1. / std::sqrt(2.);
  if (dlambda < 0.)
    dlambda_ds *= -1.;

  // This also implies the initial value of ds
  ds_current = dlambda / dlambda_ds;

  if (!quiet)
    libMesh::out << "Setting ds_current|_0=" << ds_current << std::endl;

  // Set y = -du/dlambda using finite difference approximation
  *y = *solution;
  y->add(-1., *previous_u);
  y->scale(-1./dlambda);
  y->close();
  const Real ynorm=y->l2_norm();

  // Finally, set the value of du_ds to be used in the upcoming
  // tangent calculation. du/ds = du/dlambda * dlambda/ds
  *du_ds = *y;
  du_ds->scale(-dlambda_ds);
  du_ds->close();

  // Determine additional solution normalization parameter
  // (Since we just set du/ds, it will be:  ||du||*||du/ds||)
  set_Theta();

  // The value of Theta_LOCA which makes dlambda_ds = 1/sqrt(2),
  // assuming our Theta = ||du||^2.
  // Theta_LOCA = std::abs(dlambda);

  // Assuming general Theta
  Theta_LOCA = std::sqrt(1./Theta/ynorm/ynorm);


  if (!quiet)
    {
      libMesh::out << "Setting initial Theta_LOCA = " << Theta_LOCA << std::endl;
      libMesh::out << "Theta_LOCA^2*Theta         = " << Theta_LOCA*Theta_LOCA*Theta << std::endl;
      libMesh::out << "initial d(lambda)/ds|_0    = " << dlambda_ds << std::endl;
      libMesh::out << "initial ||du_ds||_0        = " << du_ds->l2_norm() << std::endl;
    }



  // OK, we estimated the tangent at point u0.
  // Now, to estimate the tangent at point u1, we call the solve_tangent routine.

  // Set the flag which tells us the method has been initialized.
  tangent_initialized = true;

  solve_tangent();

  // Advance the solution and the parameter to the next value.
  update_solution();
}






// This is most of the "guts" of this class.  This is where we implement
// our custom Newton iterations and perform most of the solves.
void ContinuationSystem::continuation_solve()
{
  // Be sure the user has set the continuation parameter pointer
  if (!continuation_parameter)
    libmesh_error_msg("You must set the continuation_parameter pointer " \
                      << "to a member variable of the derived class, preferably in the " \
                      << "Derived class's init_data function.  This is how the ContinuationSystem " \
                      << "updates the continuation parameter.");

  // Use extra precision for all the numbers printed in this function.
  std::streamsize old_precision = libMesh::out.precision();
  libMesh::out.precision(16);
  libMesh::out.setf(std::ios_base::scientific);

  // We can't start solving the augmented PDE system unless the tangent
  // vectors have been initialized.  This only needs to occur once.
  if (!tangent_initialized)
    initialize_tangent();

  // Save the old value of -du/dlambda.  This will be used after the Newton iterations
  // to compute the angle between previous tangent vectors.  This cosine of this angle is
  //
  // tau := abs( (du/d(lambda)_i , du/d(lambda)_{i-1}) / (||du/d(lambda)_i|| * ||du/d(lambda)_{i-1}||) )
  //
  // The scaling factor tau (which should vary between 0 and 1) is used to shrink the step-size ds
  // when we are approaching a turning point.  Note that it can only shrink the step size.
  *y_old = *y;

  // Set pointer to underlying Newton solver
  if (!newton_solver)
    newton_solver = cast_ptr<NewtonSolver *> (this->time_solver->diff_solver().get());

  // A pair for catching return values from linear system solves.
  std::pair<unsigned int, Real> rval;

  // Convergence flag for the entire arcstep
  bool arcstep_converged = false;

  // Begin loop over arcstep reductions.
  for (unsigned int ns=0; ns<n_arclength_reductions; ++ns)
    {
      if (!quiet)
        {
          libMesh::out << "Current arclength stepsize, ds_current=" << ds_current << std::endl;
          libMesh::out << "Current parameter value, lambda=" << *continuation_parameter << std::endl;
        }

      // Upon exit from the nonlinear loop, the newton_converged flag
      // will tell us the convergence status of Newton's method.
      bool newton_converged = false;

      // The nonlinear residual before *any* nonlinear steps have been taken.
      Real nonlinear_residual_firststep = 0.;

      // The nonlinear residual from the current "k" Newton step, before the Newton step
      Real nonlinear_residual_beforestep = 0.;

      // The nonlinear residual from the current "k" Newton step, after the Newton step
      Real nonlinear_residual_afterstep = 0.;

      // The linear solver tolerance, can be updated dynamically at each Newton step.
      Real current_linear_tolerance = 0.;

      // The nonlinear loop
      for (newton_step=0; newton_step<newton_solver->max_nonlinear_iterations; ++newton_step)
        {
          libMesh::out << "\n === Starting Newton step " << newton_step << " ===" << std::endl;

          // Set the linear system solver tolerance
          //   // 1.) Set the current linear tolerance based as a multiple of the current residual of the system.
          //   const Real residual_multiple = 1.e-4;
          //   Real current_linear_tolerance = residual_multiple*nonlinear_residual_beforestep;

          //   // But if the current residual isn't small, don't let the solver exit with zero iterations!
          //   if (current_linear_tolerance > 1.)
          //     current_linear_tolerance = residual_multiple;

          // 2.) Set the current linear tolerance based on the method based on technique of Eisenstat & Walker.
          if (newton_step==0)
            {
              // At first step, only try reducing the residual by a small amount
              current_linear_tolerance = initial_newton_tolerance;//0.01;
            }

          else
            {
              // The new tolerance is based on the ratio of the most recent tolerances
              const Real alp=0.5*(1.+std::sqrt(5.));
              const Real gam=0.9;

              libmesh_assert_not_equal_to (nonlinear_residual_beforestep, 0.0);
              libmesh_assert_not_equal_to (nonlinear_residual_afterstep, 0.0);

              current_linear_tolerance = std::min(gam*std::pow(nonlinear_residual_afterstep/nonlinear_residual_beforestep, alp),
                                                  current_linear_tolerance*current_linear_tolerance
                                                  );

              // Don't let it get ridiculously small!!
              if (current_linear_tolerance < 1.e-12)
                current_linear_tolerance = 1.e-12;
            }

          if (!quiet)
            libMesh::out << "Using current_linear_tolerance=" << current_linear_tolerance << std::endl;


          // Assemble the residual (and Jacobian).
          rhs_mode = Residual;
          assembly(true,   // Residual
                   true); // Jacobian
          rhs->close();

          // Save the current nonlinear residual.  We don't need to recompute the residual unless
          // this is the first step, since it was already computed as part of the convergence check
          // at the end of the last loop iteration.
          if (newton_step==0)
            {
              nonlinear_residual_beforestep = rhs->l2_norm();

              // Store the residual before any steps have been taken.  This will *not*
              // be updated at each step, and can be used to see if any progress has
              // been made from the initial residual at later steps.
              nonlinear_residual_firststep = nonlinear_residual_beforestep;

              const Real old_norm_u = solution->l2_norm();
              libMesh::out << "  (before step) ||R||_{L2} = " << nonlinear_residual_beforestep << std::endl;
              libMesh::out << "  (before step) ||R||_{L2}/||u|| = " << nonlinear_residual_beforestep / old_norm_u << std::endl;

              // In rare cases (very small arcsteps), it's possible that the residual is
              // already below our absolute linear tolerance.
              if (nonlinear_residual_beforestep  < solution_tolerance)
                {
                  if (!quiet)
                    libMesh::out << "Initial guess satisfied linear tolerance, exiting with zero Newton iterations!" << std::endl;

                  // Since we go straight from here to the solve of the next tangent, we
                  // have to close the matrix before it can be assembled again.
                  matrix->close();
                  newton_converged=true;
                  break; // out of Newton iterations, with newton_converged=true
                }
            }

          else
            {
              nonlinear_residual_beforestep = nonlinear_residual_afterstep;
            }


          // Solve the linear system G_u*z = G
          // Initial guess?
          z->zero(); // It seems to be extremely important to zero z here, otherwise the solver quits early.
          z->close();

          // It's possible that we have selected the current_linear_tolerance so large that
          // a guess of z=zero yields a linear system residual |Az + R| small enough that the
          // linear solver exits in zero iterations.  If this happens, we will reduce the
          // current_linear_tolerance until the linear solver does at least 1 iteration.
          do
            {
              rval =
                linear_solver->solve(*matrix,
                                     *z,
                                     *rhs,
                                     //1.e-12,
                                     current_linear_tolerance,
                                     newton_solver->max_linear_iterations);   // max linear iterations

              if (rval.first==0)
                {
                  if (newton_step==0)
                    {
                      libMesh::out << "Repeating initial solve with smaller linear tolerance!" << std::endl;
                      current_linear_tolerance *= initial_newton_tolerance; // reduce the linear tolerance to force the solver to do some work
                    }
                  else
                    // We shouldn't get here ... it means the linear solver did no work on a Newton
                    // step other than the first one.  If this happens, we need to think more about our
                    // tolerance selection.
                    libmesh_error_msg("Linear solver did no work!");
                }

            } while (rval.first==0);


          if (!quiet)
            libMesh::out << "  G_u*z = G solver converged at step "
                         << rval.first
                         << " linear tolerance = "
                         << rval.second
                         << "."
                         << std::endl;

          // Sometimes (I am not sure why) the linear solver exits after zero iterations.
          // Perhaps it is hitting PETSc's divergence tolerance dtol???  If this occurs,
          // we should break out of the Newton iteration loop because nothing further is
          // going to happen...  Of course if the tolerance is already small enough after
          // zero iterations (how can this happen?!) we should not quit.
          if ((rval.first == 0) && (rval.second > current_linear_tolerance*nonlinear_residual_beforestep))
            {
              if (!quiet)
                libMesh::out << "Linear solver exited in zero iterations!" << std::endl;

              // Try to find out the reason for convergence/divergence
              linear_solver->print_converged_reason();

              break; // out of Newton iterations
            }

          // Note: need to scale z by -1 since our code always solves Jx=R
          // instead of Jx=-R.
          z->scale(-1.);
          z->close();






          // Assemble the G_Lambda vector, skip residual.
          rhs_mode = G_Lambda;

          // Assemble both rhs and Jacobian
          assembly(true,  // Residual
                   false); // Jacobian

          // Not sure if this is really necessary
          rhs->close();
          const Real yrhsnorm=rhs->l2_norm();
          if (yrhsnorm == 0.0)
            libmesh_error_msg("||G_Lambda|| = 0");

          // We select a tolerance for the y-system which is based on the inexact Newton
          // tolerance but scaled by an extra term proportional to the RHS (which is not -> 0 in this case)
          const Real ysystemtol=current_linear_tolerance*(nonlinear_residual_beforestep/yrhsnorm);
          if (!quiet)
            libMesh::out << "ysystemtol=" << ysystemtol << std::endl;

          // Solve G_u*y = G_{\lambda}
          // FIXME: Initial guess?  This is really a solve for -du/dlambda so we could try
          // initializing it with the latest approximation to that... du/dlambda ~ du/ds * ds/dlambda
          //*y = *solution;
          //y->add(-1., *previous_u);
          //y->scale(-1. / (*continuation_parameter - old_continuation_parameter)); // Be careful of divide by zero...
          //y->close();

          //  const unsigned int max_attempts=1;
          // unsigned int attempt=0;
          //   do
          //     {
          //       if (!quiet)
          // libMesh::out << "Trying to solve tangent system, attempt " << attempt << std::endl;

          rval =
            linear_solver->solve(*matrix,
                                 *y,
                                 *rhs,
                                 //1.e-12,
                                 ysystemtol,
                                 newton_solver->max_linear_iterations);   // max linear iterations

          if (!quiet)
            libMesh::out << "  G_u*y = G_{lambda} solver converged at step "
                         << rval.first
                         << ", linear tolerance = "
                         << rval.second
                         << "."
                         << std::endl;

          // Sometimes (I am not sure why) the linear solver exits after zero iterations.
          // Perhaps it is hitting PETSc's divergence tolerance dtol???  If this occurs,
          // we should break out of the Newton iteration loop because nothing further is
          // going to happen...
          if ((rval.first == 0) && (rval.second > ysystemtol))
            {
              if (!quiet)
                libMesh::out << "Linear solver exited in zero iterations!" << std::endl;

              break; // out of Newton iterations
            }

          //       ++attempt;
          //     } while ((attempt<max_attempts) && (rval.first==newton_solver->max_linear_iterations));





          // Compute N, the residual of the arclength constraint eqn.
          // Note 1: N(u,lambda,s) := (u-u_{old}, du_ds) + (lambda-lambda_{old}, dlambda_ds) - _ds
          // We temporarily use the delta_u vector as a temporary vector for this calculation.
          *delta_u = *solution;
          delta_u->add(-1., *previous_u);

          // First part of the arclength constraint
          const Number N1 = Theta_LOCA*Theta_LOCA*Theta*delta_u->dot(*du_ds);
          const Number N2 = ((*continuation_parameter) - old_continuation_parameter)*dlambda_ds;
          const Number N3 = ds_current;

          if (!quiet)
            {
              libMesh::out << "  N1=" << N1 << std::endl;
              libMesh::out << "  N2=" << N2 << std::endl;
              libMesh::out << "  N3=" << N3 << std::endl;
            }

          // The arclength constraint value
          const Number N = N1+N2-N3;

          if (!quiet)
            libMesh::out << "  N=" << N << std::endl;

          const Number duds_dot_z = du_ds->dot(*z);
          const Number duds_dot_y = du_ds->dot(*y);

          //libMesh::out << "duds_dot_z=" << duds_dot_z << std::endl;
          //libMesh::out << "duds_dot_y=" << duds_dot_y << std::endl;
          //libMesh::out << "dlambda_ds=" << dlambda_ds << std::endl;

          const Number delta_lambda_numerator   = -(N          + Theta_LOCA*Theta_LOCA*Theta*duds_dot_z);
          const Number delta_lambda_denominator =  (dlambda_ds - Theta_LOCA*Theta_LOCA*Theta*duds_dot_y);

          libmesh_assert_not_equal_to (delta_lambda_denominator, 0.0);

          // Now, we are ready to compute the step delta_lambda
          const Number delta_lambda_comp = delta_lambda_numerator /
            delta_lambda_denominator;
          // Lambda is real-valued
          const Real delta_lambda = libmesh_real(delta_lambda_comp);

          // Knowing delta_lambda, we are ready to update delta_u
          // delta_u = z - delta_lambda*y
          delta_u->zero();
          delta_u->add(1., *z);
          delta_u->add(-delta_lambda, *y);
          delta_u->close();

          // Update the system solution and the continuation parameter.
          solution->add(1., *delta_u);
          solution->close();
          *continuation_parameter += delta_lambda;

          // Did the Newton step actually reduce the residual?
          rhs_mode = Residual;
          assembly(true,   // Residual
                   false); // Jacobian
          rhs->close();
          nonlinear_residual_afterstep = rhs->l2_norm();


          // In a "normal" Newton step, ||du||/||R|| > 1 since the most recent
          // step is where you "just were" and the current residual is where
          // you are now.  It can occur that ||du||/||R|| < 1, but these are
          // likely not good cases to attempt backtracking (?).
          const Real norm_du_norm_R = delta_u->l2_norm() / nonlinear_residual_afterstep;
          if (!quiet)
            libMesh::out << "  norm_du_norm_R=" << norm_du_norm_R << std::endl;


          // Factor to decrease the stepsize by for backtracking
          Real newton_stepfactor = 1.;

          const bool attempt_backtracking =
            (nonlinear_residual_afterstep > solution_tolerance)
            && (nonlinear_residual_afterstep > nonlinear_residual_beforestep)
            && (n_backtrack_steps>0)
            && (norm_du_norm_R > 1.)
            ;

          // If residual is not reduced, do Newton back tracking.
          if (attempt_backtracking)
            {
              if (!quiet)
                libMesh::out << "Newton step did not reduce residual." << std::endl;

              // back off the previous step.
              solution->add(-1., *delta_u);
              solution->close();
              *continuation_parameter -= delta_lambda;

              // Backtracking: start cutting the Newton stepsize by halves until
              // the new residual is actually smaller...
              for (unsigned int backtrack_step=0; backtrack_step<n_backtrack_steps; ++backtrack_step)
                {
                  newton_stepfactor *= 0.5;

                  if (!quiet)
                    libMesh::out << "Shrinking step size by " << newton_stepfactor << std::endl;

                  // Take fractional step
                  solution->add(newton_stepfactor, *delta_u);
                  solution->close();
                  *continuation_parameter += newton_stepfactor*delta_lambda;

                  rhs_mode = Residual;
                  assembly(true,   // Residual
                           false); // Jacobian
                  rhs->close();
                  nonlinear_residual_afterstep = rhs->l2_norm();

                  if (!quiet)
                    libMesh::out << "At shrink step "
                                 << backtrack_step
                                 << ", nonlinear_residual_afterstep="
                                 << nonlinear_residual_afterstep
                                 << std::endl;

                  if (nonlinear_residual_afterstep < nonlinear_residual_beforestep)
                    {
                      if (!quiet)
                        libMesh::out << "Backtracking succeeded!" << std::endl;

                      break; // out of backtracking loop
                    }

                  else
                    {
                      // Back off that step
                      solution->add(-newton_stepfactor, *delta_u);
                      solution->close();
                      *continuation_parameter -= newton_stepfactor*delta_lambda;
                    }

                  // Save a copy of the solution from before the Newton step.
                  //std::unique_ptr<NumericVector<Number>> prior_iterate = solution->clone();
                }
            } // end if (attempte_backtracking)


          // If we tried backtracking but the residual is still not reduced, print message.
          if ((attempt_backtracking) && (nonlinear_residual_afterstep > nonlinear_residual_beforestep))
            {
              //libMesh::err << "Backtracking failed." << std::endl;
              libMesh::out << "Backtracking failed." << std::endl;

              // 1.) Quit, exit program.
              //libmesh_error_msg("Backtracking failed!");

              // 2.) Continue with last newton_stepfactor
              if (newton_step<3)
                {
                  solution->add(newton_stepfactor, *delta_u);
                  solution->close();
                  *continuation_parameter += newton_stepfactor*delta_lambda;
                  if (!quiet)
                    libMesh::out << "Backtracking could not reduce residual ... continuing anyway!" << std::endl;
                }

              // 3.) Break out of Newton iteration loop with newton_converged = false,
              //     reduce the arclength stepsize, and try again.
              else
                {
                  break; // out of Newton iteration loop, with newton_converged=false
                }
            }

          // Another type of convergence check: suppose the residual has not been reduced
          // from its initial value after half of the allowed Newton steps have occurred.
          // In our experience, this typically means that it isn't going to converge and
          // we could probably save time by dropping out of the Newton iteration loop and
          // trying a smaller arcstep.
          if (this->newton_progress_check)
            {
              if ((nonlinear_residual_afterstep > nonlinear_residual_firststep) &&
                  (newton_step+1 > static_cast<unsigned int>(0.5*newton_solver->max_nonlinear_iterations)))
                {
                  libMesh::out << "Progress check failed: the current residual: "
                               << nonlinear_residual_afterstep
                               << ", is\n"
                               << "larger than the initial residual, and half of the allowed\n"
                               << "number of Newton iterations have elapsed.\n"
                               << "Exiting Newton iterations with converged==false." << std::endl;

                  break; // out of Newton iteration loop, newton_converged = false
                }
            }

          // Safety check: Check the current continuation parameter against user-provided min-allowable parameter value
          if (*continuation_parameter < min_continuation_parameter)
            {
              libMesh::out << "Continuation parameter fell below min-allowable value." << std::endl;
              break; // out of Newton iteration loop, newton_converged = false
            }

          // Safety check: Check the current continuation parameter against user-provided max-allowable parameter value
          if ( (max_continuation_parameter != 0.0) &&
               (*continuation_parameter > max_continuation_parameter) )
            {
              libMesh::out << "Current continuation parameter value: "
                           << *continuation_parameter
                           << " exceeded max-allowable value."
                           << std::endl;
              break; // out of Newton iteration loop, newton_converged = false
            }


          // Check the convergence of the parameter and the solution.  If they are small
          // enough, we can break out of the Newton iteration loop.
          const Real norm_delta_u = delta_u->l2_norm();
          const Real norm_u = solution->l2_norm();
          libMesh::out << "  delta_lambda                   = " << delta_lambda << std::endl;
          libMesh::out << "  newton_stepfactor*delta_lambda = " << newton_stepfactor*delta_lambda << std::endl;
          libMesh::out << "  lambda_current                 = " << *continuation_parameter << std::endl;
          libMesh::out << "  ||delta_u||                    = " << norm_delta_u << std::endl;
          libMesh::out << "  ||delta_u||/||u||              = " << norm_delta_u / norm_u << std::endl;


          // Evaluate the residual at the current Newton iterate.  We don't want to detect
          // convergence due to a small Newton step when the residual is still not small.
          rhs_mode = Residual;
          assembly(true,   // Residual
                   false); // Jacobian
          rhs->close();
          const Real norm_residual = rhs->l2_norm();
          libMesh::out << "  ||R||_{L2} = " << norm_residual << std::endl;
          libMesh::out << "  ||R||_{L2}/||u|| = " << norm_residual / norm_u << std::endl;


          // FIXME: The norm_delta_u tolerance (at least) should be relative.
          // It doesn't make sense to converge a solution whose size is ~ 10^5 to
          // a tolerance of 1.e-6.  Oh, and we should also probably check the
          // (relative) size of the residual as well, instead of just the step.
          if ((std::abs(delta_lambda) < continuation_parameter_tolerance) &&
              //(norm_delta_u       < solution_tolerance)               && // This is a *very* strict criterion we can probably skip
              (norm_residual      < solution_tolerance))
            {
              if (!quiet)
                libMesh::out << "Newton iterations converged!" << std::endl;

              newton_converged = true;
              break; // out of Newton iterations
            }
        } // end nonlinear loop

      if (!newton_converged)
        {
          libMesh::out << "Newton iterations of augmented system did not converge!" << std::endl;

          // Reduce ds_current, recompute the solution and parameter, and continue to next
          // arcstep, if there is one.
          ds_current *= 0.5;

          // Go back to previous solution and parameter value.
          *solution = *previous_u;
          *continuation_parameter = old_continuation_parameter;

          // Compute new predictor with smaller ds
          apply_predictor();
        }
      else
        {
          // Set step convergence and break out
          arcstep_converged=true;
          break; // out of arclength reduction loop
        }

    } // end loop over arclength reductions

  // Check for convergence of the whole arcstep.  If not converged at this
  // point, we have no choice but to quit.
  if (!arcstep_converged)
    libmesh_error_msg("Arcstep failed to converge after max number of reductions! Exiting...");

  // Print converged solution control parameter and max value.
  libMesh::out << "lambda_current=" << *continuation_parameter << std::endl;
  //libMesh::out << "u_max=" << solution->max() << std::endl;

  // Reset old stream precision and flags.
  libMesh::out.precision(old_precision);
  libMesh::out.unsetf(std::ios_base::scientific);

  // Note: we don't want to go on to the next guess yet, since the user may
  // want to post-process this data.  It's up to the user to call advance_arcstep()
  // when they are ready to go on.
}



void ContinuationSystem::advance_arcstep()
{
  // Solve for the updated tangent du1/ds, d(lambda1)/ds
  solve_tangent();

  // Advance the solution and the parameter to the next value.
  update_solution();
}



// This function solves the tangent system:
// [ G_u                G_{lambda}        ][(du/ds)_new      ] = [  0 ]
// [ Theta*(du/ds)_old  (dlambda/ds)_old  ][(dlambda/ds)_new ]   [-N_s]
// The solution is given by:
// .) Let G_u y = G_lambda, then
// .) 2nd row yields:
//    (dlambda/ds)_new = 1.0 / ( (dlambda/ds)_old - Theta*(du/ds)_old*y )
// .) 1st row yields
//    (du_ds)_new = -(dlambda/ds)_new * y
void ContinuationSystem::solve_tangent()
{
  // We shouldn't call this unless the current tangent already makes sense.
  libmesh_assert (tangent_initialized);

  // Set pointer to underlying Newton solver
  if (!newton_solver)
    newton_solver =
      cast_ptr<NewtonSolver *> (this->time_solver->diff_solver().get());

  // Assemble the system matrix AND rhs, with rhs = G_{\lambda}
  this->rhs_mode = G_Lambda;

  // Assemble Residual and Jacobian
  this->assembly(true,   // Residual
                 true); // Jacobian

  // Not sure if this is really necessary
  rhs->close();

  // Solve G_u*y =  G_{\lambda}
  std::pair<unsigned int, Real> rval =
    linear_solver->solve(*matrix,
                         *y,
                         *rhs,
                         1.e-12, // relative linear tolerance
                         2*newton_solver->max_linear_iterations);   // max linear iterations

  // FIXME: If this doesn't converge at all, the new tangent vector is
  // going to be really bad...

  if (!quiet)
    libMesh::out << "G_u*y = G_{lambda} solver converged at step "
                 << rval.first
                 << " linear tolerance = "
                 << rval.second
                 << "."
                 << std::endl;

  // Save old solution and parameter tangents for possible use in higher-order
  // predictor schemes.
  previous_dlambda_ds = dlambda_ds;
  *previous_du_ds     = *du_ds;


  // 1.) Previous, probably wrong, technique!
  //   // Solve for the updated d(lambda)/ds
  //   // denom = N_{lambda}   - (du_ds)^t y
  //   //       = d(lambda)/ds - (du_ds)^t y
  //   Real denom = dlambda_ds - du_ds->dot(*y);

  //   //libMesh::out << "denom=" << denom << std::endl;
  //   libmesh_assert_not_equal_to (denom, 0.0);

  //   dlambda_ds = 1.0 / denom;


  //   if (!quiet)
  //     libMesh::out << "dlambda_ds=" << dlambda_ds << std::endl;

  //   // Compute the updated value of du/ds = -_dlambda_ds * y
  //   du_ds->zero();
  //   du_ds->add(-dlambda_ds, *y);
  //   du_ds->close();


  // 2.) From Brian Carnes' paper...
  // According to Carnes, y comes from solving G_u * y = -G_{\lambda}
  y->scale(-1.);
  const Real ynorm = y->l2_norm();
  dlambda_ds = 1. / std::sqrt(1. + Theta_LOCA*Theta_LOCA*Theta*ynorm*ynorm);

  // Determine the correct sign for dlambda_ds.

  // We will use delta_u to temporarily compute this sign.
  *delta_u = *solution;
  delta_u->add(-1., *previous_u);
  delta_u->close();

  const Real sgn_dlambda_ds =
    libmesh_real(Theta_LOCA*Theta_LOCA*Theta*y->dot(*delta_u) +
                 (*continuation_parameter-old_continuation_parameter));

  if (sgn_dlambda_ds < 0.)
    {
      if (!quiet)
        libMesh::out << "dlambda_ds is negative." << std::endl;

      dlambda_ds *= -1.;
    }

  // Finally, set the new tangent vector, du/ds = dlambda/ds * y.
  du_ds->zero();
  du_ds->add(dlambda_ds, *y);
  du_ds->close();

  if (!quiet)
    {
      libMesh::out << "d(lambda)/ds = " << dlambda_ds << std::endl;
      libMesh::out << "||du_ds||    = " << du_ds->l2_norm() << std::endl;
    }

  // Our next solve expects y ~ -du/dlambda, so scale it back by -1 again now.
  y->scale(-1.);
  y->close();
}



void ContinuationSystem::set_Theta()
{
  // // Use the norm of the latest solution, squared.
  //const Real normu = solution->l2_norm();
  //libmesh_assert_not_equal_to (normu, 0.0);
  //Theta = 1./normu/normu;

  // // 1.) Use the norm of du, squared
  //   *delta_u = *solution;
  //   delta_u->add(-1, *previous_u);
  //   delta_u->close();
  //   const Real normdu = delta_u->l2_norm();

  //   if (normdu < 1.) // don't divide by zero or make a huge scaling parameter.
  //     Theta = 1.;
  //   else
  //     Theta = 1./normdu/normdu;

  // 2.) Use 1.0, i.e. don't scale
  Theta=1.;

  // 3.) Use a formula which attempts to make the "solution triangle" isosceles.
  //   libmesh_assert_less (std::abs(dlambda_ds), 1.);

  //   *delta_u = *solution;
  //   delta_u->add(-1, *previous_u);
  //   delta_u->close();
  //   const Real normdu = delta_u->l2_norm();

  //   Theta = std::sqrt(1. - dlambda_ds*dlambda_ds) / normdu * tau * ds;


  //   // 4.) Use the norm of du and the norm of du/ds
  //   *delta_u = *solution;
  //   delta_u->add(-1, *previous_u);
  //   delta_u->close();
  //   const Real normdu   = delta_u->l2_norm();
  //   du_ds->close();
  //   const Real normduds = du_ds->l2_norm();

  //   if (normduds < 1.e-12)
  //     {
  //       libMesh::out << "Setting initial Theta= 1./normdu/normdu" << std::endl;
  //       libMesh::out << "normdu=" << normdu << std::endl;

  //       // Don't use this scaling if the solution delta is already O(1)
  //       if (normdu > 1.)
  // Theta = 1./normdu/normdu;
  //       else
  // Theta = 1.;
  //     }
  //   else
  //     {
  //       libMesh::out << "Setting Theta= 1./normdu/normduds" << std::endl;
  //       libMesh::out << "normdu=" << normdu << std::endl;
  //       libMesh::out << "normduds=" << normduds << std::endl;

  //       // Don't use this scaling if the solution delta is already O(1)
  //       if ((normdu>1.) || (normduds>1.))
  // Theta = 1./normdu/normduds;
  //       else
  // Theta = 1.;
  //     }

  if (!quiet)
    libMesh::out << "Setting Normalization Parameter Theta=" << Theta << std::endl;
}



void ContinuationSystem::set_Theta_LOCA()
{
  // We also recompute the LOCA normalization parameter based on the
  // most recently computed value of dlambda_ds
  // if (!quiet)
  //   libMesh::out << "(Theta_LOCA) dlambda_ds=" << dlambda_ds << std::endl;

  // Formula makes no sense if |dlambda_ds| > 1
  libmesh_assert_less (std::abs(dlambda_ds), 1.);

  // 1.) Attempt to implement the method in LOCA paper
  //   const Real g = 1./std::sqrt(2.); // "desired" dlambda_ds

  //   // According to the LOCA people, we only renormalize for
  //   // when |dlambda_ds| exceeds some pre-selected maximum (which they take to be zero, btw).
  //   if (std::abs(dlambda_ds) > .9)
  //     {
  //       // Note the *= ... This is updating the previous value of Theta_LOCA
  //       // Note: The LOCA people actually use Theta_LOCA^2 to normalize their arclength constraint.
  //       Theta_LOCA *= std::abs( (dlambda_ds/g)*std::sqrt( (1.-g*g) / (1.-dlambda_ds*dlambda_ds) ) );

  //       // Suggested max-allowable value for Theta_LOCA
  //       if (Theta_LOCA > 1.e8)
  // {
  //   Theta_LOCA = 1.e8;

  //   if (!quiet)
  //     libMesh::out << "max Theta_LOCA=" << Theta_LOCA << " has been selected." << std::endl;
  // }
  //     }
  //   else
  //     Theta_LOCA=1.0;

  // 2.) FIXME: Should we do *= or just =?  This function is of dlambda_ds is
  //  < 1,  |dlambda_ds| < 1/sqrt(2) ~~ .7071
  //  > 1,  |dlambda_ds| > 1/sqrt(2) ~~ .7071
  Theta_LOCA *= std::abs( dlambda_ds / std::sqrt( (1.-dlambda_ds*dlambda_ds) ) );

  // Suggested max-allowable value for Theta_LOCA.  I've never come close
  // to this value in my code.
  if (Theta_LOCA > 1.e8)
    {
      Theta_LOCA = 1.e8;

      if (!quiet)
        libMesh::out << "max Theta_LOCA=" << Theta_LOCA << " has been selected." << std::endl;
    }

  // 3.) Use 1.0, i.e. don't scale
  //Theta_LOCA=1.0;

  if (!quiet)
    libMesh::out << "Setting Theta_LOCA=" << Theta_LOCA << std::endl;
}



void ContinuationSystem::update_solution()
{
  // Set some stream formatting flags
  std::streamsize old_precision = libMesh::out.precision();
  libMesh::out.precision(16);
  libMesh::out.setf(std::ios_base::scientific);

  // We must have a tangent that makes sense before we can update the solution.
  libmesh_assert (tangent_initialized);

  // Compute tau, the stepsize scaling parameter which attempts to
  // reduce ds when the angle between the most recent two tangent
  // vectors becomes large.  tau is actually the (absolute value of
  // the) cosine of the angle between these two vectors... so if tau ~
  // 0 the angle is ~ 90 degrees, while if tau ~ 1 the angle is ~ 0
  // degrees.
  y_old->close();
  y->close();
  const Real yoldnorm = y_old->l2_norm();
  const Real ynorm = y->l2_norm();
  const Number yoldy = y_old->dot(*y);
  const Real yold_over_y = yoldnorm/ynorm;

  if (!quiet)
    {
      libMesh::out << "yoldnorm=" << yoldnorm << std::endl;
      libMesh::out << "ynorm="    << ynorm << std::endl;
      libMesh::out << "yoldy="    << yoldy << std::endl;
      libMesh::out << "yoldnorm/ynorm=" << yoldnorm/ynorm << std::endl;
    }

  // Save the current value of ds before updating it
  previous_ds = ds_current;

  // // 1.) Cosine method (for some reason this always predicts the angle is ~0)
  // // Don't try dividing by zero
  // if ((yoldnorm > 1.e-12) && (ynorm > 1.e-12))
  //   tau = std::abs(yoldy) / yoldnorm  / ynorm;
  // else
  //   tau = 1.;

  // // 2.) Relative size of old and new du/dlambda method with cutoff of 0.9
  // if ((yold_over_y < 0.9) && (yold_over_y > 1.e-6))
  //   tau = yold_over_y;
  // else
  //   tau = 1.;

  // 3.) Grow (or shrink) the arclength stepsize by the ratio of du/dlambda, but do not
  // exceed the user-specified value of ds.
  if (yold_over_y > 1.e-6)
    {
      // // 1.) Scale current ds by the ratio of successive tangents.
      //       ds_current *= yold_over_y;
      //       if (ds_current > ds)
      // ds_current = ds;

      // 2.) Technique 1 tends to shrink the step fairly well (and even if it doesn't
      // get very small, we still have step reduction) but it seems to grow the step
      // very slowly.  Another possible technique is step-doubling:
      //       if (yold_over_y > 1.)
      //       ds_current *= 2.;
      //       else
      // ds_current *= yold_over_y;

      // 3.) Technique 2 may over-zealous when we are also using the Newton stepgrowth
      // factor.  For technique 3 we multiply by yold_over_y unless yold_over_y > 2
      // in which case we use 2.
      //       if (yold_over_y > 2.)
      // ds_current *= 2.;
      //       else
      // ds_current *= yold_over_y;

      // 4.) Double-or-halve.  We double the arc-step if the ratio of successive tangents
      // is larger than 'double_threshold', halve it if it is less than 'halve_threshold'
      const Real double_threshold = 0.5;
      const Real halve_threshold  = 0.5;
      if (yold_over_y > double_threshold)
        ds_current *= 2.;
      else if (yold_over_y < halve_threshold)
        ds_current *= 0.5;


      // Also possibly use the number of Newton iterations required to compute the previous
      // step (relative to the maximum-allowed number of Newton iterations) to grow the step.
      if (newton_stepgrowth_aggressiveness > 0.)
        {
          libmesh_assert(newton_solver);
          const unsigned int Nmax = newton_solver->max_nonlinear_iterations;

          // // The LOCA Newton step growth technique (note: only grows step length)
          // const Real stepratio = static_cast<Real>(Nmax-(newton_step+1))/static_cast<Real>(Nmax-1.);
          // const Real newtonstep_growthfactor = 1. + newton_stepgrowth_aggressiveness*stepratio*stepratio;

          // The "Nopt/N" method, may grow or shrink the step.  Assume Nopt=Nmax/2.
          const Real newtonstep_growthfactor =
            newton_stepgrowth_aggressiveness * 0.5 *
            static_cast<Real>(Nmax) / static_cast<Real>(newton_step+1);

          if (!quiet)
            libMesh::out << "newtonstep_growthfactor=" << newtonstep_growthfactor << std::endl;

          ds_current *= newtonstep_growthfactor;
        }
    }


  // Don't let the stepsize get above the user's maximum-allowed stepsize.
  if (ds_current > ds)
    ds_current = ds;

  // Check also for a minimum allowed stepsize.
  if (ds_current < ds_min)
    {
      libMesh::out << "Enforcing minimum-allowed arclength stepsize of " << ds_min << std::endl;
      ds_current = ds_min;
    }

  if (!quiet)
    {
      libMesh::out << "Current step size: ds_current=" << ds_current << std::endl;
    }

  // Recompute scaling factor Theta for
  // the current solution before updating.
  set_Theta();

  // Also, recompute the LOCA scaling factor, which attempts to
  // maintain a reasonable value of dlambda/ds
  set_Theta_LOCA();

  libMesh::out << "Theta*Theta_LOCA^2=" << Theta*Theta_LOCA*Theta_LOCA << std::endl;

  // Based on the asymptotic singular behavior of du/dlambda near simple turning points,
  // we can compute a single parameter which may suggest that we are close to a singularity.
  *delta_u = *solution;
  delta_u->add(-1, *previous_u);
  delta_u->close();
  const Real normdu   = delta_u->l2_norm();
  const Real C = (std::log (Theta_LOCA*normdu) /
                  std::log (std::abs(*continuation_parameter-old_continuation_parameter))) - 1.0;
  if (!quiet)
    libMesh::out << "C=" << C << std::endl;

  // Save the current value of u and lambda before updating.
  save_current_solution();

  if (!quiet)
    {
      libMesh::out << "Updating the solution with the tangent guess." << std::endl;
      libMesh::out << "||u_old||=" << this->solution->l2_norm() << std::endl;
      libMesh::out << "lambda_old=" << *continuation_parameter << std::endl;
    }

  // Since we solved for the tangent vector, now we can compute an
  // initial guess for the new solution, and an initial guess for the
  // new value of lambda.
  apply_predictor();

  if (!quiet)
    {
      libMesh::out << "||u_new||=" << this->solution->l2_norm() << std::endl;
      libMesh::out << "lambda_new=" << *continuation_parameter << std::endl;
    }

  // Unset previous stream flags
  libMesh::out.precision(old_precision);
  libMesh::out.unsetf(std::ios_base::scientific);
}




void ContinuationSystem::save_current_solution()
{
  // Save the old solution vector
  *previous_u = *solution;

  // Save the old value of lambda
  old_continuation_parameter = *continuation_parameter;
}



void ContinuationSystem::apply_predictor()
{
  if (predictor == Euler)
    {
      // 1.) Euler Predictor
      // Predict next the solution
      solution->add(ds_current, *du_ds);
      solution->close();

      // Predict next parameter value
      *continuation_parameter += ds_current*dlambda_ds;
    }


  else if (predictor == AB2)
    {
      // 2.) 2nd-order explicit AB predictor
      libmesh_assert_not_equal_to (previous_ds, 0.0);
      const Real stepratio = ds_current/previous_ds;

      // Build up next solution value.
      solution->add( 0.5*ds_current*(2.+stepratio), *du_ds);
      solution->add(-0.5*ds_current*stepratio     , *previous_du_ds);
      solution->close();

      // Next parameter value
      *continuation_parameter +=
        0.5*ds_current*((2.+stepratio)*dlambda_ds -
                        stepratio*previous_dlambda_ds);
    }

  else
    libmesh_error_msg("Unknown predictor!");
}

} // namespace libMesh
