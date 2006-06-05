
#include "diff_system.h"
#include "equation_systems.h"
#include "linear_solver.h"
#include "newton_solver.h"
#include "numeric_vector.h"
#include "sparse_matrix.h"


NewtonSolver::NewtonSolver (sys_type& s) :
  Parent(s), linear_solver(LinearSolver<Number>::build())
{
}



NewtonSolver::~NewtonSolver ()
{
}



void NewtonSolver::solve()
{
// The number of steps and the stopping criterion
// for the nonlinear iterations.
const unsigned int max_nonlinear_steps = 100;
const unsigned int max_linear_iterations = 10000;

// Stopping criteria for nonlinear iterations
const Real nonlinear_abs_step_tolerance = 1.e-9;
const Real nonlinear_rel_step_tolerance = 1.e-6;
const Real nonlinear_abs_res_tolerance = 1.e-9;
const Real nonlinear_rel_res_tolerance = 1.e-6;
// Maximum amount by which to reduce Newton steps
const Real minsteplength = 0.1;
// Initial linear solver tolerance in main nonlinear solver
const Real linear_tolerance = 1.e-3;
// Amount by which nonlinear residual should exceed linear solver
// tolerance
const Real relative_tolerance = 1.e-3;
// We'll shut up eventually...
const bool verbose_convergence_chatter = true;
const bool require_residual_reduction = false;


  EquationSystems &equation_systems = _system.get_equation_systems();

  equation_systems.parameters.set<unsigned int>
                   ("linear solver maximum iterations") =
                   max_linear_iterations;

  NumericVector<Number> &solution = *(_system.solution);
  NumericVector<Number> &newton_iterate =
    _system.get_vector("_nonlinear_solution");
  NumericVector<Number> &rhs = *(_system.rhs);

  SparseMatrix<Number> &matrix = *(_system.matrix);

  // Now we begin the nonlinear loop
  newton_iterate.zero();
  newton_iterate.close();

  // Prepare to take incomplete steps
  Real last_residual=0., first_residual=0.;

  // Set starting linear tolerance
  Real current_linear_tolerance = linear_tolerance;

  for (unsigned int l=0; l<max_nonlinear_steps; ++l)
    {
      _system.assembly(true, true);
      rhs.close();
      Real current_residual = rhs.l2_norm();
      if (!l)
        first_residual = current_residual;
      last_residual = current_residual;
std::cout << "Nonlinear Residual: " << current_residual << std::endl;

      // Make sure our linear tolerance is low enough
      if (current_linear_tolerance > current_residual * relative_tolerance)
        {
          current_linear_tolerance = current_residual * relative_tolerance;
          equation_systems.parameters.set<Real> 
            ("linear solver tolerance") = current_linear_tolerance;
        }

      // At this point newton_iterate is the current guess, and
      // solution is now about to become the NEGATIVE of the next
      // Newton step

std::cout << "Linear solve starting" << std::endl;

      // Solve the linear system.  Two cases:
      const std::pair<unsigned int, Real> rval =
        (_system.have_matrix("Preconditioner")) ?
      // 1.) User-supplied preconditioner
        linear_solver->solve (matrix, _system.get_matrix("Preconditioner"),
                              solution, rhs, current_linear_tolerance,
                              max_linear_iterations) :
      // 2.) Use system matrix for the preconditioner
        linear_solver->solve (matrix, solution, rhs,
                              current_linear_tolerance, 
                              max_linear_iterations);

std::cout << "Linear solve converged, step " << rval.first
          << ", residual " << rval.second
          << ", tolerance " << current_linear_tolerance << std::endl;

      // Compute the l2 norm of the nonlinear update
      Real norm_delta = solution.l2_norm();

std::cout << "Taking full Newton step" << std::endl;
          // Take a full Newton step
      newton_iterate.add (-1., solution);
      newton_iterate.close();

      // Check residual with full Newton step
      Real steplength = 1.;
      _system.assembly(true, false);

      // backtrack if necessary
      if (require_residual_reduction)
        {
          rhs.close();
          current_residual = rhs.l2_norm();
          while (current_residual > last_residual)
            {
              // Reduce step size to 1/2, 1/4, etc.
              steplength /= 2.;
              norm_delta /= 2.;
std::cout << "Shrinking Newton step to " << steplength << std::endl;
              newton_iterate.add (steplength, solution);
              newton_iterate.close();

              // Check residual with fractional Newton step
              std::cout << "          Checking " << std::flush;
              _system.assembly (true, false);
              current_residual = rhs.l2_norm();
              std::cout << "Current Residual: " << current_residual << std::endl;

              if (steplength/2. < minsteplength && 
                  current_residual > last_residual)
                {
std::cout << "Inexact Newton step FAILED at step " << l << std::endl;

                  error();
                }
            }
        }

      // Compute the l2 norm of the whole solution
      const Real norm_total = newton_iterate.l2_norm();

      // Print out information for the 
      // nonlinear iterations.
      std::cout << "    Nonlinear step: |du|/|u| = "
                << norm_delta / norm_total
                << ", |du| = "
                << norm_delta
                << std::endl;

      // We haven't converged unless we pass a convergence test
      bool has_converged = false;

      // Is our absolute residual low enough?
      if (current_residual < nonlinear_abs_res_tolerance)
        {
          std::cout << "  Nonlinear solver converged, step "
                    << l
                    << ", residual "
                    << current_residual
                    << std::endl;
          has_converged = true;
        }
      else if (verbose_convergence_chatter)
        {
          std::cout << "  Nonlinear solver current_residual "
                    << current_residual
                    << " > "
                    << (nonlinear_abs_res_tolerance)
                    << std::endl;
        }

      // Is our relative residual low enough?
      if ((current_residual / first_residual) < nonlinear_rel_res_tolerance)
        {
          std::cout << "  Nonlinear solver converged, step "
                    << l
                    << ", residual reduction "
                    << current_residual / first_residual
                    << std::endl;
          has_converged = true;
        }
      else if (verbose_convergence_chatter)
        {
          std::cout << "  Nonlinear solver current/first residual "
                    << (current_residual / first_residual)
                    << " > "
                    << nonlinear_rel_res_tolerance
                    << std::endl;
        }

      // Is our absolute Newton step size small enough?
      if (l != 0 && (norm_delta / steplength <
          nonlinear_abs_step_tolerance))
        {
          std::cout << "  Nonlinear solver converged, step "
                    << l
                    << ", absolute step size "
                    << norm_delta / steplength
                    << std::endl;
          has_converged = true;
        }

      // Is our relative Newton step size small enough?
      if (l != 0 && (norm_delta / steplength / norm_total <
          nonlinear_rel_step_tolerance))
        {
          std::cout << "  Nonlinear solver converged, step "
                    << l
                    << ", relative step size "
                    << (norm_delta / steplength / norm_total)
                    << std::endl;
          has_converged = true;
        }

      // Terminate the solution iteration if the difference between
      // this iteration and the last is sufficiently small.
      if (has_converged)
        {
          break;
        }
      if (l >= max_nonlinear_steps - 1)
        {
          std::cout << "  Nonlinear solver DIVERGED at step "
                    << l
                    << std::endl;
          error();
          continue;
        }
    } // end nonlinear loop
}
