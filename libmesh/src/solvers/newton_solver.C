
#include "math.h"

#include "diff_system.h"
#include "libmesh_logging.h"
#include "linear_solver.h"
#include "newton_solver.h"
#include "numeric_vector.h"
#include "sparse_matrix.h"
#include "dof_map.h"


NewtonSolver::NewtonSolver (sys_type& s)
  : Parent(s),
    require_residual_reduction(false),
    minsteplength(0.0),
    linear_solver(LinearSolver<Number>::build())
{
}



NewtonSolver::~NewtonSolver ()
{
}



void NewtonSolver::reinit()
{
  Parent::reinit();

  linear_solver->clear();
}



void NewtonSolver::solve()
{
  START_LOG("solve()", "NewtonSolver");

// Amount by which nonlinear residual should exceed linear solver
// tolerance
const Real relative_tolerance = 1.e-3;


  NumericVector<Number> &solution = *(_system.solution);
  NumericVector<Number> &newton_iterate =
    _system.get_vector("_nonlinear_solution");
  newton_iterate.close();
//  solution.close();

  NumericVector<Number> &rhs = *(_system.rhs);

  SparseMatrix<Number> &matrix = *(_system.matrix);

  // Prepare to take incomplete steps
  Real last_residual=0.;

  // Set starting linear tolerance
  Real current_linear_tolerance = initial_linear_tolerance;

  // Now we begin the nonlinear loop
  for (unsigned int l=0; l<max_nonlinear_iterations; ++l)
    {
      if (!quiet)
        std::cout << "Assembling System" << std::endl;

      PAUSE_LOG("solve()", "NewtonSolver");
      _system.assembly(true, true);
      RESTART_LOG("solve()", "NewtonSolver");
      rhs.close();
      Real current_residual = rhs.l2_norm();
      last_residual = current_residual;

      if (isnan(current_residual))
        {
          std::cout << "  Nonlinear solver DIVERGED at step " << l
                    << " with norm Not-a-Number"
                    << std::endl;
          error();
          continue;
        }

      max_residual_norm = std::max (current_residual,
                                    max_residual_norm);
 
      // Compute the l2 norm of the whole solution
      Real norm_total = newton_iterate.l2_norm();

      max_solution_norm = std::max(max_solution_norm, norm_total);

      if (!quiet)
        std::cout << "Nonlinear Residual: "
                  << current_residual << std::endl;

      // Make sure our linear tolerance is low enough
      if (current_linear_tolerance > current_residual * relative_tolerance)
        {
          current_linear_tolerance = current_residual * relative_tolerance;
        }

      // But don't let it be zero
      if (current_linear_tolerance == 0.)
        {
          current_linear_tolerance = TOLERANCE * TOLERANCE;
        }

      // At this point newton_iterate is the current guess, and
      // solution is now about to become the NEGATIVE of the next
      // Newton step.

      // Our best initial guess for the solution is zero!
      solution.zero();

      if (!quiet)
        std::cout << "Linear solve starting" << std::endl;

      PAUSE_LOG("solve()", "NewtonSolver");
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
      // We may need to localize a parallel solution
      _system.update ();
      RESTART_LOG("solve()", "NewtonSolver");
      // The linear solver may not have fit our constraints exactly
      _system.get_dof_map().enforce_constraints_exactly(_system);

      if (!quiet)
        std::cout << "Linear solve finished, step " << rval.first
                  << ", residual " << rval.second
                  << ", tolerance " << current_linear_tolerance
                  << std::endl;

      // Compute the l2 norm of the nonlinear update
      Real norm_delta = solution.l2_norm();

      if (!quiet)
        std::cout << "Trying full Newton step" << std::endl;
      // Take a full Newton step
      newton_iterate.add (-1., solution);
//      newton_iterate.close();

      // Check residual with full Newton step
      Real steplength = 1.;
      PAUSE_LOG("solve()", "NewtonSolver");
      _system.assembly(true, false);
      RESTART_LOG("solve()", "NewtonSolver");

      rhs.close();
      current_residual = rhs.l2_norm();

      // backtrack if necessary
      if (require_residual_reduction)
        {
          // but don't fiddle around if we've already converged
          if (test_convergence(current_residual, norm_delta))
            {
              if (!quiet)
		print_convergence(l, current_residual, norm_delta);
              break;
            }

          while (current_residual > last_residual)
            {
              // Reduce step size to 1/2, 1/4, etc.
              steplength /= 2.;
              norm_delta /= 2.;
              if (!quiet)
                std::cout << "Shrinking Newton step to "
                          << steplength << std::endl;
              newton_iterate.add (steplength, solution);
//              newton_iterate.close();

              // Check residual with fractional Newton step
              PAUSE_LOG("solve()", "NewtonSolver");
              _system.assembly (true, false);
              RESTART_LOG("solve()", "NewtonSolver");

              rhs.close();
              current_residual = rhs.l2_norm();
              if (!quiet)
                std::cout << "Current Residual: "
                          << current_residual << std::endl;

              if (steplength/2. < minsteplength && 
                  current_residual > last_residual)
                {
                  std::cout << "Inexact Newton step FAILED at step "
                            << l << std::endl;

                  error();
                }
            }
        }

      // Compute the l2 norm of the whole solution
      norm_total = newton_iterate.l2_norm();

      max_solution_norm = std::max(max_solution_norm, norm_total);

      // Print out information for the 
      // nonlinear iterations.
      if (!quiet)
        std::cout << "  Nonlinear step: |du|/|u| = "
                  << norm_delta / norm_total
                  << ", |du| = " << norm_delta
                  << std::endl;

      // Terminate the solution iteration if the difference between
      // this iteration and the last is sufficiently small.
      if (!quiet)
        print_convergence(l, current_residual,
                          norm_delta / steplength);
      if (test_convergence(current_residual, norm_delta / steplength))
        {
          break;
        }
      if (l >= max_nonlinear_iterations - 1)
        {
          std::cout << "  Nonlinear solver FAILED TO CONVERGE by step " << l
                    << " with norm " << norm_total
                    << std::endl;
          error();
          continue;
        }
    } // end nonlinear loop

  // Copy the final nonlinear iterate into the current_solution,
  // for other libMesh functions that expect it

  solution = newton_iterate;
//  solution.close();

  // We may need to localize a parallel solution
  _system.update ();

  STOP_LOG("solve()", "NewtonSolver");
}



bool NewtonSolver::test_convergence(Real current_residual,
                                    Real step_norm)
{
  // We haven't converged unless we pass a convergence test
  bool has_converged = false;

  // Is our absolute residual low enough?
  if (current_residual < absolute_residual_tolerance)
    has_converged = true;

  // Is our relative residual low enough?
  if ((current_residual / max_residual_norm) <
      relative_residual_tolerance)
    has_converged = true;

  // Is our absolute Newton step size small enough?
  if (step_norm < absolute_step_tolerance)
    has_converged = true;

  // Is our relative Newton step size small enough?
  if (step_norm / max_solution_norm <
      relative_step_tolerance)
    has_converged = true;

  return has_converged;
}


void NewtonSolver::print_convergence(unsigned int step_num,
                                     Real current_residual,
                                     Real step_norm)
{
  // Is our absolute residual low enough?
  if (current_residual < absolute_residual_tolerance)
    {
      std::cout << "  Nonlinear solver converged, step " << step_num
                << ", residual " << current_residual
                << std::endl;
    }
  else if (absolute_residual_tolerance)
    {
      std::cout << "  Nonlinear solver current_residual "
                << current_residual << " > "
                << (absolute_residual_tolerance) << std::endl;
    }

  // Is our relative residual low enough?
  if ((current_residual / max_residual_norm) <
      relative_residual_tolerance)
    {
      std::cout << "  Nonlinear solver converged, step " << step_num
                << ", residual reduction "
                << current_residual / max_residual_norm
                << " < " << relative_residual_tolerance
                << std::endl;
    }
  else if (relative_residual_tolerance)
    {
      if (!quiet)
        std::cout << "  Nonlinear solver relative residual "
                  << (current_residual / max_residual_norm)
                  << " > " << relative_residual_tolerance
                  << std::endl;
    }

  // Is our absolute Newton step size small enough?
  if (step_norm < absolute_step_tolerance)
    {
      std::cout << "  Nonlinear solver converged, step " << step_num
                << ", absolute step size "
                << step_norm
                << " < " << absolute_step_tolerance
                << std::endl;
    }
  else if (absolute_step_tolerance)
    {
      std::cout << "  Nonlinear solver absolute step size "
                << step_norm
                << " > " << absolute_step_tolerance
                << std::endl;
    }

  // Is our relative Newton step size small enough?
  if (step_norm / max_solution_norm <
      relative_step_tolerance)
    {
      std::cout << "  Nonlinear solver converged, step " << step_num
                << ", relative step size "
                << (step_norm / max_solution_norm)
                << " < " << relative_step_tolerance
                << std::endl;
    }
  else if (relative_step_tolerance)
    {
      std::cout << "  Nonlinear solver relative step size "
                << (step_norm / max_solution_norm)
                << " > " << relative_step_tolerance
                << std::endl;
    }
}
