// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


#include "libmesh/diff_system.h"
#include "libmesh/dof_map.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/linear_solver.h"
#include "libmesh/newton_solver.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/sparse_matrix.h"

namespace libMesh
{

// SIGN from Numerical Recipes
template <typename T>
inline
T SIGN(T a, T b)
{
  return b >= 0 ? std::abs(a) : -std::abs(a);
}

Real NewtonSolver::line_search(Real tol,
                               Real last_residual,
                               Real & current_residual,
                               NumericVector<Number> & newton_iterate,
                               const NumericVector<Number> & linear_solution)
{
  // Take a full step if we got a residual reduction or if we
  // aren't substepping
  if ((current_residual < last_residual) ||
      (!require_residual_reduction &&
       (!require_finite_residual || !libmesh_isnan(current_residual))))
    return 1.;

  // The residual vector
  NumericVector<Number> & rhs = *(_system.rhs);

  Real ax = 0.;  // First abscissa, don't take negative steps
  Real cx = 1.;  // Second abscissa, don't extrapolate steps

  // Find bx, a step length that gives lower residual than ax or cx
  Real bx = 1.;

  while (libmesh_isnan(current_residual) ||
         (current_residual > last_residual &&
          require_residual_reduction))
    {
      // Reduce step size to 1/2, 1/4, etc.
      Real substepdivision;
      if (brent_line_search && !libmesh_isnan(current_residual))
        {
          substepdivision = std::min(0.5, last_residual/current_residual);
          substepdivision = std::max(substepdivision, tol*2.);
        }
      else
        substepdivision = 0.5;

      newton_iterate.add (bx * (1.-substepdivision),
                          linear_solution);
      newton_iterate.close();
      bx *= substepdivision;
      if (verbose)
        libMesh::out << "  Shrinking Newton step to "
                     << bx << std::endl;

      // Check residual with fractional Newton step
      _system.assembly (true, false);

      rhs.close();
      current_residual = rhs.l2_norm();
      if (verbose)
        libMesh::out << "  Current Residual: "
                     << current_residual << std::endl;

      if (bx/2. < minsteplength &&
          (libmesh_isnan(current_residual) ||
           (current_residual > last_residual)))
        {
          libMesh::out << "Inexact Newton step FAILED at step "
                       << _outer_iterations << std::endl;

          if (!continue_after_backtrack_failure)
            {
              libmesh_convergence_failure();
            }
          else
            {
              libMesh::out << "Continuing anyway ..." << std::endl;
              _solve_result = DiffSolver::DIVERGED_BACKTRACKING_FAILURE;
              return bx;
            }
        }
    } // end while (current_residual > last_residual)

  // Now return that reduced-residual step, or  use Brent's method to
  // find a more optimal step.

  if (!brent_line_search)
    return bx;

  // Brent's method adapted from Numerical Recipes in C, ch. 10.2
  Real e = 0.;

  Real x = bx, w = bx, v = bx;

  // Residuals at bx
  Real fx = current_residual,
    fw = current_residual,
    fv = current_residual;

  // Max iterations for Brent's method loop
  const unsigned int max_i = 20;

  // for golden ratio steps
  const Real golden_ratio = 1.-(std::sqrt(5.)-1.)/2.;

  for (unsigned int i=1; i <= max_i; i++)
    {
      Real xm = (ax+cx)*0.5;
      Real tol1 = tol * std::abs(x) + tol*tol;
      Real tol2 = 2.0 * tol1;

      // Test if we're done
      if (std::abs(x-xm) <= (tol2 - 0.5 * (cx - ax)))
        return x;

      Real d;

      // Construct a parabolic fit
      if (std::abs(e) > tol1)
        {
          Real r = (x-w)*(fx-fv);
          Real q = (x-v)*(fx-fw);
          Real p = (x-v)*q-(x-w)*r;
          q = 2. * (q-r);
          if (q > 0.)
            p = -p;
          else
            q = std::abs(q);
          if (std::abs(p) >= std::abs(0.5*q*e) ||
              p <= q * (ax-x) ||
              p >= q * (cx-x))
            {
              // Take a golden section step
              e = x >= xm ? ax-x : cx-x;
              d = golden_ratio * e;
            }
          else
            {
              // Take a parabolic fit step
              d = p/q;
              if (x+d-ax < tol2 || cx-(x+d) < tol2)
                d = SIGN(tol1, xm - x);
            }
        }
      else
        {
          // Take a golden section step
          e = x >= xm ? ax-x : cx-x;
          d = golden_ratio * e;
        }

      Real u = std::abs(d) >= tol1 ? x+d : x + SIGN(tol1,d);

      // Assemble the residual at the new steplength u
      newton_iterate.add (bx - u, linear_solution);
      newton_iterate.close();
      bx = u;
      if (verbose)
        libMesh::out << "  Shrinking Newton step to "
                     << bx << std::endl;

      _system.assembly (true, false);

      rhs.close();
      Real fu = current_residual = rhs.l2_norm();
      if (verbose)
        libMesh::out << "  Current Residual: "
                     << fu << std::endl;

      if (fu <= fx)
        {
          if (u >= x)
            ax = x;
          else
            cx = x;
          v = w;   w = x;   x = u;
          fv = fw; fw = fx; fx = fu;
        }
      else
        {
          if (u < x)
            ax = u;
          else
            cx = u;
          if (fu <= fw || w == x)
            {
              v = w;   w = u;
              fv = fw; fw = fu;
            }
          else if (fu <= fv || v == x || v == w)
            {
              v = u;
              fv = fu;
            }
        }
    }

  if (!quiet)
    libMesh::out << "Warning!  Too many iterations used in Brent line search!"
                 << std::endl;
  return bx;
}


NewtonSolver::NewtonSolver (sys_type & s)
  : Parent(s),
    require_residual_reduction(true),
    require_finite_residual(true),
    brent_line_search(true),
    track_linear_convergence(false),
    minsteplength(1e-5),
    linear_tolerance_multiplier(1e-3),
    _linear_solver(LinearSolver<Number>::build(s.comm()))
{
}



NewtonSolver::~NewtonSolver ()
{
}



void NewtonSolver::init()
{
  Parent::init();

  if (libMesh::on_command_line("--solver_system_names"))
    _linear_solver->init((_system.name()+"_").c_str());
  else
    _linear_solver->init();

  _linear_solver->init_names(_system);
}



void NewtonSolver::reinit()
{
  Parent::reinit();

  _linear_solver->clear();

  _linear_solver->init_names(_system);
}



unsigned int NewtonSolver::solve()
{
  LOG_SCOPE("solve()", "NewtonSolver");

  // Reset any prior solve result
  _solve_result = INVALID_SOLVE_RESULT;

  NumericVector<Number> & newton_iterate = *(_system.solution);

  UniquePtr<NumericVector<Number> > linear_solution_ptr = newton_iterate.zero_clone();
  NumericVector<Number> & linear_solution = *linear_solution_ptr;
  NumericVector<Number> & rhs = *(_system.rhs);

  newton_iterate.close();
  linear_solution.close();
  rhs.close();

#ifdef LIBMESH_ENABLE_CONSTRAINTS
  _system.get_dof_map().enforce_constraints_exactly(_system);
#endif

  SparseMatrix<Number> & matrix = *(_system.matrix);

  // Set starting linear tolerance
  Real current_linear_tolerance = initial_linear_tolerance;

  // Start counting our linear solver steps
  _inner_iterations = 0;

  // Now we begin the nonlinear loop
  for (_outer_iterations=0; _outer_iterations<max_nonlinear_iterations;
       ++_outer_iterations)
    {
      if (verbose)
        libMesh::out << "Assembling the System" << std::endl;

      _system.assembly(true, true);
      rhs.close();
      Real current_residual = rhs.l2_norm();

      if (libmesh_isnan(current_residual))
        {
          libMesh::out << "  Nonlinear solver DIVERGED at step "
                       << _outer_iterations
                       << " with residual Not-a-Number"
                       << std::endl;
          libmesh_convergence_failure();
          continue;
        }

      if (current_residual <= absolute_residual_tolerance)
        {
          if (verbose)
            libMesh::out << "Linear solve unnecessary; residual "
                         << current_residual
                         << " meets absolute tolerance "
                         << absolute_residual_tolerance
                         << std::endl;

          // We're not doing a solve, but other code may reuse this
          // matrix.
          matrix.close();

          _solve_result |= CONVERGED_ABSOLUTE_RESIDUAL;
          if (current_residual == 0)
            {
              if (relative_residual_tolerance > 0)
                _solve_result |= CONVERGED_RELATIVE_RESIDUAL;
              if (absolute_step_tolerance > 0)
                _solve_result |= CONVERGED_ABSOLUTE_STEP;
              if (relative_step_tolerance > 0)
                _solve_result |= CONVERGED_RELATIVE_STEP;
            }

          break;
        }

      // Prepare to take incomplete steps
      Real last_residual = current_residual;

      max_residual_norm = std::max (current_residual,
                                    max_residual_norm);

      // Compute the l2 norm of the whole solution
      Real norm_total = newton_iterate.l2_norm();

      max_solution_norm = std::max(max_solution_norm, norm_total);

      if (verbose)
        libMesh::out << "Nonlinear Residual: "
                     << current_residual << std::endl;

      // Make sure our linear tolerance is low enough
      current_linear_tolerance = std::min (current_linear_tolerance,
                                           current_residual * linear_tolerance_multiplier);

      // But don't let it be too small
      if (current_linear_tolerance < minimum_linear_tolerance)
        {
          current_linear_tolerance = minimum_linear_tolerance;
        }

      // If starting the nonlinear solve with a really good initial guess, we dont want to set an absurd linear tolerance
      current_linear_tolerance = std::max(current_linear_tolerance, absolute_residual_tolerance / current_residual / 10.0);

      // At this point newton_iterate is the current guess, and
      // linear_solution is now about to become the NEGATIVE of the next
      // Newton step.

      // Our best initial guess for the linear_solution is zero!
      linear_solution.zero();

      if (verbose)
        libMesh::out << "Linear solve starting, tolerance "
                     << current_linear_tolerance << std::endl;

      // Solve the linear system.
      const std::pair<unsigned int, Real> rval =
        _linear_solver->solve (matrix, _system.request_matrix("Preconditioner"),
                               linear_solution, rhs, current_linear_tolerance,
                               max_linear_iterations);

      if (track_linear_convergence)
        {
          LinearConvergenceReason linear_c_reason = _linear_solver->get_converged_reason();

          // Check if something went wrong during the linear solve
          if (linear_c_reason < 0)
            {
              // The linear solver failed somehow
              _solve_result |= DiffSolver::DIVERGED_LINEAR_SOLVER_FAILURE;
              // Print a message
              libMesh::out << "Linear solver failed during Newton step, dropping out."
                           << std::endl;
              break;
            }
        }

      // We may need to localize a parallel solution
      _system.update ();
      // The linear solver may not have fit our constraints exactly
#ifdef LIBMESH_ENABLE_CONSTRAINTS
      _system.get_dof_map().enforce_constraints_exactly
        (_system, &linear_solution, /* homogeneous = */ true);
#endif

      const unsigned int linear_steps = rval.first;
      libmesh_assert_less_equal (linear_steps, max_linear_iterations);
      _inner_iterations += linear_steps;

      const bool linear_solve_finished =
        !(linear_steps == max_linear_iterations);

      if (verbose)
        libMesh::out << "Linear solve finished, step " << linear_steps
                     << ", residual " << rval.second
                     << std::endl;

      // Compute the l2 norm of the nonlinear update
      Real norm_delta = linear_solution.l2_norm();

      if (verbose)
        libMesh::out << "Trying full Newton step" << std::endl;
      // Take a full Newton step
      newton_iterate.add (-1., linear_solution);
      newton_iterate.close();

      if (this->linear_solution_monitor.get())
        {
          // Compute the l2 norm of the whole solution
          norm_total = newton_iterate.l2_norm();
          rhs.close();
          (*this->linear_solution_monitor)(linear_solution, norm_delta,
                                           newton_iterate, norm_total,
                                           rhs, rhs.l2_norm(), _outer_iterations);
        }

      // Check residual with full Newton step, if that's useful for determining
      // whether to line search, whether to quit early, or whether to die after
      // hitting our max iteration count
      if (this->require_residual_reduction ||
          this->require_finite_residual ||
          _outer_iterations+1 < max_nonlinear_iterations ||
          !continue_after_max_iterations)
        {
          _system.assembly(true, false);

          rhs.close();
          current_residual = rhs.l2_norm();
          if (verbose)
            libMesh::out << "  Current Residual: "
                         << current_residual << std::endl;

          // don't fiddle around if we've already converged
          if (test_convergence(current_residual, norm_delta,
                               linear_solve_finished &&
                               current_residual <= last_residual))
            {
              if (!quiet)
                print_convergence(_outer_iterations, current_residual,
                                  norm_delta, linear_solve_finished &&
                                  current_residual <= last_residual);
              _outer_iterations++;
              break; // out of _outer_iterations for loop
            }
        }

      // since we're not converged, backtrack if necessary
      Real steplength =
        this->line_search(std::sqrt(TOLERANCE),
                          last_residual, current_residual,
                          newton_iterate, linear_solution);
      norm_delta *= steplength;

      // Check to see if backtracking failed,
      // and break out of the nonlinear loop if so...
      if (_solve_result == DiffSolver::DIVERGED_BACKTRACKING_FAILURE)
        {
          _outer_iterations++;
          break; // out of _outer_iterations for loop
        }

      if (_outer_iterations + 1 >= max_nonlinear_iterations)
        {
          libMesh::out << "  Nonlinear solver reached maximum step "
                       << max_nonlinear_iterations << ", latest evaluated residual "
                       << current_residual << std::endl;
          if (continue_after_max_iterations)
            {
              _solve_result = DiffSolver::DIVERGED_MAX_NONLINEAR_ITERATIONS;
              libMesh::out << "  Continuing..." << std::endl;
            }
          else
            {
              libmesh_convergence_failure();
            }
          continue;
        }

      // Compute the l2 norm of the whole solution
      norm_total = newton_iterate.l2_norm();

      max_solution_norm = std::max(max_solution_norm, norm_total);

      // Print out information for the
      // nonlinear iterations.
      if (verbose)
        libMesh::out << "  Nonlinear step: |du|/|u| = "
                     << norm_delta / norm_total
                     << ", |du| = " << norm_delta
                     << std::endl;

      // Terminate the solution iteration if the difference between
      // this iteration and the last is sufficiently small.
      if (test_convergence(current_residual, norm_delta / steplength,
                           linear_solve_finished))
        {
          if (!quiet)
            print_convergence(_outer_iterations, current_residual,
                              norm_delta / steplength,
                              linear_solve_finished);
          _outer_iterations++;
          break; // out of _outer_iterations for loop
        }
    } // end nonlinear loop

  // The linear solver may not have fit our constraints exactly
#ifdef LIBMESH_ENABLE_CONSTRAINTS
  _system.get_dof_map().enforce_constraints_exactly(_system);
#endif

  // We may need to localize a parallel solution
  _system.update ();

  // Make sure we are returning something sensible as the
  // _solve_result, except in the edge case where we weren't really asked to
  // solve.
  libmesh_assert (_solve_result != DiffSolver::INVALID_SOLVE_RESULT ||
                  !max_nonlinear_iterations);

  return _solve_result;
}



bool NewtonSolver::test_convergence(Real current_residual,
                                    Real step_norm,
                                    bool linear_solve_finished)
{
  // We haven't converged unless we pass a convergence test
  bool has_converged = false;

  // Is our absolute residual low enough?
  if (current_residual < absolute_residual_tolerance)
    {
      _solve_result |= CONVERGED_ABSOLUTE_RESIDUAL;
      has_converged = true;
    }

  // Is our relative residual low enough?
  if ((current_residual / max_residual_norm) <
      relative_residual_tolerance)
    {
      _solve_result |= CONVERGED_RELATIVE_RESIDUAL;
      has_converged = true;
    }

  // For incomplete linear solves, it's not safe to test step sizes
  if (!linear_solve_finished)
    {
      return has_converged;
    }

  // Is our absolute Newton step size small enough?
  if (step_norm < absolute_step_tolerance)
    {
      _solve_result |= CONVERGED_ABSOLUTE_STEP;
      has_converged = true;
    }

  // Is our relative Newton step size small enough?
  if (step_norm / max_solution_norm <
      relative_step_tolerance)
    {
      _solve_result |= CONVERGED_RELATIVE_STEP;
      has_converged = true;
    }

  return has_converged;
}


void NewtonSolver::print_convergence(unsigned int step_num,
                                     Real current_residual,
                                     Real step_norm,
                                     bool linear_solve_finished)
{
  // Is our absolute residual low enough?
  if (current_residual < absolute_residual_tolerance)
    {
      libMesh::out << "  Nonlinear solver converged, step " << step_num
                   << ", residual " << current_residual
                   << std::endl;
    }
  else if (absolute_residual_tolerance)
    {
      if (verbose)
        libMesh::out << "  Nonlinear solver current_residual "
                     << current_residual << " > "
                     << (absolute_residual_tolerance) << std::endl;
    }

  // Is our relative residual low enough?
  if ((current_residual / max_residual_norm) <
      relative_residual_tolerance)
    {
      libMesh::out << "  Nonlinear solver converged, step " << step_num
                   << ", residual reduction "
                   << current_residual / max_residual_norm
                   << " < " << relative_residual_tolerance
                   << std::endl;
    }
  else if (relative_residual_tolerance)
    {
      if (verbose)
        libMesh::out << "  Nonlinear solver relative residual "
                     << (current_residual / max_residual_norm)
                     << " > " << relative_residual_tolerance
                     << std::endl;
    }

  // For incomplete linear solves, it's not safe to test step sizes
  if (!linear_solve_finished)
    return;

  // Is our absolute Newton step size small enough?
  if (step_norm < absolute_step_tolerance)
    {
      libMesh::out << "  Nonlinear solver converged, step " << step_num
                   << ", absolute step size "
                   << step_norm
                   << " < " << absolute_step_tolerance
                   << std::endl;
    }
  else if (absolute_step_tolerance)
    {
      if (verbose)
        libMesh::out << "  Nonlinear solver absolute step size "
                     << step_norm
                     << " > " << absolute_step_tolerance
                     << std::endl;
    }

  // Is our relative Newton step size small enough?
  if (step_norm / max_solution_norm <
      relative_step_tolerance)
    {
      libMesh::out << "  Nonlinear solver converged, step " << step_num
                   << ", relative step size "
                   << (step_norm / max_solution_norm)
                   << " < " << relative_step_tolerance
                   << std::endl;
    }
  else if (relative_step_tolerance)
    {
      if (verbose)
        libMesh::out << "  Nonlinear solver relative step size "
                     << (step_norm / max_solution_norm)
                     << " > " << relative_step_tolerance
                     << std::endl;
    }
}

} // namespace libMesh
