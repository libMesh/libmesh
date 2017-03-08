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



#ifndef LIBMESH_NEWTON_SOLVER_H
#define LIBMESH_NEWTON_SOLVER_H

// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/linear_solver.h"
#include "libmesh/reference_counted_object.h"
#include "libmesh/diff_solver.h"

// C++ includes

namespace libMesh
{

/**
 * This class defines a solver which uses the default
 * libMesh linear solver in a quasiNewton method to handle a
 * DifferentiableSystem
 *
 * This class is part of the new DifferentiableSystem framework,
 * which is still experimental.  Users of this framework should
 * beware of bugs and future API changes.
 *
 * \author Roy H. Stogner
 * \date 2006
 */
class NewtonSolver : public DiffSolver
{
public:
  /**
   * Constructor. Requires a reference to the system
   * to be solved.
   */
  explicit
  NewtonSolver (sys_type & system);

  /**
   * Destructor.
   */
  virtual ~NewtonSolver ();

  typedef DiffSolver Parent;

  /**
   * The initialization function.  This method is used to
   * initialize internal data structures before a simulation begins.
   */
  virtual void init () libmesh_override;

  /**
   * The reinitialization function.  This method is used after
   * changes in the mesh.
   */
  virtual void reinit () libmesh_override;

  /**
   * This method performs a solve, using an inexact Newton-Krylov
   * method with line search.
   */
  virtual unsigned int solve () libmesh_override;

  LinearSolver<Number> & get_linear_solver()
  { libmesh_assert(_linear_solver);
    return *_linear_solver;
  }

  const LinearSolver<Number> & get_linear_solver() const
  { libmesh_assert(_linear_solver);
    return *_linear_solver;
  }

  /**
   * If this is set to true, the solver is forced to test the residual
   * after each Newton step, and to reduce the length of its steps
   * whenever necessary to avoid a residual increase.
   * It is currently set to true by default; set it to false to
   * avoid unnecessary residual assembly on well-behaved systems.
   */
  bool require_residual_reduction;

  /**
   * If this is set to true, the solver is forced to test the residual
   * after each Newton step, and to reduce the length of its steps
   * whenever necessary to avoid an infinite or NaN residual.
   * It is currently set to true by default; set it to false to
   * avoid unnecessary residual assembly on well-behaved systems.
   */
  bool require_finite_residual;

  /**
   * If require_residual_reduction is true, the solver may reduce step
   * lengths when required.  If so, brent_line_search is an option.
   * If brent_line_search is set to false, the solver reduces the
   * length of its steps by 1/2 iteratively until it finds residual
   * reduction.  If true, step lengths are first reduced by 1/2 or
   * more to find some residual reduction, then Brent's method is used
   * to find as much residual reduction as possible.
   *
   * brent_line_search is currently set to true by default.
   */
  bool brent_line_search;

  /**
   * If set to true, check for convergence of the linear solve. If no
   * convergence is acquired during the linear solve, the nonlinear solve
   * fails with DiffSolver::DIVERGED_LINEAR_SOLVER_FAILURE.
   * Enabled by default as nonlinear convergence is very difficult, if the
   * linear solver is not converged.
   */
  bool track_linear_convergence;

  /**
   * If the quasi-Newton step length must be reduced to below this
   * factor to give a residual reduction, then the Newton solver
   * dies with an error message.
   * It is currently set to 1e-5 by default.
   */
  Real minsteplength;

  /**
   * The tolerance for linear solves is kept below this multiplier (which
   * defaults to 1e-3) times the norm of the current nonlinear residual
   */
  Real linear_tolerance_multiplier;

protected:

  /**
   * The \p LinearSolver defines the interface used to
   * solve the linear_implicit system.  This class handles all the
   * details of interfacing with various linear algebra packages
   * like PETSc or LASPACK.
   */
  UniquePtr<LinearSolver<Number> > _linear_solver;

  /**
   * This does a line search in the direction opposite linear_solution
   * to try and minimize the residual of newton_iterate.
   * newton_iterate is moved to the end of the quasiNewton step, and
   * the return value is the substep size.
   */
  Real line_search(Real tol,
                   Real last_residual,
                   Real & current_residual,
                   NumericVector<Number> & newton_iterate,
                   const NumericVector<Number> & linear_solution);

  /**
   * This prints output for the convergence criteria based on
   * by the given residual and step size.
   */
  void print_convergence(unsigned int step_num,
                         Real current_residual,
                         Real step_norm,
                         bool linear_solve_finished);

  /**
   * This returns true if a convergence criterion has been passed
   * by the given residual and step size; false otherwise.
   */
  bool test_convergence(Real current_residual,
                        Real step_norm,
                        bool linear_solve_finished);
};


} // namespace libMesh


#endif // LIBMESH_NEWTON_SOLVER_H
