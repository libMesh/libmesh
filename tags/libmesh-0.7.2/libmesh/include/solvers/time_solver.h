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



#ifndef __time_solver_h__
#define __time_solver_h__

// C++ includes

// Local includes
#include "auto_ptr.h"
#include "libmesh_common.h"
#include "linear_solver.h"
#include "numeric_vector.h"
#include "reference_counted_object.h"

namespace libMesh
{

// Forward Declarations
class DiffContext;
class DiffSolver;
class DifferentiableSystem;
class ParameterVector;
class SystemNorm;
class TimeSolver;

/**
 * This is a generic class that defines a solver to handle
 * time integration of DifferentiableSystems.
 *
 * A user can define a solver by deriving from this class and 
 * implementing certain functions.
 *
 * This class is part of the new DifferentiableSystem framework,
 * which is still experimental.  Users of this framework should
 * beware of bugs and future API changes.
 *
 * @author Roy H. Stogner 2006
 */

// ------------------------------------------------------------
// Solver class definition
class TimeSolver : public ReferenceCountedObject<TimeSolver>
{
public:
  /**
   * The type of system
   */
  typedef DifferentiableSystem sys_type;
  
  /**
   * Constructor. Requires a reference to the system
   * to be solved.
   */
  TimeSolver (sys_type& s);
  
  /**
   * Destructor.
   */
  virtual ~TimeSolver ();

  /**
   * The initialization function.  This method is used to
   * initialize internal data structures before a simulation begins.
   */
  virtual void init ();

  /**
   * The reinitialization function.  This method is used after
   * changes in the mesh
   */
  virtual void reinit ();

  /**
   * This method solves for the solution at the next timestep (or
   * solves for a steady-state solution).  Usually we will only need
   * to solve one (non)linear system per timestep, but more complex
   * subclasses may override this.
   */
  virtual void solve ();

  /**
   * This method advances the solution to the next timestep, after a
   * solve() has been performed.  Often this will be done after every
   * UnsteadySolver::solve(), but adaptive mesh refinement and/or adaptive
   * time step selection may require some solve() steps to be repeated.
   */
  virtual void advance_timestep ();

  /**
   * This method advances the adjoint solution to the previous
   * timestep, after an adjoint_solve() has been performed.  This will
   * probably be done after every UnsteadySolver::adjoint_solve().
   */
  virtual void adjoint_recede_timestep ();

  /**
   * This method uses the DifferentiableSystem's
   * element_time_derivative() and element_constraint()
   * to build a full residual on an element.  What combination
   * it uses will depend on the type of solver.  See
   * the subclasses for more details.
   */
  virtual bool element_residual (bool request_jacobian,
                                 DiffContext &) = 0;

  /**
   * This method uses the DifferentiableSystem's
   * side_time_derivative() and side_constraint()
   * to build a full residual on an element's side.
   * What combination it uses will depend on the type
   * of solver.  See the subclasses for more details.
   */
  virtual bool side_residual (bool request_jacobian,
                              DiffContext &) = 0;

  /**
   * This method is for subclasses or users to override
   * to do arbitrary processing between timesteps
   */
  virtual void before_timestep () {}

  /**
   * @returns a constant reference to the system we are solving.
   */
  const sys_type & system () const { return _system; }

  /**
   * @returns a writeable reference to the system we are solving.
   */
  sys_type & system () { return _system; }
  
  /**
   * An implicit linear or nonlinear solver to use at each timestep.
   */
  virtual AutoPtr<DiffSolver> &diff_solver() { return _diff_solver; }

  /**
   * An implicit linear solver to use for adjoint and sensitivity problems.
   */
  virtual AutoPtr<LinearSolver<Number> > &linear_solver() { return _linear_solver; }

  /**
   * Print extra debugging information if quiet ==  false.
   */
  bool quiet;

  /**
   * Computes the size of ||u^{n+1} - u^{n}|| in some norm.
   * 
   * Note that, while you can always call this function, its
   * result may or may not be very meaningful.  For example, if
   * you call this function right after calling advance_timestep()
   * then you'll get a result of zero since old_nonlinear_solution
   * is set equal to nonlinear_solution in this function.
   */
  virtual Real du(const SystemNorm& norm) const = 0;

  /**
   * Is this effectively a steady-state solver?
   */
  virtual bool is_steady() const = 0;

  /**
   * This value (which defaults to zero) is the number of times the
   * TimeSolver is allowed to halve deltat and let the DiffSolver
   * repeat the latest failed solve with a reduced timestep.  Note
   * that this has no effect for SteadySolvers.  Note that you must
   * set at least one of the DiffSolver flags
   * "continue_after_max_iterations" or
   * "continue_after_backtrack_failure" to allow the TimeSolver to
   * retry the solve.
   */
  unsigned int reduce_deltat_on_diffsolver_failure;
  
protected:

  /**
   * An implicit linear or nonlinear solver to use at each timestep.
   */
  AutoPtr<DiffSolver> _diff_solver;

  /**
   * An implicit linear solver to use for adjoint problems.
   */
  AutoPtr<LinearSolver<Number> > _linear_solver;

  /**
   * A reference to the system we are solving.
   */
  sys_type& _system;

  /**
   * A bool that will be true the first time solve() is called,
   * and false thereafter
   */
  bool first_solve;

  /**
   * Serial vector of _system.get_vector("_old_nonlinear_solution")
   */
  AutoPtr<NumericVector<Number> > old_local_nonlinear_solution;
};


} // namespace libMesh


#endif // #define __time_solver_h__
