// $Id: time_solver.h 2501 2007-11-20 02:33:29Z benkirk $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2007  Benjamin S. Kirk, John W. Peterson
  
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



#ifndef __unsteady_solver_h__
#define __unsteady_solver_h__

// C++ includes

// Local includes
#include "auto_ptr.h"
#include "libmesh_common.h"
#include "numeric_vector.h"
#include "time_solver.h"

// Forward Declarations
class UnsteadySolver;

/**
 * This is a generic class that defines a solver to handle
 * time integration of DifferentiableSystems.
 *
 * A user can define a solver for unsteady problems by deriving
 * from this class and implementing certain functions.
 *
 * This class is part of the new DifferentiableSystem framework,
 * which is still experimental.  Users of this framework should
 * beware of bugs and future API changes.
 *
 * @author Roy H. Stogner 2008
 */

// ------------------------------------------------------------
// UnsteadySolver class definition
class UnsteadySolver : public TimeSolver
{
public:
  /**
   * Constructor. Requires a reference to the system
   * to be solved.
   */
  UnsteadySolver (sys_type& s);
  
  /**
   * Destructor.
   */
  virtual ~UnsteadySolver ();

  /**
   * The initialization function.  This method is used to
   * initialize internal data structures before a simulation begins.
   */
  virtual void init ();

  /**
   * This method solves for the solution at the next timestep.
   * Usually we will only need to solve one (non)linear system per timestep,
   * but more complex subclasses may override this.
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
   * This method should return the expected convergence order of the
   * (non-local) error of the time discretization scheme - e.g. 2 for the
   * O(deltat^2) Crank-Nicholson, or 1 for the O(deltat) Backward Euler.
   *
   * Useful for adaptive timestepping schemes.
   */
  virtual Real error_order () const = 0;

  /**
   * @returns the old nonlinear solution for the specified global
   * DOF.
   */
  Number old_nonlinear_solution (const unsigned int global_dof_number) const;

  /**
   * Computes the size of ||u^{n+1} - u^{n}|| in some norm.
   * Supported norms are
   * l2: norm_type=0
   * l1: norm_type=1
   * 
   * Note that, while you can always call this function, its
   * result may or may not be very meaningful.  For example, if
   * you call this function right after calling advance_timestep()
   * then you'll get a result of zero since old_nonlinear_solution
   * is set equal to nonlinear_solution in this function.
   */
  Real du(unsigned char norm_type=0);

  /**
   * This value (which defaults to zero) is the number of times the
   * UnsteadySolver is allowed to halve deltat and let the DiffSolver
   * repeat the latest failed solve with a reduced timestep.  Note
   * that this has no effect for SteadySolvers.  Note that you must
   * set at least one of the DiffSolver flags
   * "continue_after_max_iterations" or
   * "continue_after_backtrack_failure" to allow the UnsteadySolver to
   * retry the solve.
   */
  unsigned int reduce_deltat_on_diffsolver_failure;
  
protected:

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



#endif // #define __time_solver_h__
