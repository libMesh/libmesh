// $Id: time_solver.h,v 1.10 2007-02-21 21:13:37 roystgnr Exp $

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



#ifndef __time_solver_h__
#define __time_solver_h__

// C++ includes

// Local includes
#include "auto_ptr.h"
#include "libmesh_common.h"
#include "numeric_vector.h"
#include "reference_counted_object.h"

// Forward Declarations
class DiffSolver;
class TimeSolver;
class DifferentiableSystem;

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
   * This method solves for the solution at the next timestep (or solves for a
   * steady-state solution).  Usually we will only need to solve one
   * (non)linear system per timestep, but more complex subclasses
   * may override this.
   */
  virtual void solve ();

  /**
   * This method advances the solution to the next timestep, after a
   * solve() has been performed.  Often this will be
   * done after every TimeSolver::solve(), but adaptive mesh refinement
   * and/or adaptive timestep selection may require some solve() steps
   * to be repeated.
   */
  virtual void advance_timestep ();

  /**
   * This method uses the DifferentiableSystem's
   * element_time_derivative() and element_constraint()
   * to build a full residual on an element.  What combination
   * it uses will depend on the type of solver.  See
   * the subclasses for more details.
   */
  virtual bool element_residual (bool get_jacobian) = 0;

  /**
   * This method uses the DifferentiableSystem's
   * side_time_derivative() and side_constraint()
   * to build a full residual on an element's side.
   * What combination it uses will depend on the type
   * of solver.  See the subclasses for more details.
   */
  virtual bool side_residual (bool get_jacobian) = 0;

  /**
   * This method is for subclasses or users to override
   * to do arbitrary processing between timesteps
   */
  virtual void before_timestep () {}

  /**
   * @returns the old nonlinear solution for the specified global
   * DOF.
   */
  Number old_nonlinear_solution (const unsigned int global_dof_number) const;

  /**
   * @returns a constant reference to the system we are solving.
   */
  const sys_type & system () const { return _system; }

  /**
   * An implicit linear or nonlinear solver to use at each timestep.
   */
  virtual AutoPtr<DiffSolver> &diff_solver();

protected:

  /**
   * An implicit linear or nonlinear solver to use at each timestep.
   */
  AutoPtr<DiffSolver> _diff_solver;

  /**
   * @returns a writeable reference to the system we are solving.
   */
  sys_type & system () { return _system; }
  
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



#endif // #define __time_solver_h__
