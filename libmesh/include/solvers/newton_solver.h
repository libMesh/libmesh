// $Id: newton_solver.h,v 1.7 2006-10-23 23:26:39 roystgnr Exp $

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



#ifndef __newton_solver_h__
#define __newton_solver_h__

// C++ includes

// Local includes
#include "libmesh_common.h"
#include "linear_solver.h"
#include "reference_counted_object.h"
#include "diff_solver.h"

/**
 * This class defines a solver which uses the default
 * libMesh linear solver in a quasiNewton method to handle a 
 * DifferentiableSystem
 *
 * This class is part of the new DifferentiableSystem framework,
 * which is still experimental.  Users of this framework should
 * beware of bugs and future API changes.
 *
 * @author Roy H. Stogner 2006
 */

// ------------------------------------------------------------
// Solver class definition
class NewtonSolver : public DiffSolver
{
public:
  /**
   * Constructor. Requires a reference to the system
   * to be solved.
   */
  NewtonSolver (sys_type& system);
  
  /**
   * Destructor.
   */
  virtual ~NewtonSolver ();

  typedef DiffSolver Parent;

  /**
   * The reinitialization function.  This method is used after
   * changes in the mesh.
   */
  virtual void reinit ();

  /**
   * This method performs a solve.  What occurs in
   * this method will depend on the type of solver.  See
   * the subclasses for more details.
   */
  virtual void solve ();

  /**
   * If this is set to true, the solver is forced to test the residual
   * after each Newton step, and to reduce the length of its steps
   * whenever necessary to avoid a residual increase.
   * It is currently set to true by default; set it to false to
   * avoid unnecessary residual assembly on well-behaved systems.
   */
  bool require_residual_reduction;

  /**
   * If the quasi-Newton step length must be reduced to below this
   * factor to give a residual reduction, then the Newton solver
   * dies with an error()
   * It is currently set to 1e-5 by default.
   */
  Real minsteplength;

protected:

  /**
   * The \p LinearSolver defines the interface used to
   * solve the linear_implicit system.  This class handles all the
   * details of interfacing with various linear algebra packages
   * like PETSc or LASPACK.
   */
  AutoPtr<LinearSolver<Number> > linear_solver;

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



#endif // #define __newton_solver_h__
