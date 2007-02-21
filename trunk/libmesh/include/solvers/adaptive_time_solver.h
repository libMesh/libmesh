// $Id: adaptive_time_solver.h,v 1.1 2007-02-21 21:13:37 roystgnr Exp $

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



#ifndef __adaptive_time_solver_h__
#define __adaptive_time_solver_h__

// C++ includes

// Local includes
#include "time_solver.h"

/**
 * This class wraps another TimeSolver derived class, and compares the results
 * of timestepping with deltat and timestepping with 2*deltat to adjust
 * future timestep lengths.
 *
 * Currently this class only works on fully coupled Systems
 *
 * This class is part of the new DifferentiableSystem framework,
 * which is still experimental.  Users of this framework should
 * beware of bugs and future API changes.
 *
 * @author Roy H. Stogner 2007
 */

// ------------------------------------------------------------
// Solver class definition
class AdaptiveTimeSolver : public TimeSolver
{
public:
  /**
   * The type of system
   */
  typedef DifferentiableSystem sys_type;
  
  /**
   * The parent class
   */
  typedef TimeSolver Parent;
  
  /**
   * Constructor. Requires a reference to the system
   * to be solved.
   */
  AdaptiveTimeSolver (sys_type& s);
  
  /**
   * Destructor.
   */
  virtual ~AdaptiveTimeSolver ();

  virtual void init();

  virtual void reinit();

  virtual void solve();

  /**
   * This method is passed on to the core_time_solver
   */
  virtual bool element_residual (bool get_jacobian);

  /**
   * This method is passed on to the core_time_solver
   */
  virtual bool side_residual (bool get_jacobian);

  /**
   * An implicit linear or nonlinear solver to use at each timestep.
   */
  virtual AutoPtr<DiffSolver> &diff_solver();

  /**
   * This method is used to take timesteps
   */
  AutoPtr<TimeSolver> core_time_solver;

  /**
   * This tolerance is the target error between double-deltat
   * and single-deltat timesteps
   */
  Real tolerance;
};



#endif // #define __adaptive_time_solver_h__
