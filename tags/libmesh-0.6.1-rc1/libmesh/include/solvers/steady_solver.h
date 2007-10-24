
// $Id: steady_solver.h,v 1.6 2007-05-14 20:33:50 jwpeterson Exp $

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



#ifndef __steady_solver_h__
#define __steady_solver_h__

// C++ includes

// Local includes
#include "time_solver.h"

// Forward Declarations
class TimeSolver;
class DifferentiableSystem;

/**
 * This class implements a TimeSolver which does a single
 * solve of the steady state problem.
 *
 * This class is part of the new DifferentiableSystem framework,
 * which is still experimental.  Users of this framework should
 * beware of bugs and future API changes.
 *
 * @author Roy H. Stogner 2006
 */

// ------------------------------------------------------------
// Solver class definition
class SteadySolver : public TimeSolver
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
  SteadySolver (sys_type& s) : Parent(s) {}
  
  /**
   * Destructor.
   */
  virtual ~SteadySolver ();

  /**
   * error convergence order against deltat is
   * not applicable to a steady problem.
   */
  virtual Real error_order() const { return 0.; }

  /**
   * This method uses the DifferentiableSystem's
   * element_time_derivative() and element_constraint()
   * to build a full residual/jacobian on an element.
   */
  virtual bool element_residual (bool get_jacobian);

  /**
   * This method uses the DifferentiableSystem's
   * side_time_derivative() and side_constraint()
   * to build a full residual/jacobian on an element's side.
   */
  virtual bool side_residual (bool get_jacobian);

  /**
   * Nominally computes the size of the difference between
   * successive solution iterates ||u^{n+1} - u^{n}|| in some norm,
   * but for this class just returns 0.
   */
  Real du(unsigned char) { return 0.; }
};



#endif // #define __steady_solver_h__
