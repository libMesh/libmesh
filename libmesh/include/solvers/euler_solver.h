// $Id: euler_solver.h,v 1.1 2006-06-05 21:51:15 roystgnr Exp $

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



#ifndef __euler_solver_h__
#define __euler_solver_h__

// C++ includes

// Local includes
#include "time_solver.h"

/**
 * This class defines a theta-method Euler (defaulting to Backward
 * Euler with theta = 1.0) solver to handle
 * time integration of DifferentiableSystems.
 *
 * @author Roy H. Stogner 2006
 */

// ------------------------------------------------------------
// Solver class definition
class EulerSolver : public TimeSolver
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
  EulerSolver (sys_type& s);
  
  /**
   * Destructor.
   */
  virtual ~EulerSolver ();

  /**
   * Initialization.
   */
  virtual void init ();

  /**
   * This method uses the DifferentiableSystem's
   * element_time_derivative() and element_constraint()
   * to build a full residual on an element.  What combination
   * it uses will depend on theta.
   */
  virtual bool element_residual (bool get_jacobian);

  /**
   * This method uses the DifferentiableSystem's
   * side_time_derivative() and side_constraint()
   * to build a full residual on an element's side.
   * What combination it uses will depend on theta.
   */
  virtual bool side_residual (bool get_jacobian);

  /**
   * The value for the theta method to employ: 1.0 corresponds
   * to backwards Euler, 0.0 corresponds to forwards Euler,
   * 0.5 corresponds to Crank-Nicholson.
   */
  Real theta;
};



#endif // #define __euler_solver_h__
