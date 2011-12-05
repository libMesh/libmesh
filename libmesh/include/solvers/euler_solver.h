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



#ifndef __euler_solver_h__
#define __euler_solver_h__

// C++ includes

// Local includes
#include "unsteady_solver.h"

namespace libMesh
{

/**
 * This class defines a theta-method Euler (defaulting to Backward
 * Euler with theta = 1.0) solver to handle
 * time integration of DifferentiableSystems.
 *
 * This class is part of the new DifferentiableSystem framework,
 * which is still experimental.  Users of this framework should
 * beware of bugs and future API changes.
 *
 * @author Roy H. Stogner 2006
 */

// ------------------------------------------------------------
// Solver class definition
class EulerSolver : public UnsteadySolver
{
public:
  /**
   * The parent class
   */
  typedef UnsteadySolver Parent;

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
   * Error convergence order: 2 for Crank-Nicolson, 1 otherwise
   */
  virtual Real error_order() const;

  /**
   * This method uses the DifferentiableSystem's
   * element_time_derivative() and element_constraint()
   * to build a full residual on an element.  What combination
   * it uses will depend on theta.
   */
  virtual bool element_residual (bool request_jacobian,
                                 DiffContext&);

  /**
   * This method uses the DifferentiableSystem's
   * side_time_derivative() and side_constraint()
   * to build a full residual on an element's side.
   * What combination it uses will depend on theta.
   */
  virtual bool side_residual (bool request_jacobian,
                              DiffContext&);

  /**
   * The value for the theta method to employ: 1.0 corresponds
   * to backwards Euler, 0.0 corresponds to forwards Euler,
   * 0.5 corresponds to Crank-Nicolson.
   */
  Real theta;
};


} // namespace libMesh


#endif // #define __euler_solver_h__
