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



#ifndef LIBMESH_EULER2_SOLVER_H
#define LIBMESH_EULER2_SOLVER_H

// Local includes
#include "libmesh/first_order_unsteady_solver.h"

// C++ includes

namespace libMesh
{

/**
 * This class defines a theta-method (defaulting to Backward
 * Euler with theta = 1.0) solver to handle
 * time integration of DifferentiableSystems.
 * The "Euler2" solver differs from Euler in how it evaluates
 * residuals at intermediate theta values:
 * Euler solves m(u,u') = f(theta*u_new + (1-theta)*u_old),
 * Euler2 solves m(u') = theta*f(u_new) + (1-theta)*f(u_old)
 * i.e. the trapezoidal rule for theta = 0.5
 *
 * This class is part of the new DifferentiableSystem framework,
 * which is still experimental.  Users of this framework should
 * beware of bugs and future API changes.
 *
 * \author Roy H. Stogner
 * \date 2006
 */
class Euler2Solver : public FirstOrderUnsteadySolver
{
public:
  /**
   * The parent class
   */
  typedef FirstOrderUnsteadySolver Parent;

  /**
   * Constructor. Requires a reference to the system
   * to be solved.
   */
  explicit
  Euler2Solver (sys_type & s);

  /**
   * Destructor.
   */
  virtual ~Euler2Solver ();

  /**
   * Error convergence order: 2 for Crank-Nicolson, 1 otherwise
   */
  virtual Real error_order() const libmesh_override;

  /**
   * This method uses the DifferentiablePhysics'
   * element_time_derivative() and element_constraint()
   * to build a full residual on an element.  What combination
   * it uses will depend on theta.
   */
  virtual bool element_residual (bool request_jacobian,
                                 DiffContext &) libmesh_override;

  /**
   * This method uses the DifferentiablePhysics'
   * side_time_derivative() and side_constraint()
   * to build a full residual on an element's side.
   * What combination it uses will depend on theta.
   */
  virtual bool side_residual (bool request_jacobian,
                              DiffContext &) libmesh_override;

  /**
   * This method uses the DifferentiablePhysics'
   * nonlocal_time_derivative() and nonlocal_constraint()
   * to build a full residual for non-local terms.
   * What combination it uses will depend on theta.
   */
  virtual bool nonlocal_residual (bool request_jacobian,
                                  DiffContext &) libmesh_override;

  /**
   * The value for the theta method to employ: 1.0 corresponds
   * to backwards Euler, 0.0 corresponds to forwards Euler,
   * 0.5 corresponds to a Crank-Nicolson-like scheme.
   */
  Real theta;

protected:

  /**
   * This method is the underlying implementation of the public
   * residual methods.
   */
  virtual bool _general_residual (bool request_jacobian,
                                  DiffContext &,
                                  ResFuncType mass,
                                  ResFuncType damping,
                                  ResFuncType time_deriv,
                                  ResFuncType constraint,
                                  ReinitFuncType reinit,
                                  bool compute_second_order_eqns);

};


} // namespace libMesh


#endif // LIBMESH_EULER2_SOLVER_H
