// The libMesh Finite Element Library.
// Copyright (C) 2002-2016 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_STEADY_SOLVER_H
#define LIBMESH_STEADY_SOLVER_H

// Local includes
#include "libmesh/time_solver.h"

// C++ includes

namespace libMesh
{

// Forward Declarations
class DiffContext;
class DifferentiableSystem;
class TimeSolver;

/**
 * This class implements a TimeSolver which does a single
 * solve of the steady state problem.
 *
 * This class is part of the new DifferentiableSystem framework,
 * which is still experimental.  Users of this framework should
 * beware of bugs and future API changes.
 *
 * \author Roy H. Stogner
 * \date 2006
 */
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
  explicit
  SteadySolver (sys_type & s) : Parent(s) {}

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
   * This method uses the DifferentiablePhysics'
   * element_time_derivative() and element_constraint()
   * to build a full residual/jacobian on an element.
   */
  virtual bool element_residual (bool request_jacobian,
                                 DiffContext &) libmesh_override;

  /**
   * This method uses the DifferentiablePhysics'
   * side_time_derivative() and side_constraint()
   * to build a full residual/jacobian on an element's side.
   */
  virtual bool side_residual (bool request_jacobian,
                              DiffContext &) libmesh_override;

  /**
   * This method uses the DifferentiablePhysics'
   * nonlocal_time_derivative() and nonlocal_constraint()
   * to build a full residual/jacobian for non-local terms.
   */
  virtual bool nonlocal_residual (bool request_jacobian,
                                  DiffContext &) libmesh_override;

  /**
   * Nominally computes the size of the difference between
   * successive solution iterates ||u^{n+1} - u^{n}|| in some norm,
   * but for this class just returns 0.
   */
  virtual Real du(const SystemNorm &) const libmesh_override { return 0; }

  /**
   * This is a steady-state solver.
   */
  virtual bool is_steady() const libmesh_override { return true; }

protected:

  /**
   * This method is the underlying implementation of the public
   * residual methods.
   */
  virtual bool _general_residual (bool request_jacobian,
                                  DiffContext &,
                                  ResFuncType time_deriv,
                                  ResFuncType constraint);
};


} // namespace libMesh


#endif // LIBMESH_STEADY_SOLVER_H
