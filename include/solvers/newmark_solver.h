// The libMesh Finite Element Library.
// Copyright (C) 2002-2014 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_NEWMARK_SOLVER_H
#define LIBMESH_NEWMARK_SOLVER_H

// Local includes
#include "libmesh/unsteady_solver.h"

namespace libMesh
{
/**
 * This class defines a Newmark time integrator for
 * second order (in time) DifferentiableSystems.
 * There are two parameters
 * \f$\gamma\f$ and \f$\beta\f$ (defaulting to 0.5 and 0.25,
 * respectively.
 *
 * This class is part of the new DifferentiableSystem framework,
 * which is still experimental.  Users of this framework should
 * beware of bugs and future API changes.
 *
 * @author Paul T. Bauman 2015
 */

// ------------------------------------------------------------
// Solver class definition
class NewmarkSolver : public UnsteadySolver
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
  explicit
  NewmarkSolver (sys_type& s);

  /**
   * Destructor.
   */
  virtual ~NewmarkSolver ();

  /**
   * The initialization function.  This method is used to
   * initialize internal data structures before a simulation begins.
   */
  virtual void init ();

  /**
   * The data initialization function.  This method is used to
   * initialize internal data structures after the underlying System
   * has been initialized
   */
  virtual void init_data ();

  /**
   * The reinitialization function.  This method is used to
   * resize internal data vectors after a mesh change.
   */
  virtual void reinit ();

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
   * be done before every UnsteadySolver::adjoint_solve().
   */
  virtual void adjoint_advance_timestep ();

  /**
   * This method retrieves all the stored solutions at the current
   * system.time
   */
  virtual void retrieve_timestep ();

  /**
   * Error convergence order: 2 for \f$\gamma=0.5\f$, 1 otherwise
   */
  virtual Real error_order() const;

  /**
   * This method uses the DifferentiablePhysics'
   * element_time_derivative() and element_constraint()
   * to build a full residual on an element.  What combination
   * it uses will depend on theta.
   */
  virtual bool element_residual (bool request_jacobian,
                                 DiffContext&);

  /**
   * This method uses the DifferentiablePhysics'
   * side_time_derivative() and side_constraint()
   * to build a full residual on an element's side.
   * What combination it uses will depend on theta.
   */
  virtual bool side_residual (bool request_jacobian,
                              DiffContext&);

  /**
   * This method uses the DifferentiablePhysics'
   * nonlocal_time_derivative() and nonlocal_constraint()
   * to build a full residual for non-local terms.
   * What combination it uses will depend on theta.
   */
  virtual bool nonlocal_residual (bool request_jacobian,
                                  DiffContext&);

  /**
   * @returns the solution rate at the previous time step, \f$\dot{u}_n\f$,
   * for the specified global DOF.
   */
  Number old_solution_rate (const dof_id_type global_dof_number) const;

  /**
   * @returns the solution acceleration at the previous time step, \f$\ddot{u}_n\f$,
   * for the specified global DOF.
   */
  Number old_solution_accel (const dof_id_type global_dof_number) const;


protected:

  /**
   * The value for the \f$\beta\f$ to employ. Method is
   * unconditionally stable for
   * \f$ \beta \ge \frac{1}{4} \left( \gamma + \frac{1}{2}\right)^2 \f$
   */
  Real _beta;

  /**
   * The value for \f$\gamma\f$ to employ. Newmark
   * is 2nd order iff \f$\gamma=0.5\f$.
   */
  Real _gamma;

  /**
   * Serial vector of previous time step velocity \f$ \dot{u}_n \f$
   */
  UniquePtr<NumericVector<Number> > _old_local_solution_rate;

  /**
   * Serial vector of previous time step accleration \f$ \ddot{u}_n \f$
   */
  UniquePtr<NumericVector<Number> > _old_local_solution_accel;

  /**
   * This method is the underlying implementation of the public
   * residual methods.
   */
  virtual bool _general_residual (bool request_jacobian,
                                  DiffContext&,
                                  ResFuncType mass,
                                  ResFuncType damping,
                                  ResFuncType time_deriv,
                                  ResFuncType constraint);

};

} // namespace libMesh

#endif // LIBMESH_NEWMARK_SOLVER_H
