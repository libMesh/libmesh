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



#ifndef LIBMESH_NEWMARK_SOLVER_H
#define LIBMESH_NEWMARK_SOLVER_H

// Local includes
#include "libmesh/second_order_unsteady_solver.h"

namespace libMesh
{
/**
 * This class defines a Newmark time integrator for
 * second order (in time) DifferentiableSystems.
 * There are two parameters
 * \f$\gamma\f$ and \f$\beta\f$ (defaulting to 0.5 and 0.25,
 * respectively). Note that projectors are included for the
 * initial velocity and acceleration; the initial displacement
 * can be set by calling System::project_solution().
 *
 * This class is part of the new DifferentiableSystem framework,
 * which is still experimental.  Users of this framework should
 * beware of bugs and future API changes.
 *
 * \author Paul T. Bauman
 * \date 2015
 */
class NewmarkSolver : public SecondOrderUnsteadySolver
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
  NewmarkSolver (sys_type & s);

  /**
   * Destructor.
   */
  virtual ~NewmarkSolver ();

  /**
   * This method advances the solution to the next timestep, after a
   * solve() has been performed.  Often this will be done after every
   * UnsteadySolver::solve(), but adaptive mesh refinement and/or adaptive
   * time step selection may require some solve() steps to be repeated.
   */
  virtual void advance_timestep () libmesh_override;

  /**
   * This method advances the adjoint solution to the previous
   * timestep, after an adjoint_solve() has been performed.  This will
   * be done before every UnsteadySolver::adjoint_solve().
   */
  virtual void adjoint_advance_timestep () libmesh_override;

  /**
   * This method uses the specified initial displacement and velocity
   * to compute the initial acceleration \f$a_0\f$.
   */
  virtual void compute_initial_accel();

  /**
   * Specify non-zero initial acceleration. Should be called before solve().
   * This is an alternative to compute_initial_acceleration() if the
   * initial acceleration is actually known.
   * The function value f and its gradient g are user-provided cloneable functors.
   * A gradient g is only required/used for projecting onto finite element spaces
   * with continuous derivatives.
   */
  void project_initial_accel (FunctionBase<Number> * f,
                              FunctionBase<Gradient> * g = libmesh_nullptr);

  /**
   * Allow the user to (re)set whether the initial acceleration is available.
   * This is not needed if either compute_initial_accel() or project_initial_accel()
   * are called. This is useful is the user is restarting a calculation and the acceleration
   * is available from the restart.
   */
  void set_initial_accel_avail (bool initial_accel_set);

  /**
   * Error convergence order: 2 for \f$\gamma=0.5\f$, 1 otherwise
   */
  virtual Real error_order() const libmesh_override;

  /**
   * This method solves for the solution at the next timestep.
   * Usually we will only need to solve one (non)linear system per timestep,
   * but more complex subclasses may override this.
   */
  virtual void solve () libmesh_override;

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
   * Need to be able to indicate to _general_residual if we
   * are doing an acceleration solve or not.
   */
  bool _is_accel_solve;


  /**
   * This method requires an initial acceleration. So, we force the
   * user to call either compute_initial_accel or project_initial_accel
   * to set the initial acceleration.
   */
  bool _initial_accel_set;

  /**
   * This method is the underlying implementation of the public
   * residual methods.
   */
  virtual bool _general_residual (bool request_jacobian,
                                  DiffContext &,
                                  ResFuncType mass,
                                  ResFuncType damping,
                                  ResFuncType time_deriv,
                                  ResFuncType constraint);
};

} // namespace libMesh

#endif // LIBMESH_NEWMARK_SOLVER_H
