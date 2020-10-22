// The libMesh Finite Element Library.
// Copyright (C) 2002-2020 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_UNSTEADY_SOLVER_H
#define LIBMESH_UNSTEADY_SOLVER_H

// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/time_solver.h"

// C++ includes
#include <memory>

namespace libMesh
{

// Forward declarations
template <typename T> class NumericVector;

/**
 * This is a generic class that defines a solver to handle
 * time integration of DifferentiableSystems.
 *
 * A user can define a solver for unsteady problems by deriving
 * from this class and implementing certain functions.
 *
 * This class is part of the new DifferentiableSystem framework,
 * which is still experimental.  Users of this framework should
 * beware of bugs and future API changes.
 *
 * \author Roy H. Stogner
 * \date 2008
 */
class UnsteadySolver : public TimeSolver
{
public:
  /**
   * Constructor. Requires a reference to the system
   * to be solved.
   */
  explicit
  UnsteadySolver (sys_type & s);

  /**
   * Destructor.
   */
  virtual ~UnsteadySolver ();

  /**
   * The initialization function.  This method is used to
   * initialize internal data structures before a simulation begins.
   */
  virtual void init () override;

  /**
   * The data initialization function.  This method is used to
   * initialize internal data structures after the underlying System
   * has been initialized
   */
  virtual void init_data () override;

  /**
   * The reinitialization function.  This method is used to
   * resize internal data vectors after a mesh change.
   */
  virtual void reinit () override;

  /**
   * This method solves for the solution at the next timestep.
   * Usually we will only need to solve one (non)linear system per timestep,
   * but more complex subclasses may override this.
   */
  virtual void solve () override;

  /**
   * This method advances the solution to the next timestep, after a
   * solve() has been performed.  Often this will be done after every
   * UnsteadySolver::solve(), but adaptive mesh refinement and/or adaptive
   * time step selection may require some solve() steps to be repeated.
   */
  virtual void advance_timestep () override;

  void update();

  /**
   * This method solves for the adjoint solution at the next adjoint timestep
   * (or a steady state adjoint solve)
   */
  virtual std::pair<unsigned int, Real> adjoint_solve (const QoISet & qoi_indices) override;

  /**
   * This method advances the adjoint solution to the previous
   * timestep, after an adjoint_solve() has been performed.  This will
   * be done before every UnsteadySolver::adjoint_solve().
   */
  virtual void adjoint_advance_timestep () override;

  /**
   * This method retrieves all the stored solutions at the current
   * system.time
   */
  virtual void retrieve_timestep () override;

  /**
   * A method to integrate the system::QoI functionals.
   */
  virtual void integrate_qoi_timestep() override;

  /**
   * A method to integrate the adjoint sensitivity w.r.t a given parameter
   * vector. int_{tstep_start}^{tstep_end} dQ/dp dt = int_{tstep_start}^{tstep_end} (\partialQ / \partial p) - ( \partial R (u,z) / \partial p ) dt
   * The trapezoidal rule is used to numerically integrate the timestep.
   */
  virtual void integrate_adjoint_sensitivity(const QoISet & qois, const ParameterVector & parameter_vector, SensitivityData & sensitivities) override;

  /**
   * A method to compute the adjoint refinement error estimate at the current timestep.
   * int_{tstep_start}^{tstep_end} R(u^h,z) dt
   * The user provides an initialized ARefEE object.
   * Fills in an ErrorVector that contains the weighted sum of errors from all the QoIs and can be used to guide AMR.
   * CURRENTLY ONLY SUPPORTED for Backward Euler.
   */
  virtual void integrate_adjoint_refinement_error_estimate(AdjointRefinementEstimator & /*adjoint_refinement_error_estimator*/, ErrorVector & /*QoI_elementwise_error*/) override;

  /**
   * This method should return the expected convergence order of the
   * (non-local) error of the time discretization scheme - e.g. 2 for the
   * O(deltat^2) Crank-Nicholson, or 1 for the O(deltat) Backward Euler.
   *
   * Useful for adaptive timestepping schemes.
   */
  virtual Real error_order () const = 0;

  /**
   * \returns The maximum order of time derivatives for which the
   * UnsteadySolver subclass is capable of handling.
   *
   * For example, EulerSolver will have \p time_order() = 1 and
   * NewmarkSolver will have \p time_order() = 2.
   */
  virtual unsigned int time_order () const = 0;

  /**
   * \returns The old nonlinear solution for the specified global
   * DOF.
   */
  Number old_nonlinear_solution (const dof_id_type global_dof_number) const;

  /**
   * Serial vector of _system.get_vector("_old_nonlinear_solution")
   * This is a shared_ptr so that it can be shared between different
   * derived class instances, as in e.g. AdaptiveTimeSolver.
   */
  std::shared_ptr<NumericVector<Number>> old_local_nonlinear_solution;

  /**
   * Computes the size of ||u^{n+1} - u^{n}|| in some norm.
   *
   * \note While you can always call this function, its
   * result may or may not be very meaningful.  For example, if
   * you call this function right after calling advance_timestep()
   * then you'll get a result of zero since old_nonlinear_solution
   * is set equal to nonlinear_solution in this function.
   */
  virtual Real du(const SystemNorm & norm) const override;

  /**
   * This is not a steady-state solver.
   */
  virtual bool is_steady() const override { return false; }

  /**
   * A setter for the first_adjoint_step boolean. Needed for nested time solvers.
   */
  void set_first_adjoint_step(bool first_adjoint_step_setting)
  {
    first_adjoint_step = first_adjoint_step_setting;
  }

protected:

  /**
   * A bool that will be true the first time solve() is called,
   * and false thereafter
   */
  bool first_solve;

  /**
   * A bool that will be true the first time adjoint_advance_timestep() is called,
   * (when the primal solution is to be used to set adjoint boundary conditions) and false thereafter
   */
  bool first_adjoint_step;

  /**
   * A vector of pointers to vectors holding the adjoint solution at the last time step
   */
  std::vector< std::unique_ptr<NumericVector<Number>> > old_adjoints;

  /**
   * We will need to move the system.time around to ensure that residuals
   * are built with the right deltat and the right time.
   */
  Real last_step_deltat;
  Real next_step_deltat;
};



} // namespace libMesh


#endif // LIBMESH_UNSTEADY_SOLVER_H
