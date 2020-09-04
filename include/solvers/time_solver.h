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



#ifndef LIBMESH_TIME_SOLVER_H
#define LIBMESH_TIME_SOLVER_H

// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/reference_counted_object.h"

// C++ includes
#include <memory>

namespace libMesh
{

// Forward Declarations
class DiffContext;
class DiffSolver;
class DifferentiablePhysics;
class DifferentiableSystem;
class ParameterVector;
class SensitivityData;
class SolutionHistory;
class SystemNorm;
class QoISet;

template <typename T>
class LinearSolver;

/**
 * This is a generic class that defines a solver to handle
 * time integration of DifferentiableSystems.
 *
 * A user can define a solver by deriving from this class and
 * implementing certain functions.
 *
 * This class is part of the new DifferentiableSystem framework,
 * which is still experimental.  Users of this framework should
 * beware of bugs and future API changes.
 *
 * \author Roy H. Stogner
 * \date 2006
 */
class TimeSolver : public ReferenceCountedObject<TimeSolver>
{
public:
  /**
   * The type of system
   */
  typedef DifferentiableSystem sys_type;

  /**
   * Constructor. Requires a reference to the system
   * to be solved.
   */
  explicit
  TimeSolver (sys_type & s);

  /**
   * Destructor.
   */
  virtual ~TimeSolver ();

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
   * The reinitialization function.  This method is used after
   * changes in the mesh
   */
  virtual void reinit ();

  /**
   * This method solves for the solution at the next timestep (or
   * solves for a steady-state solution).  Usually we will only need
   * to solve one (non)linear system per timestep, but more complex
   * subclasses may override this.
   */
  virtual void solve ();

  /**
   * This method advances the solution to the next timestep, after a
   * solve() has been performed.  Often this will be done after every
   * UnsteadySolver::solve(), but adaptive mesh refinement and/or adaptive
   * time step selection may require some solve() steps to be repeated.
   */
  virtual void advance_timestep ();

  /**
   * This method solves for the adjoint solution at the next adjoint timestep
   * (or a steady state adjoint solve)
   */
  virtual std::pair<unsigned int, Real> adjoint_solve (const QoISet & qoi_indices);

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
  virtual void retrieve_timestep();

  /**
   * A method to integrate the adjoint sensitivity w.r.t a given parameter
   * vector. int_{tstep_start}^{tstep_end} dQ/dp dt = int_{tstep_start}^{tstep_end} (\partialQ / \partial p) - ( \partial R (u,z) / \partial p ) dt
   */
  virtual void integrate_adjoint_sensitivity(const QoISet & qois, const ParameterVector & parameter_vector, SensitivityData & sensitivities);

  /**
   * This method uses the DifferentiablePhysics
   * element_time_derivative(), element_constraint(), and
   * mass_residual() to build a full residual on an element.  What
   * combination
   *
   * it uses will depend on the type of solver.  See
   * the subclasses for more details.
   */
  virtual bool element_residual (bool request_jacobian,
                                 DiffContext &) = 0;

  /**
   * This method uses the DifferentiablePhysics
   * side_time_derivative(), side_constraint(), and
   * side_mass_residual() to build a full residual on an element's
   * side.
   *
   * What combination it uses will depend on the type
   * of solver.  See the subclasses for more details.
   */
  virtual bool side_residual (bool request_jacobian,
                              DiffContext &) = 0;

  /**
   * This method uses the DifferentiablePhysics
   * nonlocal_time_derivative(), nonlocal_constraint(), and
   * nonlocal_mass_residual() to build a full residual of non-local
   * terms.
   *
   * What combination it uses will depend on the type
   * of solver.  See the subclasses for more details.
   */
  virtual bool nonlocal_residual (bool request_jacobian,
                                  DiffContext &) = 0;

  /**
   * This method is for subclasses or users to override
   * to do arbitrary processing between timesteps
   */
  virtual void before_timestep () {}

  /**
   * \returns A constant reference to the system we are solving.
   */
  const sys_type & system () const { return _system; }

  /**
   * \returns A writable reference to the system we are solving.
   */
  sys_type & system () { return _system; }

  /**
   * An implicit linear or nonlinear solver to use at each timestep.
   */
  virtual std::unique_ptr<DiffSolver> & diff_solver() { return _diff_solver; }

  /**
   * An implicit linear solver to use for adjoint and sensitivity problems.
   */
  virtual std::unique_ptr<LinearSolver<Number>> & linear_solver() { return _linear_solver; }

  /**
   * Print extra debugging information if quiet ==  false.
   */
  bool quiet;

  /**
   * Computes the size of ||u^{n+1} - u^{n}|| in some norm.
   *
   * \note While you can always call this function, its
   * result may or may not be very meaningful.  For example, if
   * you call this function right after calling advance_timestep()
   * then you'll get a result of zero since old_nonlinear_solution
   * is set equal to nonlinear_solution in this function.
   */
  virtual Real du(const SystemNorm & norm) const = 0;

  /**
   * Is this effectively a steady-state solver?
   */
  virtual bool is_steady() const = 0;

  /**
   * This value (which defaults to zero) is the number of times the
   * TimeSolver is allowed to halve deltat and let the DiffSolver
   * repeat the latest failed solve with a reduced timestep.
   *
   * \note This has no effect for SteadySolvers.
   * \note You must set at least one of the DiffSolver flags
   * "continue_after_max_iterations" or
   * "continue_after_backtrack_failure" to allow the TimeSolver to
   * retry the solve.
   */
  unsigned int reduce_deltat_on_diffsolver_failure;

  /**
   * A setter function users will employ if they need to do something
   * other than save no solution history
   */
  void set_solution_history(const SolutionHistory & _solution_history);

  /**
￼   * A getter function that returns a reference to the solution history
￼   * object owned by TimeSolver
￼   * */
  SolutionHistory & get_solution_history();

  /**
   * Accessor for querying whether we need to do a primal
   * or adjoint solve
   */
  bool is_adjoint() const
  { return _is_adjoint; }

  /**
   * Accessor for setting whether we need to do a primal
   * or adjoint solve
   */
  void set_is_adjoint(bool _is_adjoint_value)
  { _is_adjoint = _is_adjoint_value; }

  /**
   * Function to return 'last_deltat()', returns system.deltat if
   * fixed timestep solver is used, completed_deltat if the adaptive
   * time solver is used.
   */
  virtual Real last_complete_deltat();


protected:

  /**
   * An implicit linear or nonlinear solver to use at each timestep.
   */
  std::unique_ptr<DiffSolver> _diff_solver;

  /**
   * An implicit linear solver to use for adjoint problems.
   */
  std::unique_ptr<LinearSolver<Number>> _linear_solver;

  /**
   * A reference to the system we are solving.
   */
  sys_type & _system;

  /**
   * A std::unique_ptr to a SolutionHistory object. Default is
   * NoSolutionHistory, which the user can override by declaring a
   * different kind of SolutionHistory in the application
   */
  std::unique_ptr<SolutionHistory> solution_history;

  /**
   * Definitions of argument types for use in refactoring subclasses.
   */

  typedef bool (DifferentiablePhysics::*ResFuncType) (bool, DiffContext &);

  typedef void (DiffContext::*ReinitFuncType) (Real);

private:

  /**
   * This boolean tells the TimeSolver whether we are solving a primal or
   * adjoint problem
   */
  bool _is_adjoint;

};


} // namespace libMesh


#endif // LIBMESH_TIME_SOLVER_H
