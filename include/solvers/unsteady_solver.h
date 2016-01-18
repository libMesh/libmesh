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



#ifndef LIBMESH_UNSTEADY_SOLVER_H
#define LIBMESH_UNSTEADY_SOLVER_H

// Local includes
#include "libmesh/auto_ptr.h"
#include "libmesh/libmesh_common.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/time_solver.h"

// C++ includes

namespace libMesh
{

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
  virtual void init () libmesh_override;

  /**
   * The data initialization function.  This method is used to
   * initialize internal data structures after the underlying System
   * has been initialized
   */
  virtual void init_data () libmesh_override;

  /**
   * The reinitialization function.  This method is used to
   * resize internal data vectors after a mesh change.
   */
  virtual void reinit () libmesh_override;

  /**
   * This method solves for the solution at the next timestep.
   * Usually we will only need to solve one (non)linear system per timestep,
   * but more complex subclasses may override this.
   */
  virtual void solve () libmesh_override;

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
   * This method retrieves all the stored solutions at the current
   * system.time
   */
  virtual void retrieve_timestep () libmesh_override;

  /**
   * This method should return the expected convergence order of the
   * (non-local) error of the time discretization scheme - e.g. 2 for the
   * O(deltat^2) Crank-Nicholson, or 1 for the O(deltat) Backward Euler.
   *
   * Useful for adaptive timestepping schemes.
   */
  virtual Real error_order () const = 0;

  /**
   * Returns the maximum order of time derivatives for which the
   * UnsteadySolver subclass is capable of handling. E.g. EulerSolver
   * will have time_order = 1 and NewmarkSolver will have time_order = 2
   */
  virtual unsigned int time_order () const = 0;

  /**
   * @returns the old nonlinear solution for the specified global
   * DOF.
   */
  Number old_nonlinear_solution (const dof_id_type global_dof_number) const;

  /**
   * Serial vector of _system.get_vector("_old_nonlinear_solution")
   */
  UniquePtr<NumericVector<Number> > old_local_nonlinear_solution;

  /**
   * Computes the size of ||u^{n+1} - u^{n}|| in some norm.
   *
   * Note that, while you can always call this function, its
   * result may or may not be very meaningful.  For example, if
   * you call this function right after calling advance_timestep()
   * then you'll get a result of zero since old_nonlinear_solution
   * is set equal to nonlinear_solution in this function.
   */
  virtual Real du(const SystemNorm & norm) const libmesh_override;

  /**
   * This is not a steady-state solver.
   */
  virtual bool is_steady() const libmesh_override { return false; }

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
};



} // namespace libMesh


#endif // LIBMESH_UNSTEADY_SOLVER_H
