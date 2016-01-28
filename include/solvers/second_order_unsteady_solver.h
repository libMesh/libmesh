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

#ifndef LIBMESH_SECOND_ORDER_UNSTEADY_SOLVER_H
#define LIBMESH_SECOND_ORDER_UNSTEADY_SOLVER_H

#include "libmesh/unsteady_solver.h"

namespace libMesh
{
/**
 * Generic class from which second order UnsteadySolvers should subclass.
 *
 * Subclasses of this class are meant to solve problems of the form
 * \f[ M(u)\ddot{u} + C(u)\dot{u} + F(u) = 0 \f]
 *
 * This class is part of the new DifferentiableSystem framework,
 * which is still experimental.  Users of this framework should
 * beware of bugs and future API changes.
 *
 * \author Paul T. Bauman
 * \date 2015
 */
class SecondOrderUnsteadySolver : public UnsteadySolver
{
public:
  /**
   * Constructor. Requires a reference to the system
   * to be solved.
   */
  explicit
  SecondOrderUnsteadySolver (sys_type & s);

  /**
   * Destructor.
   */
  virtual ~SecondOrderUnsteadySolver ();

  virtual unsigned int time_order() const libmesh_override
  { return 2; }

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
   * This method retrieves all the stored solutions at the current
   * system.time
   */
  virtual void retrieve_timestep () libmesh_override;

  /**
   * Specify non-zero initial velocity. Should be called before solve().
   * The function value f and its gradient g are user-provided cloneable functors.
   * A gradient g is only required/used for projecting onto finite element spaces
   * with continuous derivatives.
   */
  void project_initial_rate(FunctionBase<Number> * f,
                            FunctionBase<Gradient> * g = libmesh_nullptr);

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
   * Serial vector of previous time step velocity \f$ \dot{u}_n \f$
   */
  UniquePtr<NumericVector<Number> > _old_local_solution_rate;

  /**
   * Serial vector of previous time step accleration \f$ \ddot{u}_n \f$
   */
  UniquePtr<NumericVector<Number> > _old_local_solution_accel;
};

} // end namespace libMesh

# endif // LIBMESH_SECOND_ORDER_UNSTEADY_SOLVER_H
