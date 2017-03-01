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

#ifndef LIBMESH_FIRST_ORDER_UNSTEADY_SOLVER_H
#define LIBMESH_FIRST_ORDER_UNSTEADY_SOLVER_H

#include "libmesh/unsteady_solver.h"

namespace libMesh
{
/**
 * Generic class from which first order UnsteadySolvers should subclass.
 *
 * Subclasses of this class are meant to solve problems of the form
 * \f[ M(u)\dot{u} = F(u)\f]
 *
 * There is also infrastructure for subclasses to support solving second
 * order-in-time systems (or both first order and second order).
 * In particular, consider the second system
 * \f[ M(u)\ddot{u} + C(u)\dot{u} + F(u) = 0 \f]
 * If we introduce the equation \f$\dot{u} = v\f$, then, we have the
 * following first order system:
 * \f{eqnarray}{ \dot{u} &=& v \\ M(u)\dot{v} + C(u)v + F(u) &=& 0\f}
 * Such systems are supported by the user specifying that the time order
 * of the variable \f$u\f$ is 2 when adding the variable using
 * FEMSystem::time_evolving. This class will then add the \f$v\f$ variables
 * to the FEMSystem (named "dot_u" if the second order variable name is "u").
 * Sublasses will then need to make sure to use the function
 * FirstOrderUnsteadySolver::prepare_accel() to prepare data structures for the
 * users to use from FEMContext. Furthermore, subclasses will need to appropriately
 * call damping_residual functions. Finally, subclasses will call
 * FirstOrderUnsteadySolver::compute_second_order_eqns() to actually assemble
 * residual and Jacobians for the \f$ \dot{u} = v\f$ equations. These aspects
 * should be invisible ot users during their element assembly.
 *
 * Unfortunately, complete usage of the feature of treating second order
 * equations as a system of first order equations is not completely invisible
 * to the user. The user must assemble their residual in the correct equation.
 * In particular, the residual must go in the \f$v\f$ residual equation, i.e.
 * \f$ M(u)\dot{v} + C(u)v + F(u) = 0 \f$, and,
 * subsequently, they must also be careful to populate to the correct Jacobian
 * blocks. The function to help facilitate this is get_second_order_dot_var.
 * If you have a second order variable, you pass that variable index and
 * get_second_order_dot_var will return the index of the corresponding velocity
 * variable, \f$v\f$ in the notation above. Then, the user knows what blocks to
 * populate. Note that the API is designed so that if the user codes to this
 * paradigm, the TimeSolver will be interchangeable for those element kernels.
 * That is, they'll be able to use either a FirstOrderUnsteadySolver or a
 * SecondOrderUnsteadySolver.
 *
 * This class is part of the new DifferentiableSystem framework,
 * which is still experimental.  Users of this framework should
 * beware of bugs and future API changes.
 *
 * \author Paul T. Bauman
 * \date 2015
 */
class FirstOrderUnsteadySolver : public UnsteadySolver
{
public:
  /**
   * Constructor. Requires a reference to the system
   * to be solved.
   */
  explicit
  FirstOrderUnsteadySolver (sys_type & s)
    : UnsteadySolver(s) {}

  /**
   * Destructor.
   */
  virtual ~FirstOrderUnsteadySolver (){}

  virtual unsigned int time_order() const libmesh_override
  { return 1; }

protected:

  /**
   * If there are second order variables in the system,
   * then we also prepare the accel for those variables
   * so the user can treat them as such.
   */
  void prepare_accel(DiffContext & context);

  /**
   * If there are second order variables, then we need to compute their residual equations
   * and corresponding Jacobian. The residual equation will simply be
   * \f$ \dot{u} - v = 0 \f$, where \f$ u \f$ is the second order variable add
   * by the user and \f$ v \f$ is the variable added by the time-solver as the
   * "velocity" variable.
   */
  bool compute_second_order_eqns(bool compute_jacobian, DiffContext & c);

};

} // end namespace libMesh

#endif // LIBMESH_FIRST_ORDER_UNSTEADY_SOLVER_H
