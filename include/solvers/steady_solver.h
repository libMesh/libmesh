// The libMesh Finite Element Library.
// Copyright (C) 2002-2021 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
                                 DiffContext &) override;

  /**
   * This method uses the DifferentiablePhysics'
   * side_time_derivative() and side_constraint()
   * to build a full residual/jacobian on an element's side.
   */
  virtual bool side_residual (bool request_jacobian,
                              DiffContext &) override;

  /**
   * This method uses the DifferentiablePhysics'
   * nonlocal_time_derivative() and nonlocal_constraint()
   * to build a full residual/jacobian for non-local terms.
   */
  virtual bool nonlocal_residual (bool request_jacobian,
                                  DiffContext &) override;

  /**
   * \returns 0, but derived classes should override this function to
   * compute the size of the difference between successive solution
   * iterates ||u^{n+1} - u^{n}|| in some norm.
   */
  virtual Real du(const SystemNorm &) const override { return 0; }

  /**
   * This is a steady-state solver.
   */
  virtual bool is_steady() const override { return true; }

  /**
   * A method to integrate the system::QoI functionals.
   */
  virtual void integrate_qoi_timestep() override;

  /**
   * A method to integrate the adjoint sensitivity w.r.t a given parameter
   * vector. int_{tstep_start}^{tstep_end} dQ/dp dt = int_{tstep_start}^{tstep_end} (\partialQ / \partial p) - ( \partial R (u,z) / \partial p ) dt
   */
  virtual void integrate_adjoint_sensitivity(const QoISet & qois, const ParameterVector & parameter_vector, SensitivityData & sensitivities) override;

#ifdef LIBMESH_ENABLE_AMR
  /**
   * A method to compute the adjoint refinement error estimate at the current timestep.
   * int_{tstep_start}^{tstep_end} R(u^h,z) dt
   * The user provides an initialized ARefEE object.
   * Fills in an ErrorVector that contains the weighted sum of errors from all the QoIs and can be used to guide AMR.
   * CURRENTLY ONLY SUPPORTED for Backward Euler.
   */
  virtual void integrate_adjoint_refinement_error_estimate(AdjointRefinementEstimator & adjoint_refinement_error_estimator, ErrorVector & QoI_elementwise_error) override;
#endif // LIBMESH_ENABLE_AMR

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
