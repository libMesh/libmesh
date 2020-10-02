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


#include "libmesh/diff_solver.h"
#include "libmesh/diff_system.h"
#include "libmesh/steady_solver.h"

#include "libmesh/adjoint_refinement_estimator.h"
#include "libmesh/error_vector.h"

namespace libMesh
{



SteadySolver::~SteadySolver ()
{
}



bool SteadySolver::element_residual(bool request_jacobian,
                                    DiffContext & context)
{
  return this->_general_residual(request_jacobian,
                                 context,
                                 &DifferentiablePhysics::element_time_derivative,
                                 &DifferentiablePhysics::element_constraint);
}



bool SteadySolver::side_residual(bool request_jacobian,
                                 DiffContext & context)
{
  return this->_general_residual(request_jacobian,
                                 context,
                                 &DifferentiablePhysics::side_time_derivative,
                                 &DifferentiablePhysics::side_constraint);
}



bool SteadySolver::nonlocal_residual(bool request_jacobian,
                                     DiffContext & context)
{
  return this->_general_residual(request_jacobian,
                                 context,
                                 &DifferentiablePhysics::nonlocal_time_derivative,
                                 &DifferentiablePhysics::nonlocal_constraint);
}



bool SteadySolver::_general_residual(bool request_jacobian,
                                     DiffContext & context,
                                     ResFuncType time_deriv,
                                     ResFuncType constraint)
{
  // If a fixed solution is requested, it will just be the current
  // solution
  if (_system.use_fixed_solution)
    {
      context.get_elem_fixed_solution() = context.get_elem_solution();
      context.fixed_solution_derivative = 1.0;
    }

  bool jacobian_computed =
    (_system.get_physics()->*time_deriv)(request_jacobian, context);

  // The user shouldn't compute a jacobian unless requested
  libmesh_assert (request_jacobian || !jacobian_computed);

  bool jacobian_computed2 =
    (_system.get_physics()->*constraint)(jacobian_computed, context);

  // The user shouldn't compute a jacobian unless requested
  libmesh_assert (jacobian_computed || !jacobian_computed2);

  return jacobian_computed2;
}

void SteadySolver::integrate_qoi_timestep()
{
  this->_system.assemble_qoi();

  return;
}

void SteadySolver::integrate_adjoint_sensitivity(const QoISet & qois, const ParameterVector & parameter_vector, SensitivityData & sensitivities)
{
  this->_system.ImplicitSystem::adjoint_qoi_parameter_sensitivity(qois, parameter_vector, sensitivities);

  return;
}

void SteadySolver::integrate_adjoint_refinement_error_estimate
  (AdjointRefinementEstimator & adjoint_refinement_error_estimator,
   ErrorVector & QoI_elementwise_error,
   std::vector<Real *> /* QoI_time_instant */)
{
  // Make sure the system::qoi_error_estimates vector is of the same size as system::qoi
  if(_system.qoi_error_estimates.size() != _system.qoi.size())
      _system.qoi_error_estimates.resize(_system.qoi.size());

  // Base class assumes a direct steady state error estimate
  adjoint_refinement_error_estimator.estimate_error(_system, QoI_elementwise_error);

  // Also get the spatially integrated errors for all the QoIs in the QoI set
  for (auto j : make_range(_system.n_qois()))
  {
    // Skip this QoI if not in the QoI Set
    if (adjoint_refinement_error_estimator.qoi_set().has_index(j))
    {
      (_system.qoi_error_estimates)[j] = adjoint_refinement_error_estimator.get_global_QoI_error_estimate(j);
    }
  }

  return;
}

} // namespace libMesh
