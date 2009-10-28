// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
  
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


#include "diff_system.h"
#include "dof_map.h"
#include "numeric_vector.h"
#include "parameter_vector.h"
#include "qoi_set.h"
#include "sensitivity_data.h"
#include "time_solver.h"



DifferentiableSystem::DifferentiableSystem
                      (EquationSystems& es,
                       const std::string& name,
                       const unsigned int number) :
  Parent      (es, name, number),
  compute_internal_sides(false),
  postprocess_sides(false),
  time_solver (NULL),
  time(0.),
  deltat(1.),
  print_solution_norms(false),
  print_solutions(false),
  print_residual_norms(false),
  print_residuals(false),
  print_jacobian_norms(false),
  print_jacobians(false),
  print_element_jacobians(false),
  use_fixed_solution(false)
{
}



DifferentiableSystem::~DifferentiableSystem ()
{
  this->clear();
}



void DifferentiableSystem::clear ()
{
  _time_evolving.resize(0);

  use_fixed_solution = false;
}



void DifferentiableSystem::reinit ()
{
  Parent::reinit();

  time_solver->reinit();
}



void DifferentiableSystem::init_data ()
{
  // give us flags for every variable that might be time evolving
  _time_evolving.resize(this->n_vars(), false);

  // Do any initialization our solvers need
  libmesh_assert (time_solver.get() != NULL);
  time_solver->init();

  // Next initialize ImplicitSystem data
  Parent::init_data();
}



AutoPtr<DiffContext> DifferentiableSystem::build_context ()
{
  return AutoPtr<DiffContext>(new DiffContext(*this));
}



void DifferentiableSystem::assemble ()
{
  this->assembly(true, true);
}



void DifferentiableSystem::assemble_residual_derivatives(const ParameterVector& parameters)
{
  const unsigned int Np = parameters.size();
  Real deltap = TOLERANCE;

  for (unsigned int p=0; p != Np; ++p)
    {
      NumericVector<Number> &sensitivity_rhs = this->add_sensitivity_rhs(p);

      // Approximate -(partial R / partial p) by
      // (R(p-dp) - R(p+dp)) / (2*dp)

      Number old_parameter = *parameters[p];
      *parameters[p] -= deltap;

      this->assembly(true, false);
      this->rhs->close();
      sensitivity_rhs = *this->rhs;

      *parameters[p] = old_parameter + deltap;

      this->assembly(true, false);
      this->rhs->close();

      sensitivity_rhs -= *this->rhs;
      sensitivity_rhs /= (2*deltap);
      sensitivity_rhs.close();

      *parameters[p] = old_parameter;
    }
}



void DifferentiableSystem::solve ()
{
  time_solver->solve();
}



void DifferentiableSystem::sensitivity_solve (const ParameterVector& parameters)
{
  time_solver->sensitivity_solve(parameters);
}



void DifferentiableSystem::adjoint_solve (const QoISet& qoi_indices)
{
  time_solver->adjoint_solve(qoi_indices);
}



void DifferentiableSystem::adjoint_qoi_parameter_sensitivity
  (const QoISet& qoi_indices,
   const ParameterVector& parameters,
   SensitivityData& sensitivities)
{
  const unsigned int Np = parameters.size();
  const unsigned int Nq = qoi.size();

  // An introduction to the problem:
  //
  // Residual R(u(p),p) = 0
  // partial R / partial u = J = system matrix
  //
  // This implies that:
  // d/dp(R) = 0
  // (partial R / partial p) + 
  // (partial R / partial u) * (partial u / partial p) = 0

  // We first do an adjoint solve:
  // J^T * z = (partial q / partial u)

  this->adjoint_solve(qoi_indices);

  // Get ready to fill in senstivities:
  sensitivities.allocate_data(qoi_indices, *this, parameters);

  // We use the identities:
  // dq/dp = (partial q / partial p) + (partial q / partial u) *
  //         (partial u / partial p)
  // dq/dp = (partial q / partial p) + (J^T * z) *
  //         (partial u / partial p)
  // dq/dp = (partial q / partial p) + z * J *
  //         (partial u / partial p)
 
  // Leading to our final formula:
  // dq/dp = (partial q / partial p) - z * (partial R / partial p)

  for (unsigned int j=0; j != Np; ++j)
    {
      // We currently get partial derivatives via central differencing
      Number delta_p = 1e-6;

      // (partial q / partial p) ~= (q(p+dp)-q(p-dp))/(2*dp)
      // (partial R / partial p) ~= (rhs(p+dp) - rhs(p-dp))/(2*dp)

      Number old_parameter = *parameters[j];
      // Number old_qoi = this->qoi;

      *parameters[j] = old_parameter - delta_p;
      this->assemble_qoi(qoi_indices);
      std::vector<Number> qoi_minus = this->qoi;

      this->assembly(true, false);
      this->rhs->close();
      AutoPtr<NumericVector<Number> > partialR_partialp = this->rhs->clone();
      *partialR_partialp *= -1;

      *parameters[j] = old_parameter + delta_p;
      this->assemble_qoi(qoi_indices);
      std::vector<Number>& qoi_plus = this->qoi;

      std::vector<Number> partialq_partialp(Nq, 0);
      for (unsigned int i=0; i != Nq; ++i)
        if (qoi_indices.has_index(i))
          partialq_partialp[i] = (qoi_plus[i] - qoi_minus[i]) / (2.*delta_p);

      this->assembly(true, false);
      this->rhs->close();
      *partialR_partialp += *this->rhs;
      *partialR_partialp /= (2.*delta_p);

      // Don't leave the parameter changed
      *parameters[j] = old_parameter;

      for (unsigned int i=0; i != Nq; ++i)
        if (qoi_indices.has_index(i))
          sensitivities[i][j] = partialq_partialp[i] -
            partialR_partialp->dot(this->get_adjoint_solution(i));
    }

  // All parameters have been reset.
  // Don't leave the qoi or system changed - principle of least
  // surprise.
  this->assembly(true, false);
  this->rhs->close();
  this->assemble_qoi();
}



void DifferentiableSystem::forward_qoi_parameter_sensitivity
  (const QoISet& qoi_indices,
   const ParameterVector& parameters,
   SensitivityData& sensitivities)
{
  const unsigned int Np = parameters.size();
  const unsigned int Nq = qoi.size();

  // An introduction to the problem:
  //
  // Residual R(u(p),p) = 0
  // partial R / partial u = J = system matrix
  //
  // This implies that:
  // d/dp(R) = 0
  // (partial R / partial p) + 
  // (partial R / partial u) * (partial u / partial p) = 0

  // We first solve for (partial u / partial p) for each parameter:
  // J * (partial u / partial p) = - (partial R / partial p)

  this->sensitivity_solve(parameters);

  // Get ready to fill in senstivities:
  sensitivities.allocate_data(qoi_indices, *this, parameters);

  // We use the identity
  // dq/dp = (partial q / partial p) + (partial q / partial u) *
  //         (partial u / partial p)

  // We get (partial q / partial u) from the user
  this->assemble_qoi_derivative(qoi_indices);
 
  for (unsigned int j=0; j != Np; ++j)
    {
      // We currently get partial derivatives via central differencing
      Number delta_p = 1e-6;

      // (partial q / partial p) ~= (q(p+dp)-q(p-dp))/(2*dp)

      Number old_parameter = *parameters[j];

      *parameters[j] = old_parameter - delta_p;
      this->assemble_qoi(qoi_indices);
      std::vector<Number> qoi_minus = this->qoi;

      *parameters[j] = old_parameter + delta_p;
      this->assemble_qoi(qoi_indices);
      std::vector<Number>& qoi_plus = this->qoi;

      std::vector<Number> partialq_partialp(Nq, 0);
      for (unsigned int i=0; i != Nq; ++i)
        if (qoi_indices.has_index(i))
          partialq_partialp[i] = (qoi_plus[i] - qoi_minus[i]) / (2.*delta_p);

      // Don't leave the parameter changed
      *parameters[j] = old_parameter;

      for (unsigned int i=0; i != Nq; ++i)
        if (qoi_indices.has_index(i))
          sensitivities[i][j] = partialq_partialp[i] +
            this->get_adjoint_rhs(i).dot(this->get_sensitivity_solution(i));
    }

  // All parameters have been reset.
  // Don't leave the qoi or system changed - principle of least
  // surprise.
  this->assembly(true, false);
  this->rhs->close();
  this->assemble_qoi(qoi_indices);
}
