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
  untested();
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



void DifferentiableSystem::solve ()
{
  time_solver->solve();
}



void DifferentiableSystem::adjoint_solve ()
{
  time_solver->adjoint_solve();
}



void DifferentiableSystem::qoi_parameter_sensitivity
  (std::vector<Number *>& parameters,
   std::vector<Number>&   sensitivities)
{
  // Get ready to fill in senstivities:
  sensitivities.clear();
  sensitivities.resize(parameters.size(), 0);

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

  this->adjoint_solve();

  // We use the identities:
  // dq/dp = (partial q / partial p) + (partial q / partial u) *
  //         (partial u / partial p)
  // dq/dp = (partial q / partial p) + (J^T * z) *
  //         (partial u / partial p)
  // dq/dp = (partial q / partial p) + z * J *
  //         (partial u / partial p)
 
  // Leading to our final formula:
  // dq/dp = (partial q / partial p) - z * (partial R / partial p)

  for (unsigned int i=0; i != parameters.size(); ++i)
    {
      // We currently get partial derivatives via central differencing
      Number delta_p = 1e-6;

      // (partial q / partial p) ~= (q(p+dp)-q(p-dp))/(2*dp)
      // (partial R / partial p) ~= (rhs(p+dp) - rhs(p-dp))/(2*dp)

      Number old_parameter = *parameters[i];
      // Number old_qoi = this->qoi;

      *parameters[i] = old_parameter - delta_p;
      this->assemble_qoi();
      Number qoi_minus = this->qoi;

      this->assembly(true, false);
      this->rhs->close();
      AutoPtr<NumericVector<Number> > partialR_partialp = this->rhs->clone();
      *partialR_partialp *= -1;

      *parameters[i] = old_parameter + delta_p;
      this->assemble_qoi();
      Number qoi_plus = this->qoi;
      Number partialq_partialp = (qoi_plus - qoi_minus) / (2.*delta_p);

      this->assembly(true, false);
      this->rhs->close();
      *partialR_partialp += *this->rhs;
      *partialR_partialp /= (2.*delta_p);

      // Don't leave the parameter changed
      *parameters[i] = old_parameter;

      sensitivities[i] = partialq_partialp -
			 partialR_partialp->dot(this->get_adjoint_solution());
    }

  // All parameters have been reset.
  // Don't leave the qoi or system changed - principle of least
  // surprise.
  this->assembly(true, false);
  this->rhs->close();
  this->assemble_qoi();
}
