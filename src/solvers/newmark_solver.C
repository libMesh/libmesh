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

#include "libmesh/diff_system.h"
#include "libmesh/newmark_solver.h"
#include "libmesh/dof_map.h"
#include "libmesh/diff_solver.h"

namespace libMesh
{
NewmarkSolver::NewmarkSolver (sys_type & s)
  : SecondOrderUnsteadySolver(s),
    _beta(0.25),
    _gamma(0.5),
    _is_accel_solve(false),
    _initial_accel_set(false)
{}

NewmarkSolver::~NewmarkSolver ()
{}

Real NewmarkSolver::error_order() const
{
  if (_gamma == 0.5)
    return 2.;
  return 1.;
}

void NewmarkSolver::advance_timestep ()
{
  // We need to update velocity and acceleration before
  // we update the nonlinear solution (displacement) and
  // delta_t

  NumericVector<Number> & old_solution_rate =
    _system.get_vector("_old_solution_rate");

  NumericVector<Number> & old_solution_accel =
    _system.get_vector("_old_solution_accel");

  if (!first_solve)
    {
      NumericVector<Number> & old_nonlinear_soln =
        _system.get_vector("_old_nonlinear_solution");

      NumericVector<Number> & nonlinear_solution =
        *(_system.solution);

      // We need to cache the new solution_rate before updating the old_solution_rate
      // so we can update acceleration with the proper old_solution_rate
      // v_{n+1} = gamma/(beta*Delta t)*(x_{n+1}-x_n)
      //         - ((gamma/beta)-1)*v_n
      //         - (gamma/(2*beta)-1)*(Delta t)*a_n
      std::unique_ptr<NumericVector<Number>> new_solution_rate = nonlinear_solution.clone();
      (*new_solution_rate) -= old_nonlinear_soln;
      (*new_solution_rate) *= (_gamma/(_beta*_system.deltat));
      new_solution_rate->add( (1.0-_gamma/_beta), old_solution_rate );
      new_solution_rate->add( (1.0-_gamma/(2.0*_beta))*_system.deltat, old_solution_accel );

      // a_{n+1} = (1/(beta*(Delta t)^2))*(x_{n+1}-x_n)
      //         - 1/(beta*Delta t)*v_n
      //         - (1-1/(2*beta))*a_n
      std::unique_ptr<NumericVector<Number>> new_solution_accel = old_solution_accel.clone();
      (*new_solution_accel) *=  -(1.0/(2.0*_beta)-1.0);
      new_solution_accel->add( -1.0/(_beta*_system.deltat), old_solution_rate );
      new_solution_accel->add( 1.0/(_beta*_system.deltat*_system.deltat), nonlinear_solution );
      new_solution_accel->add( -1.0/(_beta*_system.deltat*_system.deltat), old_nonlinear_soln );

      // Now update old_solution_rate
      old_solution_rate = (*new_solution_rate);
      old_solution_accel = (*new_solution_accel);
    }

  // Localize updated vectors
  old_solution_rate.localize
    (*_old_local_solution_rate,
     _system.get_dof_map().get_send_list());

  old_solution_accel.localize
    (*_old_local_solution_accel,
     _system.get_dof_map().get_send_list());

  // Now we can finish advancing the timestep
  UnsteadySolver::advance_timestep();
}

void NewmarkSolver::adjoint_advance_timestep ()
{
  libmesh_not_implemented();
}

void NewmarkSolver::compute_initial_accel()
{
  // We need to compute the initial acceleration based off of
  // the initial position and velocity and, thus, acceleration
  // is the unknown in diff_solver and not the displacement. So,
  // We swap solution and acceleration. NewmarkSolver::_general_residual
  // will check _is_accel_solve and provide the correct
  // values to the FEMContext assuming this swap was made.
  this->_is_accel_solve = true;

  //solution->accel, accel->solution
  _system.solution->swap(_system.get_vector("_old_solution_accel"));
  _system.update();

  this->_diff_solver->solve();

  // solution->solution, accel->accel
  _system.solution->swap(_system.get_vector("_old_solution_accel"));
  _system.update();

  // We're done, so no longer doing an acceleration solve
  this->_is_accel_solve = false;

  this->set_initial_accel_avail(true);
}

void NewmarkSolver::project_initial_accel(FunctionBase<Number> * f,
                                          FunctionBase<Gradient> * g)
{
  NumericVector<Number> & old_solution_accel =
    _system.get_vector("_old_solution_accel");

  _system.project_vector( old_solution_accel, f, g );

  this->set_initial_accel_avail(true);
}

void NewmarkSolver::set_initial_accel_avail( bool initial_accel_set )
{
  _initial_accel_set = initial_accel_set;
}

void NewmarkSolver::solve ()
{
  // First, check that the initial accel was set one way or another
  if (!_initial_accel_set)
    {
      std::string error = "ERROR: Must first set initial acceleration using one of:\n";
      error += "NewmarkSolver::compute_initial_accel()\n";
      error += "NewmarkSolver::project_initial_accel()\n";
      libmesh_error_msg(error);
    }

  // That satisfied, now we can solve
  UnsteadySolver::solve();
}

bool NewmarkSolver::element_residual (bool request_jacobian,
                                      DiffContext & context)
{
  return this->_general_residual(request_jacobian,
                                 context,
                                 &DifferentiablePhysics::mass_residual,
                                 &DifferentiablePhysics::damping_residual,
                                 &DifferentiablePhysics::_eulerian_time_deriv,
                                 &DifferentiablePhysics::element_constraint,
                                 &DiffContext::elem_reinit);
}



bool NewmarkSolver::side_residual (bool request_jacobian,
                                   DiffContext & context)
{
  return this->_general_residual(request_jacobian,
                                 context,
                                 &DifferentiablePhysics::side_mass_residual,
                                 &DifferentiablePhysics::side_damping_residual,
                                 &DifferentiablePhysics::side_time_derivative,
                                 &DifferentiablePhysics::side_constraint,
                                 &DiffContext::elem_side_reinit);
}



bool NewmarkSolver::nonlocal_residual (bool request_jacobian,
                                       DiffContext & context)
{
  return this->_general_residual(request_jacobian,
                                 context,
                                 &DifferentiablePhysics::nonlocal_mass_residual,
                                 &DifferentiablePhysics::nonlocal_damping_residual,
                                 &DifferentiablePhysics::nonlocal_time_derivative,
                                 &DifferentiablePhysics::nonlocal_constraint,
                                 &DiffContext::nonlocal_reinit);
}



bool NewmarkSolver::_general_residual (bool request_jacobian,
                                       DiffContext & context,
                                       ResFuncType mass,
                                       ResFuncType damping,
                                       ResFuncType time_deriv,
                                       ResFuncType constraint,
                                       ReinitFuncType reinit_func)
{
  unsigned int n_dofs = context.get_elem_solution().size();

  // We might need to save the old jacobian in case one of our physics
  // terms later is unable to update it analytically.
  DenseMatrix<Number> old_elem_jacobian(n_dofs, n_dofs);

  // Local velocity at old time step
  DenseVector<Number> old_elem_solution_rate(n_dofs);
  for (unsigned int i=0; i != n_dofs; ++i)
    old_elem_solution_rate(i) =
      old_solution_rate(context.get_dof_indices()[i]);

  // The user is computing the initial acceleration
  // So upstream we've swapped _system.solution and _old_local_solution_accel
  // So we need to give the context the correct entries since we're solving for
  // acceleration here.
  if (_is_accel_solve)
    {
      // System._solution is actually the acceleration right now so we need
      // to reset the elem_solution to the right thing, which in this case
      // is the initial guess for displacement, which got swapped into
      // _old_solution_accel vector
      DenseVector<Number> old_elem_solution(n_dofs);
      for (unsigned int i=0; i != n_dofs; ++i)
        old_elem_solution(i) =
          old_solution_accel(context.get_dof_indices()[i]);

      context.elem_solution_derivative = 0.0;
      context.elem_solution_rate_derivative = 0.0;
      context.elem_solution_accel_derivative = 1.0;

      // Acceleration is currently the unknown so it's already sitting
      // in elem_solution() thanks to FEMContext::pre_fe_reinit
      context.get_elem_solution_accel() = context.get_elem_solution();

      // Now reset elem_solution() to what the user is expecting
      context.get_elem_solution() = old_elem_solution;

      context.get_elem_solution_rate() = old_elem_solution_rate;

      // The user's Jacobians will be targeting derivatives w.r.t. u_{n+1}.
      // Although the vast majority of cases will have the correct analytic
      // Jacobians in this iteration, since we reset elem_solution_derivative*,
      // if there are coupled/overlapping problems, there could be
      // mismatches in the Jacobian. So we force finite differencing for
      // the first iteration.
      request_jacobian = false;
    }
  // Otherwise, the unknowns are the displacements and everything is straight
  // forward and is what you think it is
  else
    {
      if (request_jacobian)
        old_elem_jacobian.swap(context.get_elem_jacobian());

      // Local displacement at old timestep
      DenseVector<Number> old_elem_solution(n_dofs);
      for (unsigned int i=0; i != n_dofs; ++i)
        old_elem_solution(i) =
          old_nonlinear_solution(context.get_dof_indices()[i]);

      // Local acceleration at old time step
      DenseVector<Number> old_elem_solution_accel(n_dofs);
      for (unsigned int i=0; i != n_dofs; ++i)
        old_elem_solution_accel(i) =
          old_solution_accel(context.get_dof_indices()[i]);

      // Convenience
      libMesh::Real dt = _system.deltat;

      context.elem_solution_derivative = 1.0;

      // Local velocity at current time step
      // v_{n+1} = gamma/(beta*Delta t)*(x_{n+1}-x_n)
      //         + (1-(gamma/beta))*v_n
      //         + (1-gamma/(2*beta))*(Delta t)*a_n
      context.elem_solution_rate_derivative = (_gamma/(_beta*dt));

      context.get_elem_solution_rate()  = context.get_elem_solution();
      context.get_elem_solution_rate() -= old_elem_solution;
      context.get_elem_solution_rate() *= context.elem_solution_rate_derivative;
      context.get_elem_solution_rate().add( (1.0-_gamma/_beta), old_elem_solution_rate);
      context.get_elem_solution_rate().add( (1.0-_gamma/(2.0*_beta))*dt, old_elem_solution_accel);



      // Local acceleration at current time step
      // a_{n+1} = (1/(beta*(Delta t)^2))*(x_{n+1}-x_n)
      //         - 1/(beta*Delta t)*v_n
      //         - (1/(2*beta)-1)*a_n
      context.elem_solution_accel_derivative = 1.0/(_beta*dt*dt);

      context.get_elem_solution_accel()  = context.get_elem_solution();
      context.get_elem_solution_accel() -= old_elem_solution;
      context.get_elem_solution_accel() *= context.elem_solution_accel_derivative;
      context.get_elem_solution_accel().add(-1.0/(_beta*dt), old_elem_solution_rate);
      context.get_elem_solution_accel().add(-(1.0/(2.0*_beta)-1.0), old_elem_solution_accel);

      // Move the mesh into place first if necessary, set t = t_{n+1}
      (context.*reinit_func)(1.);
    }

  // If a fixed solution is requested, we'll use x_{n+1}
  if (_system.use_fixed_solution)
    context.get_elem_fixed_solution() = context.get_elem_solution();

  // Get the time derivative at t_{n+1}, F(u_{n+1})
  bool jacobian_computed = (_system.get_physics()->*time_deriv)(request_jacobian, context);

  // Damping at t_{n+1}, C(u_{n+1})
  jacobian_computed = (_system.get_physics()->*damping)(jacobian_computed, context) &&
    jacobian_computed;

  // Mass at t_{n+1}, M(u_{n+1})
  jacobian_computed = (_system.get_physics()->*mass)(jacobian_computed, context) &&
    jacobian_computed;

  // Add the constraint term
  jacobian_computed = (_system.get_physics()->*constraint)(jacobian_computed, context) &&
    jacobian_computed;

  // Add back (or restore) the old jacobian
  if (request_jacobian)
    {
      if (jacobian_computed)
        context.get_elem_jacobian() += old_elem_jacobian;
      else
        context.get_elem_jacobian().swap(old_elem_jacobian);
    }

  return jacobian_computed;
}

}
