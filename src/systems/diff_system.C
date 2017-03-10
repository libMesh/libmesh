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


#include "libmesh/diff_solver.h"
#include "libmesh/diff_system.h"
#include "libmesh/time_solver.h"
#include "libmesh/unsteady_solver.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/dof_map.h"
#include "libmesh/zero_function.h"

namespace libMesh
{



DifferentiableSystem::DifferentiableSystem(EquationSystems & es,
                                           const std::string & name_in,
                                           const unsigned int number_in) :
  Parent      (es, name_in, number_in),
  time_solver (),
  deltat(1.),
  print_solution_norms(false),
  print_solutions(false),
  print_residual_norms(false),
  print_residuals(false),
  print_jacobian_norms(false),
  print_jacobians(false),
  print_element_solutions(false),
  print_element_residuals(false),
  print_element_jacobians(false),
  _diff_physics(this),
  diff_qoi(this)
{
}



DifferentiableSystem::~DifferentiableSystem ()
{
  // If we had an attached Physics object, delete it.
  if (this->_diff_physics != this)
    delete this->_diff_physics;

  // If we had an attached QoI object, delete it.
  if (this->diff_qoi != this)
    delete this->diff_qoi;
}



void DifferentiableSystem::clear ()
{
  // If we had an attached Physics object, delete it.
  if (this->_diff_physics != this)
    {
      delete this->_diff_physics;
      this->_diff_physics = this;
    }
  // If we had no attached Physics object, clear our own Physics data
  else
    this->clear_physics();

  // If we had an attached QoI object, delete it.
  if (this->diff_qoi != this)
    {
      delete this->diff_qoi;
      this->diff_qoi = this;
    }
  // If we had no attached QoI object, clear our own QoI data
  else
    this->clear_qoi();

  use_fixed_solution = false;
}



void DifferentiableSystem::reinit ()
{
  Parent::reinit();

  libmesh_assert(time_solver.get());
  libmesh_assert_equal_to (&(time_solver->system()), this);

  time_solver->reinit();
}



void DifferentiableSystem::init_data ()
{
  // If it isn't a separate initialized-upon-attachment object, do any
  // initialization our physics needs.
  if (this->_diff_physics == this)
    this->init_physics(*this);

  // Do any initialization our solvers need
  libmesh_assert(time_solver.get());
  libmesh_assert_equal_to (&(time_solver->system()), this);

  // Now check for second order variables and add their velocities to the System.
  if( !time_solver->is_steady() )
    {
      const UnsteadySolver & unsteady_solver =
        cast_ref<const UnsteadySolver &>(*(time_solver.get()));

      if( unsteady_solver.time_order() == 1 )
        this->add_second_order_dot_vars();
    }

  time_solver->init();

  // Next initialize ImplicitSystem data
  Parent::init_data();

  time_solver->init_data();
}

UniquePtr<DiffContext> DifferentiableSystem::build_context ()
{
  DiffContext * context = new DiffContext(*this);
  context->set_deltat_pointer( &this->deltat );
  return UniquePtr<DiffContext>(context);
}


void DifferentiableSystem::assemble ()
{
  this->assembly(true, true);
}



void DifferentiableSystem::solve ()
{
  // Get the time solver object associated with the system, and tell it that
  // we are not solving the adjoint problem
  this->get_time_solver().set_is_adjoint(false);

  libmesh_assert_equal_to (&(time_solver->system()), this);
  time_solver->solve();
}



std::pair<unsigned int, Real> DifferentiableSystem::adjoint_solve (const QoISet & qoi_indices)
{
  // Get the time solver object associated with the system, and tell it that
  // we are solving the adjoint problem
  this->get_time_solver().set_is_adjoint(true);

  return this->ImplicitSystem::adjoint_solve(qoi_indices);
}



LinearSolver<Number> * DifferentiableSystem::get_linear_solver() const
{
  libmesh_assert(time_solver.get());
  libmesh_assert_equal_to (&(time_solver->system()), this);
  return this->time_solver->linear_solver().get();
}



std::pair<unsigned int, Real> DifferentiableSystem::get_linear_solve_parameters() const
{
  libmesh_assert(time_solver.get());
  libmesh_assert_equal_to (&(time_solver->system()), this);
  return std::make_pair(this->time_solver->diff_solver()->max_linear_iterations,
                        this->time_solver->diff_solver()->relative_residual_tolerance);
}



void DifferentiableSystem::release_linear_solver(LinearSolver<Number> *) const
{
}

void DifferentiableSystem::add_second_order_dot_vars()
{
  const std::set<unsigned int> & second_order_vars = this->get_second_order_vars();
  if( !second_order_vars.empty() )
    {
      for( std::set<unsigned int>::const_iterator var_it = second_order_vars.begin();
           var_it != second_order_vars.end(); ++var_it )
        {
          const Variable & var = this->variable(*var_it);
          std::string new_var_name = std::string("dot_")+var.name();

          unsigned int v_var_idx;

          if( var.active_subdomains().empty() )
            v_var_idx = this->add_variable( new_var_name, var.type() );
          else
            v_var_idx = this->add_variable( new_var_name, var.type(), &var.active_subdomains() );

          _second_order_dot_vars.insert( std::pair<unsigned int,unsigned int>(*var_it,v_var_idx) );

          // The new velocities are time evolving variables of first order
          this->time_evolving( v_var_idx, 1 );

          // And if there are any boundary conditions set on the second order
          // variable, we also need to set it on its velocity variable.
          this->add_dot_var_dirichlet_bcs( *var_it, v_var_idx );
        }
    }
}

void DifferentiableSystem::add_dot_var_dirichlet_bcs( unsigned int var_idx,
                                                      unsigned int dot_var_idx )
{
  // We're assuming that there could be a lot more variables than
  // boundary conditions, so we search each of the boundary conditions
  // for this variable rather than looping over boundary conditions
  // in a separate loop and searching through all the variables.
  const DirichletBoundaries * all_dbcs =
    this->get_dof_map().get_dirichlet_boundaries();

  if( all_dbcs )
    {
      // We need to cache the DBCs to be added so that we add them
      // after looping over the existing DBCs. Otherwise, we're polluting
      // the thing we're looping over.
      std::vector<DirichletBoundary*> new_dbcs;

      DirichletBoundaries::const_iterator dbc_it = all_dbcs->begin();
      for( ; dbc_it != all_dbcs->end(); ++dbc_it )
        {
          libmesh_assert(*dbc_it);
          DirichletBoundary & dbc = *(*dbc_it);

          // Look for second order variable in the current
          // DirichletBoundary object
          std::vector<unsigned int>::const_iterator dbc_var_it =
            std::find( dbc.variables.begin(), dbc.variables.end(), var_idx );

          // If we found it, then we also need to add it's corresponding
          // "dot" variable to a DirichletBoundary
          std::vector<unsigned int> vars_to_add;
          if( dbc_var_it != dbc.variables.end() )
            vars_to_add.push_back(dot_var_idx);

          if( !vars_to_add.empty() )
            {
              // We need to check if the boundary condition is time-dependent.
              // Currently, we cannot automatically differentiate w.r.t. time
              // so if the user supplies a time-dependent Dirichlet BC, then
              // we can't automatically support the Dirichlet BC for the
              // "velocity" boundary condition, so we error. Otherwise,
              // the "velocity boundary condition will just be zero.
              bool is_time_evolving_bc = false;
              if( dbc.f )
                is_time_evolving_bc = dbc.f->is_time_dependent();
              else if( dbc.f_fem )
                // We it's a FEMFunctionBase object, it will be implictly
                // time-dependent since it is assumed to depend on the solution.
                is_time_evolving_bc = true;
              else
                libmesh_error_msg("Could not find valid boundary function!");

              if( is_time_evolving_bc )
                libmesh_error_msg("Cannot currently support time-dependent Dirichlet BC for dot variables!");


              DirichletBoundary * new_dbc;

              if( dbc.f )
                {
                  ZeroFunction<Number> zf;

                  new_dbc = new DirichletBoundary(dbc.b, vars_to_add, zf);
                }
              else
                libmesh_error();

              new_dbcs.push_back(new_dbc);
            }
        }

      // Now add the new DBCs for the "dot" vars to the DofMap
      std::vector<DirichletBoundary*>::iterator new_dbc_it =
        new_dbcs.begin();

      for( ; new_dbc_it != new_dbcs.end(); ++new_dbc_it )
        {
          const DirichletBoundary & dbc = *(*new_dbc_it);
          this->get_dof_map().add_dirichlet_boundary(dbc);
          delete *new_dbc_it;
        }

    } // if(all_dbcs)
}

unsigned int DifferentiableSystem::get_second_order_dot_var( unsigned int var ) const
{
  // For SteadySolver or SecondOrderUnsteadySolvers, we just give back var
  unsigned int dot_var = var;

  if( !time_solver->is_steady() )
    {
      const UnsteadySolver & unsteady_solver =
        cast_ref<const UnsteadySolver &>(*(time_solver.get()));

      if( unsteady_solver.time_order() == 1 )
        dot_var = this->_second_order_dot_vars.find(var)->second;
    }

  return dot_var;
}

bool DifferentiableSystem::have_first_order_scalar_vars() const
{
  bool have_first_order_scalar_vars = false;

  if(this->have_first_order_vars())
    {
      for( std::set<unsigned int>::const_iterator var_it = this->get_first_order_vars().begin();
           var_it != this->get_first_order_vars().end();
           ++var_it )
        {
          if( this->variable(*var_it).type().family == SCALAR )
            have_first_order_scalar_vars = true;
        }
    }

  return have_first_order_scalar_vars;
}

bool DifferentiableSystem::have_second_order_scalar_vars() const
{
  bool have_second_order_scalar_vars = false;

  if(this->have_second_order_vars())
    {
      for( std::set<unsigned int>::const_iterator var_it = this->get_second_order_vars().begin();
           var_it != this->get_second_order_vars().end();
           ++var_it )
        {
          if( this->variable(*var_it).type().family == SCALAR )
            have_second_order_scalar_vars = true;
        }
    }

  return have_second_order_scalar_vars;
}

} // namespace libMesh
