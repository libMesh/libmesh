// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

namespace libMesh
{



DifferentiableSystem::DifferentiableSystem
                      (EquationSystems& es,
                       const std::string& name,
                       const unsigned int number) :
  Parent      (es, name, number),
  time_solver (NULL),
  deltat(1.),
  print_solution_norms(false),
  print_solutions(false),
  print_residual_norms(false),
  print_residuals(false),
  print_jacobian_norms(false),
  print_jacobians(false),
  print_element_jacobians(false),
  _diff_physics(this),
  diff_qoi(this)
{
}



DifferentiableSystem::~DifferentiableSystem ()
{
  this->clear();
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
  time_solver->init();

  // Next initialize ImplicitSystem data
  Parent::init_data();
}



AutoPtr<DiffContext> DifferentiableSystem::build_context ()
{
  AutoPtr<DiffContext> ap(new DiffContext(*this));

  ap->set_deltat_pointer( &this->deltat );

  return ap;
}


void DifferentiableSystem::assemble ()
{
  this->assembly(true, true);
}



void DifferentiableSystem::solve ()
{
  libmesh_assert(time_solver.get());
  libmesh_assert_equal_to (&(time_solver->system()), this);
  time_solver->solve();
}



LinearSolver<Number>* DifferentiableSystem::get_linear_solver() const
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



void DifferentiableSystem::release_linear_solver(LinearSolver<Number>*) const
{
}

} // namespace libMesh
