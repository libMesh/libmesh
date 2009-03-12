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
