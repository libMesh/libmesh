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



// C++ includes

// Local includes
#include "nonlinear_implicit_system.h"
#include "equation_systems.h"
#include "libmesh_logging.h"
#include "nonlinear_solver.h"


// ------------------------------------------------------------
// NonlinearImplicitSystem implementation
NonlinearImplicitSystem::NonlinearImplicitSystem (EquationSystems& es,
						  const std::string& name,
						  const unsigned int number) :
  
  Parent                    (es, name, number),
  nonlinear_solver          (NonlinearSolver<Number>::build(*this)),
  _n_nonlinear_iterations   (0),
  _final_nonlinear_residual (1.e20)
{}



NonlinearImplicitSystem::~NonlinearImplicitSystem ()
{
  // Clear data
  this->clear();
}



void NonlinearImplicitSystem::clear ()
{
  // clear the nonlinear solver
  nonlinear_solver->clear();
  
  // clear the parent data
  Parent::clear();
}



void NonlinearImplicitSystem::reinit ()
{
  // re-initialize the nonlinear solver interface
  nonlinear_solver->clear();
  
  // initialize parent data
  Parent::reinit();  
}



void NonlinearImplicitSystem::solve ()
{
  if (this->assemble_before_solve)
    // Assemble the nonlinear system
    this->assemble (); 

  // Log how long the nonlinear solve takes.
  START_LOG("solve()", "System");
  
  // Get a reference to the EquationSystems
  const EquationSystems& es =
    this->get_equation_systems();  
  
  // Get the user-specifiied nonlinear solver tolerance
  const Real tol            =
    es.parameters.get<Real>("nonlinear solver tolerance");

  // Get the user-specified maximum # of nonlinear solver iterations
  const unsigned int maxits =
    es.parameters.get<unsigned int>("nonlinear solver maximum iterations");

  // Solve the nonlinear system.  Two cases:
  const std::pair<unsigned int, Real> rval =
    nonlinear_solver->solve (*matrix, *solution, *rhs, tol, maxits);

  // Store the number of nonlinear iterations required to
  // solve and the final residual.
  _n_nonlinear_iterations   = rval.first;
  _final_nonlinear_residual = rval.second;
    
  // Stop logging the nonlinear solve
  STOP_LOG("solve()", "System");

  // Update the system after the solve
  this->update();  
}
