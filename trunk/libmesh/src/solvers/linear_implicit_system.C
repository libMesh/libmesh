// $Id: linear_implicit_system.C,v 1.2 2005-02-22 22:17:43 jwpeterson Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2005  Benjamin S. Kirk, John W. Peterson
  
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
#include "linear_implicit_system.h"
#include "equation_systems.h"
#include "sparse_matrix.h"
#include "linear_solver.h"
#include "numeric_vector.h"
#include "utility.h"
#include "libmesh_logging.h"


// ------------------------------------------------------------
// LinearImplicitSystem implementation
LinearImplicitSystem::LinearImplicitSystem (EquationSystems& es,
					    const std::string& name,
					    const unsigned int number) :
  
  Parent                 (es, name, number),
  linear_solver          (LinearSolver<Number>::build()),
  _n_linear_iterations   (0),
  _final_linear_residual (1.e20)
{
}



LinearImplicitSystem::~LinearImplicitSystem ()
{
  // Clear data
  this->clear();
}



void LinearImplicitSystem::clear ()
{
  // clear the linear solver
  linear_solver->clear();
  
  // clear the parent data
  Parent::clear();
}



void LinearImplicitSystem::reinit ()
{
  // re-initialize the linear solver interface
  linear_solver->clear();
  
  // initialize parent data
  Parent::reinit();  
}



void LinearImplicitSystem::solve ()
{
  // Assemble the linear system
  this->assemble (); 

  // Log how long the linear solve takes.
  START_LOG("solve()", "System");
  
  // Get the user-specifiied linear solver tolerance
  const Real tol            =
    _equation_systems.parameters.get<Real>("linear solver tolerance");

  // Get the user-specified maximum # of linear solver iterations
  const unsigned int maxits =
    _equation_systems.parameters.get<unsigned int>("linear solver maximum iterations");

  // Solve the linear system.  Two cases:
  const std::pair<unsigned int, Real> rval =
    (this->have_matrix("Preconditioner")) ?
    // 1.) User-supplied preconditioner
    linear_solver->solve (*matrix, this->get_matrix("Preconditioner"), *solution, *rhs, tol, maxits) :
    // 2.) Use system matrix for the preconditioner
    linear_solver->solve (*matrix, *solution, *rhs, tol, maxits);

  // Store the number of linear iterations required to
  // solve and the final residual.
  _n_linear_iterations   = rval.first;
  _final_linear_residual = rval.second;
    
  // Stop logging the linear solve
  STOP_LOG("solve()", "System");

  // Update the system after the solve
  this->update();  
}
