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
#include "linear_implicit_system.h"
#include "linear_solver.h"
#include "equation_systems.h"
#include "libmesh_logging.h"
#include "sparse_matrix.h" // for get_transpose


// ------------------------------------------------------------
// LinearImplicitSystem implementation
LinearImplicitSystem::LinearImplicitSystem (EquationSystems& es,
					    const std::string& name,
					    const unsigned int number) :
  
  Parent                 (es, name, number),
  linear_solver          (LinearSolver<Number>::build()),
  _n_linear_iterations   (0),
  _final_linear_residual (1.e20),
  _shell_matrix(NULL)
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
  if (this->assemble_before_solve)
    // Assemble the linear system
    this->assemble (); 

  // Log how long the linear solve takes.
  // This gets done by the LinearSolver classes now [RHS]
  // START_LOG("solve()", "System");

  // Get a reference to the EquationSystems
  const EquationSystems& es =
    this->get_equation_systems();
  
  // Get the user-specifiied linear solver tolerance
  const Real tol            =
    es.parameters.get<Real>("linear solver tolerance");

  // Get the user-specified maximum # of linear solver iterations
  const unsigned int maxits =
    es.parameters.get<unsigned int>("linear solver maximum iterations");

  // Solve the linear system.  Several cases:
  std::pair<unsigned int, Real> rval = std::make_pair(0,0.0);
  if(_shell_matrix!=NULL)
    {
      if(this->have_matrix("Preconditioner"))
	{
	  // 1.) Shell matrix plus user-supplied preconditioner.
	  rval = linear_solver->solve(*_shell_matrix, this->get_matrix("Preconditioner"), *solution, *rhs, tol, maxits);
	}
      else
	{
	  // 2.) Shell matrix without user-supplied preconditioner.
	  rval = linear_solver->solve(*_shell_matrix, *solution, *rhs, tol, maxits);
	}
    }
  else
    {
      if(this->have_matrix("Preconditioner"))
	{
	  // 3.) No shell matrix, but with user-supplied preconditioner
	  rval = linear_solver->solve (*matrix, this->get_matrix("Preconditioner"), *solution, *rhs, tol, maxits);
	}
      else
	{
	  // 4.) No shell matrix, and use system matrix for the preconditioner
	  rval = linear_solver->solve (*matrix, *solution, *rhs, tol, maxits);
	}
    }

  // Store the number of linear iterations required to
  // solve and the final residual.
  _n_linear_iterations   = rval.first;
  _final_linear_residual = rval.second;
    
  // Stop logging the linear solve
  // This gets done by the LinearSolver classes now [RHS]
  // STOP_LOG("solve()", "System");

  // Update the system after the solve
  this->update();  
}



void LinearImplicitSystem::adjoint_solve ()
{
  // We currently don't support adjoint solves of shell matrices
  // FIXME - we should let shell matrices support
  // vector_transpose_mult so that we can use them here.
  if(_shell_matrix!=NULL)
    libmesh_not_implemented();

  // Adding an adjoint_solution vector, allocate an adjoint_solution if it doesn't already exist 

  NumericVector<Number> & adjoint_solution = System::add_vector("adjoint_solution");     

  if (this->assemble_before_solve)
    {
      // Assemble the linear system
      this->assemble (); 

      // And take the adjoint
      matrix->get_transpose(*matrix);

      // Including of any separate preconditioner
      if(this->have_matrix("Preconditioner"))
        {
          SparseMatrix<Number> &pre = this->get_matrix("Preconditioner");
          pre.get_transpose(pre);
        }

      // But now replace the right hand side with the quantity of
      // interest functional
      this->assemble_qoi();
    }

  // Get a reference to the EquationSystems
  const EquationSystems& es =
    this->get_equation_systems();
  
  // Get the user-specifiied linear solver tolerance
  const Real tol            =
    es.parameters.get<Real>("adjoint solver tolerance");

  // Get the user-specified maximum # of linear solver iterations
  const unsigned int maxits =
    es.parameters.get<unsigned int>("adjoint solver maximum iterations");

  // Solve the linear system.  Several cases:
  std::pair<unsigned int, Real> rval = std::make_pair(0,0.0);
  if(this->have_matrix("Preconditioner"))
    {
      // 3.) No shell matrix, but with user-supplied preconditioner
      rval = linear_solver->solve (*matrix, this->get_matrix("Preconditioner"), adjoint_solution, *rhs, tol, maxits);
    }
  else
    {
      // 4.) No shell matrix, and use system matrix for the preconditioner
      rval = linear_solver->solve (*matrix, adjoint_solution, *rhs, tol, maxits);
    }

  // Store the number of linear iterations required to
  // solve and the final residual.
  _n_linear_iterations   = rval.first;
  _final_linear_residual = rval.second;
    
  // Stop logging the linear solve
  // This gets done by the LinearSolver classes now [RHS]
  // STOP_LOG("solve()", "System");

  // Update the system after the solve
  this->update();  
}



void LinearImplicitSystem::attach_shell_matrix (ShellMatrix<Number>* shell_matrix)
{
  _shell_matrix = shell_matrix;
}



