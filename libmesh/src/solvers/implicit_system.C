// $Id: implicit_system.C,v 1.4 2004-10-12 19:46:58 benkirk Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2004  Benjamin S. Kirk, John W. Peterson
  
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
#include "implicit_system.h"
#include "equation_systems.h"
#include "sparse_matrix.h"
#include "linear_solver_interface.h"
#include "numeric_vector.h"
#include "utility.h"
#include "libmesh_logging.h"

// typedefs
typedef std::map<std::string, SparseMatrix<Number>* >::iterator        matrices_iterator;
typedef std::map<std::string, SparseMatrix<Number>* >::const_iterator  matrices_const_iterator;


// ------------------------------------------------------------
// ImplicitSystem implementation
ImplicitSystem::ImplicitSystem (EquationSystems& es,
				const std::string& name,
				const unsigned int number) :
  
  ExplicitSystem          (es, name, number),
  matrix                  (NULL),
  linear_solver_interface (LinearSolverInterface<Number>::build()),
  _can_add_matrices       (true),
  _n_linear_iterations    (0),
  _final_linear_residual         (1.e20)
{
}



ImplicitSystem::~ImplicitSystem ()
{
  // Clear data
  this->clear();
}



void ImplicitSystem::clear ()
{
  // clear the parent data
  ExplicitSystem::clear();

  // clear any user-added matrices
  {
    for (matrices_iterator pos = _matrices.begin();
	 pos != _matrices.end(); ++pos)
      {
        pos->second->clear ();
	delete pos->second;
	pos->second = NULL;
      }

    _matrices.clear();
    _can_add_matrices = true;
  }

  // NULL the matrix.
  matrix = NULL;
}



void ImplicitSystem::init_data ()
{
  // initialize parent data
  ExplicitSystem::init_data();

  // Add the system matrix.
  this->add_system_matrix ();
  
  // Initialize the matrices for the system
  this->init_matrices ();
}



void ImplicitSystem::init_matrices ()
{
  assert (matrix != NULL);

  // Check for quick return in case the system matrix
  // (and by extension all the matrices) has already
  // been initialized
  if (matrix->initialized())
    return;
  
  // no chance to add other matrices
  _can_add_matrices = false;
  
  // Tell the matrices about the dof map, and vice versa
  for (matrices_iterator pos = _matrices.begin();
       pos != _matrices.end(); ++pos)
    {
      assert (!pos->second->initialized());
      _dof_map.attach_matrix (*(pos->second));
    }
  
  // Compute the sparsity pattern for the current
  // mesh and DOF distribution.  This also updates
  // additional matrices, \p DofMap now knows them
  _dof_map.compute_sparsity (this->get_mesh());
  
  // Initialize matrices
  for (matrices_iterator pos = _matrices.begin(); 
       pos != _matrices.end(); ++pos)
    pos->second->init ();
  
  // Set the additional matrices to 0.
  for (matrices_iterator pos = _matrices.begin(); 
       pos != _matrices.end(); ++pos)
    pos->second->zero ();
}



void ImplicitSystem::reinit ()
{
  // initialize parent data
  ExplicitSystem::reinit();  

  // Clear the matrices
  for (matrices_iterator pos = _matrices.begin();
       pos != _matrices.end(); ++pos)
    pos->second->clear();

  // Clear the linear solver interface
  linear_solver_interface->clear();

  // Re-initialize the matrices
  this->init_matrices ();
}



void ImplicitSystem::assemble ()
{
  assert (matrix != NULL);
  assert (matrix->initialized());
  assert (rhs    != NULL);
  assert (rhs->initialized());

  // Zero the matrix and RHS
  matrix->zero ();
  rhs->zero ();

  // Call the base class assemble function
  ExplicitSystem::assemble ();
}



void ImplicitSystem::solve ()
{
  // Assemble the linear system
  this->assemble (); 

  // Log how long the linear solve takes.
  START_LOG("solve()", "System");
  
  // Get the user-specifiied linear solver tolerance
  const Real tol            =
    _equation_systems.parameter("linear solver tolerance");

  // Get the user-specified maximum # of linear solver iterations
  const unsigned int maxits =
    static_cast<unsigned int>(_equation_systems.parameter("linear solver maximum iterations"));

  // Solve the linear system.  Two cases:
  const std::pair<unsigned int, Real> rval =
    (this->have_matrix("Preconditioner")) ?
    // 1.) User-supplied preconditioner
    linear_solver_interface->solve (*matrix, this->get_matrix("Preconditioner"), *solution, *rhs, tol, maxits) :
    // 2.) Use system matrix for the preconditioner
    linear_solver_interface->solve (*matrix, *solution, *rhs, tol, maxits);

  // Store the number of linear iterations required to
  // solve and the final residual.
  _n_linear_iterations   = rval.first;
  _final_linear_residual = rval.second;
    
  // Stop logging the linear solve
  STOP_LOG("solve()", "System");

  // Update the system after the solve
  this->update();  
}



SparseMatrix<Number> & ImplicitSystem::add_matrix (const std::string& mat_name)
{
  // only add matrices before initializing...
  if (!_can_add_matrices)
    {
      std::cerr << "ERROR: Too late.  Cannot add matrices to the system after initialization"
		<< std::endl
		<< " any more.  You should have done this earlier."
		<< std::endl;
      error();
    }

  // Return the matrix if it is already there.
  if (this->have_matrix(mat_name))
    return *(_matrices[mat_name]);

  // Otherwise build the matrix and return it.
  SparseMatrix<Number>* buf = SparseMatrix<Number>::build().release();
  _matrices.insert (std::make_pair (mat_name, buf));

  return *buf;
}



const SparseMatrix<Number> & ImplicitSystem::get_matrix (const std::string& mat_name) const
{
  // Make sure the matrix exists
  matrices_const_iterator pos = _matrices.find (mat_name);
  
  if (pos == _matrices.end())
    {
      std::cerr << "ERROR: matrix "
		<< mat_name
		<< " does not exist in this system!"
		<< std::endl;      
      error();
    }
  
  return *(pos->second);
}



SparseMatrix<Number> & ImplicitSystem::get_matrix (const std::string& mat_name)
{
  // Make sure the matrix exists
  matrices_iterator pos = _matrices.find (mat_name);
  
  if (pos == _matrices.end())
    {
      std::cerr << "ERROR: matrix "
		<< mat_name
		<< " does not exist in this system!"
		<< std::endl;      
      error();
    }
  
  return *(pos->second);
}



void ImplicitSystem::add_system_matrix ()
{
  // Possible that we cleared the _matrices but
  // forgot to NULL-out the matrix?
  if (_matrices.empty()) matrix = NULL;


  // Only need to add the matrix if it isn't there
  // already!
  if (matrix == NULL)
    matrix = &(this->add_matrix ("System Matrix"));

  assert (matrix != NULL);
}
