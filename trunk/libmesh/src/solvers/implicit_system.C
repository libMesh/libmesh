// $Id: implicit_system.C,v 1.8 2005-02-22 22:17:43 jwpeterson Exp $

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
#include "implicit_system.h"
#include "equation_systems.h"
#include "sparse_matrix.h"
#include "utility.h"
#include "libmesh_logging.h"



// ------------------------------------------------------------
// ImplicitSystem implementation
ImplicitSystem::ImplicitSystem (EquationSystems& es,
				const std::string& name,
				const unsigned int number) :
  
  Parent            (es, name, number),
  matrix            (NULL),
  _can_add_matrices (true)
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
  Parent::clear();

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
  Parent::init_data();

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
  Parent::reinit();  

  // Clear the matrices
  for (matrices_iterator pos = _matrices.begin();
       pos != _matrices.end(); ++pos)
    pos->second->clear();

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
  Parent::assemble ();
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
  const_matrices_iterator pos = _matrices.find (mat_name);
  
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
