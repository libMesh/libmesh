// $Id: eigen_system.C,v 1.1 2005-05-02 13:12:29 spetersen Exp $

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


#include "libmesh_config.h"

// Currently, the EigenSystem should only be available
// if SLEPc support is enabled. 
#if defined(HAVE_SLEPC)

// C++ includes

// Local includes
#include "eigen_system.h"
#include "equation_systems.h"
#include "sparse_matrix.h"
#include "eigen_solver.h"


// ------------------------------------------------------------
// EigenSystem implementation
EigenSystem::EigenSystem (EquationSystems& es,
			  const std::string& name,
			  const unsigned int number) :
  Parent           (es, name, number),
  matrix           (NULL),
  eigen_solver     (EigenSolver<Number>::build()),
  _n_converged_eigenpairs (0),
  _n_iterations           (0)
{
}



EigenSystem::~EigenSystem ()
{
  // clear data
  this->clear();
}



void EigenSystem::clear ()
{
  // Clear the parent data
  Parent::clear();

  // delete the matrix
  delete matrix;

  // NULL-out the matrix.
  matrix = NULL;

  // clear the solver
  eigen_solver->clear();

}



void EigenSystem::init_data ()
{
 
  // initialize parent data
  Parent::init_data();

  // build the system matrix
  matrix = SparseMatrix<Number>::build().release();

  // add matrix to the _dof_map
  // and compute the sparsity
  DofMap& dof_map = this->get_dof_map();

  dof_map.attach_matrix(*matrix);
  dof_map.compute_sparsity(this->get_mesh());

  // initialize and zero system matrix
  matrix->init();
  matrix->zero();

}



void EigenSystem::reinit ()
{
  // initialize parent data
  Parent::reinit();  
  
}



void EigenSystem::solve ()
{

  // A reference to the EquationSystems
  EquationSystems& es = this->get_equation_systems();

  // check that necessary parameters have been set
  assert (es.parameters.have_parameter<unsigned int>("eigenpairs"));
  assert (es.parameters.have_parameter<unsigned int>("basis vectors"));

  // Assemble the linear system
  this->assemble ();

  // Get the tolerance for the solver and the maximum
  // number of iterations. Here, we simply adopt the linear solver
  // specific parameters.
  const Real tol            =
    es.parameters.get<Real>("linear solver tolerance");

  const unsigned int maxits =
    es.parameters.get<unsigned int>("linear solver maximum iterations");

  const unsigned int nev    =
    es.parameters.get<unsigned int>("eigenpairs");

  const unsigned int ncv    =
    es.parameters.get<unsigned int>("basis vectors");

  // call the solver
  const std::pair<unsigned int, unsigned int> solve_data =
  eigen_solver->solve(*matrix, nev, ncv, tol, maxits);

  this->_n_converged_eigenpairs = solve_data.first;
  this->_n_iterations           = solve_data.second;

}



void EigenSystem::assemble ()
{

  // Assemble the linear system
  Parent::assemble (); 

}


std::pair<Real, Real> EigenSystem::get_eigenpair (unsigned int i)
{
  // call the eigen_solver get_eigenpair method
  return eigen_solver->get_eigenpair (i, *solution);
}

#endif // HAVE_SLEPC
