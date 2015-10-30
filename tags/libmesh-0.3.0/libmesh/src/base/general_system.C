// $Id: general_system.C,v 1.5 2003-02-20 04:59:58 benkirk Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002  Benjamin S. Kirk, John W. Peterson
  
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
#include <math.h>
#include <algorithm>

// Local includes
#include "general_system.h"
#include "libmesh.h"



// ------------------------------------------------------------
// GeneralSystem implementation
GeneralSystem::GeneralSystem (EquationSystems&    es,
			      const std::string&  name,
			      const unsigned int  number,
			      const SolverPackage solver_package) :
  SystemBase             (es.get_mesh(), name, number, solver_package),
  current_local_solution (NumericVector<Number>::build(solver_package)),
  old_local_solution     (NumericVector<Number>::build(solver_package)),
  older_local_solution   (NumericVector<Number>::build(solver_package)),
  init_system_fptr       (NULL),
  assemble_fptr          (NULL),
  equation_systems       (es)
{
}



GeneralSystem::~GeneralSystem ()
{
  //init_system_fptr = assemble_fptr = NULL;

  current_local_solution->clear ();

  old_local_solution->clear ();
  
  older_local_solution->clear ();
}



void GeneralSystem::clear ()
{
  SystemBase::clear();

  //init_system_fptr = assemble_fptr = NULL;

  current_local_solution->clear ();

  old_local_solution->clear ();
  
  older_local_solution->clear ();
}




void GeneralSystem::init ()
{
  // initialize parent data
  SystemBase::init();

  current_local_solution->init (n_dofs ());

  old_local_solution->init     (n_dofs ());
  
  older_local_solution->init   (n_dofs ());

  // Possibly call a user-supplied initialization
  // method.
  if (init_system_fptr != NULL)
    {
      init_system_fptr (equation_systems, name());

      update();

      *old_local_solution   = *current_local_solution;
      *older_local_solution = *current_local_solution; 
    }
}



void GeneralSystem::update ()
{
  assert (current_local_solution->local_size() == solution->size());

  const std::vector<unsigned int>& send_list = _dof_map.get_send_list ();

  assert (!send_list.empty());

  assert (send_list.size() <= solution->size());
  
  solution->localize (*current_local_solution, send_list);
}



void GeneralSystem::assemble ()
{
  assert (assemble_fptr != NULL);

  // prepare matrix with the help of the _dof_map, 
  // fill with sparsity pattern
  SystemBase::assemble();

  // Log how long the user's matrix assembly code takes
  libMesh::log.start_event("assemble()");
  
  // Call the user-specified matrix assembly function
  assemble_fptr (equation_systems, name());

  // Stop logging the user code
  libMesh::log.stop_event("assemble()");
}



std::pair<unsigned int, Real>
GeneralSystem::solve ()
{
  // Assemble the linear system
  assemble (); 

  // Log how long the linear solve takes.
  libMesh::log.start_event("solve()");
  
  // Get the user-specifiied linear solver tolerance
  const Real tol            =
    equation_systems.parameter("linear solver tolerance");

  // Get the user-specified maximum # of linear solver iterations
  const unsigned int maxits =
    static_cast<unsigned int>(equation_systems.parameter("linear solver maximum iterations"));

  // Solve the linear system
  const std::pair<unsigned int, Real> rval = 
    linear_solver_interface->solve (*matrix, *solution, *rhs, tol, maxits);

  // Update the local solution to reflect the new values
  update ();

  // Stop logging the linear solve
  libMesh::log.stop_event("solve()");

  return rval; 
}



void GeneralSystem::attach_init_function(void fptr(EquationSystems& es,
						   const std::string& name))
{
  assert (fptr != NULL);
  
  init_system_fptr = fptr;
}



void GeneralSystem::attach_assemble_function(void fptr(EquationSystems& es,
						       const std::string& name))
{
  assert (fptr != NULL);
  
  assemble_fptr = fptr;  
}
