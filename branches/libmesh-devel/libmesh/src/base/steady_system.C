// $Id: steady_system.C,v 1.1.2.3 2003-05-06 21:53:35 benkirk Exp $

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

// Local includes
#include "steady_system.h"
#include "equation_systems.h"
#include "mesh.h"
#include "linear_solver_interface.h"
#include "mesh_logging.h"




// ------------------------------------------------------------
// SteadySystem implementation
SteadySystem::SteadySystem (EquationSystems& es,
			    const std::string& name,
			    const unsigned int number) :
  
  SystemBase             (es, name, number),
  current_local_solution (NumericVector<Number>::build())
{
}



SteadySystem::~SteadySystem ()
{
}



void SteadySystem::clear ()
{
  SystemBase::clear();

  //init_system_fptr = assemble_fptr = NULL;
  
  // clear the current local solution  
  current_local_solution->clear();
}




void SteadySystem::init ()
{
  assert (_mesh.is_prepared());
  
  // initialize parent data
  SystemBase::init();

  // Possibly call a user-supplied initialization
  // method.
  if (init_system != NULL)
    this->init_system (_equation_systems, this->name());
}



void SteadySystem::reinit ()
{
  assert (_mesh.is_prepared());
   
  // initialize parent data
  SystemBase::reinit();
   
  // Project the solution to the new mesh
  {
    AutoPtr<NumericVector<Number> >
      old_solution (current_local_solution->clone());
    
    // Form the L2 projectiom
    this->project_vector (old_solution.get(),
			  solution.get());
    
    
    // Initialize the current local solution to
    // be the right size
    current_local_solution->init (this->n_dofs ());
    
    
    // Clear the linear solver interface
    linear_solver_interface->clear();
    
    // Get the user-specifiied linear solver tolerance
    const Real tol =
      _equation_systems.parameter("linear solver tolerance");
    
    // Get the user-specified maximum # of linear solver iterations
    const unsigned int maxits =
      static_cast<unsigned int>(_equation_systems.parameter("linear solver maximum iterations"));
    
    // Solve the linear system
    linear_solver_interface->solve (*matrix, *solution, *rhs, tol, maxits);
    
    // Update the local solution to reflect the new values
    this->update ();
  }
}



void SteadySystem::update ()
{
  const std::vector<unsigned int>& send_list = _dof_map.get_send_list ();

  // Check sizes
  assert (current_local_solution->local_size() == solution->size());
  assert (!send_list.empty());
  assert (send_list.size() <= solution->size());

  // Create current_local_solution from solution.  This will
  // put a local copy of solution into current_local_solution.
  // Only the necessary values (specified by the send_list)
  // are copied to minimize communication
  solution->localize (*current_local_solution, send_list); 
}



void SteadySystem::assemble ()
{
  assert (assemble_system != NULL);

  // prepare matrix with the help of the _dof_map, 
  // fill with sparsity pattern
  SystemBase::assemble();

  // Log how long the user's matrix assembly code takes
  START_LOG("assemble()", "SteadySystem");
  
  // Call the user-specified matrix assembly function
  if (this->assemble_system != NULL)
    this->assemble_system (_equation_systems, this->name());

  // Stop logging the user code
  STOP_LOG("assemble()", "SteadySystem");
}



std::pair<unsigned int, Real>
SteadySystem::solve ()
{
  // Assemble the linear system
  this->assemble (); 

  // Log how long the linear solve takes.
  START_LOG("solve()", "SteadySystem");
  
  // Get the user-specifiied linear solver tolerance
  const Real tol            =
    _equation_systems.parameter("linear solver tolerance");

  // Get the user-specified maximum # of linear solver iterations
  const unsigned int maxits =
    static_cast<unsigned int>(_equation_systems.parameter("linear solver maximum iterations"));

  // Solve the linear system
  const std::pair<unsigned int, Real> rval = 
    linear_solver_interface->solve (*matrix, *solution, *rhs, tol, maxits);

  // Stop logging the linear solve
  STOP_LOG("solve()", "SteadySystem");

  return rval; 
}
