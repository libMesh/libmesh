// $Id: transient_system.C,v 1.1.2.2 2003-05-07 20:47:15 benkirk Exp $

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
#include "transient_system.h"
#include "equation_systems.h"
#include "linear_solver_interface.h"
#include "mesh_logging.h"




// ------------------------------------------------------------
// TransientSystem implementation
TransientSystem::TransientSystem (EquationSystems& es,
			    const std::string& name,
			    const unsigned int number) :
  
  SteadySystem         (es, name, number),
  old_local_solution   (NumericVector<Number>::build()),
  older_local_solution (NumericVector<Number>::build())
{
}



TransientSystem::~TransientSystem ()
{
  // Null-out the function pointers.  Since this
  // class is getting destructed it is pointless,
  // but a good habit.
  init_system = assemble_system = NULL;
}



void TransientSystem::clear ()
{
  SteadySystem::clear();

  // clear the old & older local solutions
  old_local_solution->clear();  
  older_local_solution->clear();  
}




void TransientSystem::init_data ()
{
  // initialize parent data
  SteadySystem::init_data();

  // Initialize the old & older solutions
  old_local_solution->init   (this->n_dofs());
  older_local_solution->init (this->n_dofs());
}



void TransientSystem::reinit ()
{
  // initialize parent data
  SteadySystem::reinit();
   
   
  // Clear the linear solver interface
  linear_solver_interface->clear();
  
  // Project the older solution to the new mesh
  {
    // Form the L2 projectiom
    this->project_vector (older_local_solution.get(),
			  solution.get());
      
    // Initialize the older local solution to
    // be the right size
    older_local_solution->init (this->n_dofs ());
    
    // Get the user-specifiied linear solver tolerance
    const Real tol =
      _equation_systems.parameter("linear solver tolerance");
    
    // Get the user-specified maximum # of linear solver iterations
    const unsigned int maxits =
      static_cast<unsigned int>(_equation_systems.parameter("linear solver maximum iterations"));
    
    // Solve the linear system
    linear_solver_interface->solve (*matrix, *solution, *rhs, tol, maxits);
    
    // Update the local solution to reflect the new values
    {
      const std::vector<unsigned int>& send_list = _dof_map.get_send_list ();
      
      // Check sizes
      assert (older_local_solution->local_size() == solution->size());
      assert (!send_list.empty());
      assert (send_list.size() <= solution->size());
      
      // Create current_local_solution from solution.  This will
      // put a local copy of solution into current_local_solution.
      // Only the necessary values (specified by the send_list)
      // are copied to minimize communication
      solution->localize (*older_local_solution, send_list); 
    }    
  }


  
  // Project the old solution to the new mesh
  {
    // Form the L2 projectiom
    this->project_vector (old_local_solution.get(),
			  solution.get());
      
    // Initialize the old local solution to
    // be the right size
    old_local_solution->init (this->n_dofs ());
    
    // Get the user-specifiied linear solver tolerance
    const Real tol =
      _equation_systems.parameter("linear solver tolerance");
    
    // Get the user-specified maximum # of linear solver iterations
    const unsigned int maxits =
      static_cast<unsigned int>(_equation_systems.parameter("linear solver maximum iterations"));
    
    // Solve the linear system
    linear_solver_interface->solve (*matrix, *solution, *rhs, tol, maxits);
    
    // Update the local solution to reflect the new values
    {
      const std::vector<unsigned int>& send_list = _dof_map.get_send_list ();
      
      // Check sizes
      assert (old_local_solution->local_size() == solution->size());
      assert (!send_list.empty());
      assert (send_list.size() <= solution->size());
      
      // Create current_local_solution from solution.  This will
      // put a local copy of solution into current_local_solution.
      // Only the necessary values (specified by the send_list)
      // are copied to minimize communication
      solution->localize (*old_local_solution, send_list); 
    }    
  }
}
