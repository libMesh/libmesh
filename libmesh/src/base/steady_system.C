// $Id: steady_system.C,v 1.3 2003-05-16 14:37:49 benkirk Exp $

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
#include "sparse_matrix.h"
#include "linear_solver_interface.h"
#include "numeric_vector.h"



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
  
  // clear the current local solution  
  current_local_solution->clear();
}




void SteadySystem::init_data ()
{
  // initialize parent data
  SystemBase::init_data();

  // Initialize the local copy of the solution vector
  current_local_solution->init (this->n_dofs());
}



void SteadySystem::reinit ()
{
  // initialize parent data
  SystemBase::reinit();
   
  // Project the solution to the new mesh
  {
    // Project the current local solution
    this->project_vector (*current_local_solution);

    // Update the solution based on the projected
    // current_local_solution.
    {
      solution->init (this->n_dofs(), this->n_local_dofs());

      const unsigned int first_local_dof = solution->first_local_index();
      const unsigned int local_size      = solution->local_size();
      
      for (unsigned int i=0; i<local_size; i++)
	solution->set(i+first_local_dof,
		      (*current_local_solution)(i+first_local_dof));
    }
    
    // Update the local solution to reflect the new values
    //this->update ();
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
