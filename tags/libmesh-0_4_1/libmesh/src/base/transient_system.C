// $Id: transient_system.C,v 1.7 2003-11-05 22:26:44 benkirk Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002-2003  Benjamin S. Kirk, John W. Peterson
  
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
#include "libmesh_logging.h"
#include "utility.h"



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

  this->clear();
}



void TransientSystem::clear ()
{
  // clear the parent data
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
    
  // Project the old & older vectors to the new mesh
  this->project_vector (*old_local_solution);
  this->project_vector (*older_local_solution);
}



void TransientSystem::update ()
{
  // Update the parent system
  SteadySystem::update ();
}


void TransientSystem::re_update ()
{
  // re_update the parent system
  SteadySystem::re_update ();
  
  //const std::vector<unsigned int>& send_list = _dof_map.get_send_list ();

  // Explicitly build a send_list
  std::vector<unsigned int> send_list(solution->size());
  Utility::iota (send_list.begin(), send_list.end(), 0);
  
  const unsigned int first_local_dof = _dof_map.first_dof();
  const unsigned int last_local_dof  = _dof_map.last_dof();

  // Check sizes
  assert (last_local_dof > first_local_dof);
  assert (send_list.size() >= (last_local_dof - first_local_dof + 1));
  assert (older_local_solution->size() >= send_list.size());
  assert (old_local_solution->size()   >= send_list.size());
  
  // Update the old & older solutions with the send_list,
  // which may have changed since their last update.
  older_local_solution->localize (first_local_dof,
				  last_local_dof,
				  send_list);
  
  old_local_solution->localize (first_local_dof,
				last_local_dof,
				send_list);  
}
