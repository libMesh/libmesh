// $Id: transient_system.C,v 1.5 2005-02-22 22:17:43 jwpeterson Exp $

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
#include "transient_system.h"
#include "implicit_system.h"
#include "explicit_system.h"
#include "libmesh_logging.h"
#include "utility.h"



// ------------------------------------------------------------
// TransientSystem implementation
template <class Base>
TransientSystem<Base>::TransientSystem (EquationSystems& es,
					const std::string& name,
					const unsigned int number) :
  
  Base                 (es, name, number),
  old_local_solution   (NumericVector<Number>::build()),
  older_local_solution (NumericVector<Number>::build())
{
}



template <class Base>
TransientSystem<Base>::~TransientSystem ()
{
  this->clear();
}



template <class Base>
void TransientSystem<Base>::clear ()
{
  // clear the parent data
  Base::clear();

  // clear the old & older local solutions
  old_local_solution->clear();
  older_local_solution->clear();  
}




template <class Base>
void TransientSystem<Base>::init_data ()
{
  // initialize parent data
  Base::init_data();

  // Initialize the old & older solutions
  old_local_solution->init   (this->n_dofs());
  older_local_solution->init (this->n_dofs());
}



template <class Base>
void TransientSystem<Base>::reinit ()
{
  // initialize parent data
  Base::reinit();
    
  // Project the old & older vectors to the new mesh
  this->project_vector (*old_local_solution);
  this->project_vector (*older_local_solution);
}



template <class Base>
void TransientSystem<Base>::re_update ()
{
  // re_update the parent system
  Base::re_update ();
  
  //const std::vector<unsigned int>& send_list = Base::_dof_map.get_send_list ();

  // Explicitly build a send_list
  std::vector<unsigned int> send_list(Base::solution->size());
  Utility::iota (send_list.begin(), send_list.end(), 0);
  
  const unsigned int first_local_dof = Base::_dof_map.first_dof();
  const unsigned int last_local_dof  = Base::_dof_map.last_dof();

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



// ------------------------------------------------------------
// TransientSystem instantiations
template class TransientSystem<LinearImplicitSystem>;
template class TransientSystem<NonlinearImplicitSystem>;
template class TransientSystem<ExplicitSystem>;
