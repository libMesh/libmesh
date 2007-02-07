// $Id: transient_system.C,v 1.10 2006-03-29 20:56:45 roystgnr Exp $

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
#include "explicit_system.h"
#include "linear_implicit_system.h"
#include "nonlinear_implicit_system.h"
#include "utility.h"
#include "dof_map.h"
#include "numeric_vector.h"


// ------------------------------------------------------------
// TransientSystem implementation
template <class Base>
TransientSystem<Base>::TransientSystem (EquationSystems& es,
					const std::string& name,
					const unsigned int number) :
  
  Base                 (es, name, number)
{
  old_local_solution =
    AutoPtr<NumericVector<Number> >
      (&(this->add_vector("_transient_old_local_solution")));
  older_local_solution =
    AutoPtr<NumericVector<Number> >
      (&(this->add_vector("_transient_older_local_solution")));
}



template <class Base>
TransientSystem<Base>::~TransientSystem ()
{
  this->clear();

  // We still have AutoPtrs for API compatibility, but
  // now that we're System::add_vector()ing these, we can trust
  // the base class to handle memory management
  old_local_solution.release();  
  older_local_solution.release();  
}



template <class Base>
void TransientSystem<Base>::clear ()
{
  // clear the parent data
  Base::clear();

  // the old & older local solutions
  // are now deleted by System!
  // old_local_solution->clear();
  // older_local_solution->clear();  

  // FIXME: This preserves maximum backwards compatibility,
  // but is probably grossly unnecessary:
  old_local_solution.release();  
  older_local_solution.release();  

  old_local_solution =
    AutoPtr<NumericVector<Number> >
      (&(this->add_vector("_transient_old_local_solution")));
  older_local_solution =
    AutoPtr<NumericVector<Number> >
      (&(this->add_vector("_transient_older_local_solution")));
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
  // The System::reinit handles this now
  // this->project_vector (*old_local_solution);
  // this->project_vector (*older_local_solution);
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
  
  const unsigned int first_local_dof = Base::get_dof_map().first_dof();
  const unsigned int last_local_dof  = Base::get_dof_map().last_dof();

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




template <class Base>
Number TransientSystem<Base>::old_solution (const unsigned int global_dof_number) const
{
  // Check the sizes
  assert (global_dof_number < this->get_dof_map().n_dofs());
  assert (global_dof_number < old_local_solution->size());
   
  return (*old_local_solution)(global_dof_number);
}



template <class Base>
Number TransientSystem<Base>::older_solution (const unsigned int global_dof_number) const
{
  // Check the sizes
  assert (global_dof_number < this->get_dof_map().n_dofs());
  assert (global_dof_number < older_local_solution->size());
   
  return (*older_local_solution)(global_dof_number);
}




// ------------------------------------------------------------
// TransientSystem instantiations
template class TransientSystem<LinearImplicitSystem>;
template class TransientSystem<NonlinearImplicitSystem>;
template class TransientSystem<ExplicitSystem>;
