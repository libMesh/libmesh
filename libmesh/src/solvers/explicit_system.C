// $Id: explicit_system.C,v 1.4 2005-01-06 21:55:04 benkirk Exp $

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
#include "explicit_system.h"


// ------------------------------------------------------------
// ExplicitSystem implementation
ExplicitSystem::ExplicitSystem (EquationSystems& es,
				const std::string& name,
				const unsigned int number) :
  Parent (es, name, number),
  rhs(NULL)

{
  //rhs = &(this->add_vector ("RHS Vector"));
}



ExplicitSystem::~ExplicitSystem ()
{
  // clear data
  this->clear();
}



void ExplicitSystem::clear ()
{
  // Clear the parent data
  Parent::clear();

  // NULL-out the matrix.  Note that
  // System::clear() actually deleted it.
  rhs = NULL;
}



void ExplicitSystem::init_data ()
{
  // Add the system RHS.
  // (We must do this before initializing the System data,
  //  then we lose the opportunity to add vectors).
  this->add_system_rhs ();
  
  // initialize parent data
  Parent::init_data();
}



void ExplicitSystem::reinit ()
{
  // initialize parent data
  Parent::reinit();  
  
  // Resize the RHS conformal to the current mesh
  rhs->init (this->n_dofs(), this->n_local_dofs());
}



void ExplicitSystem::solve ()
{
  // Assemble the linear system
  this->assemble (); 

  // Update the system after the solve
  this->update();  
}



void ExplicitSystem::add_system_rhs ()
{
  // Possible that we cleared the _vectors but
  // forgot to NULL-out the rhs?
  if (_vectors.empty()) rhs = NULL;


  // Only need to add the rhs if it isn't there
  // already!
  if (rhs == NULL)
    rhs = &(this->add_vector ("RHS Vector"));

  assert (rhs != NULL);
}
