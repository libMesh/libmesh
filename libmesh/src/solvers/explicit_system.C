// $Id: explicit_system.C,v 1.1 2004-01-03 15:37:44 benkirk Exp $

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
  
  System (es, name, number)
{
}



ExplicitSystem::~ExplicitSystem ()
{
}



void ExplicitSystem::solve ()
{
  // Assemble the linear system
  this->assemble (); 

  // Update the system after the solve
  this->update();  
}
