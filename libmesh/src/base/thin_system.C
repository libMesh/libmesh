// $Id: thin_system.C,v 1.2 2003-04-30 21:09:29 benkirk Exp $

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
#include "thin_system.h"



// ------------------------------------------------------------
// ThinSystem implementation
ThinSystem::ThinSystem (EquationSystems<ThinSystem>& es,
			const std::string&           name,
			const unsigned int           number,
			const SolverPackage          solver_package) :
  SystemBase       (es.get_mesh(), name, number, solver_package),
  init_system_fptr (NULL),
  assemble_fptr    (NULL),
  _equation_systems (es)
{
}



ThinSystem::~ThinSystem ()
{
  //init_system_fptr = assemble_fptr = NULL;
}



void ThinSystem::clear ()
{
  SystemBase::clear();

  //init_system_fptr = assemble_fptr = NULL;
}




void ThinSystem::init ()
{
  assert (_mesh.is_prepared());
  
  // initialize parent data
  SystemBase::init();

  // Possibly call a user-supplied initialization
  // method.
  if (init_system_fptr != NULL)
    this->init_system_fptr (_equation_systems, this->name());
}



void ThinSystem::reinit ()
{
  assert (_mesh.is_prepared());
  
  // initialize parent data
  SystemBase::reinit();

  error();
}



void ThinSystem::assemble ()
{
  assert (assemble_fptr != NULL);

  // prepare matrix with the help of the _dof_map, 
  // fill with sparsity pattern
  SystemBase::assemble();

  // Log how long the user's matrix assembly code takes
  START_LOG("assemble()", "ThinSystem");
  
  // Call the user-specified matrix assembly function
  this->assemble_fptr (_equation_systems, this->name());

  // Stop logging the user code
  STOP_LOG("assemble()", "ThinSystem");
}



std::pair<unsigned int, Real>
ThinSystem::solve ()
{
  // Assemble the linear system
  this->assemble (); 

  // Log how long the linear solve takes.
  START_LOG("solve()", "ThinSystem");
  
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
  STOP_LOG("solve()", "ThinSystem");

  return rval; 
}




bool ThinSystem::compare (const ThinSystem& other_system, 
			  const Real threshold,
			  const bool verbose) const
{
  // we have no additional data to compare,
  // let SystemBase do the job
  const SystemBase& other_system_base = static_cast<const SystemBase&>(other_system);
  return SystemBase::compare (other_system_base, threshold, verbose);
}




void ThinSystem::attach_init_function(void fptr(EquationSystems<ThinSystem>& es,
						const std::string& name))
{
  assert (fptr != NULL);
  
  init_system_fptr = fptr;
}



void ThinSystem::attach_assemble_function(void fptr(EquationSystems<ThinSystem>& es,
						    const std::string& name))
{
  assert (fptr != NULL);
  
  assemble_fptr = fptr;  
}
