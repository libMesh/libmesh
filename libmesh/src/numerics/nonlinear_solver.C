// $Id: nonlinear_solver.C,v 1.1 2005-01-03 20:14:41 benkirk Exp $

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

// Local Includes
#include "nonlinear_solver.h"
#include "petsc_nonlinear_solver.h"
#include "auto_ptr.h"


//------------------------------------------------------------------
// NonlinearSolver members
template <typename T>
AutoPtr<NonlinearSolver<T> >
NonlinearSolver<T>::build(const SolverPackage solver_package)
{
  // Build the appropriate solver
  switch (solver_package)
    {

#ifdef HAVE_PETSC
    case PETSC_SOLVERS:
      {
	AutoPtr<NonlinearSolver<T> > ap(new PetscNonlinearSolver<T>);
	return ap;
      }
#endif

    default:
      std::cerr << "ERROR:  Unrecognized solver package: "
		<< solver_package
		<< std::endl;
      error();
    }
    
  AutoPtr<NonlinearSolver<T> > ap(NULL);
  return ap;    
}



//------------------------------------------------------------------
// Explicit instantiations
template class NonlinearSolver<Number>;



