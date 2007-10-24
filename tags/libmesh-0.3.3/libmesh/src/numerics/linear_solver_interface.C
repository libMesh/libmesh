// $Id: linear_solver_interface.C,v 1.6 2003-03-17 01:27:08 benkirk Exp $

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
#include <math.h>

// Local Includes
#include "linear_solver_interface.h"
#include "laspack_interface.h"
#include "petsc_interface.h"



//------------------------------------------------------------------
// LinearSolverInterface members
template <typename T>
AutoPtr<LinearSolverInterface<T> >
LinearSolverInterface<T>::build(const SolverPackage solver_package)
{

  switch (solver_package)
    {


#if defined(HAVE_LASPACK)
    case LASPACK_SOLVERS:
      {
	AutoPtr<LinearSolverInterface<T> > ap(new LaspackInterface<T>);
	return ap;
      }
#endif


#if defined(HAVE_PETSC)
    case PETSC_SOLVERS:
      {
	AutoPtr<LinearSolverInterface<T> > ap(new PetscInterface<T>);
	return ap;
      }
#endif

    default:
      std::cerr << "ERROR:  Unrecognized solver package: "
		<< solver_package
		<< std::endl;
      error();
    }
    
  AutoPtr<LinearSolverInterface<T> > ap(NULL);
  return ap;    
}



//------------------------------------------------------------------
// Explicit instantiations
template class LinearSolverInterface<Number>;



