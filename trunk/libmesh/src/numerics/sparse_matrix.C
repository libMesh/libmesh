// $Id: sparse_matrix.C,v 1.7 2003-04-09 01:20:25 benkirk Exp $

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
#include "sparse_matrix.h"
#include "laspack_matrix.h"
#include "petsc_matrix.h"



//------------------------------------------------------------------
// SparseMatrix Methods

// Full specialization for Real datatypes
template <typename T>
AutoPtr<SparseMatrix<T> >
SparseMatrix<T>::build(const SolverPackage solver_package_in)
{
  // Possibly overload the solver package based on
  // command-line arguments
  SolverPackage solver_package = solver_package_in;
  
#ifdef HAVE_PETSC
  if (libMesh::on_command_line ("--use-petsc"))
    solver_package = PETSC_SOLVERS;
#endif
  
#ifdef HAVE_LASPACK
  if (libMesh::on_command_line("--use-laspack"))
    solver_package = LASPACK_SOLVERS;
#endif




  // Build the appropriate vector
  switch (solver_package)
    {


#ifdef HAVE_LASPACK
    case LASPACK_SOLVERS:
      {
	AutoPtr<SparseMatrix<T> > ap(new LaspackMatrix<T>);
	return ap;
      }
#endif


#ifdef HAVE_PETSC
    case PETSC_SOLVERS:
      {
	AutoPtr<SparseMatrix<T> > ap(new PetscMatrix<T>);
	return ap;
      }
#endif

    default:
      std::cerr << "ERROR:  Unrecognized solver package: "
		<< solver_package
		<< std::endl;
      error();
    }

  AutoPtr<SparseMatrix<T> > ap(NULL);
  return ap;    
}



//------------------------------------------------------------------
// Explicit instantiations
template class SparseMatrix<Number>;
