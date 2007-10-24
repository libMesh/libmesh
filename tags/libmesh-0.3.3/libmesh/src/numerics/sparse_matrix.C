// $Id: sparse_matrix.C,v 1.5 2003-03-14 09:56:41 ddreyer Exp $

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
template <>
AutoPtr<SparseMatrix<Real> >
SparseMatrix<Real>::build(const SolverPackage solver_package)
{

  switch (solver_package)
    {


#if defined(HAVE_LASPACK) && defined(USE_REAL_NUMBERS)
    case LASPACK_SOLVERS:
      {
	AutoPtr<SparseMatrix<Real> > ap(new LaspackMatrix<Real>);
	return ap;
      }
#endif


#if defined(HAVE_PETSC) && defined(USE_REAL_NUMBERS)
    case PETSC_SOLVERS:
      {
	AutoPtr<SparseMatrix<Real> > ap(new PetscMatrix<Real>);
	return ap;
      }
#endif

    default:
      std::cerr << "ERROR:  Unrecognized solver package: "
		<< solver_package
		<< std::endl;
      error();
    }

  AutoPtr<SparseMatrix<Real> > ap(NULL);
  return ap;    
}




// Full specialization for Complex datatypes
template <>
AutoPtr<SparseMatrix<Complex> >
SparseMatrix<Complex>::build(const SolverPackage solver_package)
{

  switch (solver_package)
    {


#if defined(HAVE_LASPACK) && defined(USE_COMPLEX_NUMBERS)
    case LASPACK_SOLVERS:
      {
	AutoPtr<SparseMatrix<Complex> > ap(new LaspackMatrix<Complex>);
	return ap;
      }
#endif


#if defined(HAVE_PETSC) && defined(USE_COMPLEX_NUMBERS)
    case PETSC_SOLVERS:
      {
	AutoPtr<SparseMatrix<Complex> > ap(new PetscMatrix<Complex>);
	return ap;
      }
#endif

    default:
      std::cerr << "ERROR:  Unrecognized solver package: "
		<< solver_package
		<< std::endl;
      error();
    }

  AutoPtr<SparseMatrix<Complex> > ap(NULL);
  return ap;    
}


//------------------------------------------------------------------
// Explicit instantiations
template class SparseMatrix<Number>;
