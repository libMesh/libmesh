// $Id: numeric_vector.C,v 1.2 2003-02-11 00:58:24 benkirk Exp $

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
#include "numeric_vector.h"
#include "laspack_vector.h"
#include "petsc_vector.h"



//------------------------------------------------------------------
//
AutoPtr<NumericVector>
NumericVector::build(const SolverPackage solver_package)
{

  switch (solver_package)
    {


#ifdef HAVE_LASPACK
    case LASPACK_SOLVERS:
      {
	AutoPtr<NumericVector> ap(new LaspackVector);
	return ap;
      };
#endif


#ifdef HAVE_PETSC
    case PETSC_SOLVERS:
      {
	AutoPtr<NumericVector> ap(new PetscVector);
	return ap;
      };
#endif

    default:
      std::cerr << "ERROR:  Unrecognized solver package: "
		<< solver_package
		<< std::endl;
      error();
    };
    
  AutoPtr<NumericVector> ap(NULL);
  return ap;    
};
