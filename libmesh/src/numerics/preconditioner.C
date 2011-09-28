// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
  
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
#include "auto_ptr.h"
#include "preconditioner.h"
#include "petsc_preconditioner.h"
#include "trilinos_preconditioner.h"

namespace libMesh
{

//------------------------------------------------------------------
// Preconditioner members
template <typename T>
Preconditioner<T> *
Preconditioner<T>::build(const SolverPackage solver_package)
{
  // Build the appropriate solver
  switch (solver_package)
    {

/*
#ifdef LIBMESH_HAVE_LASPACK
    case LASPACK_SOLVERS:
      {
	AutoPtr<Preconditioner<T> > ap(new LaspackPreconditioner<T>);
	return ap;
      }
#endif
*/

#ifdef LIBMESH_HAVE_PETSC
    case PETSC_SOLVERS:
      {
	return new PetscPreconditioner<T>();
      }
#endif

#ifdef LIBMESH_HAVE_TRILINOS
    case TRILINOS_SOLVERS:
      return new TrilinosPreconditioner<T>();
#endif

    default:
      libMesh::err << "ERROR:  Unrecognized solver package: "
		    << solver_package
		    << std::endl;
      libmesh_error();
    }
    
  return NULL;    
}



//------------------------------------------------------------------
// Explicit instantiations
template class Preconditioner<Number>;

} // namespace libMesh



