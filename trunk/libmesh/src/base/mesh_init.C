// $Id: mesh_init.C,v 1.3 2003-02-13 22:56:08 benkirk Exp $

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




// Local includes
#include "mesh_init.h"


#if defined(HAVE_PETSC)
#  if  !defined(USE_COMPLEX_NUMBERS)
/**
 * Petsc include files.  PETSc with complex numbers 
 * is actually C++.
 */
namespace Petsc
{
  extern "C"
  {
#   include <petsc.h>
  }
}
using namespace Petsc;
// for easy switching between Petsc 2.1.0/2.1.1
// typedef Scalar PetscScalar;
#  else
#    include <petsc.h>
#  endif
#elif defined(HAVE_MPI)
namespace Mpi
{
  extern "C"
  {
#    include <mpi.h>
  }
}
using namespace Mpi;
#endif


// ------------------------------------------------------------
//
void libMesh::init (int argc, char** argv)
{

#if defined(HAVE_PETSC)

  PetscInitialize (&argc, &argv, NULL, NULL);

#elif defined(HAVE_MPI)

  MPI_Init (&argc, &argv);
  
#endif
  
}



void libMesh::close ()
{

#if defined(HAVE_PETSC)

  PetscFinalize ();

#elif defined(HAVE_MPI)

  MPI_Finalize();
  
#endif
  
}
