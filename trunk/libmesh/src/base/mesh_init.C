// $Id: mesh_init.C,v 1.1 2003-02-10 03:55:51 benkirk Exp $

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



#ifdef HAVE_PETSC

#ifndef USE_COMPLEX_NUMBERS
/**
 * Petsc include files.  PETSc with complex numbers 
 * is actually C++.
 */
namespace Petsc {
extern "C" {
#include <petsc.h>
}
}
using namespace Petsc;
// for easy switching between Petsc 2.1.0/2.1.1
// typedef Scalar PetscScalar;
#else
#include <petsc.h>
#endif 
#endif


// ------------------------------------------------------------
//
void libMesh::init (int argc, char** argv)
{

#ifdef HAVE_PETSC

  PetscInitialize (&argc, &argv, NULL, NULL);
  
#endif

  
};



void libMesh::close ()
{

#ifdef HAVE_PETSC

  PetscFinalize ();
  
#endif
  
};
