// $Id: petsc_macro.h,v 1.1 2007-07-16 15:12:36 jwpeterson Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2005  Benjamin S. Kirk, John W. Peterson
  
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

#ifndef __petsc_macro_h__
#define __petsc_macro_h__

// Local includes
#include "libmesh_config.h"

#ifdef HAVE_PETSC

// Petsc include files
//#include <petscversion.h> 

// Can't include petscversion unless petsc.h is also included
#ifndef USE_COMPLEX_NUMBERS
extern "C" {
#include <petsc.h>
#include <petscversion.h>
}
#else
#include <petsc.h>
#include <petscversion.h>
#endif

// A convenient macro for comparing PETSc versions.
// Returns 1 if the current PETSc version is < major.minor.subminor
// and zero otherwise.
#define PETSC_VERSION_LESS_THAN(major,minor,subminor)			            \
  ((PETSC_VERSION_MAJOR < (major) ||						    \
    (PETSC_VERSION_MAJOR == (major) && (PETSC_VERSION_MINOR < (minor) ||	    \
				  (PETSC_VERSION_MINOR == (minor) &&		    \
				   PETSC_VERSION_SUBMINOR < (subminor))))) ? 1 : 0)


#endif // HAVE_PETSC

#endif // __petsc_macro_h__
