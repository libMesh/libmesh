// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2007  Benjamin S. Kirk, John W. Peterson
  
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

// A convenient macro for comparing PETSc versions.
// Returns 1 if the current PETSc version is < major.minor.subminor
// and zero otherwise.
#define PETSC_VERSION_LESS_THAN(major,minor,subminor)			            \
  ((PETSC_VERSION_MAJOR < (major) ||						    \
    (PETSC_VERSION_MAJOR == (major) && (PETSC_VERSION_MINOR < (minor) ||	    \
				  (PETSC_VERSION_MINOR == (minor) &&		    \
				   PETSC_VERSION_SUBMINOR < (subminor))))) ? 1 : 0)

// We have to jump through some hoops here: in Petsc 2.3.1 you cannot
// include petscversion.h without including petsc.h beforehand.  This
// was because petscversion.h relied on the existence of
// PETSC_EXTERN_CXX_BEGIN/END in 2.3.1.  The problem is, we don't know
// how to include petsc.h (wrapped or not wrapped in extern "C") until
// we know the version!

// First figure out if we need to define PETSC_EXTERN_CXX_BEGIN/END to
// make petscversion.h happy
#ifndef PETSC_EXTERN_CXX_BEGIN
#  define PETSC_EXTERN_CXX_BEGIN
#  define PETSC_EXTERN_CXX_END
#  define LIBMESH_PETSC_EXTERN_C_WORKAROUND
#endif

// Now actually include it
#include <petscversion.h>

// And finally, get rid of our bogus-definitions.  <petsc.h> will set these itself.
#ifdef LIBMESH_PETSC_EXTERN_C_WORKAROUND
#  undef PETSC_EXTERN_CXX_BEGIN
#  undef PETSC_EXTERN_CXX_END
#endif

// Make up for missing extern "C" in old PETSc versions
#if !defined(USE_COMPLEX_NUMBERS) && PETSC_VERSION_LESS_THAN(2,3,0)
#  define EXTERN_C_FOR_PETSC_BEGIN extern "C" {
#  define EXTERN_C_FOR_PETSC_END }
#else
#  define EXTERN_C_FOR_PETSC_BEGIN 
#  define EXTERN_C_FOR_PETSC_END
#endif


// Petsc include files
EXTERN_C_FOR_PETSC_BEGIN
#include <petsc.h>
EXTERN_C_FOR_PETSC_END


#endif // HAVE_PETSC

#endif // __petsc_macro_h__
