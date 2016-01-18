// The libMesh Finite Element Library.
// Copyright (C) 2002-2016 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_SLEPC_MACRO_H
#define LIBMESH_SLEPC_MACRO_H

// C++ includes

// Local includes

/**
 * SLEPc include files. SLEPc can only be used
 * together with PETSc.
 */
#ifdef LIBMESH_HAVE_SLEPC

# include <slepcversion.h>

// A convenient macro for comparing SLEPc versions.
// Returns 1 if the current SLEPc version is < major.minor.subminor
// and zero otherwise.
#define SLEPC_VERSION_LESS_THAN(major,minor,subminor)                   \
  ((SLEPC_VERSION_MAJOR < (major) ||                                    \
    (SLEPC_VERSION_MAJOR == (major) && (SLEPC_VERSION_MINOR < (minor) || \
                                        (SLEPC_VERSION_MINOR == (minor) && \
                                         SLEPC_VERSION_SUBMINOR < (subminor))))) ? 1 : 0)

// Make up for missing extern "C" in old SLEPc versions
#if !defined(LIBMESH_USE_COMPLEX_NUMBERS) && SLEPC_VERSION_LESS_THAN(3,0,0)
#  define EXTERN_C_FOR_SLEPC_BEGIN extern "C" {
#  define EXTERN_C_FOR_SLEPC_END }
#else
#  define EXTERN_C_FOR_SLEPC_BEGIN
#  define EXTERN_C_FOR_SLEPC_END
#endif

#if SLEPC_VERSION_RELEASE && SLEPC_VERSION_LESS_THAN(3,1,1)
#  define LibMeshEPSDestroy(x)         EPSDestroy(*(x))
#else
#  define LibMeshEPSDestroy(x)         EPSDestroy(x)
#endif

#endif // #if LIBMESH_HAVE_SLEPC
#endif // LIBMESH_SLEPC_MACRO_H
