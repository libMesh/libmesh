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

#ifndef LIBMESH_PETSC_MACRO_H
#define LIBMESH_PETSC_MACRO_H

// Local includes
#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_HAVE_PETSC

// A convenient macro for comparing PETSc versions.  Returns 1 if the
// current PETSc version is < major.minor.subminor and zero otherwise.
//
// This macro does not require petscversion.h to be included for it to work correctly.
// It instead relies on the PETSc version numbers detected during configure.  Note that if
// LIBMESH_HAVE_PETSC is not defined, none of the LIBMESH_DETECTED_PETSC_VERSION_* variables will
// be defined either.
#define PETSC_VERSION_LESS_THAN(major,minor,subminor)                   \
  ((LIBMESH_DETECTED_PETSC_VERSION_MAJOR < (major) ||                   \
    (LIBMESH_DETECTED_PETSC_VERSION_MAJOR == (major) && (LIBMESH_DETECTED_PETSC_VERSION_MINOR < (minor) || \
                                                         (LIBMESH_DETECTED_PETSC_VERSION_MINOR == (minor) && \
                                                          LIBMESH_DETECTED_PETSC_VERSION_SUBMINOR < (subminor))))) ? 1 : 0)

// The PETSC_VERSION_RELEASE constant was introduced just prior to 2.3.0 (ca. Apr 22 2005),
// so fall back to using PETSC_VERSION_LESS_THAN in case it doesn't exist.
#ifdef LIBMESH_DETECTED_PETSC_VERSION_RELEASE

#define PETSC_RELEASE_LESS_THAN(major,minor,subminor)                   \
  (PETSC_VERSION_LESS_THAN(major,minor,subminor) && LIBMESH_DETECTED_PETSC_VERSION_RELEASE)

#else

#define PETSC_RELEASE_LESS_THAN(major,minor,subminor)   \
  (PETSC_VERSION_LESS_THAN(major,minor,subminor))

#endif

// We used to have workarounds for missing extern "C" in old PETSc
// versions.  We no longer support PETSc versions so old, but we do
// still support libMesh applications old enough to have used these
// macros.
#define EXTERN_C_FOR_PETSC_BEGIN
#define EXTERN_C_FOR_PETSC_END

// Petsc include files
#include <petsc.h>

#if PETSC_RELEASE_LESS_THAN(3,1,1)
typedef PetscTruth PetscBool;
#endif

#if PETSC_RELEASE_LESS_THAN(3,1,1)
#  define LibMeshVecDestroy(x)         VecDestroy(*(x))
#  define LibMeshVecScatterDestroy(x)  VecScatterDestroy(*(x))
#  define LibMeshMatDestroy(x)         MatDestroy(*(x))
#  define LibMeshISDestroy(x)          ISDestroy(*(x))
#  define LibMeshKSPDestroy(x)         KSPDestroy(*(x))
#  define LibMeshSNESDestroy(x)        SNESDestroy(*(x))
#  define LibMeshPetscViewerDestroy(x) PetscViewerDestroy(*(x))
#  define LibMeshPCDestroy(x)          PCDestroy(*(x))
#else
#  define LibMeshVecDestroy(x)         VecDestroy(x)
#  define LibMeshVecScatterDestroy(x)  VecScatterDestroy(x)
#  define LibMeshMatDestroy(x)         MatDestroy(x)
#  define LibMeshISDestroy(x)          ISDestroy(x)
#  define LibMeshKSPDestroy(x)         KSPDestroy(x)
#  define LibMeshSNESDestroy(x)        SNESDestroy(x)
#  define LibMeshPetscViewerDestroy(x) PetscViewerDestroy(x)
#  define LibMeshPCDestroy(x)          PCDestroy(x)
#endif

#if PETSC_RELEASE_LESS_THAN(3,1,1)
typedef enum { PETSC_COPY_VALUES, PETSC_OWN_POINTER, PETSC_USE_POINTER} PetscCopyMode;
#  define ISCreateLibMesh(comm,n,idx,mode,is)           \
  ((mode) == PETSC_USE_POINTER                          \
   ? ISCreateGeneralWithArray((comm),(n),(idx),(is))    \
   : ((mode) == PETSC_OWN_POINTER                       \
      ? ISCreateGeneralNC((comm),(n),(idx),(is))        \
      : ISCreateGeneral((comm),(n),(idx),(is))))
#else
#  define ISCreateLibMesh(comm,n,idx,mode,is) ISCreateGeneral((comm),(n),(idx),(mode),(is))
#endif

#else // LIBMESH_HAVE_PETSC

#define PETSC_VERSION_LESS_THAN(major,minor,subminor) 1
#define PETSC_RELEASE_LESS_THAN(major,minor,subminor) 1

#endif // LIBMESH_HAVE_PETSC

#endif // LIBMESH_PETSC_MACRO_H
