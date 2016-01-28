// The libMesh Finite Element Library.
// Copyright (C) 2002-2014 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#ifndef LIBMESH_PETSC_SOLVER_EXCEPTION_H
#define LIBMESH_PETSC_SOLVER_EXCEPTION_H

#include "libmesh/libmesh_common.h"

#ifdef LIBMESH_HAVE_PETSC

#include "libmesh_exceptions.h"
#include <petscsys.h>

namespace libMesh
{

/**
 * A specialization of the SolverException class for PETSc.
 */
class PetscSolverException : public SolverException
{
public:
  PetscSolverException(int error_code_in) :
    SolverException(error_code_in)
  {
    const char * text;
    // This is one scenario where we don't catch the error code
    // returned by a PETSc function :)
    PetscErrorMessage(error_code, &text, libmesh_nullptr);
    what_message = text;
  }
};



// Macro which we call after every PETSc function that returns an error code.
#ifdef LIBMESH_ENABLE_EXCEPTIONS

#define LIBMESH_CHKERR(ierr)                    \
  do {                                          \
    if (ierr != 0) {                            \
      throw PetscSolverException(ierr);         \
    } } while(0)

#else

// If we don't have exceptions enabled, just fall back on calling
// PETSc's CHKERRABORT macro.
#define LIBMESH_CHKERR(ierr) CHKERRABORT(this->comm().get(), ierr);

// Let's also be backwards-compatible with the old macro name.
#define LIBMESH_CHKERRABORT(ierr) LIBMESH_CHKERR(ierr)

#endif

} // namespace libMesh

#endif // LIBMESH_HAVE_PETSC
#endif // LIBMESH_PETSC_SOLVER_EXCEPTION_H
