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
  {}
};


// Macro which we call after every PETSc function that returns an error code.
#define LIBMESH_CHKERRABORT(ierr)               \
  do {                                          \
    if (ierr != 0) {                            \
      throw PetscSolverException(ierr);         \
    } } while(0)

} // namespace libMesh

#endif // LIBMESH_HAVE_PETSC
#endif // LIBMESH_PETSC_SOLVER_EXCEPTION_H
