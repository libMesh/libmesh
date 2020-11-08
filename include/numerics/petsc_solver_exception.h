// The libMesh Finite Element Library.
// Copyright (C) 2002-2020 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#include "libmesh/libmesh_exceptions.h"

#ifdef I
# define LIBMESH_SAW_I
#endif

#include "libmesh/ignore_warnings.h"
#include <petscsys.h>
#include "libmesh/restore_warnings.h"

#ifndef LIBMESH_SAW_I
# undef I // Avoid complex.h contamination
#endif

namespace libMesh
{

// Forward declarations
namespace Parallel
{
class Communicator;
}

// The SolverException class is only defined when exceptions are enabled.
#ifdef LIBMESH_ENABLE_EXCEPTIONS

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
    char * specific;
    // This is one scenario where we don't catch the error code
    // returned by a PETSc function :)
    PetscErrorMessage(error_code, &text, &specific);

    // Usually the "specific" error message string is more useful than
    // the generic text corresponding to the error_code, since many
    // SETERRQ calls just use error_code == 1
    if (specific)
      what_message = std::string(specific);
    else if (text)
      what_message = std::string(text);
  }
};



// Macro which we call after every PETSc function that returns an error code.
#define LIBMESH_CHKERR(ierr)                    \
  do {                                          \
    if (ierr != 0) {                            \
      throw PetscSolverException(ierr);         \
    } } while (0)

// Two-argument CHKERR macro that takes both a comm and an error
// code. When exceptions are enabled, the comm is not used for
// anything, so we libmesh_ignore() it.  This macro is useful when you
// need to call LIBMESH_CHKERR but you're not in a ParallelObject.
#define LIBMESH_CHKERR2(comm, ierr)             \
  do {                                          \
    libmesh_ignore(comm);                       \
    LIBMESH_CHKERR(ierr);                       \
  } while (0)

#else

// If we don't have exceptions enabled, just fall back on calling
// PETSc's CHKERRABORT macro.
#define LIBMESH_CHKERR(ierr) CHKERRABORT(this->comm().get(), ierr);

// Two argument version of the function above where you pass in the comm
// instead of relying on it being available from the "this" pointer.
#define LIBMESH_CHKERR2(comm, ierr) CHKERRABORT(comm.get(), ierr);

// Let's also be backwards-compatible with the old macro name.
#define LIBMESH_CHKERRABORT(ierr) LIBMESH_CHKERR(ierr)

#endif

#define PETSC_BEGIN_END(Function)                                       \
  template<class ...Args>                                               \
  inline                                                                \
  void Function ## BeginEnd(const Parallel::Communicator & comm, const Args&... args) \
  {                                                                     \
    PetscErrorCode ierr = 0;                                            \
    ierr = Function ## Begin(args...);                                  \
    LIBMESH_CHKERR2(comm, ierr);                                        \
    ierr = Function ## End(args...);                                    \
    LIBMESH_CHKERR2(comm, ierr);                                        \
  }

PETSC_BEGIN_END(VecScatter) // VecScatterBeginEnd
PETSC_BEGIN_END(MatAssembly) // MatAssemblyBeginEnd
PETSC_BEGIN_END(VecAssembly) // VecAssemblyBeginEnd
PETSC_BEGIN_END(VecGhostUpdate) // VecGhostUpdateBeginEnd

} // namespace libMesh

#endif // LIBMESH_HAVE_PETSC
#endif // LIBMESH_PETSC_SOLVER_EXCEPTION_H
