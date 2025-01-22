// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/petsc_macro.h"
#include "timpi/communicator.h"

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
    auto ierr = PetscErrorMessage(cast_int<PetscErrorCode>(error_code), &text, &specific);
    // This is one scenario where we ignore the error code
    // returned by a PETSc function :)
    libmesh_ignore(ierr);

    // Usually the "specific" error message string is more useful than
    // the generic text corresponding to the error_code, since many
    // SETERRQ calls just use error_code == 1
    if (specific)
      what_message = std::string(specific);
    else if (text)
      what_message = std::string(text);
  }
};



// Macro which we call in functions returning a datatype convertible
// to PetscErrorCode after every PETSc function that returns an error code.
#define LIBMESH_CHKERRQ(ierr)                    \
  do {                                           \
    if (ierr != 0)                               \
      throw libMesh::PetscSolverException(ierr); \
  } while (0)

// Two-argument CHKERR macro that takes both a comm and an error
// code. When exceptions are enabled, the comm is not used for
// anything, so we libmesh_ignore() it.
#define LIBMESH_CHKERRA(comm, ierr)             \
  do {                                          \
    libmesh_ignore(comm);                       \
    LIBMESH_CHKERRQ(ierr);                      \
  } while (0)

// Remove me: for backward compatibility with MOOSE only
#define LIBMESH_CHKERR(ierr) LIBMESH_CHKERRQ(ierr)
#define LIBMESH_CHKERR2(comm, ierr) LIBMESH_CHKERRA(comm, ierr)
#define LibmeshPetscCall(...)                                           \
  do                                                                    \
  {                                                                     \
    PetscErrorCode libmesh_petsc_call_ierr;                             \
    libmesh_petsc_call_ierr = __VA_ARGS__;                              \
    LIBMESH_CHKERRQ(libmesh_petsc_call_ierr);                           \
  } while (0)

#else

// If we don't have exceptions enabled, just fall back on calling
// PETSc's CHKERRQ or CHKERRABORT macros.
#define LIBMESH_CHKERRQ(ierr) CHKERRQ(ierr);
#define LIBMESH_CHKERRA(comm, ierr) CHKERRABORT(comm, ierr);

// Remove me: for backward compatibility with MOOSE only
#define LibmeshPetscCall(...)                                           \
  do                                                                    \
  {                                                                     \
    PetscErrorCode libmesh_petsc_call_ierr;                             \
    libmesh_petsc_call_ierr = __VA_ARGS__;                              \
    LIBMESH_CHKERRA(this->comm().get(), libmesh_petsc_call_ierr);       \
  } while (0)

#endif

#define PETSC_BEGIN_END(Function)                                       \
  template<class ...Args>                                               \
  inline                                                                \
  void Function ## BeginEnd(const Parallel::Communicator & comm, const Args&... args) \
  {                                                                     \
    PetscErrorCode ierr = LIBMESH_PETSC_SUCCESS;                        \
    ierr = Function ## Begin(args...);                                  \
    LIBMESH_CHKERRA(comm.get(), ierr);                                  \
    ierr = Function ## End(args...);                                    \
    LIBMESH_CHKERRA(comm.get(), ierr);                                  \
  }

PETSC_BEGIN_END(VecScatter) // VecScatterBeginEnd
PETSC_BEGIN_END(MatAssembly) // MatAssemblyBeginEnd
PETSC_BEGIN_END(VecAssembly) // VecAssemblyBeginEnd
PETSC_BEGIN_END(VecGhostUpdate) // VecGhostUpdateBeginEnd

// To use instead of PETSc idioms that'd call CHKERRQ or PetscCall().
#define LibmeshPetscCallQ(...)                                          \
  do                                                                    \
  {                                                                     \
    PetscErrorCode libmesh_petsc_call_ierr;                             \
    libmesh_petsc_call_ierr = __VA_ARGS__;                              \
    LIBMESH_CHKERRQ(libmesh_petsc_call_ierr);                           \
  } while (0)

// To use instead of PETSc idioms that'd call CHKERRABORT or PetscCallAbort().
#define LibmeshPetscCallA(comm, ...)                                    \
  do                                                                    \
  {                                                                     \
    PetscErrorCode libmesh_petsc_call_ierr;                             \
    libmesh_petsc_call_ierr = __VA_ARGS__;                              \
    LIBMESH_CHKERRA(comm, libmesh_petsc_call_ierr);                     \
  } while (0)

// Shortcut for LibmeshPetscCallA for use within a ParallelObject, i.e. when
// we can rely on the communicator being available from the "this" pointer.
/*
#define LibmeshPetscCall(...)                                           \
  LibmeshPetscCallA(this->comm().get(), __VA_ARGS__)
*/

// Shortcut for LibmeshPetscCallA for use when we have a Parallel::Communicator
// available instead of just a bare MPI communicator.
#define LibmeshPetscCall2(comm, ...)                                    \
  LibmeshPetscCallA(comm.get(), __VA_ARGS__)

#ifdef LIBMESH_ENABLE_EXCEPTIONS
#define LibmeshPetscCallExternal(func, ...)                             \
  do {                                                                  \
    const auto libmesh_petsc_call_external_ierr = cast_int<int>(func(__VA_ARGS__)); \
    if (libmesh_petsc_call_external_ierr != 0)                          \
      throw PetscSolverException(libmesh_petsc_call_external_ierr);     \
  } while (0)
#else
#define LibmeshPetscCallExternal(func, ...)                             \
  do {                                                                  \
    const auto libmesh_petsc_call_external_ierr = cast_int<int>(func(__VA_ARGS__)); \
    if (libmesh_petsc_call_external_ierr != 0)                          \
      libmesh_terminate_handler();                                      \
  } while (0)
#endif


} // namespace libMesh

#endif // LIBMESH_HAVE_PETSC
#endif // LIBMESH_PETSC_SOLVER_EXCEPTION_H
