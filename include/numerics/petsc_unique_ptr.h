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

#ifndef LIBMESH_PETSC_UNIQUE_PTR_H
#define LIBMESH_PETSC_UNIQUE_PTR_H

#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_HAVE_PETSC

// libMesh includes
#include "libmesh/petsc_solver_exception.h"

// PETSc includes
#ifdef I
# define LIBMESH_SAW_I
#endif

#include <petscksp.h>
#include <petscvec.h>
#include <petscis.h>
#include <petscmat.h>

#ifndef LIBMESH_SAW_I
# undef I // Avoid complex.h contamination
#endif

// C++ includes
#include <memory>

namespace libMesh
{

/**
 * Generic class which defines an operator() and can be used as a
 * custom deleter object for a std::unique_ptr
 */
template <typename U>
struct PetscDeleter
{
  // To be specialized for different PETSc types
  void operator()(U * obj);
};

/**
 * Macro that defines a full specialization of the PetscDeleter class
 * for a specific Type.  Explicit specializations have to appear at
 * namespace scope.
 */
#define PETSC_DELETER(Type)                             \
  template <>                                           \
  struct PetscDeleter<Type>                             \
  {                                                     \
    void operator()(Type * obj)                         \
    {                                                   \
      PetscErrorCode ierr = Type ## Destroy(obj);       \
      LIBMESH_CHKERR(ierr);                             \
    }                                                   \
  }

PETSC_DELETER(KSP);
PETSC_DELETER(IS);
PETSC_DELETER(Vec);
PETSC_DELETER(Mat);
PETSC_DELETER(VecScatter);

// Template alias for a unique_ptr with a custom deleter
template <typename T>
using petsc_unique_ptr = std::unique_ptr<T, PetscDeleter<T>>;

} // namespace libMesh

#endif // LIBMESH_HAVE_PETSC

#endif
