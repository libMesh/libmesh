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



// Local Includes
#include "libmesh/preconditioner.h"
#include "libmesh/auto_ptr.h"
#include "libmesh/eigen_preconditioner.h"
#include "libmesh/petsc_preconditioner.h"
#include "libmesh/trilinos_preconditioner.h"
#include "libmesh/enum_solver_package.h"
#include "libmesh/enum_preconditioner_type.h"

namespace libMesh
{

template <typename T>
inline
Preconditioner<T>::Preconditioner (const libMesh::Parallel::Communicator & comm_in) :
  ParallelObject(comm_in),
  _matrix(nullptr),
  _preconditioner_type (ILU_PRECOND),
  _is_initialized      (false)
{
}



template <typename T>
std::unique_ptr<Preconditioner<T>>
Preconditioner<T>::build_preconditioner(const libMesh::Parallel::Communicator & comm,
                                        const SolverPackage solver_package)
{
  // Avoid unused parameter warnings when no solver packages are enabled.
  libmesh_ignore(comm);

  // Build and return the appropriate Preconditioner object.
  switch (solver_package)
    {

#ifdef LIBMESH_HAVE_PETSC
    case PETSC_SOLVERS:
      {
        return libmesh_make_unique<PetscPreconditioner<T>>(comm);
      }
#endif

#ifdef LIBMESH_TRILINOS_HAVE_EPETRA
    case TRILINOS_SOLVERS:
      return libmesh_make_unique<TrilinosPreconditioner<T>>(comm);
#endif

#ifdef LIBMESH_HAVE_EIGEN
    case EIGEN_SOLVERS:
      return libmesh_make_unique<EigenPreconditioner<T>>(comm);
#endif

    default:
      libmesh_error_msg("ERROR:  Unrecognized solver package: " << solver_package);
    }
}



//------------------------------------------------------------------
// Explicit instantiations
template class Preconditioner<Number>;

} // namespace libMesh
