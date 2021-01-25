// The libMesh Finite Element Library.
// Copyright (C) 2002-2021 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



// C++ includes

// Local Includes
#include "libmesh/petsc_matrix.h"
#include "libmesh/shell_matrix.h"
#include "libmesh/petsc_shell_matrix.h"

namespace libMesh
{

  // Full specialization for Real datatypes
  template <typename T>
  std::unique_ptr<ShellMatrix<T>>
  ShellMatrix<T>::build(const Parallel::Communicator & comm,
                           const SolverPackage solver_package)
  {
    // Avoid unused parameter warnings when no solver packages are enabled.
    libmesh_ignore(comm);

    // Build the appropriate vector
    switch (solver_package)
      {
  #ifdef LIBMESH_HAVE_PETSC
      case PETSC_SOLVERS:
        return libmesh_make_unique<PetscShellMatrix<T>>(comm);
  #endif

      default:
        libmesh_error_msg("ERROR:  Unrecognized solver package: " << solver_package);
      }
  }

//------------------------------------------------------------------
// Explicit instantiations
template class ShellMatrix<Number>;

} // namespace libMesh
