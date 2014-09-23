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


#include "libmesh/libmesh_config.h"
#ifdef LIBMESH_HAVE_SLEPC

// C++ includes

// Local Includes
#include "libmesh/eigen_solver.h"
#include "libmesh/slepc_eigen_solver.h"

namespace libMesh
{


//------------------------------------------------------------------
// EigenSolver members
template <typename T>
UniquePtr<EigenSolver<T> >
EigenSolver<T>::build(const Parallel::Communicator &comm,
                      const SolverPackage solver_package)
{
  // Build the appropriate solver
  switch (solver_package)
    {

#ifdef LIBMESH_HAVE_SLEPC
    case SLEPC_SOLVERS:
      return UniquePtr<EigenSolver<T> >(new SlepcEigenSolver<T>(comm));
#endif

    default:
      libmesh_error_msg("ERROR:  Unrecognized eigen solver package: " << solver_package);
    }

  return UniquePtr<EigenSolver<T> >();
}



//------------------------------------------------------------------
// Explicit instantiations
template class EigenSolver<Number>;

} // namespace libMesh


#endif // LIBMESH_HAVE_SLEPC
