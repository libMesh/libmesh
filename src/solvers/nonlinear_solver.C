// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/auto_ptr.h" // libmesh_make_unique
#include "libmesh/nonlinear_solver.h"
#include "libmesh/petsc_nonlinear_solver.h"
#include "libmesh/trilinos_nox_nonlinear_solver.h"
#include "libmesh/solver_configuration.h"
#include "libmesh/enum_solver_package.h"

namespace libMesh
{


//------------------------------------------------------------------
// NonlinearSolver members
#if defined(LIBMESH_HAVE_PETSC) || defined(LIBMESH_TRILINOS_HAVE_NOX)
template <typename T>
std::unique_ptr<NonlinearSolver<T>>
NonlinearSolver<T>::build(sys_type & s, const SolverPackage solver_package)
{
  // Build the appropriate solver
  switch (solver_package)
    {

#ifdef LIBMESH_HAVE_PETSC
    case PETSC_SOLVERS:
      return libmesh_make_unique<PetscNonlinearSolver<T>>(s);
#endif // LIBMESH_HAVE_PETSC

#if defined(LIBMESH_TRILINOS_HAVE_NOX) && defined(LIBMESH_TRILINOS_HAVE_EPETRA)
    case TRILINOS_SOLVERS:
      return libmesh_make_unique<NoxNonlinearSolver<T>>(s);
#endif

    default:
      libmesh_error_msg("ERROR:  Unrecognized solver package: " << solver_package);
    }
}

#else // LIBMESH_HAVE_PETSC || LIBMESH_TRILINOS_HAVE_NOX

template <typename T>
std::unique_ptr<NonlinearSolver<T>>
NonlinearSolver<T>::build(sys_type &, const SolverPackage)
{
  libmesh_not_implemented_msg("ERROR: libMesh was compiled without nonlinear solver support");
}
#endif


template <typename T>
void
NonlinearSolver<T>::attach_preconditioner(Preconditioner<T> * preconditioner)
{
  if (this->_is_initialized)
    libmesh_error_msg("Preconditioner must be attached before the solver is initialized!");

  _preconditioner = preconditioner;
}

template <typename T>
void NonlinearSolver<T>::set_solver_configuration(SolverConfiguration & solver_configuration)
{
  _solver_configuration = &solver_configuration;
}


//------------------------------------------------------------------
// Explicit instantiations
template class NonlinearSolver<Number>;

} // namespace libMesh
