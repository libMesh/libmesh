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



// C++ includes

// Local Includes
#include "libmesh/nonlinear_solver.h"
#include "libmesh/petsc_nonlinear_solver.h"
#include "libmesh/petsc_dm_nonlinear_solver.h"
#include "libmesh/trilinos_nox_nonlinear_solver.h"
#include "libmesh/auto_ptr.h"

namespace libMesh
{


//------------------------------------------------------------------
// NonlinearSolver members
#if defined(LIBMESH_HAVE_PETSC) || defined(LIBMESH_HAVE_NOX)
template <typename T>
AutoPtr<NonlinearSolver<T> >
NonlinearSolver<T>::build(sys_type& s, const SolverPackage solver_package)
{
  // Build the appropriate solver
  switch (solver_package)
    {

#ifdef LIBMESH_HAVE_PETSC
    case PETSC_SOLVERS:
#if PETSC_VERSION_LESS_THAN(3,3,0)
      return AutoPtr<NonlinearSolver<T> >(new PetscNonlinearSolver<T>(s));
#else
      if (libMesh::on_command_line ("--use-petsc-dm"))
        return AutoPtr<NonlinearSolver<T> >(new PetscDMNonlinearSolver<T>(s));

      else
        return AutoPtr<NonlinearSolver<T> >(new PetscNonlinearSolver<T>(s));

#endif
#endif // LIBMESH_HAVE_PETSC

#ifdef LIBMESH_HAVE_NOX
    case TRILINOS_SOLVERS:
      return AutoPtr<NonlinearSolver<T> >(new NoxNonlinearSolver<T>(s));
#endif

    default:
      libmesh_error_msg("ERROR:  Unrecognized solver package: " << solver_package);
    }

  libmesh_error_msg("We'll never get here!");
  return AutoPtr<NonlinearSolver<T> >();
}

#else // LIBMESH_HAVE_PETSC || LIBMESH_HAVE_NOX

template <typename T>
AutoPtr<NonlinearSolver<T> >
NonlinearSolver<T>::build(sys_type&, const SolverPackage)
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


//------------------------------------------------------------------
// Explicit instantiations
template class NonlinearSolver<Number>;

} // namespace libMesh
