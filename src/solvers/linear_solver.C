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
#include "libmesh/libmesh_logging.h"
#include "libmesh/auto_ptr.h"
#include "libmesh/linear_solver.h"
#include "libmesh/laspack_linear_solver.h"
#include "libmesh/eigen_sparse_linear_solver.h"
#include "libmesh/petsc_linear_solver.h"
#include "libmesh/trilinos_aztec_linear_solver.h"
#include "libmesh/preconditioner.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/string_to_enum.h"

namespace libMesh
{

//------------------------------------------------------------------
// LinearSolver members
template <typename T>
UniquePtr<LinearSolver<T> >
LinearSolver<T>::build(const libMesh::Parallel::Communicator &comm,
                       const SolverPackage solver_package)
{
  // Build the appropriate solver
  switch (solver_package)
    {


#ifdef LIBMESH_HAVE_LASPACK
    case LASPACK_SOLVERS:
      return UniquePtr<LinearSolver<T> >(new LaspackLinearSolver<T>(comm));
#endif


#ifdef LIBMESH_HAVE_PETSC
    case PETSC_SOLVERS:
      return UniquePtr<LinearSolver<T> >(new PetscLinearSolver<T>(comm));
#endif


#ifdef LIBMESH_HAVE_TRILINOS
    case TRILINOS_SOLVERS:
      return UniquePtr<LinearSolver<T> >(new AztecLinearSolver<T>(comm));
#endif


#ifdef LIBMESH_HAVE_EIGEN
    case EIGEN_SOLVERS:
      return UniquePtr<LinearSolver<T> >(new EigenSparseLinearSolver<T>(comm));
#endif

    default:
      libmesh_error_msg("ERROR:  Unrecognized solver package: " << solver_package);
    }

  return UniquePtr<LinearSolver<T> >();
}

template <typename T>
PreconditionerType
LinearSolver<T>::preconditioner_type () const
{
  if(_preconditioner)
    return _preconditioner->type();

  return _preconditioner_type;
}

template <typename T>
void
LinearSolver<T>::set_preconditioner_type (const PreconditionerType pct)
{
  if(_preconditioner)
    _preconditioner->set_type(pct);
  else
    _preconditioner_type = pct;
}

template <typename T>
void
LinearSolver<T>::attach_preconditioner(Preconditioner<T> * preconditioner)
{
  if (this->_is_initialized)
    libmesh_error_msg("Preconditioner must be attached before the solver is initialized!");

  _preconditioner_type = SHELL_PRECOND;
  _preconditioner = preconditioner;
}

template <typename T>
void
LinearSolver<T>::reuse_preconditioner(bool reuse_flag)
{
  same_preconditioner = reuse_flag;
}

template <typename T>
void
LinearSolver<T>::restrict_solve_to(const std::vector<unsigned int>* const dofs,
                                   const SubsetSolveMode /*subset_solve_mode*/)
{
  if(dofs!=NULL)
    {
      libmesh_not_implemented();
    }
}


template <typename T>
std::pair<unsigned int, Real> LinearSolver<T>::adjoint_solve (SparseMatrix<T> & mat,
                                                              NumericVector<T>& sol,
                                                              NumericVector<T>& rhs,
                                                              const double tol,
                                                              const unsigned int n_iter)
{
  // Log how long the linear solve takes.
  START_LOG("adjoint_solve()", "LinearSolver");

  // Take the discrete adjoint
  mat.close();
  mat.get_transpose(mat);

  // Call the solve function for the relevant linear algebra library and
  // solve the transpose matrix
  const std::pair<unsigned int, Real> totalrval =  this->solve (mat, sol, rhs, tol, n_iter);

  // Now transpose back and restore the original matrix
  // by taking the discrete adjoint
  mat.get_transpose(mat);

  // Stop logging the nonlinear solve
  STOP_LOG("adjoint_solve()", "LinearSolver");

  return totalrval;

}

template <typename T>
void LinearSolver<T>::print_converged_reason() const
{
  LinearConvergenceReason reason = this->get_converged_reason();
  libMesh::out << "Linear solver convergence/divergence reason: " << Utility::enum_to_string(reason) << std::endl;
}

//------------------------------------------------------------------
// Explicit instantiations
template class LinearSolver<Number>;



} // namespace libMesh
