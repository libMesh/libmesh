// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh_logging.h"
#include "auto_ptr.h"
#include "linear_solver.h"
#include "laspack_linear_solver.h"
#include "petsc_linear_solver.h"
#include "trilinos_aztec_linear_solver.h"
#include "preconditioner.h"
#include "sparse_matrix.h"

namespace libMesh
{

//------------------------------------------------------------------
// LinearSolver members
template <typename T>
AutoPtr<LinearSolver<T> >
LinearSolver<T>::build(const SolverPackage solver_package)
{
  // Build the appropriate solver
  switch (solver_package)
    {


#ifdef LIBMESH_HAVE_LASPACK
    case LASPACK_SOLVERS:
      {
	AutoPtr<LinearSolver<T> > ap(new LaspackLinearSolver<T>);
	return ap;
      }
#endif


#ifdef LIBMESH_HAVE_PETSC
    case PETSC_SOLVERS:
      {
	AutoPtr<LinearSolver<T> > ap(new PetscLinearSolver<T>);
	return ap;
      }
#endif


#ifdef LIBMESH_HAVE_TRILINOS
    case TRILINOS_SOLVERS:
      {
	AutoPtr<LinearSolver<T> > ap(new AztecLinearSolver<T>);
	return ap;
      }
#endif

    default:
      libMesh::err << "ERROR:  Unrecognized solver package: "
		    << solver_package
		    << std::endl;
      libmesh_error();
    }

  AutoPtr<LinearSolver<T> > ap(NULL);
  return ap;
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
  if(this->_is_initialized)
  {
    libMesh::err<<"Preconditioner must be attached before the solver is initialized!"<<std::endl;
    libmesh_error();
  }

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


//------------------------------------------------------------------
// Explicit instantiations
template class LinearSolver<Number>;



} // namespace libMesh



