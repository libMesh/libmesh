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



#include "libmesh/libmesh_common.h"

#if defined(LIBMESH_HAVE_SLEPC) && defined(LIBMESH_HAVE_PETSC)

// Local Includes
#include "libmesh/libmesh_logging.h"
#include "libmesh/petsc_matrix_base.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/slepc_eigen_solver.h"
#include "libmesh/shell_matrix.h"
#include "libmesh/enum_to_string.h"
#include "libmesh/solver_configuration.h"
#include "libmesh/enum_eigen_solver_type.h"
#include "libmesh/petsc_shell_matrix.h"

// PETSc 3.15 includes a non-release SLEPc 3.14.2 that has already
// deprecated STPrecondSetMatForPC but that hasn't upgraded its
// version number to let us switch to the non-deprecated version.  If
// we're in that situation we need to ignore the deprecated warning.
#if SLEPC_VERSION_LESS_THAN(3,15,0) && !SLEPC_VERSION_LESS_THAN(3,14,2)
#include "libmesh/ignore_warnings.h"
#endif

namespace libMesh
{



template <typename T>
SlepcEigenSolver<T>::SlepcEigenSolver (const Parallel::Communicator & comm_in) :
  EigenSolver<T>(comm_in),
  _initial_space(nullptr)
{
  this->_eigen_solver_type  = ARNOLDI;
  this->_eigen_problem_type = NHEP;
}



template <typename T>
SlepcEigenSolver<T>::~SlepcEigenSolver ()
{
  this->SlepcEigenSolver::clear ();
}



template <typename T>
void SlepcEigenSolver<T>::clear () noexcept
{
  if (this->initialized())
    {
      this->_is_initialized = false;

      PetscErrorCode ierr = LibMeshEPSDestroy(&_eps);
      if (ierr)
        libmesh_warning("Warning: EPSDestroy returned a non-zero error code which we ignored.");

      // SLEPc default eigenproblem solver
      this->_eigen_solver_type = KRYLOVSCHUR;
    }
}



template <typename T>
void SlepcEigenSolver<T>::init ()
{
  // Initialize the data structures if not done so already.
  if (!this->initialized())
    {
      this->_is_initialized = true;

      // Create the eigenproblem solver context
      LibmeshPetscCall(EPSCreate (this->comm().get(), &_eps));

      // Set user-specified  solver
      set_slepc_solver_type();
    }
}



template <typename T>
std::pair<unsigned int, unsigned int>
SlepcEigenSolver<T>::solve_standard (SparseMatrix<T> & matrix_A_in,
                                     int nev,                  // number of requested eigenpairs
                                     int ncv,                  // number of basis vectors
                                     const double tol,         // solver tolerance
                                     const unsigned int m_its) // maximum number of iterations
{
  LOG_SCOPE("solve_standard()", "SlepcEigenSolver");

  this->clear ();

  this->init ();

  // Make sure the SparseMatrix passed in is really a PetscMatrix
  auto * const matrix_A = cast_ptr<PetscMatrixBase<T> *>(&matrix_A_in);

  libmesh_error_msg_if(!matrix_A, "Error: input matrix to solve_standard() must be a PetscMatrixBase.");

  // Close the matrix and vectors in case this wasn't already done.
  if (this->_close_matrix_before_solve)
    matrix_A->close ();

  return _solve_standard_helper(matrix_A->mat(), nullptr, nev, ncv, tol, m_its);
}


template <typename T>
std::pair<unsigned int, unsigned int>
SlepcEigenSolver<T>::solve_standard (ShellMatrix<T> & shell_matrix,
                                     int nev,                  // number of requested eigenpairs
                                     int ncv,                  // number of basis vectors
                                     const double tol,         // solver tolerance
                                     const unsigned int m_its) // maximum number of iterations
{
  this->clear ();

  this->init ();

  // Prepare the matrix.  Note that the const_cast is only necessary
  // because PETSc does not accept a const void *.  Inside the member
  // function _petsc_shell_matrix() below, the pointer is casted back
  // to a const ShellMatrix<T> *.
  Mat mat;
  LibmeshPetscCall(MatCreateShell(this->comm().get(),
                                  shell_matrix.m(), // Specify the number of local rows
                                  shell_matrix.n(), // Specify the number of local columns
                                  PETSC_DETERMINE,
                                  PETSC_DETERMINE,
                                  const_cast<void *>(static_cast<const void *>(&shell_matrix)),
                                  &mat));

  LibmeshPetscCall(MatShellSetOperation(mat,MATOP_MULT,reinterpret_cast<void(*)(void)>(_petsc_shell_matrix_mult)));
  LibmeshPetscCall(MatShellSetOperation(mat,MATOP_GET_DIAGONAL,reinterpret_cast<void(*)(void)>(_petsc_shell_matrix_get_diagonal)));

  return _solve_standard_helper(mat, nullptr, nev, ncv, tol, m_its);
}

template <typename T>
std::pair<unsigned int, unsigned int>
SlepcEigenSolver<T>::solve_standard (ShellMatrix<T> & shell_matrix,
                                     SparseMatrix<T> & precond_in,
                                     int nev,                  // number of requested eigenpairs
                                     int ncv,                  // number of basis vectors
                                     const double tol,         // solver tolerance
                                     const unsigned int m_its) // maximum number of iterations
{
  this->clear ();

  this->init ();

  // Make sure the SparseMatrix passed in is really a PetscMatrix
  auto * const precond = cast_ptr<PetscMatrixBase<T> *>(&precond_in);

  libmesh_error_msg_if(!precond, "Error: input preconditioning matrix to solve_standard() must be a PetscMatrixBase.");

  auto * const matrix = cast_ptr<PetscShellMatrix<T> *> (&shell_matrix);

  libmesh_error_msg_if(!matrix, "Error: input operator matrix to solve_standard() must be a PetscShellMatrix.");

  return _solve_standard_helper(matrix->mat(), precond->mat(), nev, ncv, tol, m_its);
}

template <typename T>
std::pair<unsigned int, unsigned int>
SlepcEigenSolver<T>::solve_standard (ShellMatrix<T> & shell_matrix,
                                     ShellMatrix<T> & precond_in,
                                     int nev,                  // number of requested eigenpairs
                                     int ncv,                  // number of basis vectors
                                     const double tol,         // solver tolerance
                                     const unsigned int m_its) // maximum number of iterations
{
  this->clear ();

  this->init ();

  // Make sure the SparseMatrix passed in is really a PetscMatrix
  auto * const precond = cast_ptr<PetscShellMatrix<T> *>(&precond_in);

  libmesh_error_msg_if(!precond, "Error: input preconditioning matrix to solve_standard() must be a PetscShellMatrix.");

  auto * const matrix = cast_ptr<PetscShellMatrix<T> *> (&shell_matrix);

  libmesh_error_msg_if(!matrix, "Error: input operator matrix to solve_standard() must be a PetscShellMatrix.");

  return _solve_standard_helper(matrix->mat(), precond->mat(), nev, ncv, tol, m_its);
}



template <typename T>
std::pair<unsigned int, unsigned int>
SlepcEigenSolver<T>::_solve_standard_helper(Mat mat,
                                            Mat precond,
                                            int nev,                  // number of requested eigenpairs
                                            int ncv,                  // number of basis vectors
                                            const double tol,         // solver tolerance
                                            const unsigned int m_its) // maximum number of iterations
{
  LOG_SCOPE("solve_standard()", "SlepcEigenSolver");

  // Set operators.
  LibmeshPetscCall(EPSSetOperators (_eps, mat, LIBMESH_PETSC_NULLPTR));

  return this->_solve_helper(precond, nev, ncv, tol, m_its);
}



template <typename T>
std::pair<unsigned int, unsigned int>
SlepcEigenSolver<T>::solve_generalized (SparseMatrix<T> & matrix_A_in,
                                        SparseMatrix<T> & matrix_B_in,
                                        int nev,                  // number of requested eigenpairs
                                        int ncv,                  // number of basis vectors
                                        const double tol,         // solver tolerance
                                        const unsigned int m_its) // maximum number of iterations
{
  this->clear ();

  this->init ();

  // Make sure the data passed in are really of Petsc types
  auto * const matrix_A = cast_ptr<PetscMatrixBase<T> *>(&matrix_A_in);
  auto * const matrix_B = cast_ptr<PetscMatrixBase<T> *>(&matrix_B_in);

  libmesh_error_msg_if(!matrix_A || !matrix_B,
                       "Error: inputs to solve_generalized() must be of type PetscMatrixBase.");

  // Close the matrix and vectors in case this wasn't already done.
  if (this->_close_matrix_before_solve)
    {
      matrix_A->close ();
      matrix_B->close ();
    }

  return _solve_generalized_helper (matrix_A->mat(), matrix_B->mat(), nullptr, nev, ncv, tol, m_its);
}

template <typename T>
std::pair<unsigned int, unsigned int>
SlepcEigenSolver<T>::solve_generalized (ShellMatrix<T> & shell_matrix_A,
                                        SparseMatrix<T> & matrix_B_in,
                                        int nev,                  // number of requested eigenpairs
                                        int ncv,                  // number of basis vectors
                                        const double tol,         // solver tolerance
                                        const unsigned int m_its) // maximum number of iterations
{
  this->clear();

  this->init ();

  // Prepare the matrix. Note that the const_cast is only necessary
  // because PETSc does not accept a const void *.  Inside the member
  // function _petsc_shell_matrix() below, the pointer is casted back
  // to a const ShellMatrix<T> *.
  Mat mat_A;
  LibmeshPetscCall(MatCreateShell(this->comm().get(),
                                  shell_matrix_A.m(), // Specify the number of local rows
                                  shell_matrix_A.n(), // Specify the number of local columns
                                  PETSC_DETERMINE,
                                  PETSC_DETERMINE,
                                  const_cast<void *>(static_cast<const void *>(&shell_matrix_A)),
                                  &mat_A));

  auto * const matrix_B = cast_ptr<PetscMatrixBase<T> *>(&matrix_B_in);

  libmesh_error_msg_if(!matrix_B, "Error: inputs to solve_generalized() must be of type PetscMatrixBase.");

  // Close the matrix and vectors in case this wasn't already done.
  if (this->_close_matrix_before_solve)
    matrix_B->close ();

  LibmeshPetscCall(MatShellSetOperation(mat_A,MATOP_MULT,reinterpret_cast<void(*)(void)>(_petsc_shell_matrix_mult)));
  LibmeshPetscCall(MatShellSetOperation(mat_A,MATOP_GET_DIAGONAL,reinterpret_cast<void(*)(void)>(_petsc_shell_matrix_get_diagonal)));

  return _solve_generalized_helper (mat_A, matrix_B->mat(), nullptr, nev, ncv, tol, m_its);
}

template <typename T>
std::pair<unsigned int, unsigned int>
SlepcEigenSolver<T>::solve_generalized (SparseMatrix<T> & matrix_A_in,
                                        ShellMatrix<T> & shell_matrix_B,
                                        int nev,                  // number of requested eigenpairs
                                        int ncv,                  // number of basis vectors
                                        const double tol,         // solver tolerance
                                        const unsigned int m_its) // maximum number of iterations
{
  this->clear();

  this->init ();

  auto * const matrix_A = cast_ptr<PetscMatrixBase<T> *>(&matrix_A_in);

  libmesh_error_msg_if(!matrix_A, "Error: inputs to solve_generalized() must be of type PetscMatrixBase.");

  // Close the matrix and vectors in case this wasn't already done.
  if (this->_close_matrix_before_solve)
    matrix_A->close ();

  // Prepare the matrix.  Note that the const_cast is only necessary
  // because PETSc does not accept a const void *.  Inside the member
  // function _petsc_shell_matrix() below, the pointer is casted back
  // to a const ShellMatrix<T> *.
  Mat mat_B;
  LibmeshPetscCall(MatCreateShell(this->comm().get(),
                                  shell_matrix_B.m(), // Specify the number of local rows
                                  shell_matrix_B.n(), // Specify the number of local columns
                                  PETSC_DETERMINE,
                                  PETSC_DETERMINE,
                                  const_cast<void *>(static_cast<const void *>(&shell_matrix_B)),
                                  &mat_B));

  LibmeshPetscCall(MatShellSetOperation(mat_B,MATOP_MULT,reinterpret_cast<void(*)(void)>(_petsc_shell_matrix_mult)));
  LibmeshPetscCall(MatShellSetOperation(mat_B,MATOP_GET_DIAGONAL,reinterpret_cast<void(*)(void)>(_petsc_shell_matrix_get_diagonal)));

  return _solve_generalized_helper (matrix_A->mat(), mat_B, nullptr, nev, ncv, tol, m_its);
}

template <typename T>
std::pair<unsigned int, unsigned int>
SlepcEigenSolver<T>::solve_generalized (ShellMatrix<T> & shell_matrix_A,
                                        ShellMatrix<T> & shell_matrix_B,
                                        int nev,                  // number of requested eigenpairs
                                        int ncv,                  // number of basis vectors
                                        const double tol,         // solver tolerance
                                        const unsigned int m_its) // maximum number of iterations
{
  this->clear();

  this->init ();

  // Prepare the matrices.  Note that the const_casts are only
  // necessary because PETSc does not accept a const void *.  Inside
  // the member function _petsc_shell_matrix() below, the pointer is
  // casted back to a const ShellMatrix<T> *.
  Mat mat_A;
  LibmeshPetscCall(MatCreateShell(this->comm().get(),
                                  shell_matrix_A.m(), // Specify the number of local rows
                                  shell_matrix_A.n(), // Specify the number of local columns
                                  PETSC_DETERMINE,
                                  PETSC_DETERMINE,
                                  const_cast<void *>(static_cast<const void *>(&shell_matrix_A)),
                                  &mat_A));

  Mat mat_B;
  LibmeshPetscCall(MatCreateShell(this->comm().get(),
                                  shell_matrix_B.m(), // Specify the number of local rows
                                  shell_matrix_B.n(), // Specify the number of local columns
                                  PETSC_DETERMINE,
                                  PETSC_DETERMINE,
                                  const_cast<void *>(static_cast<const void *>(&shell_matrix_B)),
                                  &mat_B));

  LibmeshPetscCall(MatShellSetOperation(mat_A,MATOP_MULT,reinterpret_cast<void(*)(void)>(_petsc_shell_matrix_mult)));
  LibmeshPetscCall(MatShellSetOperation(mat_A,MATOP_GET_DIAGONAL,reinterpret_cast<void(*)(void)>(_petsc_shell_matrix_get_diagonal)));

  LibmeshPetscCall(MatShellSetOperation(mat_B,MATOP_MULT,reinterpret_cast<void(*)(void)>(_petsc_shell_matrix_mult)));
  LibmeshPetscCall(MatShellSetOperation(mat_B,MATOP_GET_DIAGONAL,reinterpret_cast<void(*)(void)>(_petsc_shell_matrix_get_diagonal)));

  return _solve_generalized_helper (mat_A, mat_B, nullptr, nev, ncv, tol, m_its);
}

template <typename T>
std::pair<unsigned int, unsigned int>
SlepcEigenSolver<T>::solve_generalized (ShellMatrix<T> & shell_matrix_A,
                                        ShellMatrix<T> & shell_matrix_B,
                                        SparseMatrix<T> & precond_in,
                                        int nev,                  // number of requested eigenpairs
                                        int ncv,                  // number of basis vectors
                                        const double tol,         // solver tolerance
                                        const unsigned int m_its) // maximum number of iterations
{
  this->clear();

  this->init ();

  // Make sure the SparseMatrix passed in is really a PetscMatrix
  auto * const precond = cast_ptr<PetscMatrixBase<T> *>(&precond_in);

  libmesh_error_msg_if(!precond, "Error: input preconditioning matrix to solve_generalized() must be of type PetscMatrixBase.");

  auto * const matrix_A = cast_ptr<PetscShellMatrix<T> *> (&shell_matrix_A);

  libmesh_error_msg_if(!matrix_A, "Error: input operator A to solve_generalized() must be of type PetscShellMatrix.");

  auto * const matrix_B = cast_ptr<PetscShellMatrix<T> *> (&shell_matrix_B);

  libmesh_error_msg_if(!matrix_B, "Error: input operator B to solve_generalized() must be of type PetscShellMatrix.");

  return _solve_generalized_helper (matrix_A->mat(), matrix_B->mat(), precond->mat(), nev, ncv, tol, m_its);
}

template <typename T>
std::pair<unsigned int, unsigned int>
SlepcEigenSolver<T>::solve_generalized (ShellMatrix<T> & shell_matrix_A,
                                        ShellMatrix<T> & shell_matrix_B,
                                        ShellMatrix<T> & precond_in,
                                        int nev,                  // number of requested eigenpairs
                                        int ncv,                  // number of basis vectors
                                        const double tol,         // solver tolerance
                                        const unsigned int m_its) // maximum number of iterations
{
  this->clear();

  this->init ();

  // Make sure the ShellMatrix passed in is really a PetscShellMatrix
  auto * const precond = cast_ptr<PetscShellMatrix<T> *>(&precond_in);

  libmesh_error_msg_if(!precond, "Error: input preconditioning matrix to solve_generalized() must be of type PetscShellMatrix.");

  auto * const matrix_A = cast_ptr<PetscShellMatrix<T> *> (&shell_matrix_A);

  libmesh_error_msg_if(!matrix_A, "Error: input operator A to solve_generalized() must be of type PetscShellMatrix.");

  auto * const matrix_B = cast_ptr<PetscShellMatrix<T> *> (&shell_matrix_B);

  libmesh_error_msg_if(!matrix_B, "Error: input operator B to solve_generalized() must be of type PetscShellMatrix.");

  return _solve_generalized_helper (matrix_A->mat(), matrix_B->mat(), precond->mat(), nev, ncv, tol, m_its);
}

template <typename T>
std::pair<unsigned int, unsigned int>
SlepcEigenSolver<T>::_solve_generalized_helper (Mat mat_A,
                                                Mat mat_B,
                                                Mat precond,
                                                int nev,                  // number of requested eigenpairs
                                                int ncv,                  // number of basis vectors
                                                const double tol,         // solver tolerance
                                                const unsigned int m_its) // maximum number of iterations
{
  LOG_SCOPE("solve_generalized()", "SlepcEigenSolver");

  // Set operators.
  LibmeshPetscCall(EPSSetOperators (_eps, mat_A, mat_B));

  return this->_solve_helper(precond, nev, ncv, tol, m_its);
}



template <typename T>
std::pair<unsigned int, unsigned int>
SlepcEigenSolver<T>::_solve_helper(Mat precond,
                                   int nev,                  // number of requested eigenpairs
                                   int ncv,                  // number of basis vectors
                                   const double tol,         // solver tolerance
                                   const unsigned int m_its) // maximum number of iterations
{
  // converged eigen pairs and number of iterations
  PetscInt nconv=0;
  PetscInt its=0;
  ST st=nullptr;

  //set the problem type and the position of the spectrum
  set_slepc_problem_type();
  set_slepc_position_of_spectrum();

  // Set eigenvalues to be computed.
#if SLEPC_VERSION_LESS_THAN(3,0,0)
  LibmeshPetscCall(EPSSetDimensions (_eps, nev, ncv));
#else
  LibmeshPetscCall(EPSSetDimensions (_eps, nev, ncv, PETSC_DECIDE));
#endif
  // Set the tolerance and maximum iterations.
  LibmeshPetscCall(EPSSetTolerances (_eps, tol, m_its));

  // Set runtime options, e.g.,
  //      -eps_type <type>, -eps_nev <nev>, -eps_ncv <ncv>
  // Similar to PETSc, these options will override those specified
  // above as long as EPSSetFromOptions() is called _after_ any
  // other customization routines.
  LibmeshPetscCall(EPSSetFromOptions (_eps));

  // Set a preconditioning matrix to ST
  if (precond) {
    LibmeshPetscCall(EPSGetST(_eps,&st));
#if SLEPC_VERSION_LESS_THAN(3,15,0)
    LibmeshPetscCall(STPrecondSetMatForPC(st, precond));
#else
    LibmeshPetscCall(STSetPreconditionerMat(st, precond));
#endif
  }

  // If the SolverConfiguration object is provided, use it to override
  // solver options.
  if (this->_solver_configuration)
    {
      this->_solver_configuration->configure_solver();
    }

  // If an initial space is provided, let us attach it to EPS
  if (_initial_space) {
    // Get a handle for the underlying Vec.
    Vec initial_vector = _initial_space->vec();

    LibmeshPetscCall(EPSSetInitialSpace(_eps, 1, &initial_vector));
  }

  // Solve the eigenproblem.
  LibmeshPetscCall(EPSSolve (_eps));

  // Get the number of iterations.
  LibmeshPetscCall(EPSGetIterationNumber (_eps, &its));

  // Get number of converged eigenpairs.
  LibmeshPetscCall(EPSGetConverged(_eps,&nconv));

  // return the number of converged eigenpairs
  // and the number of iterations
  return std::make_pair(nconv, its);
}


template <typename T>
void SlepcEigenSolver<T>::print_eigenvalues() const
{
  // converged eigen pairs and number of iterations
  PetscInt nconv=0;

  // The relative error.
  PetscReal error, re, im;

  // Pointer to vectors of the real parts, imaginary parts.
  PetscScalar kr, ki;

  // Get number of converged eigenpairs.
  LibmeshPetscCall(EPSGetConverged(_eps,&nconv));

  // Display eigenvalues and relative errors.
  LibmeshPetscCall(PetscPrintf(this->comm().get(),
                               "           k           ||Ax-kx||/|kx|\n"
                               "   ----------------- -----------------\n" ));

  for (PetscInt i=0; i<nconv; i++ )
    {
      LibmeshPetscCall(EPSGetEigenpair(_eps, i, &kr, &ki, LIBMESH_PETSC_NULLPTR,
                                       LIBMESH_PETSC_NULLPTR));

#if SLEPC_VERSION_LESS_THAN(3,6,0)
      LibmeshPetscCall(EPSComputeRelativeError(_eps, i, &error));
#else
      LibmeshPetscCall(EPSComputeError(_eps, i, EPS_ERROR_RELATIVE, &error));
#endif

#ifdef LIBMESH_USE_COMPLEX_NUMBERS
      re = PetscRealPart(kr);
      im = PetscImaginaryPart(kr);
#else
      re = kr;
      im = ki;
#endif

      if (im != .0)
        LibmeshPetscCall(PetscPrintf(this->comm().get()," %9f%+9f i %12f\n", re, im, error));
      else
        LibmeshPetscCall(PetscPrintf(this->comm().get(),"   %12f       %12f\n", re, error));
    }

  LibmeshPetscCall(PetscPrintf(this->comm().get(),"\n" ));
}


template <typename T>
void SlepcEigenSolver<T>::set_slepc_solver_type()
{
  switch (this->_eigen_solver_type)
    {
    case POWER:
      LibmeshPetscCall(EPSSetType (_eps, EPSPOWER)); return;
    case SUBSPACE:
      LibmeshPetscCall(EPSSetType (_eps, EPSSUBSPACE)); return;
    case LAPACK:
      LibmeshPetscCall(EPSSetType (_eps, EPSLAPACK)); return;
    case ARNOLDI:
      LibmeshPetscCall(EPSSetType (_eps, EPSARNOLDI)); return;
    case LANCZOS:
      LibmeshPetscCall(EPSSetType (_eps, EPSLANCZOS)); return;
    case KRYLOVSCHUR:
      LibmeshPetscCall(EPSSetType (_eps, EPSKRYLOVSCHUR)); return;
      // case ARPACK:
      // LibmeshPetscCall(EPSSetType (_eps, (char *) EPSARPACK)); return;

    default:
      libMesh::err << "ERROR:  Unsupported SLEPc Eigen Solver: "
                   << Utility::enum_to_string(this->_eigen_solver_type) << std::endl
                   << "Continuing with SLEPc defaults" << std::endl;
    }
}




template <typename T>
void SlepcEigenSolver<T>:: set_slepc_problem_type()
{
  switch (this->_eigen_problem_type)
    {
    case NHEP:
      LibmeshPetscCall(EPSSetProblemType (_eps, EPS_NHEP)); return;
    case GNHEP:
      LibmeshPetscCall(EPSSetProblemType (_eps, EPS_GNHEP)); return;
    case HEP:
      LibmeshPetscCall(EPSSetProblemType (_eps, EPS_HEP)); return;
    case GHEP:
      LibmeshPetscCall(EPSSetProblemType (_eps, EPS_GHEP)); return;
#if !SLEPC_VERSION_LESS_THAN(3,3,0)
      // EPS_GHIEP added in 3.3.0
    case GHIEP:
      LibmeshPetscCall(EPSSetProblemType (_eps, EPS_GHIEP)); return;
#endif

    default:
      libMesh::err << "ERROR:  Unsupported SLEPc Eigen Problem: "
                   << this->_eigen_problem_type        << std::endl
                   << "Continuing with SLEPc defaults" << std::endl;
    }
}



template <typename T>
void SlepcEigenSolver<T>:: set_slepc_position_of_spectrum()
{
  switch (this->_position_of_spectrum)
    {
    case LARGEST_MAGNITUDE:
      {
        LibmeshPetscCall(EPSSetWhichEigenpairs (_eps, EPS_LARGEST_MAGNITUDE));
        return;
      }
    case SMALLEST_MAGNITUDE:
      {
        LibmeshPetscCall(EPSSetWhichEigenpairs (_eps, EPS_SMALLEST_MAGNITUDE));
        return;
      }
    case LARGEST_REAL:
      {
        LibmeshPetscCall(EPSSetWhichEigenpairs (_eps, EPS_LARGEST_REAL));
        return;
      }
    case SMALLEST_REAL:
      {
        LibmeshPetscCall(EPSSetWhichEigenpairs (_eps, EPS_SMALLEST_REAL));
        return;
      }
    case LARGEST_IMAGINARY:
      {
        LibmeshPetscCall(EPSSetWhichEigenpairs (_eps, EPS_LARGEST_IMAGINARY));
        return;
      }
    case SMALLEST_IMAGINARY:
      {
        LibmeshPetscCall(EPSSetWhichEigenpairs (_eps, EPS_SMALLEST_IMAGINARY));
        return;
      }

      // The EPS_TARGET_XXX enums were added in SLEPc 3.1
#if !SLEPC_VERSION_LESS_THAN(3,1,0)
    case TARGET_MAGNITUDE:
      {
        LibmeshPetscCall(EPSSetTarget(_eps, PS(this->_target_val)));
        LibmeshPetscCall(EPSSetWhichEigenpairs (_eps, EPS_TARGET_MAGNITUDE));
        return;
      }
    case TARGET_REAL:
      {
        LibmeshPetscCall(EPSSetTarget(_eps, PS(this->_target_val)));
        LibmeshPetscCall(EPSSetWhichEigenpairs (_eps, EPS_TARGET_REAL));
        return;
      }
    case TARGET_IMAGINARY:
      {
        LibmeshPetscCall(EPSSetTarget(_eps, PS(this->_target_val)));
        LibmeshPetscCall(EPSSetWhichEigenpairs (_eps, EPS_TARGET_IMAGINARY));
        return;
      }
#endif

    default:
      libmesh_error_msg("ERROR:  Unsupported SLEPc position of spectrum: " << this->_position_of_spectrum);
    }
}






template <typename T>
std::pair<Real, Real> SlepcEigenSolver<T>::get_eigenpair(dof_id_type i,
                                                         NumericVector<T> & solution_in)
{
  PetscReal re, im;

  // Make sure the NumericVector passed in is really a PetscVector
  PetscVector<T> * solution = dynamic_cast<PetscVector<T> *>(&solution_in);

  libmesh_error_msg_if(!solution, "Error getting eigenvector: input vector must be a PetscVector.");

  // real and imaginary part of the ith eigenvalue.
  PetscScalar kr, ki;

  solution->close();

  LibmeshPetscCall(EPSGetEigenpair(_eps, i, &kr, &ki, solution->vec(),
                                   LIBMESH_PETSC_NULLPTR));

#ifdef LIBMESH_USE_COMPLEX_NUMBERS
  re = PetscRealPart(kr);
  im = PetscImaginaryPart(kr);
#else
  re = kr;
  im = ki;
#endif

  return std::make_pair(re, im);
}


template <typename T>
std::pair<Real, Real> SlepcEigenSolver<T>::get_eigenvalue(dof_id_type i)
{
  PetscReal re, im;

  // real and imaginary part of the ith eigenvalue.
  PetscScalar kr, ki;

  LibmeshPetscCall(EPSGetEigenvalue(_eps, i, &kr, &ki));

#ifdef LIBMESH_USE_COMPLEX_NUMBERS
  re = PetscRealPart(kr);
  im = PetscImaginaryPart(kr);
#else
  re = kr;
  im = ki;
#endif

  return std::make_pair(re, im);
}


template <typename T>
Real SlepcEigenSolver<T>::get_relative_error(unsigned int i)
{
  PetscReal error;

#if SLEPC_VERSION_LESS_THAN(3,6,0)
  LibmeshPetscCall(EPSComputeRelativeError(_eps, i, &error));
#else
  LibmeshPetscCall(EPSComputeError(_eps, i, EPS_ERROR_RELATIVE, &error));
#endif

  return error;
}


template <typename T>
void SlepcEigenSolver<T>::attach_deflation_space(NumericVector<T> & deflation_vector_in)
{
  this->init();

  // Make sure the input vector is actually a PetscVector
  PetscVector<T> * deflation_vector_petsc_vec =
    dynamic_cast<PetscVector<T> *>(&deflation_vector_in);

  libmesh_error_msg_if(!deflation_vector_petsc_vec, "Error attaching deflation space: input vector must be a PetscVector.");

  // Get a handle for the underlying Vec.
  Vec deflation_vector = deflation_vector_petsc_vec->vec();

#if SLEPC_VERSION_LESS_THAN(3,1,0)
  LibmeshPetscCall(EPSAttachDeflationSpace(_eps, 1, &deflation_vector, PETSC_FALSE));
#else
  LibmeshPetscCall(EPSSetDeflationSpace(_eps, 1, &deflation_vector));
#endif
}

template <typename T>
void SlepcEigenSolver<T>::set_initial_space(NumericVector<T> & initial_space_in)
{
#if SLEPC_VERSION_LESS_THAN(3,1,0)
  libmesh_error_msg("SLEPc 3.1 is required to call EigenSolver::set_initial_space()");
#else
  // Make sure the input vector (which is still owned by caller) is
  // actually a PetscVector
  _initial_space = cast_ptr<PetscVector<T> *>(&initial_space_in);
#endif
}

template <typename T>
PetscErrorCode SlepcEigenSolver<T>::_petsc_shell_matrix_mult(Mat mat, Vec arg, Vec dest)
{
  PetscFunctionBegin;

  // Get the matrix context.
  void * ctx;
  LibmeshPetscCallQ(MatShellGetContext(mat,&ctx));

  // Get user shell matrix object.
  const ShellMatrix<T> & shell_matrix = *static_cast<const ShellMatrix<T> *>(ctx);

  // Make \p NumericVector instances around the vectors.
  PetscVector<T> arg_global(arg,   shell_matrix.comm());
  PetscVector<T> dest_global(dest, shell_matrix.comm());

  // Call the user function.
  shell_matrix.vector_mult(dest_global,arg_global);

  PetscFunctionReturn(LIBMESH_PETSC_SUCCESS);
}

template <typename T>
PetscErrorCode SlepcEigenSolver<T>::_petsc_shell_matrix_get_diagonal(Mat mat, Vec dest)
{
  PetscFunctionBegin;

  // Get the matrix context.
  void * ctx;
  LibmeshPetscCallQ(MatShellGetContext(mat,&ctx));

  // Get user shell matrix object.
  const ShellMatrix<T> & shell_matrix = *static_cast<const ShellMatrix<T> *>(ctx);

  // Make \p NumericVector instances around the vector.
  PetscVector<T> dest_global(dest, shell_matrix.comm());

  // Call the user function.
  shell_matrix.get_diagonal(dest_global);

  PetscFunctionReturn(LIBMESH_PETSC_SUCCESS);
}

//------------------------------------------------------------------
// Explicit instantiations
template class LIBMESH_EXPORT SlepcEigenSolver<Number>;

} // namespace libMesh


// In case we do unity builds someday
#if SLEPC_VERSION_LESS_THAN(3,15,0) && !SLEPC_VERSION_LESS_THAN(3,14,2)
#include "libmesh/restore_warnings.h"
#endif


#endif // #ifdef LIBMESH_HAVE_SLEPC
