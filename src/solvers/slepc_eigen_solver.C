// The libMesh Finite Element Library.
// Copyright (C) 2002-2024 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/petsc_matrix.h"
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

  PetscErrorCode ierr = LIBMESH_PETSC_SUCCESS;

  // Initialize the data structures if not done so already.
  if (!this->initialized())
    {
      this->_is_initialized = true;

      // Create the eigenproblem solver context
      ierr = EPSCreate (this->comm().get(), &_eps);
      LIBMESH_CHKERR(ierr);

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
  PetscMatrix<T> * matrix_A = dynamic_cast<PetscMatrix<T> *>(&matrix_A_in);

  libmesh_error_msg_if(!matrix_A, "Error: input matrix to solve_standard() must be a PetscMatrix.");

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

  PetscErrorCode ierr = LIBMESH_PETSC_SUCCESS;

  // Prepare the matrix.  Note that the const_cast is only necessary
  // because PETSc does not accept a const void *.  Inside the member
  // function _petsc_shell_matrix() below, the pointer is casted back
  // to a const ShellMatrix<T> *.
  Mat mat;
  ierr = MatCreateShell(this->comm().get(),
                        shell_matrix.m(), // Specify the number of local rows
                        shell_matrix.n(), // Specify the number of local columns
                        PETSC_DETERMINE,
                        PETSC_DETERMINE,
                        const_cast<void *>(static_cast<const void *>(&shell_matrix)),
                        &mat);
  LIBMESH_CHKERR(ierr);

  ierr = MatShellSetOperation(mat,MATOP_MULT,reinterpret_cast<void(*)(void)>(_petsc_shell_matrix_mult));
  LIBMESH_CHKERR(ierr);
  ierr = MatShellSetOperation(mat,MATOP_GET_DIAGONAL,reinterpret_cast<void(*)(void)>(_petsc_shell_matrix_get_diagonal));
  LIBMESH_CHKERR(ierr);

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
  PetscMatrix<T> * precond = dynamic_cast<PetscMatrix<T> *>(&precond_in);

  libmesh_error_msg_if(!precond, "Error: input preconditioning matrix to solve_standard() must be a PetscMatrix.");

  PetscShellMatrix<T> * matrix = dynamic_cast<PetscShellMatrix<T> *> (&shell_matrix);

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
  PetscShellMatrix<T> * precond = dynamic_cast<PetscShellMatrix<T> *>(&precond_in);

  libmesh_error_msg_if(!precond, "Error: input preconditioning matrix to solve_standard() must be a PetscShellMatrix.");

  PetscShellMatrix<T> * matrix = dynamic_cast<PetscShellMatrix<T> *> (&shell_matrix);

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
  PetscErrorCode ierr = EPSSetOperators (_eps, mat, LIBMESH_PETSC_NULLPTR);
  LIBMESH_CHKERR(ierr);

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
  PetscMatrix<T> * matrix_A = dynamic_cast<PetscMatrix<T> *>(&matrix_A_in);
  PetscMatrix<T> * matrix_B = dynamic_cast<PetscMatrix<T> *>(&matrix_B_in);

  libmesh_error_msg_if(!matrix_A || !matrix_B,
                       "Error: inputs to solve_generalized() must be of type PetscMatrix.");

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

  PetscErrorCode ierr = LIBMESH_PETSC_SUCCESS;

  // Prepare the matrix. Note that the const_cast is only necessary
  // because PETSc does not accept a const void *.  Inside the member
  // function _petsc_shell_matrix() below, the pointer is casted back
  // to a const ShellMatrix<T> *.
  Mat mat_A;
  ierr = MatCreateShell(this->comm().get(),
                        shell_matrix_A.m(), // Specify the number of local rows
                        shell_matrix_A.n(), // Specify the number of local columns
                        PETSC_DETERMINE,
                        PETSC_DETERMINE,
                        const_cast<void *>(static_cast<const void *>(&shell_matrix_A)),
                        &mat_A);
  LIBMESH_CHKERR(ierr);

  PetscMatrix<T> * matrix_B = dynamic_cast<PetscMatrix<T> *>(&matrix_B_in);

  libmesh_error_msg_if(!matrix_B, "Error: inputs to solve_generalized() must be of type PetscMatrix.");

  // Close the matrix and vectors in case this wasn't already done.
  if (this->_close_matrix_before_solve)
    matrix_B->close ();

  ierr = MatShellSetOperation(mat_A,MATOP_MULT,reinterpret_cast<void(*)(void)>(_petsc_shell_matrix_mult));
  LIBMESH_CHKERR(ierr);
  ierr = MatShellSetOperation(mat_A,MATOP_GET_DIAGONAL,reinterpret_cast<void(*)(void)>(_petsc_shell_matrix_get_diagonal));
  LIBMESH_CHKERR(ierr);

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

  PetscErrorCode ierr = LIBMESH_PETSC_SUCCESS;

  PetscMatrix<T> * matrix_A = dynamic_cast<PetscMatrix<T> *>(&matrix_A_in);

  libmesh_error_msg_if(!matrix_A, "Error: inputs to solve_generalized() must be of type PetscMatrix.");

  // Close the matrix and vectors in case this wasn't already done.
  if (this->_close_matrix_before_solve)
    matrix_A->close ();

  // Prepare the matrix.  Note that the const_cast is only necessary
  // because PETSc does not accept a const void *.  Inside the member
  // function _petsc_shell_matrix() below, the pointer is casted back
  // to a const ShellMatrix<T> *.
  Mat mat_B;
  ierr = MatCreateShell(this->comm().get(),
                        shell_matrix_B.m(), // Specify the number of local rows
                        shell_matrix_B.n(), // Specify the number of local columns
                        PETSC_DETERMINE,
                        PETSC_DETERMINE,
                        const_cast<void *>(static_cast<const void *>(&shell_matrix_B)),
                        &mat_B);
  LIBMESH_CHKERR(ierr);


  ierr = MatShellSetOperation(mat_B,MATOP_MULT,reinterpret_cast<void(*)(void)>(_petsc_shell_matrix_mult));
  LIBMESH_CHKERR(ierr);
  ierr = MatShellSetOperation(mat_B,MATOP_GET_DIAGONAL,reinterpret_cast<void(*)(void)>(_petsc_shell_matrix_get_diagonal));
  LIBMESH_CHKERR(ierr);

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

  PetscErrorCode ierr = LIBMESH_PETSC_SUCCESS;

  // Prepare the matrices.  Note that the const_casts are only
  // necessary because PETSc does not accept a const void *.  Inside
  // the member function _petsc_shell_matrix() below, the pointer is
  // casted back to a const ShellMatrix<T> *.
  Mat mat_A;
  ierr = MatCreateShell(this->comm().get(),
                        shell_matrix_A.m(), // Specify the number of local rows
                        shell_matrix_A.n(), // Specify the number of local columns
                        PETSC_DETERMINE,
                        PETSC_DETERMINE,
                        const_cast<void *>(static_cast<const void *>(&shell_matrix_A)),
                        &mat_A);
  LIBMESH_CHKERR(ierr);

  Mat mat_B;
  ierr = MatCreateShell(this->comm().get(),
                        shell_matrix_B.m(), // Specify the number of local rows
                        shell_matrix_B.n(), // Specify the number of local columns
                        PETSC_DETERMINE,
                        PETSC_DETERMINE,
                        const_cast<void *>(static_cast<const void *>(&shell_matrix_B)),
                        &mat_B);
  LIBMESH_CHKERR(ierr);

  ierr = MatShellSetOperation(mat_A,MATOP_MULT,reinterpret_cast<void(*)(void)>(_petsc_shell_matrix_mult));
  LIBMESH_CHKERR(ierr);
  ierr = MatShellSetOperation(mat_A,MATOP_GET_DIAGONAL,reinterpret_cast<void(*)(void)>(_petsc_shell_matrix_get_diagonal));
  LIBMESH_CHKERR(ierr);

  ierr = MatShellSetOperation(mat_B,MATOP_MULT,reinterpret_cast<void(*)(void)>(_petsc_shell_matrix_mult));
  LIBMESH_CHKERR(ierr);
  ierr = MatShellSetOperation(mat_B,MATOP_GET_DIAGONAL,reinterpret_cast<void(*)(void)>(_petsc_shell_matrix_get_diagonal));
  LIBMESH_CHKERR(ierr);

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
  PetscMatrix<T> * precond = dynamic_cast<PetscMatrix<T> *>(&precond_in);

  libmesh_error_msg_if(!precond, "Error: input preconditioning matrix to solve_generalized() must be of type PetscMatrix.");

  PetscShellMatrix<T> * matrix_A = dynamic_cast<PetscShellMatrix<T> *> (&shell_matrix_A);

  libmesh_error_msg_if(!matrix_A, "Error: input operator A to solve_generalized() must be of type PetscShellMatrix.");

  PetscShellMatrix<T> * matrix_B = dynamic_cast<PetscShellMatrix<T> *> (&shell_matrix_B);

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
  PetscShellMatrix<T> * precond = dynamic_cast<PetscShellMatrix<T> *>(&precond_in);

  libmesh_error_msg_if(!precond, "Error: input preconditioning matrix to solve_generalized() must be of type PetscShellMatrix.");

  PetscShellMatrix<T> * matrix_A = dynamic_cast<PetscShellMatrix<T> *> (&shell_matrix_A);

  libmesh_error_msg_if(!matrix_A, "Error: input operator A to solve_generalized() must be of type PetscShellMatrix.");

  PetscShellMatrix<T> * matrix_B = dynamic_cast<PetscShellMatrix<T> *> (&shell_matrix_B);

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
  PetscErrorCode ierr = EPSSetOperators (_eps, mat_A, mat_B);
  LIBMESH_CHKERR(ierr);

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
  PetscErrorCode ierr = LIBMESH_PETSC_SUCCESS;

  // converged eigen pairs and number of iterations
  PetscInt nconv=0;
  PetscInt its=0;
  ST st=nullptr;

  //set the problem type and the position of the spectrum
  set_slepc_problem_type();
  set_slepc_position_of_spectrum();

  // Set eigenvalues to be computed.
#if SLEPC_VERSION_LESS_THAN(3,0,0)
  ierr = EPSSetDimensions (_eps, nev, ncv);
#else
  ierr = EPSSetDimensions (_eps, nev, ncv, PETSC_DECIDE);
#endif
  LIBMESH_CHKERR(ierr);
  // Set the tolerance and maximum iterations.
  ierr = EPSSetTolerances (_eps, tol, m_its);
  LIBMESH_CHKERR(ierr);

  // Set runtime options, e.g.,
  //      -eps_type <type>, -eps_nev <nev>, -eps_ncv <ncv>
  // Similar to PETSc, these options will override those specified
  // above as long as EPSSetFromOptions() is called _after_ any
  // other customization routines.
  ierr = EPSSetFromOptions (_eps);
  LIBMESH_CHKERR(ierr);

  // Set a preconditioning matrix to ST
  if (precond) {
    ierr = EPSGetST(_eps,&st);LIBMESH_CHKERR(ierr);
#if SLEPC_VERSION_LESS_THAN(3,15,0)
    ierr = STPrecondSetMatForPC(st, precond);
#else
    ierr = STSetPreconditionerMat(st, precond);
#endif
    LIBMESH_CHKERR(ierr);
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

    ierr = EPSSetInitialSpace(_eps, 1, &initial_vector);
    LIBMESH_CHKERR(ierr);
  }

  // Solve the eigenproblem.
  ierr = EPSSolve (_eps);
  LIBMESH_CHKERR(ierr);

  // Get the number of iterations.
  ierr = EPSGetIterationNumber (_eps, &its);
  LIBMESH_CHKERR(ierr);

  // Get number of converged eigenpairs.
  ierr = EPSGetConverged(_eps,&nconv);
  LIBMESH_CHKERR(ierr);

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
  auto ierr = EPSGetConverged(_eps,&nconv);
  LIBMESH_CHKERR(ierr);

  // Display eigenvalues and relative errors.
  ierr = PetscPrintf(this->comm().get(),
                     "           k           ||Ax-kx||/|kx|\n"
                     "   ----------------- -----------------\n" );
  LIBMESH_CHKERR(ierr);

  for (PetscInt i=0; i<nconv; i++ )
    {
      ierr = EPSGetEigenpair(_eps, i, &kr, &ki, LIBMESH_PETSC_NULLPTR,
                             LIBMESH_PETSC_NULLPTR);
      LIBMESH_CHKERR(ierr);

#if SLEPC_VERSION_LESS_THAN(3,6,0)
      ierr = EPSComputeRelativeError(_eps, i, &error);
#else
      ierr = EPSComputeError(_eps, i, EPS_ERROR_RELATIVE, &error);
#endif
      LIBMESH_CHKERR(ierr);

#ifdef LIBMESH_USE_COMPLEX_NUMBERS
      re = PetscRealPart(kr);
      im = PetscImaginaryPart(kr);
#else
      re = kr;
      im = ki;
#endif

      if (im != .0)
        {
          ierr = PetscPrintf(this->comm().get()," %9f%+9f i %12f\n", re, im, error);
          LIBMESH_CHKERR(ierr);
        }
      else
        {
          ierr = PetscPrintf(this->comm().get(),"   %12f       %12f\n", re, error);
          LIBMESH_CHKERR(ierr);
        }
    }

  ierr = PetscPrintf(this->comm().get(),"\n" );
  LIBMESH_CHKERR(ierr);
}


template <typename T>
void SlepcEigenSolver<T>::set_slepc_solver_type()
{
  PetscErrorCode ierr = LIBMESH_PETSC_SUCCESS;

  switch (this->_eigen_solver_type)
    {
    case POWER:
      ierr = EPSSetType (_eps, EPSPOWER);    LIBMESH_CHKERR(ierr); return;
    case SUBSPACE:
      ierr = EPSSetType (_eps, EPSSUBSPACE); LIBMESH_CHKERR(ierr); return;
    case LAPACK:
      ierr = EPSSetType (_eps, EPSLAPACK);   LIBMESH_CHKERR(ierr); return;
    case ARNOLDI:
      ierr = EPSSetType (_eps, EPSARNOLDI);  LIBMESH_CHKERR(ierr); return;
    case LANCZOS:
      ierr = EPSSetType (_eps, EPSLANCZOS);  LIBMESH_CHKERR(ierr); return;
    case KRYLOVSCHUR:
      ierr = EPSSetType (_eps, EPSKRYLOVSCHUR);  LIBMESH_CHKERR(ierr); return;
      // case ARPACK:
      // ierr = EPSSetType (_eps, (char *) EPSARPACK);   LIBMESH_CHKERR(ierr); return;

    default:
      libMesh::err << "ERROR:  Unsupported SLEPc Eigen Solver: "
                   << Utility::enum_to_string(this->_eigen_solver_type) << std::endl
                   << "Continuing with SLEPc defaults" << std::endl;
    }
}




template <typename T>
void SlepcEigenSolver<T>:: set_slepc_problem_type()
{
  PetscErrorCode ierr = LIBMESH_PETSC_SUCCESS;

  switch (this->_eigen_problem_type)
    {
    case NHEP:
      ierr = EPSSetProblemType (_eps, EPS_NHEP);  LIBMESH_CHKERR(ierr); return;
    case GNHEP:
      ierr = EPSSetProblemType (_eps, EPS_GNHEP); LIBMESH_CHKERR(ierr); return;
    case HEP:
      ierr = EPSSetProblemType (_eps, EPS_HEP);   LIBMESH_CHKERR(ierr); return;
    case GHEP:
      ierr = EPSSetProblemType (_eps, EPS_GHEP);  LIBMESH_CHKERR(ierr); return;
#if !SLEPC_VERSION_LESS_THAN(3,3,0)
      // EPS_GHIEP added in 3.3.0
    case GHIEP:
      ierr = EPSSetProblemType (_eps, EPS_GHIEP);  LIBMESH_CHKERR(ierr); return;
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
  PetscErrorCode ierr = LIBMESH_PETSC_SUCCESS;

  switch (this->_position_of_spectrum)
    {
    case LARGEST_MAGNITUDE:
      {
        ierr = EPSSetWhichEigenpairs (_eps, EPS_LARGEST_MAGNITUDE);
        LIBMESH_CHKERR(ierr);
        return;
      }
    case SMALLEST_MAGNITUDE:
      {
        ierr = EPSSetWhichEigenpairs (_eps, EPS_SMALLEST_MAGNITUDE);
        LIBMESH_CHKERR(ierr);
        return;
      }
    case LARGEST_REAL:
      {
        ierr = EPSSetWhichEigenpairs (_eps, EPS_LARGEST_REAL);
        LIBMESH_CHKERR(ierr);
        return;
      }
    case SMALLEST_REAL:
      {
        ierr = EPSSetWhichEigenpairs (_eps, EPS_SMALLEST_REAL);
        LIBMESH_CHKERR(ierr);
        return;
      }
    case LARGEST_IMAGINARY:
      {
        ierr = EPSSetWhichEigenpairs (_eps, EPS_LARGEST_IMAGINARY);
        LIBMESH_CHKERR(ierr);
        return;
      }
    case SMALLEST_IMAGINARY:
      {
        ierr = EPSSetWhichEigenpairs (_eps, EPS_SMALLEST_IMAGINARY);
        LIBMESH_CHKERR(ierr);
        return;
      }

      // The EPS_TARGET_XXX enums were added in SLEPc 3.1
#if !SLEPC_VERSION_LESS_THAN(3,1,0)
    case TARGET_MAGNITUDE:
      {
        ierr = EPSSetTarget(_eps, PS(this->_target_val));
        LIBMESH_CHKERR(ierr);
        ierr = EPSSetWhichEigenpairs (_eps, EPS_TARGET_MAGNITUDE);
        LIBMESH_CHKERR(ierr);
        return;
      }
    case TARGET_REAL:
      {
        ierr = EPSSetTarget(_eps, PS(this->_target_val));
        LIBMESH_CHKERR(ierr);
        ierr = EPSSetWhichEigenpairs (_eps, EPS_TARGET_REAL);
        LIBMESH_CHKERR(ierr);
        return;
      }
    case TARGET_IMAGINARY:
      {
        ierr = EPSSetTarget(_eps, PS(this->_target_val));
        LIBMESH_CHKERR(ierr);
        ierr = EPSSetWhichEigenpairs (_eps, EPS_TARGET_IMAGINARY);
        LIBMESH_CHKERR(ierr);
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
  PetscErrorCode ierr = LIBMESH_PETSC_SUCCESS;

  PetscReal re, im;

  // Make sure the NumericVector passed in is really a PetscVector
  PetscVector<T> * solution = dynamic_cast<PetscVector<T> *>(&solution_in);

  libmesh_error_msg_if(!solution, "Error getting eigenvector: input vector must be a PetscVector.");

  // real and imaginary part of the ith eigenvalue.
  PetscScalar kr, ki;

  solution->close();

  ierr = EPSGetEigenpair(_eps, i, &kr, &ki, solution->vec(),
                         LIBMESH_PETSC_NULLPTR);
  LIBMESH_CHKERR(ierr);

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
  PetscErrorCode ierr = LIBMESH_PETSC_SUCCESS;

  PetscReal re, im;

  // real and imaginary part of the ith eigenvalue.
  PetscScalar kr, ki;

  ierr = EPSGetEigenvalue(_eps, i, &kr, &ki);
  LIBMESH_CHKERR(ierr);

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
  PetscErrorCode ierr = LIBMESH_PETSC_SUCCESS;
  PetscReal error;

#if SLEPC_VERSION_LESS_THAN(3,6,0)
  ierr = EPSComputeRelativeError(_eps, i, &error);
#else
  ierr = EPSComputeError(_eps, i, EPS_ERROR_RELATIVE, &error);
#endif
  LIBMESH_CHKERR(ierr);

  return error;
}


template <typename T>
void SlepcEigenSolver<T>::attach_deflation_space(NumericVector<T> & deflation_vector_in)
{
  this->init();

  PetscErrorCode ierr = LIBMESH_PETSC_SUCCESS;

  // Make sure the input vector is actually a PetscVector
  PetscVector<T> * deflation_vector_petsc_vec =
    dynamic_cast<PetscVector<T> *>(&deflation_vector_in);

  libmesh_error_msg_if(!deflation_vector_petsc_vec, "Error attaching deflation space: input vector must be a PetscVector.");

  // Get a handle for the underlying Vec.
  Vec deflation_vector = deflation_vector_petsc_vec->vec();

#if SLEPC_VERSION_LESS_THAN(3,1,0)
  ierr = EPSAttachDeflationSpace(_eps, 1, &deflation_vector, PETSC_FALSE);
#else
  ierr = EPSSetDeflationSpace(_eps, 1, &deflation_vector);
#endif
  LIBMESH_CHKERR(ierr);
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
  PetscErrorCode ierr = LIBMESH_PETSC_SUCCESS;
  void * ctx;
  ierr = MatShellGetContext(mat,&ctx);
  CHKERRQ(ierr);

  // Get user shell matrix object.
  const ShellMatrix<T> & shell_matrix = *static_cast<const ShellMatrix<T> *>(ctx);

  // Make \p NumericVector instances around the vectors.
  PetscVector<T> arg_global(arg,   shell_matrix.comm());
  PetscVector<T> dest_global(dest, shell_matrix.comm());

  // Call the user function.
  shell_matrix.vector_mult(dest_global,arg_global);

  PetscFunctionReturn(ierr);
}

template <typename T>
PetscErrorCode SlepcEigenSolver<T>::_petsc_shell_matrix_get_diagonal(Mat mat, Vec dest)
{
  PetscFunctionBegin;

  // Get the matrix context.
  PetscErrorCode ierr = LIBMESH_PETSC_SUCCESS;
  void * ctx;
  ierr = MatShellGetContext(mat,&ctx);
  CHKERRQ(ierr);

  // Get user shell matrix object.
  const ShellMatrix<T> & shell_matrix = *static_cast<const ShellMatrix<T> *>(ctx);

  // Make \p NumericVector instances around the vector.
  PetscVector<T> dest_global(dest, shell_matrix.comm());

  // Call the user function.
  shell_matrix.get_diagonal(dest_global);

  PetscFunctionReturn(ierr);
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
