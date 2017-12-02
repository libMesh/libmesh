// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


// C++ includes

// Local Includes
#include "libmesh/libmesh_logging.h"
#include "libmesh/petsc_matrix.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/slepc_eigen_solver.h"
#include "libmesh/shell_matrix.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/solver_configuration.h"

namespace libMesh
{

template <typename T>
void SlepcEigenSolver<T>::clear ()
{
  if (this->initialized())
    {
      this->_is_initialized = false;

      PetscErrorCode ierr=0;

      ierr = LibMeshEPSDestroy(&_eps);
      LIBMESH_CHKERR(ierr);

      // SLEPc default eigenproblem solver
      this->_eigen_solver_type = KRYLOVSCHUR;
    }
}



template <typename T>
void SlepcEigenSolver<T>::init ()
{

  PetscErrorCode ierr=0;

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

  this->init ();

  // Make sure the SparseMatrix passed in is really a PetscMatrix
  PetscMatrix<T> * matrix_A = dynamic_cast<PetscMatrix<T> *>(&matrix_A_in);

  if (!matrix_A)
    libmesh_error_msg("Error: input matrix to solve_standard() must be a PetscMatrix.");

  // Close the matrix and vectors in case this wasn't already done.
  if (this->_close_matrix_before_solve)
    matrix_A->close ();

  return _solve_standard_helper(matrix_A->mat(), nev, ncv, tol, m_its);
}


template <typename T>
std::pair<unsigned int, unsigned int>
SlepcEigenSolver<T>::solve_standard (ShellMatrix<T> & shell_matrix,
                                     int nev,                  // number of requested eigenpairs
                                     int ncv,                  // number of basis vectors
                                     const double tol,         // solver tolerance
                                     const unsigned int m_its) // maximum number of iterations
{
  this->init ();

  PetscErrorCode ierr=0;

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

  return _solve_standard_helper(mat, nev, ncv, tol, m_its);
}

template <typename T>
std::pair<unsigned int, unsigned int>
SlepcEigenSolver<T>::_solve_standard_helper(Mat mat,
                                            int nev,                  // number of requested eigenpairs
                                            int ncv,                  // number of basis vectors
                                            const double tol,         // solver tolerance
                                            const unsigned int m_its) // maximum number of iterations
{
  LOG_SCOPE("solve_standard()", "SlepcEigenSolver");

  PetscErrorCode ierr=0;

  // converged eigen pairs and number of iterations
  PetscInt nconv=0;
  PetscInt its=0;

#ifdef  DEBUG
  // The relative error.
  PetscReal error, re, im;

  // Pointer to vectors of the real parts, imaginary parts.
  PetscScalar kr, ki;
#endif

  // Set operators.
  ierr = EPSSetOperators (_eps, mat, PETSC_NULL);
  LIBMESH_CHKERR(ierr);

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

  // If the SolverConfiguration object is provided, use it to override
  // solver options.
  if (this->_solver_configuration)
    {
      this->_solver_configuration->configure_solver();
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


#ifdef DEBUG
  // ierr = PetscPrintf(this->comm().get(),
  //         "\n Number of iterations: %d\n"
  //         " Number of converged eigenpairs: %d\n\n", its, nconv);

  // Display eigenvalues and relative errors.
  ierr = PetscPrintf(this->comm().get(),
                     "           k           ||Ax-kx||/|kx|\n"
                     "   ----------------- -----------------\n" );
  LIBMESH_CHKERR(ierr);

  for (PetscInt i=0; i<nconv; i++ )
    {
      ierr = EPSGetEigenpair(_eps, i, &kr, &ki, PETSC_NULL, PETSC_NULL);
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
#endif // DEBUG

  // return the number of converged eigenpairs
  // and the number of iterations
  return std::make_pair(nconv, its);
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
  this->init ();

  // Make sure the data passed in are really of Petsc types
  PetscMatrix<T> * matrix_A = dynamic_cast<PetscMatrix<T> *>(&matrix_A_in);
  PetscMatrix<T> * matrix_B = dynamic_cast<PetscMatrix<T> *>(&matrix_B_in);

  if (!matrix_A || !matrix_B)
    libmesh_error_msg("Error: inputs to solve_generalized() must be of type PetscMatrix.");

  // Close the matrix and vectors in case this wasn't already done.
  if (this->_close_matrix_before_solve)
    {
      matrix_A->close ();
      matrix_B->close ();
    }

  return _solve_generalized_helper (matrix_A->mat(), matrix_B->mat(), nev, ncv, tol, m_its);
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
  this->init ();

  PetscErrorCode ierr=0;

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

  if (!matrix_B)
    libmesh_error_msg("Error: inputs to solve_generalized() must be of type PetscMatrix.");

  // Close the matrix and vectors in case this wasn't already done.
  if (this->_close_matrix_before_solve)
    matrix_B->close ();

  ierr = MatShellSetOperation(mat_A,MATOP_MULT,reinterpret_cast<void(*)(void)>(_petsc_shell_matrix_mult));
  LIBMESH_CHKERR(ierr);
  ierr = MatShellSetOperation(mat_A,MATOP_GET_DIAGONAL,reinterpret_cast<void(*)(void)>(_petsc_shell_matrix_get_diagonal));
  LIBMESH_CHKERR(ierr);

  return _solve_generalized_helper (mat_A, matrix_B->mat(), nev, ncv, tol, m_its);
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
  this->init ();

  PetscErrorCode ierr=0;

  PetscMatrix<T> * matrix_A = dynamic_cast<PetscMatrix<T> *>(&matrix_A_in);

  if (!matrix_A)
    libmesh_error_msg("Error: inputs to solve_generalized() must be of type PetscMatrix.");

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

  return _solve_generalized_helper (matrix_A->mat(), mat_B, nev, ncv, tol, m_its);
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
  this->init ();

  PetscErrorCode ierr=0;

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

  return _solve_generalized_helper (mat_A, mat_B, nev, ncv, tol, m_its);
}



template <typename T>
std::pair<unsigned int, unsigned int>
SlepcEigenSolver<T>::_solve_generalized_helper (Mat mat_A,
                                                Mat mat_B,
                                                int nev,                  // number of requested eigenpairs
                                                int ncv,                  // number of basis vectors
                                                const double tol,         // solver tolerance
                                                const unsigned int m_its) // maximum number of iterations
{
  LOG_SCOPE("solve_generalized()", "SlepcEigenSolver");

  PetscErrorCode ierr=0;

  // converged eigen pairs and number of iterations
  PetscInt nconv=0;
  PetscInt its=0;

#ifdef  DEBUG
  // The relative error.
  PetscReal error, re, im;

  // Pointer to vectors of the real parts, imaginary parts.
  PetscScalar kr, ki;
#endif

  // Set operators.
  ierr = EPSSetOperators (_eps, mat_A, mat_B);
  LIBMESH_CHKERR(ierr);

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

  // If the SolverConfiguration object is provided, use it to override
  // solver options.
  if (this->_solver_configuration)
    {
      this->_solver_configuration->configure_solver();
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


#ifdef DEBUG
  // ierr = PetscPrintf(this->comm().get(),
  //         "\n Number of iterations: %d\n"
  //         " Number of converged eigenpairs: %d\n\n", its, nconv);

  // Display eigenvalues and relative errors.
  ierr = PetscPrintf(this->comm().get(),
                     "           k           ||Ax-kx||/|kx|\n"
                     "   ----------------- -----------------\n" );
  LIBMESH_CHKERR(ierr);

  for (PetscInt i=0; i<nconv; i++ )
    {
      ierr = EPSGetEigenpair(_eps, i, &kr, &ki, PETSC_NULL, PETSC_NULL);
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
#endif // DEBUG

  // return the number of converged eigenpairs
  // and the number of iterations
  return std::make_pair(nconv, its);
}











template <typename T>
void SlepcEigenSolver<T>::set_slepc_solver_type()
{
  PetscErrorCode ierr = 0;

  switch (this->_eigen_solver_type)
    {
    case POWER:
      ierr = EPSSetType (_eps, (char *) EPSPOWER);    LIBMESH_CHKERR(ierr); return;
    case SUBSPACE:
      ierr = EPSSetType (_eps, (char *) EPSSUBSPACE); LIBMESH_CHKERR(ierr); return;
    case LAPACK:
      ierr = EPSSetType (_eps, (char *) EPSLAPACK);   LIBMESH_CHKERR(ierr); return;
    case ARNOLDI:
      ierr = EPSSetType (_eps, (char *) EPSARNOLDI);  LIBMESH_CHKERR(ierr); return;
    case LANCZOS:
      ierr = EPSSetType (_eps, (char *) EPSLANCZOS);  LIBMESH_CHKERR(ierr); return;
    case KRYLOVSCHUR:
      ierr = EPSSetType (_eps, (char *) EPSKRYLOVSCHUR);  LIBMESH_CHKERR(ierr); return;
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
  PetscErrorCode ierr = 0;

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
  PetscErrorCode ierr = 0;

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
        ierr = EPSSetTarget(_eps, this->_target_val);
        LIBMESH_CHKERR(ierr);
        ierr = EPSSetWhichEigenpairs (_eps, EPS_TARGET_MAGNITUDE);
        LIBMESH_CHKERR(ierr);
        return;
      }
    case TARGET_REAL:
      {
        ierr = EPSSetTarget(_eps, this->_target_val);
        LIBMESH_CHKERR(ierr);
        ierr = EPSSetWhichEigenpairs (_eps, EPS_TARGET_REAL);
        LIBMESH_CHKERR(ierr);
        return;
      }
    case TARGET_IMAGINARY:
      {
        ierr = EPSSetTarget(_eps, this->_target_val);
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
  PetscErrorCode ierr=0;

  PetscReal re, im;

  // Make sure the NumericVector passed in is really a PetscVector
  PetscVector<T> * solution = dynamic_cast<PetscVector<T> *>(&solution_in);

  if (!solution)
    libmesh_error_msg("Error getting eigenvector: input vector must be a PetscVector.");

  // real and imaginary part of the ith eigenvalue.
  PetscScalar kr, ki;

  solution->close();

  ierr = EPSGetEigenpair(_eps, i, &kr, &ki, solution->vec(), PETSC_NULL);
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
  PetscErrorCode ierr=0;

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
  PetscErrorCode ierr=0;
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

  PetscErrorCode ierr = 0;

  // Make sure the input vector is actually a PetscVector
  PetscVector<T> * deflation_vector_petsc_vec =
    dynamic_cast<PetscVector<T> *>(&deflation_vector_in);

  if (!deflation_vector_petsc_vec)
    libmesh_error_msg("Error attaching deflation space: input vector must be a PetscVector.");

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
  this->init();

  PetscErrorCode ierr = 0;

  // Make sure the input vector is actually a PetscVector
  PetscVector<T> * initial_space_petsc_vec =
    dynamic_cast<PetscVector<T> *>(&initial_space_in);

  if (!initial_space_petsc_vec)
    libmesh_error_msg("Error attaching initial space: input vector must be a PetscVector.");

  // Get a handle for the underlying Vec.
  Vec initial_vector = initial_space_petsc_vec->vec();

  ierr = EPSSetInitialSpace(_eps, 1, &initial_vector);
  LIBMESH_CHKERR(ierr);
#endif
}

template <typename T>
PetscErrorCode SlepcEigenSolver<T>::_petsc_shell_matrix_mult(Mat mat, Vec arg, Vec dest)
{
  // Get the matrix context.
  PetscErrorCode ierr=0;
  void * ctx;
  ierr = MatShellGetContext(mat,&ctx);

  Parallel::communicator comm;
  PetscObjectGetComm((PetscObject)mat,&comm);
  CHKERRABORT(comm,ierr);

  // Get user shell matrix object.
  const ShellMatrix<T> & shell_matrix = *static_cast<const ShellMatrix<T> *>(ctx);

  // Make \p NumericVector instances around the vectors.
  PetscVector<T> arg_global(arg,   shell_matrix.comm());
  PetscVector<T> dest_global(dest, shell_matrix.comm());

  // Call the user function.
  shell_matrix.vector_mult(dest_global,arg_global);

  return ierr;
}

template <typename T>
PetscErrorCode SlepcEigenSolver<T>::_petsc_shell_matrix_get_diagonal(Mat mat, Vec dest)
{
  // Get the matrix context.
  PetscErrorCode ierr=0;
  void * ctx;
  ierr = MatShellGetContext(mat,&ctx);

  Parallel::communicator comm;
  PetscObjectGetComm((PetscObject)mat,&comm);
  CHKERRABORT(comm,ierr);

  // Get user shell matrix object.
  const ShellMatrix<T> & shell_matrix = *static_cast<const ShellMatrix<T> *>(ctx);

  // Make \p NumericVector instances around the vector.
  PetscVector<T> dest_global(dest, shell_matrix.comm());

  // Call the user function.
  shell_matrix.get_diagonal(dest_global);

  return ierr;
}

//------------------------------------------------------------------
// Explicit instantiations
template class SlepcEigenSolver<Number>;

} // namespace libMesh



#endif // #ifdef LIBMESH_HAVE_SLEPC
