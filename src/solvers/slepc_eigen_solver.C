// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

namespace libMesh
{


/*----------------------- functions ----------------------------------*/
template <typename T>
void SlepcEigenSolver<T>::clear ()
{
  if (this->initialized())
    {
      this->_is_initialized = false;

      int ierr=0;

      ierr = LibMeshEPSDestroy(&_eps);
             CHKERRABORT(libMesh::COMM_WORLD,ierr);

      // SLEPc default eigenproblem solver
#if SLEPC_VERSION_LESS_THAN(2,3,2)
    this->_eigen_solver_type = ARNOLDI;
#else
    // Krylov-Schur showed up as of Slepc 2.3.2
    this->_eigen_solver_type = KRYLOVSCHUR;
#endif
    }
}



template <typename T>
void SlepcEigenSolver<T>::init ()
{

  int ierr=0;

  // Initialize the data structures if not done so already.
  if (!this->initialized())
    {
      this->_is_initialized = true;

      // Create the eigenproblem solver context
      ierr = EPSCreate (libMesh::COMM_WORLD, &_eps);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

      // Set user-specified  solver
      set_slepc_solver_type();
    }
}



template <typename T>
std::pair<unsigned int, unsigned int>
SlepcEigenSolver<T>::solve_standard (SparseMatrix<T> &matrix_A_in,
				     int nev,                  // number of requested eigenpairs
				     int ncv,                  // number of basis vectors
				     const double tol,         // solver tolerance
				     const unsigned int m_its) // maximum number of iterations
{
//   START_LOG("solve_standard()", "SlepcEigenSolver");

  this->init ();

  // Make sure the SparseMatrix passed in is really a PetscMatrix
  PetscMatrix<T>* matrix_A   = libmesh_cast_ptr<PetscMatrix<T>*>(&matrix_A_in);

  // Close the matrix and vectors in case this wasn't already done.
  matrix_A->close ();

  // just for debugging, remove this
//   char mat_file[] = "matA.petsc";
//   PetscViewer petsc_viewer;
//   ierr = PetscViewerBinaryOpen(libMesh::COMM_WORLD, mat_file, PETSC_FILE_CREATE, &petsc_viewer);
//          CHKERRABORT(libMesh::COMM_WORLD,ierr);
//   ierr = MatView(matrix_A->mat(),petsc_viewer);
//          CHKERRABORT(libMesh::COMM_WORLD,ierr);

  return _solve_standard_helper(matrix_A->mat(), nev, ncv, tol, m_its);
}


template <typename T>
std::pair<unsigned int, unsigned int>
SlepcEigenSolver<T>::solve_standard (ShellMatrix<T> &shell_matrix,
				     int nev,                  // number of requested eigenpairs
				     int ncv,                  // number of basis vectors
				     const double tol,         // solver tolerance
				     const unsigned int m_its) // maximum number of iterations
{
  this->init ();

  int ierr=0;

  // Prepare the matrix.
  Mat mat;
  ierr = MatCreateShell(libMesh::COMM_WORLD,
            shell_matrix.m(), // Specify the number of local rows
            shell_matrix.n(), // Specify the number of local columns
            PETSC_DETERMINE,
            PETSC_DETERMINE,
            const_cast<void*>(static_cast<const void*>(&shell_matrix)),
            &mat);

  /* Note that the const_cast above is only necessary because PETSc
     does not accept a const void*.  Inside the member function
     _petsc_shell_matrix() below, the pointer is casted back to a
     const ShellMatrix<T>*.  */

  CHKERRABORT(libMesh::COMM_WORLD,ierr);
  ierr = MatShellSetOperation(mat,MATOP_MULT,reinterpret_cast<void(*)(void)>(_petsc_shell_matrix_mult));
  ierr = MatShellSetOperation(mat,MATOP_GET_DIAGONAL,reinterpret_cast<void(*)(void)>(_petsc_shell_matrix_get_diagonal));
  CHKERRABORT(libMesh::COMM_WORLD,ierr);


  return _solve_standard_helper(mat, nev, ncv, tol, m_its);
}

template <typename T>
std::pair<unsigned int, unsigned int>
SlepcEigenSolver<T>::_solve_standard_helper
				    (Mat mat,
				     int nev,                  // number of requested eigenpairs
				     int ncv,                  // number of basis vectors
				     const double tol,         // solver tolerance
				     const unsigned int m_its) // maximum number of iterations
{
  START_LOG("solve_standard()", "SlepcEigenSolver");

  int ierr=0;

  // converged eigen pairs and number of iterations
  int nconv=0;
  int its=0;

#ifdef  DEBUG
  // The relative error.
  PetscReal error, re, im;

  // Pointer to vectors of the real parts, imaginary parts.
  PetscScalar kr, ki;
#endif

  // Set operators.
  ierr = EPSSetOperators (_eps, mat, PETSC_NULL);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

  //set the problem type and the position of the spectrum
  set_slepc_problem_type();
  set_slepc_position_of_spectrum();

  // Set eigenvalues to be computed.
#if SLEPC_VERSION_LESS_THAN(3,0,0)
  ierr = EPSSetDimensions (_eps, nev, ncv);
#else
  ierr = EPSSetDimensions (_eps, nev, ncv, PETSC_DECIDE);
#endif
         CHKERRABORT(libMesh::COMM_WORLD,ierr);
  // Set the tolerance and maximum iterations.
  ierr = EPSSetTolerances (_eps, tol, m_its);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

  // Set runtime options, e.g.,
  //      -eps_type <type>, -eps_nev <nev>, -eps_ncv <ncv>
  // Similar to PETSc, these options will override those specified
  // above as long as EPSSetFromOptions() is called _after_ any
  // other customization routines.
  ierr = EPSSetFromOptions (_eps);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

  // Solve the eigenproblem.
  ierr = EPSSolve (_eps);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

  // Get the number of iterations.
  ierr = EPSGetIterationNumber (_eps, &its);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

  // Get number of converged eigenpairs.
  ierr = EPSGetConverged(_eps,&nconv);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);


#ifdef DEBUG
	 // ierr = PetscPrintf(libMesh::COMM_WORLD,
	 //         "\n Number of iterations: %d\n"
	 //         " Number of converged eigenpairs: %d\n\n", its, nconv);

  // Display eigenvalues and relative errors.
  ierr = PetscPrintf(libMesh::COMM_WORLD,
		     "           k           ||Ax-kx||/|kx|\n"
		     "   ----------------- -----------------\n" );
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

  for(int i=0; i<nconv; i++ )
    {
      ierr = EPSGetEigenpair(_eps, i, &kr, &ki, PETSC_NULL, PETSC_NULL);
             CHKERRABORT(libMesh::COMM_WORLD,ierr);

      ierr = EPSComputeRelativeError(_eps, i, &error);
             CHKERRABORT(libMesh::COMM_WORLD,ierr);

#ifdef LIBMESH_USE_COMPLEX_NUMBERS
      re = PetscRealPart(kr);
      im = PetscImaginaryPart(kr);
#else
      re = kr;
      im = ki;
#endif

      if (im != .0)
	{
	  ierr = PetscPrintf(libMesh::COMM_WORLD," %9f%+9f i %12f\n", re, im, error);
	         CHKERRABORT(libMesh::COMM_WORLD,ierr);
	}
      else
	{
	  ierr = PetscPrintf(libMesh::COMM_WORLD,"   %12f       %12f\n", re, error);
	         CHKERRABORT(libMesh::COMM_WORLD,ierr);
	}
    }

  ierr = PetscPrintf(libMesh::COMM_WORLD,"\n" );
         CHKERRABORT(libMesh::COMM_WORLD,ierr);
#endif // DEBUG


  STOP_LOG("solve_standard()", "SlepcEigenSolver");

  // return the number of converged eigenpairs
  // and the number of iterations
  return std::make_pair(nconv, its);

}





template <typename T>
std::pair<unsigned int, unsigned int>
SlepcEigenSolver<T>::solve_generalized (SparseMatrix<T> &matrix_A_in,
					SparseMatrix<T> &matrix_B_in,
					int nev,                  // number of requested eigenpairs
					int ncv,                  // number of basis vectors
					const double tol,         // solver tolerance
					const unsigned int m_its) // maximum number of iterations
{
  this->init ();

  // Make sure the data passed in are really of Petsc types
  PetscMatrix<T>* matrix_A   = libmesh_cast_ptr<PetscMatrix<T>*>(&matrix_A_in);
  PetscMatrix<T>* matrix_B   = libmesh_cast_ptr<PetscMatrix<T>*>(&matrix_B_in);

  // Close the matrix and vectors in case this wasn't already done.
  matrix_A->close ();
  matrix_B->close ();


  return _solve_generalized_helper (matrix_A->mat(), matrix_B->mat(), nev, ncv, tol, m_its);
}

template <typename T>
std::pair<unsigned int, unsigned int>
SlepcEigenSolver<T>::solve_generalized (ShellMatrix<T> &shell_matrix_A,
					SparseMatrix<T> &matrix_B_in,
					int nev,                  // number of requested eigenpairs
					int ncv,                  // number of basis vectors
					const double tol,         // solver tolerance
					const unsigned int m_its) // maximum number of iterations
{
  this->init ();

  int ierr=0;

  // Prepare the matrix.
  Mat mat_A;
  ierr = MatCreateShell(libMesh::COMM_WORLD,
            shell_matrix_A.m(), // Specify the number of local rows
            shell_matrix_A.n(), // Specify the number of local columns
            PETSC_DETERMINE,
            PETSC_DETERMINE,
            const_cast<void*>(static_cast<const void*>(&shell_matrix_A)),
            &mat_A);

  PetscMatrix<T>* matrix_B   = libmesh_cast_ptr<PetscMatrix<T>*>(&matrix_B_in);

  // Close the matrix and vectors in case this wasn't already done.
  matrix_B->close ();

  /* Note that the const_cast above is only necessary because PETSc
     does not accept a const void*.  Inside the member function
     _petsc_shell_matrix() below, the pointer is casted back to a
     const ShellMatrix<T>*.  */

  CHKERRABORT(libMesh::COMM_WORLD,ierr);
  ierr = MatShellSetOperation(mat_A,MATOP_MULT,reinterpret_cast<void(*)(void)>(_petsc_shell_matrix_mult));
  ierr = MatShellSetOperation(mat_A,MATOP_GET_DIAGONAL,reinterpret_cast<void(*)(void)>(_petsc_shell_matrix_get_diagonal));
  CHKERRABORT(libMesh::COMM_WORLD,ierr);

  return _solve_generalized_helper (mat_A, matrix_B->mat(), nev, ncv, tol, m_its);
}

template <typename T>
std::pair<unsigned int, unsigned int>
SlepcEigenSolver<T>::solve_generalized (SparseMatrix<T> &matrix_A_in,
					ShellMatrix<T> &shell_matrix_B,
					int nev,                  // number of requested eigenpairs
					int ncv,                  // number of basis vectors
					const double tol,         // solver tolerance
					const unsigned int m_its) // maximum number of iterations
{
  this->init ();

  int ierr=0;

  PetscMatrix<T>* matrix_A = libmesh_cast_ptr<PetscMatrix<T>*>(&matrix_A_in);

  // Close the matrix and vectors in case this wasn't already done.
  matrix_A->close ();

  // Prepare the matrix.
  Mat mat_B;
  ierr = MatCreateShell(libMesh::COMM_WORLD,
            shell_matrix_B.m(), // Specify the number of local rows
            shell_matrix_B.n(), // Specify the number of local columns
            PETSC_DETERMINE,
            PETSC_DETERMINE,
            const_cast<void*>(static_cast<const void*>(&shell_matrix_B)),
            &mat_B);


  /* Note that the const_cast above is only necessary because PETSc
     does not accept a const void*.  Inside the member function
     _petsc_shell_matrix() below, the pointer is casted back to a
     const ShellMatrix<T>*.  */

  CHKERRABORT(libMesh::COMM_WORLD,ierr);
  ierr = MatShellSetOperation(mat_B,MATOP_MULT,reinterpret_cast<void(*)(void)>(_petsc_shell_matrix_mult));
  ierr = MatShellSetOperation(mat_B,MATOP_GET_DIAGONAL,reinterpret_cast<void(*)(void)>(_petsc_shell_matrix_get_diagonal));
  CHKERRABORT(libMesh::COMM_WORLD,ierr);

  return _solve_generalized_helper (matrix_A->mat(), mat_B, nev, ncv, tol, m_its);
}

template <typename T>
std::pair<unsigned int, unsigned int>
SlepcEigenSolver<T>::solve_generalized (ShellMatrix<T> &shell_matrix_A,
					ShellMatrix<T> &shell_matrix_B,
					int nev,                  // number of requested eigenpairs
					int ncv,                  // number of basis vectors
					const double tol,         // solver tolerance
					const unsigned int m_its) // maximum number of iterations
{
  this->init ();

  int ierr=0;

  // Prepare the matrix.
  Mat mat_A;
  ierr = MatCreateShell(libMesh::COMM_WORLD,
            shell_matrix_A.m(), // Specify the number of local rows
            shell_matrix_A.n(), // Specify the number of local columns
            PETSC_DETERMINE,
            PETSC_DETERMINE,
            const_cast<void*>(static_cast<const void*>(&shell_matrix_A)),
            &mat_A);

  Mat mat_B;
  ierr = MatCreateShell(libMesh::COMM_WORLD,
            shell_matrix_B.m(), // Specify the number of local rows
            shell_matrix_B.n(), // Specify the number of local columns
            PETSC_DETERMINE,
            PETSC_DETERMINE,
            const_cast<void*>(static_cast<const void*>(&shell_matrix_B)),
            &mat_B);

  /* Note that the const_cast above is only necessary because PETSc
     does not accept a const void*.  Inside the member function
     _petsc_shell_matrix() below, the pointer is casted back to a
     const ShellMatrix<T>*.  */

  CHKERRABORT(libMesh::COMM_WORLD,ierr);
  ierr = MatShellSetOperation(mat_A,MATOP_MULT,reinterpret_cast<void(*)(void)>(_petsc_shell_matrix_mult));
  ierr = MatShellSetOperation(mat_A,MATOP_GET_DIAGONAL,reinterpret_cast<void(*)(void)>(_petsc_shell_matrix_get_diagonal));
  CHKERRABORT(libMesh::COMM_WORLD,ierr);

  CHKERRABORT(libMesh::COMM_WORLD,ierr);
  ierr = MatShellSetOperation(mat_B,MATOP_MULT,reinterpret_cast<void(*)(void)>(_petsc_shell_matrix_mult));
  ierr = MatShellSetOperation(mat_B,MATOP_GET_DIAGONAL,reinterpret_cast<void(*)(void)>(_petsc_shell_matrix_get_diagonal));
  CHKERRABORT(libMesh::COMM_WORLD,ierr);

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
  START_LOG("solve_generalized()", "SlepcEigenSolver");

  int ierr=0;

  // converged eigen pairs and number of iterations
  int nconv=0;
  int its=0;

#ifdef  DEBUG
  // The relative error.
  PetscReal error, re, im;

  // Pointer to vectors of the real parts, imaginary parts.
  PetscScalar kr, ki;
#endif

  // Set operators.
  ierr = EPSSetOperators (_eps, mat_A, mat_B);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

  //set the problem type and the position of the spectrum
  set_slepc_problem_type();
  set_slepc_position_of_spectrum();

  // Set eigenvalues to be computed.
#if SLEPC_VERSION_LESS_THAN(3,0,0)
  ierr = EPSSetDimensions (_eps, nev, ncv);
#else
  ierr = EPSSetDimensions (_eps, nev, ncv, PETSC_DECIDE);
#endif
         CHKERRABORT(libMesh::COMM_WORLD,ierr);


  // Set the tolerance and maximum iterations.
  ierr = EPSSetTolerances (_eps, tol, m_its);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

  // Set runtime options, e.g.,
  //      -eps_type <type>, -eps_nev <nev>, -eps_ncv <ncv>
  // Similar to PETSc, these options will override those specified
  // above as long as EPSSetFromOptions() is called _after_ any
  // other customization routines.
  ierr = EPSSetFromOptions (_eps);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

  // Solve the eigenproblem.
  ierr = EPSSolve (_eps);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

  // Get the number of iterations.
  ierr = EPSGetIterationNumber (_eps, &its);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

  // Get number of converged eigenpairs.
  ierr = EPSGetConverged(_eps,&nconv);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);


#ifdef DEBUG
	 // ierr = PetscPrintf(libMesh::COMM_WORLD,
	 //         "\n Number of iterations: %d\n"
	 //         " Number of converged eigenpairs: %d\n\n", its, nconv);

  // Display eigenvalues and relative errors.
  ierr = PetscPrintf(libMesh::COMM_WORLD,
		     "           k           ||Ax-kx||/|kx|\n"
		     "   ----------------- -----------------\n" );
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

  for(int i=0; i<nconv; i++ )
    {
      ierr = EPSGetEigenpair(_eps, i, &kr, &ki, PETSC_NULL, PETSC_NULL);
             CHKERRABORT(libMesh::COMM_WORLD,ierr);

      ierr = EPSComputeRelativeError(_eps, i, &error);
             CHKERRABORT(libMesh::COMM_WORLD,ierr);

#ifdef LIBMESH_USE_COMPLEX_NUMBERS
      re = PetscRealPart(kr);
      im = PetscImaginaryPart(kr);
#else
      re = kr;
      im = ki;
#endif

      if (im != .0)
	{
	  ierr = PetscPrintf(libMesh::COMM_WORLD," %9f%+9f i %12f\n", re, im, error);
	         CHKERRABORT(libMesh::COMM_WORLD,ierr);
	}
      else
	{
	  ierr = PetscPrintf(libMesh::COMM_WORLD,"   %12f       %12f\n", re, error);
	         CHKERRABORT(libMesh::COMM_WORLD,ierr);
	}
    }

  ierr = PetscPrintf(libMesh::COMM_WORLD,"\n" );
         CHKERRABORT(libMesh::COMM_WORLD,ierr);
#endif // DEBUG

  STOP_LOG("solve_generalized()", "SlepcEigenSolver");

  // return the number of converged eigenpairs
  // and the number of iterations
  return std::make_pair(nconv, its);

}











template <typename T>
void SlepcEigenSolver<T>::set_slepc_solver_type()
{
  int ierr = 0;

  switch (this->_eigen_solver_type)
    {
    case POWER:
      ierr = EPSSetType (_eps, (char*) EPSPOWER);    CHKERRABORT(libMesh::COMM_WORLD,ierr); return;
    case SUBSPACE:
      ierr = EPSSetType (_eps, (char*) EPSSUBSPACE); CHKERRABORT(libMesh::COMM_WORLD,ierr); return;
    case LAPACK:
      ierr = EPSSetType (_eps, (char*) EPSLAPACK);   CHKERRABORT(libMesh::COMM_WORLD,ierr); return;
    case ARNOLDI:
      ierr = EPSSetType (_eps, (char*) EPSARNOLDI);  CHKERRABORT(libMesh::COMM_WORLD,ierr); return;
    case LANCZOS:
      ierr = EPSSetType (_eps, (char*) EPSLANCZOS);  CHKERRABORT(libMesh::COMM_WORLD,ierr); return;
#if !SLEPC_VERSION_LESS_THAN(2,3,2)
      // EPSKRYLOVSCHUR added in 2.3.2
    case KRYLOVSCHUR:
      ierr = EPSSetType (_eps, (char*) EPSKRYLOVSCHUR);  CHKERRABORT(libMesh::COMM_WORLD,ierr); return;
#endif
      // case ARPACK:
      // ierr = EPSSetType (_eps, (char*) EPSARPACK);   CHKERRABORT(libMesh::COMM_WORLD,ierr); return;

    default:
      libMesh::err << "ERROR:  Unsupported SLEPc Eigen Solver: "
		    << this->_eigen_solver_type         << std::endl
		    << "Continuing with SLEPc defaults" << std::endl;
    }
}




template <typename T>
void SlepcEigenSolver<T>:: set_slepc_problem_type()
{
  int ierr = 0;

  switch (this->_eigen_problem_type)
    {
    case NHEP:
      ierr = EPSSetProblemType (_eps, EPS_NHEP);  CHKERRABORT(libMesh::COMM_WORLD,ierr); return;
    case GNHEP:
      ierr = EPSSetProblemType (_eps, EPS_GNHEP); CHKERRABORT(libMesh::COMM_WORLD,ierr); return;
    case HEP:
      ierr = EPSSetProblemType (_eps, EPS_HEP);   CHKERRABORT(libMesh::COMM_WORLD,ierr); return;
    case GHEP:
      ierr = EPSSetProblemType (_eps, EPS_GHEP);  CHKERRABORT(libMesh::COMM_WORLD,ierr); return;

    default:
      libMesh::err << "ERROR:  Unsupported SLEPc Eigen Problem: "
		    << this->_eigen_problem_type        << std::endl
		    << "Continuing with SLEPc defaults" << std::endl;
    }
}



template <typename T>
void SlepcEigenSolver<T>:: set_slepc_position_of_spectrum()
{
  int ierr = 0;

  switch (this->_position_of_spectrum)
    {
    case LARGEST_MAGNITUDE:
      ierr = EPSSetWhichEigenpairs (_eps, EPS_LARGEST_MAGNITUDE);  CHKERRABORT(libMesh::COMM_WORLD,ierr); return;
    case SMALLEST_MAGNITUDE:
      ierr = EPSSetWhichEigenpairs (_eps, EPS_SMALLEST_MAGNITUDE); CHKERRABORT(libMesh::COMM_WORLD,ierr); return;
    case LARGEST_REAL:
      ierr = EPSSetWhichEigenpairs (_eps, EPS_LARGEST_REAL);       CHKERRABORT(libMesh::COMM_WORLD,ierr); return;
    case SMALLEST_REAL:
      ierr = EPSSetWhichEigenpairs (_eps, EPS_SMALLEST_REAL);      CHKERRABORT(libMesh::COMM_WORLD,ierr); return;
    case LARGEST_IMAGINARY:
      ierr = EPSSetWhichEigenpairs (_eps, EPS_LARGEST_IMAGINARY);  CHKERRABORT(libMesh::COMM_WORLD,ierr); return;
    case SMALLEST_IMAGINARY:
      ierr = EPSSetWhichEigenpairs (_eps, EPS_SMALLEST_IMAGINARY); CHKERRABORT(libMesh::COMM_WORLD,ierr); return;


    default:
      libMesh::err << "ERROR:  Unsupported SLEPc position of spectrum: "
		    << this->_position_of_spectrum        << std::endl;
      libmesh_error();
    }
}






template <typename T>
std::pair<Real, Real> SlepcEigenSolver<T>::get_eigenpair(unsigned int i,
							 NumericVector<T> &solution_in)
{
  int ierr=0;

  PetscReal re, im;

  // Make sure the NumericVector passed in is really a PetscVector
  PetscVector<T>* solution = libmesh_cast_ptr<PetscVector<T>*>(&solution_in);

  // real and imaginary part of the ith eigenvalue.
  PetscScalar kr, ki;

  solution->close();

  ierr = EPSGetEigenpair(_eps, i, &kr, &ki, solution->vec(), PETSC_NULL);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

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
std::pair<Real, Real> SlepcEigenSolver<T>::get_eigenvalue(unsigned int i)
{
  int ierr=0;

  PetscReal re, im;

  // real and imaginary part of the ith eigenvalue.
  PetscScalar kr, ki;

  ierr = EPSGetEigenvalue(_eps, i, &kr, &ki);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

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
  int ierr=0;
  PetscReal error;

  ierr = EPSComputeRelativeError(_eps, i, &error);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

  return error;
}


template <typename T>
void SlepcEigenSolver<T>::attach_deflation_space(NumericVector<T>& deflation_vector_in)
{
  this->init();

  int ierr = 0;
  Vec deflation_vector = (libmesh_cast_ptr<PetscVector<T>*>(&deflation_vector_in))->vec();
  Vec* deflation_space = &deflation_vector;
#if SLEPC_VERSION_LESS_THAN(3,1,0)
  ierr = EPSAttachDeflationSpace(_eps, 1, deflation_space, PETSC_FALSE);
#else
  ierr = EPSSetDeflationSpace(_eps, 1, deflation_space);
#endif
  CHKERRABORT(libMesh::COMM_WORLD,ierr);
}

template <typename T>
PetscErrorCode SlepcEigenSolver<T>::_petsc_shell_matrix_mult(Mat mat, Vec arg, Vec dest)
{
  /* Get the matrix context.  */
  int ierr=0;
  void* ctx;
  ierr = MatShellGetContext(mat,&ctx);
  CHKERRABORT(libMesh::COMM_WORLD,ierr);

  /* Get user shell matrix object.  */
  const ShellMatrix<T>& shell_matrix = *static_cast<const ShellMatrix<T>*>(ctx);

  /* Make \p NumericVector instances around the vectors.  */
  PetscVector<T> arg_global(arg, shell_matrix.communicator());
  PetscVector<T> dest_global(dest, shell_matrix.communicator());

  /* Call the user function.  */
  shell_matrix.vector_mult(dest_global,arg_global);

  return ierr;
}

template <typename T>
PetscErrorCode SlepcEigenSolver<T>::_petsc_shell_matrix_get_diagonal(Mat mat, Vec dest)
{
  /* Get the matrix context.  */
  int ierr=0;
  void* ctx;
  ierr = MatShellGetContext(mat,&ctx);
  CHKERRABORT(libMesh::COMM_WORLD,ierr);

  /* Get user shell matrix object.  */
  const ShellMatrix<T>& shell_matrix = *static_cast<const ShellMatrix<T>*>(ctx);

  /* Make \p NumericVector instances around the vector.  */
  PetscVector<T> dest_global(dest, shell_matrix.communicator());

  /* Call the user function.  */
  shell_matrix.get_diagonal(dest_global);

  return ierr;
}

//------------------------------------------------------------------
// Explicit instantiations
template class SlepcEigenSolver<Number>;

} // namespace libMesh



#endif // #ifdef LIBMESH_HAVE_SLEPC
