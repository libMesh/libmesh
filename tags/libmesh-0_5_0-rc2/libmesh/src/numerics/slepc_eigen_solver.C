// $Id: slepc_eigen_solver.C,v 1.2 2005-05-11 23:12:10 benkirk Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2005  Benjamin S. Kirk, John W. Peterson
  
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



#include "libmesh_common.h"

#if defined(HAVE_SLEPC) && defined(HAVE_PETSC)


// C++ includes

// Local Includes
#include "slepc_eigen_solver.h"



/*----------------------- functions ----------------------------------*/
template <typename T>
void SlepcEigenSolver<T>::clear ()
{
  if (this->initialized())
    {
      this->_is_initialized = false;

      int ierr=0;

      ierr = EPSDestroy(_eps);
             CHKERRABORT(libMesh::COMM_WORLD,ierr);

      // SLEPc default eigenproblem solver
      this->_eigen_solver_type = ARNOLDI;
 
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

      // Set modified Gram-Schmidt orthogonalization as default
      // and leave other parameters unchanged
      EPSOrthogonalizationRefinementType refinement;
      PetscReal eta;
      ierr = EPSGetOrthogonalization (_eps, PETSC_NULL, &refinement, &eta);
      ierr = EPSSetOrthogonalization (_eps, EPS_MGS_ORTH, refinement, eta);
             CHKERRABORT(libMesh::COMM_WORLD,ierr);

      // Set user-specified  solver
      set_slepc_solver_type();
      
    } 
}    



template <typename T>
std::pair<unsigned int, unsigned int> 
SlepcEigenSolver<T>::solve (SparseMatrix<T> &matrix_A_in,
			    int nev,                  // number of requested eigenpairs
			    int ncv,                  // number of basis vectors
			    const double tol,         // solver tolerance
			    const unsigned int m_its) // maximum number of iterations
{

  this->init ();
  
  PetscMatrix<T>* matrix_A   = dynamic_cast<PetscMatrix<T>*>(&matrix_A_in);

  int ierr=0;

  // converged eigen pairs and number of iterations
  int nconv=0;
  int its=0;

  // The relative error.
  PetscReal error, re, im;

  // Pointer to vectors of the real parts, imaginary parts.
  PetscScalar kr, ki;

  // Close the matrix and vectors in case this wasn't already done.
  matrix_A->close ();

  // just for debugging, remove this 
//   char mat_file[] = "matA.petsc";
//   PetscViewer petsc_viewer;
//   ierr = PetscViewerBinaryOpen(libMesh::COMM_WORLD, mat_file, PETSC_FILE_CREATE, &petsc_viewer);
//          CHKERRABORT(libMesh::COMM_WORLD,ierr);
//   ierr = MatView(matrix_A->mat(),petsc_viewer);
//          CHKERRABORT(libMesh::COMM_WORLD,ierr);

  // Set operators.
  ierr = EPSSetOperators (_eps, matrix_A->mat(), PETSC_NULL);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

  // Set eigenvalues to be computed.
  ierr = EPSSetDimensions (_eps, nev, ncv);
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

#ifdef USE_COMPLEX_NUMBERS
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
      
    default:
      std::cerr << "ERROR:  Unsupported SLEPc Eigen Solver: "
		<< this->_eigen_solver_type         << std::endl
		<< "Continuing with SLEPc defaults" << std::endl;
    }
}

template <typename T>
std::pair<Real, Real> SlepcEigenSolver<T>::get_eigenpair(unsigned int i,
							 NumericVector<T> &solution_in)
{
  int ierr=0;

  PetscReal re, im;

  PetscVector<T>* solution = dynamic_cast<PetscVector<T>*>(&solution_in);

  // real and imaginary part of the ith eigenvalue.
  PetscScalar kr, ki;

  solution->close();

  ierr = EPSGetEigenpair(_eps, i, &kr, &ki, solution->vec(), PETSC_NULL);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

#ifdef USE_COMPLEX_NUMBERS
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

//------------------------------------------------------------------
// Explicit instantiations
template class SlepcEigenSolver<Number>;
 


#endif // #ifdef HAVE_SLEPC
