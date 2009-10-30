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



#ifndef __slepc_eigen_solver_h__
#define __slepc_eigen_solver_h__

// C++ includes

// Local includes
#include "eigen_solver.h"
#include "petsc_vector.h"
#include "petsc_matrix.h"


/**
 * SLEPc include files. SLEPs can only be used
 * together with PETSc.
 */
#if defined(LIBMESH_HAVE_SLEPC) && defined(LIBMESH_HAVE_PETSC)

#ifndef LIBMESH_USE_COMPLEX_NUMBERS
extern "C"
{
# include <slepceps.h>
# include <slepcversion.h>
}
#else
# include <slepceps.h>
# include <slepcversion.h>
#endif


// A convenient macro for comparing SLEPc versions.
// Returns 1 if the current SLEPc version is < major.minor.subminor
// and zero otherwise.
#define SLEPC_VERSION_LESS_THAN(major,minor,subminor)			            \
  ((SLEPC_VERSION_MAJOR < (major) ||						    \
    (SLEPC_VERSION_MAJOR == (major) && (SLEPC_VERSION_MINOR < (minor) ||	    \
				  (SLEPC_VERSION_MINOR == (minor) &&		    \
				   SLEPC_VERSION_SUBMINOR < (subminor))))) ? 1 : 0)



/**
 * This class provides an interface to the SLEPc
 * eigenvalue solver library \p www.grycap.upv.es/slepc/.
 */

template <typename T>
class SlepcEigenSolver : public EigenSolver<T>
{

public:

  /**
   *  Constructor. Initializes Petsc data structures
   */
  SlepcEigenSolver();
    

  /**
   * Destructor.
   */
  ~SlepcEigenSolver();

  
  /**
   * Release all memory and clear data structures.
   */
  void clear();


  /**
   * Initialize data structures if not done so already.
   */
  void init();

       
  /**
   * This function calls the SLEPc solver to compute
   * the eigenpairs of the SparseMatrix matrix_A. \p nev is
   * the number of eigenpairs to be computed and
   * \p ncv is the number of basis vectors to be
   * used in the solution procedure. Return values
   * are the number of converged eigen values and the
   * number of the iterations carried out by the eigen
   * solver.
   */    
  std::pair<unsigned int, unsigned int>  solve_standard (SparseMatrix<T> &matrix_A,
							 int nev,
							 int ncv,
							 const double tol,
							 const unsigned int m_its);

  /**
   * Same as above except that matrix_A is a ShellMatrix
   * in this case.
   */
  std::pair<unsigned int, unsigned int>  solve_standard (ShellMatrix<T> &shell_matrix,
							 int nev,
							 int ncv,
							 const double tol,
							 const unsigned int m_its);
  

 /**
   * This function calls the SLEPc solver to compute
   * the eigenpairs for the generalized eigenproblem
   * defined by the matrix_A and matrix_B,
   * which are of type SparseMatrix. The argument 
   * \p nev is the number of eigenpairs to be computed
   * and \p ncv is the number of basis vectors to be
   * used in the solution procedure. Return values
   * are the number of converged eigen values and the
   * number of the iterations carried out by the eigen
   * solver.
   */
  std::pair<unsigned int, unsigned int>  solve_generalized(SparseMatrix<T> &matrix_A,
							   SparseMatrix<T> &matrix_B,
							   int nev,
							   int ncv,
							   const double tol,
							   const unsigned int m_its);

  /**
   * Solve generalized eigenproblem when matrix_A is of
   * type ShellMatrix, matrix_B is of type SparseMatrix.
   */
    std::pair<unsigned int, unsigned int>  solve_generalized(ShellMatrix<T> &matrix_A,
  							   SparseMatrix<T> &matrix_B,
  							   int nev,
  							   int ncv,
  							   const double tol,
  							   const unsigned int m_its);

  /**
   * Solve generalized eigenproblem when matrix_A is of
   * type SparseMatrix, matrix_B is of type ShellMatrix.
   * When using this function, one should use the
   * command line options:
   * -st_ksp_type gmres -st_pc_type none
   * or
   * -st_ksp_type gmres -st_pc_type jacobi
   * or similar.
   */
    std::pair<unsigned int, unsigned int>  solve_generalized(SparseMatrix<T> &matrix_A,
  							   ShellMatrix<T> &matrix_B,
  							   int nev,
  							   int ncv,
  							   const double tol,
  							   const unsigned int m_its);

  /**
   * Solve generalized eigenproblem when both matrix_A and
   * matrix_B are of type ShellMatrix.
   * When using this function, one should use the
   * command line options:
   * -st_ksp_type gmres -st_pc_type none
   * or
   * -st_ksp_type gmres -st_pc_type jacobi
   * or similar.
   */
    std::pair<unsigned int, unsigned int>  solve_generalized(ShellMatrix<T> &matrix_A,
  							   ShellMatrix<T> &matrix_B,
  							   int nev,
  							   int ncv,
  							   const double tol,
  							   const unsigned int m_its);

 
 
  /**
   * This function returns the real and imaginary part of the
   * ith eigenvalue and copies the respective eigenvector to the
   * solution vector. Note that also in case of purely real matrix
   * entries the eigenpair may be complex values.
   */                
  std::pair<Real, Real> get_eigenpair (unsigned int i,
				       NumericVector<T> &solution_in);

  /**
   * @computes and returns the relative error ||A*x-lambda*x||/|lambda*x|
   * of the ith eigenpair. (or the equivalent for a general eigenvalue problem)
   */                
  Real get_relative_error (unsigned int i);

  /**
   * Attach a deflation space defined by a single vector.
   */
  void attach_deflation_space(NumericVector<T>& deflation_vector);

private:

  /**
   * Helper function that actually performs the standard eigensolve.
   */
  std::pair<unsigned int, unsigned int>  _solve_standard_helper (Mat mat,
							 int nev,
							 int ncv,
							 const double tol,
							 const unsigned int m_its);

  /**
   * Helper function that actually performs the generalized eigensolve.
   */
  std::pair<unsigned int, unsigned int>  _solve_generalized_helper (Mat mat_A,
							         Mat mat_B,
							         int nev,
							         int ncv,
							         const double tol,
							         const unsigned int m_its);

  /**
   * Tells Slepc to use the user-specified solver stored in
   * \p _eigen_solver_type
   */
  void set_slepc_solver_type ();

  /**
   * Tells Slepc to deal with the type of problem stored in
   * \p _eigen_problem_type
   */
  void set_slepc_problem_type ();

  /**
   * Tells Slepc to compute the spectrum at the position
   * stored in \p _position_of_spectrum
   */
  void set_slepc_position_of_spectrum();

  /**
   * Internal function if shell matrix mode is used, this just
   * calls the shell matrix's matrix multiplication function.
   * See PetscLinearSolver for a similar implementation.
   */
  static PetscErrorCode _petsc_shell_matrix_mult(Mat mat, Vec arg, Vec dest);

  /**
   * Internal function if shell matrix mode is used, this just
   * calls the shell matrix's get_diagonal function.
   * Required in order to use Jacobi preconditioning.
   */
  static PetscErrorCode _petsc_shell_matrix_get_diagonal(Mat mat, Vec dest);

  /**
   * Eigenproblem solver context
   */
  EPS _eps;

};


/*----------------------- inline functions ----------------------------------*/
template <typename T>
inline
SlepcEigenSolver<T>::SlepcEigenSolver ()
{
  this->_eigen_solver_type  = ARNOLDI;
  this->_eigen_problem_type = NHEP;
}



template <typename T>
inline
SlepcEigenSolver<T>::~SlepcEigenSolver ()
{
  this->clear ();
}


#endif // #ifdef LIBMESH_HAVE_SLEPC
#endif // #ifdef __slepc_linear_solver_h__
