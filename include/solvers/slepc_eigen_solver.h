// The libMesh Finite Element Library.
// Copyright (C) 2002-2016 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_SLEPC_EIGEN_SOLVER_H
#define LIBMESH_SLEPC_EIGEN_SOLVER_H

#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_HAVE_SLEPC

// Local includes
#include "libmesh/eigen_solver.h"
#include "libmesh/slepc_macro.h"

/**
 * SLEPc include files.
 */
EXTERN_C_FOR_SLEPC_BEGIN
# include <slepceps.h>
EXTERN_C_FOR_SLEPC_END

// C++ includes


namespace libMesh
{


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
  SlepcEigenSolver(const Parallel::Communicator & comm_in
                   LIBMESH_CAN_DEFAULT_TO_COMMWORLD);


  /**
   * Destructor.
   */
  ~SlepcEigenSolver();


  /**
   * Release all memory and clear data structures.
   */
  virtual void clear() libmesh_override;


  /**
   * Initialize data structures if not done so already.
   */
  virtual void init() libmesh_override;


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
  virtual std::pair<unsigned int, unsigned int>
  solve_standard (SparseMatrix<T> & matrix_A,
                  int nev,
                  int ncv,
                  const double tol,
                  const unsigned int m_its) libmesh_override;

  /**
   * Same as above except that matrix_A is a ShellMatrix
   * in this case.
   */
  virtual std::pair<unsigned int, unsigned int>
  solve_standard (ShellMatrix<T> & shell_matrix,
                  int nev,
                  int ncv,
                  const double tol,
                  const unsigned int m_its) libmesh_override;


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
  virtual std::pair<unsigned int, unsigned int>
  solve_generalized(SparseMatrix<T> & matrix_A,
                    SparseMatrix<T> & matrix_B,
                    int nev,
                    int ncv,
                    const double tol,
                    const unsigned int m_its) libmesh_override;

  /**
   * Solve generalized eigenproblem when matrix_A is of
   * type ShellMatrix, matrix_B is of type SparseMatrix.
   */
  virtual std::pair<unsigned int, unsigned int>
  solve_generalized(ShellMatrix<T> & matrix_A,
                    SparseMatrix<T> & matrix_B,
                    int nev,
                    int ncv,
                    const double tol,
                    const unsigned int m_its) libmesh_override;

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
  virtual std::pair<unsigned int, unsigned int>
  solve_generalized(SparseMatrix<T> & matrix_A,
                    ShellMatrix<T> & matrix_B,
                    int nev,
                    int ncv,
                    const double tol,
                    const unsigned int m_its) libmesh_override;

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
  virtual std::pair<unsigned int, unsigned int>
  solve_generalized(ShellMatrix<T> & matrix_A,
                    ShellMatrix<T> & matrix_B,
                    int nev,
                    int ncv,
                    const double tol,
                    const unsigned int m_its) libmesh_override;



  /**
   * This function returns the real and imaginary part of the
   * ith eigenvalue and copies the respective eigenvector to the
   * solution vector. Note that also in case of purely real matrix
   * entries the eigenpair may be complex values.
   */
  virtual std::pair<Real, Real>
  get_eigenpair (dof_id_type i,
                 NumericVector<T> & solution_in) libmesh_override;

  /**
   * Same as above, but does not copy the eigenvector.
   */
  virtual std::pair<Real, Real>
  get_eigenvalue (dof_id_type i) libmesh_override;

  /**
   * @returns the relative error ||A*x-lambda*x||/|lambda*x|
   * of the ith eigenpair. (or the equivalent for a general eigenvalue problem)
   */
  Real get_relative_error (unsigned int i);

  /**
   * Attach a deflation space defined by a single vector.
   */
  virtual void attach_deflation_space(NumericVector<T> & deflation_vector) libmesh_override;

  /**
   * Provide one basis vector for the initial guess
   */
  virtual void
  set_initial_space(NumericVector<T> & initial_space_in) libmesh_override;

  /**
   * Set the spectral transformation to be used
   */
  virtual void 
  set_spectral_transform(SpectralTransform st) libmesh_override;

  /**
   * Returns the raw SLEPc eps context pointer.
   */
  EPS eps() { this->init(); return _eps; }

private:

  /**
   * Helper function that actually performs the standard eigensolve.
   */
  std::pair<unsigned int, unsigned int> _solve_standard_helper (Mat mat,
                                                                int nev,
                                                                int ncv,
                                                                const double tol,
                                                                const unsigned int m_its);

  /**
   * Helper function that actually performs the generalized eigensolve.
   */
  std::pair<unsigned int, unsigned int> _solve_generalized_helper (Mat mat_A,
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
SlepcEigenSolver<T>::SlepcEigenSolver (const Parallel::Communicator & comm_in) :
  EigenSolver<T>(comm_in)
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

} // namespace libMesh


#endif // #ifdef LIBMESH_HAVE_SLEPC
#endif // LIBMESH_SLEPC_EIGEN_SOLVER_H
