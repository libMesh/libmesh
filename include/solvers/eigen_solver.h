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



#ifndef LIBMESH_EIGEN_SOLVER_H
#define LIBMESH_EIGEN_SOLVER_H


#include "libmesh/libmesh_config.h"
#ifdef LIBMESH_HAVE_SLEPC

// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/enum_solver_package.h"
#include "libmesh/enum_eigen_solver_type.h"
#include "libmesh/reference_counted_object.h"
#include "libmesh/libmesh.h"
#include "libmesh/parallel_object.h"
#include "libmesh/auto_ptr.h"

// C++ includes

namespace libMesh
{

// forward declarations
template <typename T> class SparseMatrix;
template <typename T> class ShellMatrix;
template <typename T> class NumericVector;
class SolverConfiguration;

/**
 * This class provides an interface to solvers for eigenvalue
 * problems.
 */

template <typename T>
class EigenSolver : public ReferenceCountedObject<EigenSolver<T> >,
                    public ParallelObject
{
public:

  /**
   *  Constructor. Initializes Solver data structures
   */
  EigenSolver (const Parallel::Communicator & comm_in
               LIBMESH_CAN_DEFAULT_TO_COMMWORLD);

  /**
   * Destructor.
   */
  virtual ~EigenSolver ();

  /**
   * Builds an \p EigenSolver using the linear solver package specified by
   * \p solver_package
   */
  static UniquePtr<EigenSolver<T> > build(const Parallel::Communicator & comm_in
                                          LIBMESH_CAN_DEFAULT_TO_COMMWORLD,
                                          const SolverPackage solver_package = SLEPC_SOLVERS);

  /**
   * @returns true if the data structures are
   * initialized, false otherwise.
   */
  bool initialized () const { return _is_initialized; }


  /**
   * Release all memory and clear data structures.
   */
  virtual void clear () {}

  /**
   * Initialize data structures if not done so already.
   */
  virtual void init () = 0;

  /**
   * Returns the type of eigensolver to use.
   */
  EigenSolverType eigen_solver_type () const { return _eigen_solver_type; }

  /**
   * Returns the type of the eigen problem.
   */
  EigenProblemType eigen_problem_type () const { return _eigen_problem_type;}

  /**
   * Returns the position of the spectrum to compute.
   */
  PositionOfSpectrum position_of_spectrum () const
  { return _position_of_spectrum;}

  /**
   * Sets the type of eigensolver to use.
   */
  void set_eigensolver_type (const EigenSolverType est)
  { _eigen_solver_type = est; }

  /**
   * Sets the type of the eigenproblem.
   */
  void set_eigenproblem_type ( EigenProblemType ept)
  {_eigen_problem_type = ept;}

  /**
   * Sets the position of the spectrum.
   */
  void set_position_of_spectrum (PositionOfSpectrum pos)
  {_position_of_spectrum= pos;}

  /**
   * Solves the standard eigen problem when matrix_A is a
   * \p SparseMatrix, and returns the number of converged
   * eigenpairs and the number of iterations.
   */
  virtual std::pair<unsigned int, unsigned int> solve_standard (SparseMatrix<T> & matrix_A,
                                                                int nev,
                                                                int ncv,
                                                                const double tol,
                                                                const unsigned int m_its) = 0;

  /**
   * Solves the standard eigen problem when matrix_A is a
   * \p ShellMatrix, and returns the number of converged
   * eigenpairs and the number of iterations.
   */
  virtual std::pair<unsigned int, unsigned int> solve_standard (ShellMatrix<T> & matrix_A,
                                                                int nev,
                                                                int ncv,
                                                                const double tol,
                                                                const unsigned int m_its) = 0;


  /**
   * Solves the generalized eigen problem when both matrix_A
   * and matrix_B are of type \p SparseMatrix and returns the
   * number of converged eigenpairs and the number
   * of iterations.
   */
  virtual std::pair<unsigned int, unsigned int> solve_generalized (SparseMatrix<T> & matrix_A,
                                                                   SparseMatrix<T> & matrix_B,
                                                                   int nev,
                                                                   int ncv,
                                                                   const double tol,
                                                                   const unsigned int m_its) = 0;

  /**
   * Solves the generalized eigen problem when matrix_A is
   * a ShellMatrix and matrix_B is a SparseMatrix.
   */
  virtual std::pair<unsigned int, unsigned int> solve_generalized (ShellMatrix<T> & matrix_A,
                                                                   SparseMatrix<T> & matrix_B,
                                                                   int nev,
                                                                   int ncv,
                                                                   const double tol,
                                                                   const unsigned int m_its) = 0;

  /**
   * Solves the generalized eigen problem when matrix_A is
   * a SparseMatrix and matrix_B is a ShellMatrix.
   */
  virtual std::pair<unsigned int, unsigned int> solve_generalized (SparseMatrix<T> & matrix_A,
                                                                   ShellMatrix<T> & matrix_B,
                                                                   int nev,
                                                                   int ncv,
                                                                   const double tol,
                                                                   const unsigned int m_its) = 0;

  /**
   * Solves the generalized eigen problem when both matrix_A
   * and matrix_B are of type ShellMatrix.
   */
  virtual std::pair<unsigned int, unsigned int> solve_generalized (ShellMatrix<T> & matrix_A,
                                                                   ShellMatrix<T> & matrix_B,
                                                                   int nev,
                                                                   int ncv,
                                                                   const double tol,
                                                                   const unsigned int m_its) = 0;


  /**
   * Returns the \p ith eigenvalue (real and imaginary part),
   * and copies the \ ith eigen vector to the solution vector.
   */
  virtual std::pair<Real, Real> get_eigenpair (unsigned int i,
                                               NumericVector<T> & solution) = 0;

  /**
   * Returns the \p ith eigenvalue (real and imaginary part).
   * Same as above function, except it does copy the eigenvector.
   */
  virtual std::pair<Real, Real> get_eigenvalue (unsigned int i) = 0;

  /**
   * Attach a deflation space defined by a single vector.
   */
  virtual void attach_deflation_space(NumericVector<T> & deflation_vector) = 0;

  /**
   * Set the solver configuration object.
   */
  void set_solver_configuration(SolverConfiguration & solver_configuration);

protected:

  /**
   * Enum stating which type of eigensolver to use.
   */
  EigenSolverType _eigen_solver_type;

  /**
   * Enum stating which type of eigen problem we deal with.
   */
  EigenProblemType _eigen_problem_type;

  /**
   * Enum stating where to evaluate the spectrum.
   */
  PositionOfSpectrum _position_of_spectrum;

  /**
   * Flag indicating if the data structures have been initialized.
   */
  bool _is_initialized;

  /**
   * Optionally store a SolverOptions object that can be used
   * to set parameters like solver type, tolerances and iteration limits.
   */
  SolverConfiguration * _solver_configuration;

};




/*----------------------- inline functions ----------------------------------*/
template <typename T>
inline
EigenSolver<T>::EigenSolver (const Parallel::Communicator & comm_in) :
  ParallelObject(comm_in),
  _eigen_solver_type    (ARNOLDI),
  _eigen_problem_type   (NHEP),
  _position_of_spectrum (LARGEST_MAGNITUDE),
  _is_initialized       (false),
  _solver_configuration(libmesh_nullptr)
{
}



template <typename T>
inline
EigenSolver<T>::~EigenSolver ()
{
  this->clear ();
}

} // namespace libMesh

#endif // LIBMESH_HAVE_SLEPC

#endif // LIBMESH_EIGEN_SOLVER_H
