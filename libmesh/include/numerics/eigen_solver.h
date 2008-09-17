// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
  
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
  
// This library is distributd in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA



#ifndef __eigen_solver_h__
#define __eigen_solver_h__


#include "libmesh_config.h"
#ifdef LIBMESH_HAVE_SLEPC

// C++ includes

// Local includes
#include "libmesh_common.h"
#include "enum_solver_package.h"
#include "enum_eigen_solver_type.h"
#include "reference_counted_object.h"
#include "libmesh.h"

// forward declarations
template <typename T> class AutoPtr;
template <typename T> class SparseMatrix;
template <typename T> class NumericVector;


/**
 * This class provides an interface to solvers for eigenvalue
 * problems.
 */

template <typename T>
class EigenSolver : public ReferenceCountedObject<EigenSolver<T> >
{
public:
  
  /**
   *  Constructor. Initializes Solver data structures
   */
  EigenSolver ();
    
  /**
   * Destructor.
   */
  virtual ~EigenSolver ();
  
  /**
   * Builds a \p EigenSolver using the linear solver package specified by
   * \p solver_package
   */
  static AutoPtr<EigenSolver<T> > build(const SolverPackage solver_package =
					SLEPC_SOLVERS);
  
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
  PositionOfSpectrum postition_of_spectrum () const
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
   * Solves the standard eigen problem and returns the
   * number of converged eigenpairs and the number
   * of iterations.
   */
  virtual std::pair<unsigned int, unsigned int> solve_standard (SparseMatrix<T> &matrix_A,  
								int nev,
								int ncv,
								const double tol,
								const unsigned int m_its) = 0;


  
  /**
   * Solves the generalized eigen problem and returns the
   * number of converged eigenpairs and the number
   * of iterations.
   */
  virtual std::pair<unsigned int, unsigned int> solve_generalized (SparseMatrix<T> &matrix_A, 
								   SparseMatrix<T> &matrix_B,  
								   int nev,
								   int ncv,
								   const double tol,
								   const unsigned int m_its) = 0;



  /**
   * Returns the \p ith eigenvalue (real and imaginary part),
   * and copies the \ ith eigen vector to the solution vector.
   */
  virtual std::pair<Real, Real> get_eigenpair (unsigned int i,
					       NumericVector<T> &solution) = 0;

  
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

};




/*----------------------- inline functions ----------------------------------*/
template <typename T>
inline
EigenSolver<T>::EigenSolver () :
  
  _eigen_solver_type    (ARNOLDI),
  _eigen_problem_type   (NHEP),
  _position_of_spectrum (LARGEST_MAGNITUDE),
  _is_initialized       (false)
{
}



template <typename T>
inline
EigenSolver<T>::~EigenSolver ()
{
  this->clear ();
}

#endif // LIBMESH_HAVE_SLEPC

#endif // #ifdef __eigen_solver_h__
