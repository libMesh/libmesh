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



#ifndef __nonlinear_solver_h__
#define __nonlinear_solver_h__


// C++ includes

// Local includes
#include "libmesh_common.h"
#include "enum_solver_package.h"
#include "reference_counted_object.h"
#include "libmesh.h"

// forward declarations
template <typename T> class AutoPtr;
template <typename T> class SparseMatrix;
template <typename T> class NumericVector;
class NonlinearImplicitSystem;





/**
 * This class provides a uniform interface for nonlinear solvers.  This base
 * class is overloaded to provide nonlinear solvers from different packages
 * like PETSC.
 *
 * @author Benjamin Kirk, 2005
 */

template <typename T>
class NonlinearSolver : public ReferenceCountedObject<NonlinearSolver<T> >
{
public:
  /**
   * The type of system
   */
  typedef NonlinearImplicitSystem sys_type;

  /**
   *  Constructor. Initializes Solver data structures
   */
  NonlinearSolver (sys_type& s);
    
  /**
   * Destructor.
   */
  virtual ~NonlinearSolver ();
  
  /**
   * Builds a \p NonlinearSolver using the nonlinear solver package specified by
   * \p solver_package
   */
  static AutoPtr<NonlinearSolver<T> > build(sys_type& s, const SolverPackage solver_package =
					    libMesh::default_solver_package());
  
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
   * Solves the nonlinear system.
   */
  virtual std::pair<unsigned int, Real> solve (SparseMatrix<T>&,  // System Jacobian Matrix
					       NumericVector<T>&, // Solution vector
					       NumericVector<T>&, // Residual vector
					       const double,      // Stopping tolerance
					       const unsigned int) = 0; // N. Iterations

  /**
   * Function that computes the residual \p R(X) of the nonlinear system
   * at the input iterate \p X.
   */
  void (* residual) (const NumericVector<Number>& X,
		     NumericVector<Number>& R);

  /**
   * Function that computes the Jacobian \p J(X) of the nonlinear system
   * at the input iterate \p X.
   */
  void (* jacobian) (const NumericVector<Number>& X,
		     SparseMatrix<Number>& J);

  /**
   * Function that computes either the residual \f$ R(X) \f$ or the
   * Jacobian \f$ J(X) \f$ of the nonlinear system at the input
   * iterate \f$ X \f$.  Note that either \p R or \p J could be
   * \p XSNULL.
   */
  void (* matvec) (const NumericVector<Number>& X,
		   NumericVector<Number>* R,
		   SparseMatrix<Number>*  J);

  /**
   * @returns a constant reference to the system we are solving.
   */
  const sys_type & system () const { return _system; }

  /**
   * @returns a writeable reference to the system we are solving.
   */
  sys_type & system () { return _system; }


  /**
   * Maximum number of non-linear iterations.
   */
  unsigned int max_nonlinear_iterations;

  /**
   * Maximum number of function evaluations.
   */
  unsigned int max_function_evaluations;
  
  /**
   * The NonlinearSolver should exit after the residual is
   * reduced to either less than absolute_residual_tolerance
   * or less than relative_residual_tolerance times the
   * initial residual.
   *
   * Users should increase any of these tolerances that they want to use for a
   * stopping condition.
   *
   */
  Real absolute_residual_tolerance;
  Real relative_residual_tolerance;

  /**
   * The NonlinearSolver should exit after the full nonlinear step norm is
   * reduced to either less than absolute_step_tolerance
   * or less than relative_step_tolerance times the largest
   * nonlinear solution which has been seen so far.
   *
   * Users should increase any of these tolerances that they want to use for a
   * stopping condition.
   *
   * Note that not all NonlinearSolvers support relative_step_tolerance!
   */
  Real absolute_step_tolerance;
  Real relative_step_tolerance;

  /**
   * Each linear solver step should exit after \p max_linear_iterations
   * is exceeded.
   */
  unsigned int max_linear_iterations;

  /**
   * Any required linear solves will at first be done with this tolerance;
   * the NonlinearSolver may tighten the tolerance for later solves.
   */
  Real initial_linear_tolerance;

  /**
   * The tolerance for linear solves is kept above this minimum
   */
  Real minimum_linear_tolerance;

protected:
  /**
   * A reference to the system we are solving.
   */
  sys_type& _system;
  
  /**
   * Flag indicating if the data structures have been initialized.
   */
  bool _is_initialized;
};




/*----------------------- inline functions ----------------------------------*/
template <typename T>
inline
NonlinearSolver<T>::NonlinearSolver (sys_type& s) :
  residual        (NULL),
  jacobian        (NULL),
  matvec          (NULL),
  _system(s),
  _is_initialized (false),
  max_nonlinear_iterations(0),
  max_function_evaluations(0),
  absolute_residual_tolerance(0),
  relative_residual_tolerance(0),
  absolute_step_tolerance(0),
  relative_step_tolerance(0),
  max_linear_iterations(0),
  initial_linear_tolerance(0),
  minimum_linear_tolerance(0)
{
}



template <typename T>
inline
NonlinearSolver<T>::~NonlinearSolver ()
{
  this->clear ();
}



#endif // #ifdef __solver_h__
