// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2007  Benjamin S. Kirk, John W. Peterson
  
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



#ifndef __linear_solver_h__
#define __linear_solver_h__


// C++ includes

// Local includes
#include "libmesh_common.h"
#include "enum_solver_package.h"
#include "enum_solver_type.h"
#include "enum_preconditioner_type.h"
#include "reference_counted_object.h"
#include "libmesh.h"

// forward declarations
template <typename T> class AutoPtr;
template <typename T> class SparseMatrix;
template <typename T> class NumericVector;





/**
 * This class provides a uniform interface for linear solvers.  This base
 * class is overloaded to provide linear solvers from different packages
 * like PETSC or LASPACK.
 *
 * @author Benjamin Kirk, 2003
 */

template <typename T>
class LinearSolver : public ReferenceCountedObject<LinearSolver<T> >
{
public:
  
  /**
   *  Constructor. Initializes Solver data structures
   */
  LinearSolver ();
    
  /**
   * Destructor.
   */
  virtual ~LinearSolver ();
  
  /**
   * Builds a \p LinearSolver using the linear solver package specified by
   * \p solver_package
   */
  static AutoPtr<LinearSolver<T> > build(const SolverPackage solver_package =
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
   * Returns the type of solver to use.
   */
  SolverType solver_type () const { return _solver_type; }

  /**
   * Sets the type of solver to use.
   */
  void set_solver_type (const SolverType st)
  { _solver_type = st; }

  /**
   * Returns the type of preconditioner to use.
   */
  PreconditionerType preconditioner_type () const
  { return _preconditioner_type; }

  /**
   * Sets the type of preconditioner to use.
   */
  void set_preconditioner_type (const PreconditionerType pct)
  { _preconditioner_type = pct; }
  


  /**
   * This function calls the solver
   * "_solver_type" preconditioned with the
   * "_preconditioner_type" preconditioner.  Note that this method
   * will compute the preconditioner from the system matrix.
   */
  virtual std::pair<unsigned int, Real> solve (SparseMatrix<T>&,  // System Matrix
					       NumericVector<T>&, // Solution vector
					       NumericVector<T>&, // RHS vector
					       const double,      // Stopping tolerance
					       const unsigned int) = 0; // N. Iterations
  


  /**
   * This function calls the solver
   * "_solver_type" preconditioned with the
   * "_preconditioner_type" preconditioner.  Note that this method
   * will compute the preconditioner from the system matrix.
   */
  virtual std::pair<unsigned int, Real> solve (SparseMatrix<T>&,  // System Matrix
					       SparseMatrix<T>&,  // Preconditioning Matrix
					       NumericVector<T>&, // Solution vector
					       NumericVector<T>&, // RHS vector
					       const double,      // Stopping tolerance
					       const unsigned int) = 0; // N. Iterations
  
  /**
   * Prints a useful message about why the latest linear solve
   * con(di)verged.
   */
  virtual void print_converged_reason() = 0;
  
protected:

  
  /**
   * Enum stating which type of iterative solver to use.
   */
  SolverType _solver_type;

  /**
   * Enum statitng with type of preconditioner to use.
   */
  PreconditionerType _preconditioner_type;
  
  /**
   * Flag indicating if the data structures have been initialized.
   */
  bool _is_initialized;
};




/*----------------------- inline functions ----------------------------------*/
template <typename T>
inline
LinearSolver<T>::LinearSolver () :
  
  _solver_type         (GMRES),
  _preconditioner_type (ILU_PRECOND),
  _is_initialized      (false)
{
}



template <typename T>
inline
LinearSolver<T>::~LinearSolver ()
{
  this->clear ();
}



#endif // #ifdef __solver_h__
