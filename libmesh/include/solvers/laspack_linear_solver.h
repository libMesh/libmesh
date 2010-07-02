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



#ifndef __laspack_linear_solver_h__
#define __laspack_linear_solver_h__

#include "libmesh_common.h"

#if defined(LIBMESH_HAVE_LASPACK)
//#if defined(LIBMESH_HAVE_LASPACK) && !defined(LIBMESH_USE_COMPLEX_NUMBERS)


// C++ includes

// Local includes
#include "linear_solver.h"
#include "laspack_vector.h"
#include "laspack_matrix.h"


#include <itersolv.h>
#include <rtc.h>
#include <errhandl.h>

namespace libMesh
{



/**
 * This class provides an interface to Laspack
 * iterative solvers that is compatible with the \p libMesh
 * \p LinearSolver<>
 * 
 * @author Benjamin Kirk, 2002-2007
 */
template <typename T>
class LaspackLinearSolver : public LinearSolver<T>
{
 public:
  /**
   *  Constructor. Initializes Laspack data structures
   */
  LaspackLinearSolver ();
    
  /**
   * Destructor.
   */
  ~LaspackLinearSolver ();
  
  /**
   * Release all memory and clear data structures.
   */
  void clear ();

  /**
   * Initialize data structures if not done so already.
   */
  void init ();
  

  /**
   * Call the Laspack solver
   */    
  std::pair<unsigned int, Real> 
    solve (SparseMatrix<T>  &matrix,
	   NumericVector<T> &solution,
	   NumericVector<T> &rhs,
	   const double tol,
	   const unsigned int m_its);
  /**
   * Call the Laspack solver
   */    
  std::pair<unsigned int, Real> 
    solve (SparseMatrix<T>  &matrix,
	   SparseMatrix<T>  &pc,
	   NumericVector<T> &solution,
	   NumericVector<T> &rhs,
	   const double tol,
	   const unsigned int m_its);
   
  /**
   * This function solves a system whose matrix is a shell matrix.
   */
  std::pair<unsigned int, Real>
    solve (const ShellMatrix<T>& shell_matrix,
	   NumericVector<T>& solution_in,
	   NumericVector<T>& rhs_in,
	   const double tol,
	   const unsigned int m_its);
  
  /**
   * This function solves a system whose matrix is a shell matrix, but
   * a sparse matrix is used as preconditioning matrix, this allowing
   * other preconditioners than JACOBI.
   */
  virtual std::pair<unsigned int, Real>
    solve (const ShellMatrix<T>& shell_matrix,
	   const SparseMatrix<T>& precond_matrix,
	   NumericVector<T>& solution_in,
	   NumericVector<T>& rhs_in,
	   const double tol,
	   const unsigned int m_its);
  
  /**
   * Prints a useful message about why the latest linear solve
   * con(di)verged.
   */
  virtual void print_converged_reason();
  
 private:
  
  /**
   * Tells LASPACK to use the user-specified preconditioner stored in
   * \p _preconditioner_type
   */
  void set_laspack_preconditioner_type ();

  /**
   * Preconditioner type
   */
  PrecondProcType _precond_type;
};


/*----------------------- functions ----------------------------------*/
template <typename T>
inline
LaspackLinearSolver<T>::LaspackLinearSolver () :
  _precond_type (ILUPrecond)
{
}



template <typename T>
inline
LaspackLinearSolver<T>::~LaspackLinearSolver ()
{
  this->clear ();
}



template <typename T>
inline
std::pair<unsigned int, Real>
LaspackLinearSolver<T>::solve (SparseMatrix<T>&,
			       SparseMatrix<T>&,
			       NumericVector<T>&,
			       NumericVector<T>&,
			       const double,
			       const unsigned int)
{
  libMesh::err << "ERROR: LASPACK does not support a user-supplied preconditioner!"
	        << std::endl;
  libmesh_error();

  std::pair<unsigned int, Real> p;
  return p;
}

} // namespace libMesh

#endif // #ifdef LIBMESH_HAVE_LASPACK
#endif // #ifndef __laspack_linear_solver_h__
