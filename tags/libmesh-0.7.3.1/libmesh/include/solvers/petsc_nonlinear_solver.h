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



#ifndef __petsc_nonlinear_solver_h__
#define __petsc_nonlinear_solver_h__

// C++ includes

// Local includes
#include "libmesh_config.h"


// Petsc include files.
#ifdef LIBMESH_HAVE_PETSC

#include "nonlinear_solver.h"
#include "petsc_macro.h"

EXTERN_C_FOR_PETSC_BEGIN
# include <petscsnes.h>
EXTERN_C_FOR_PETSC_END

namespace libMesh
{
  // Allow users access to these functions in case they want to reuse them.  Note that users shouldn't
  // need access to these most of the time as they are used internally by this object.
  extern "C"
  {
    PetscErrorCode __libmesh_petsc_snes_monitor (SNES, PetscInt its, PetscReal fnorm, void *);
    PetscErrorCode __libmesh_petsc_snes_residual (SNES, Vec x, Vec r, void *ctx);
    PetscErrorCode __libmesh_petsc_snes_jacobian (SNES, Vec x, Mat *jac, Mat *pc, MatStructure *msflag, void *ctx);
  }

/**
 * This class provides an interface to PETSc
 * iterative solvers that is compatible with the \p libMesh
 * \p NonlinearSolver<>
 *
 * @author Benjamin Kirk, 2002-2007
 */

template <typename T>
class PetscNonlinearSolver : public NonlinearSolver<T>
{
public:
  /**
   * The type of system
   */
  typedef NonlinearImplicitSystem sys_type;

  /**
   *  Constructor. Initializes Petsc data structures
   */
  PetscNonlinearSolver (sys_type& system);

  /**
   * Destructor.
   */
  ~PetscNonlinearSolver ();

  /**
   * Release all memory and clear data structures.
   */
  virtual void clear ();

  /**
   * Initialize data structures if not done so already.
   */
  virtual void init ();

  /**
   * Returns the raw PETSc snes context pointer.
   */
  SNES snes() { this->init(); return _snes; }

  /**
   * Call the Petsc solver.  It calls the method below, using the
   * same matrix for the system and preconditioner matrices.
   */
  virtual std::pair<unsigned int, Real> solve (SparseMatrix<T>&,    // System Jacobian Matrix
					       NumericVector<T>&,   // Solution vector
					       NumericVector<T>&,   // Residual vector
					       const double,        // Stopping tolerance
					       const unsigned int); // N. Iterations

  /**
   * Prints a useful message about why the latest nonlinear solve
   * con(di)verged.
   */
  virtual void print_converged_reason();

  /**
   * Returns the currently-available (or most recently obtained, if the SNES object has
   * been destroyed) convergence reason.  Refer to PETSc docs for the meaning of different
   * SNESConvergedReasons.
   */
  SNESConvergedReason get_converged_reason();

  /**
   * Get the total number of linear iterations done in the last solve
   */
  virtual int get_total_linear_iterations();

  /**
   * If called *during* the solve(), for example by the user-specified
   * residual or Jacobian function, returns the current nonlinear iteration
   * number.
   */
  virtual unsigned get_current_nonlinear_iteration_number() const { return _current_nonlinear_iteration_number; }

  /**
   * This public setter is necessary since the value is computed in the
   * __libmesh_petsc_snes_residual()/jacobian() function and must be stored
   * somehow.
   */
  void set_current_nonlinear_iteration_number(unsigned num) { _current_nonlinear_iteration_number = num; }

private:

  /**
   * Nonlinear solver context
   */
  SNES _snes;

  /**
   * Store the reason for SNES convergence/divergence for use even after the _snes
   * has been cleared.  Note that print_converged_reason() will always *try* to
   * get the current reason with SNESGetConvergedReason(), but if the SNES object
   * has already been cleared, it will fall back on this stored value.  Note that
   * this value is therefore necessarily *not* cleared by the clear() function.
   */
  SNESConvergedReason _reason;

  /**
   * Stores the total number of linear iterations from the last solve.
   */
  int _n_linear_iterations;

  /**
   * Stores the current nonlinear iteration number
   */
  unsigned _current_nonlinear_iteration_number;
};



} // namespace libMesh


#endif // #ifdef LIBMESH_HAVE_PETSC
#endif // #ifdef __petsc_nonlinear_solver_h__
