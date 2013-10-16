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



#ifndef LIBMESH_PETSC_NONLINEAR_SOLVER_H
#define LIBMESH_PETSC_NONLINEAR_SOLVER_H

#include "libmesh/libmesh_config.h"

// Petsc include files.
#ifdef LIBMESH_HAVE_PETSC

// Local includes
#include "libmesh/nonlinear_solver.h"
#include "libmesh/petsc_macro.h"

// PETSc includes
EXTERN_C_FOR_PETSC_BEGIN
# include <petscsnes.h>
EXTERN_C_FOR_PETSC_END

// C++ includes

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
  explicit
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
   * Set if the residual should be zeroed out in the callback
   */
  void set_residual_zero_out(bool state) { _zero_out_residual = state; }

  /**
   * Set if the jacobian should be zeroed out in the callback
   */
  void set_jacobian_zero_out(bool state) { _zero_out_jacobian = state; }

  /**
   * Set to true to use the libMeash's default monitor, set to false to use your own
   */
  void use_default_monitor(bool state) { _default_monitor = state; }

protected:
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
  PetscInt _n_linear_iterations;

  /**
   * Stores the current nonlinear iteration number
   */
  unsigned _current_nonlinear_iteration_number;

  /**
   * true to zero out residual before going into application level call-back, otherwise false
   */
  bool _zero_out_residual;

  /**
   * true to zero out jacobian before going into application level call-back, otherwise false
   */
  bool _zero_out_jacobian;

  /**
   * true if we want the default monitor to be set, false for no monitor (i.e. user code can use their own)
   */
  bool _default_monitor;

 private:
#if !PETSC_VERSION_LESS_THAN(3,3,0)
  void build_mat_null_space(NonlinearImplicitSystem::ComputeVectorSubspace* computeSubspaceObject,
                            void (*)(std::vector<NumericVector<Number>*>&, sys_type&),
                            MatNullSpace*);
#endif

  friend PetscErrorCode __libmesh_petsc_snes_residual (SNES snes, Vec x, Vec r, void *ctx);
  friend PetscErrorCode __libmesh_petsc_snes_jacobian (SNES snes, Vec x, Mat *jac, Mat *pc, MatStructure *msflag, void *ctx);
};



} // namespace libMesh


#endif // #ifdef LIBMESH_HAVE_PETSC
#endif // LIBMESH_PETSC_NONLINEAR_SOLVER_H
