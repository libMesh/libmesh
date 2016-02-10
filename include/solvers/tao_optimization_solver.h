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



#ifndef LIBMESH_TAO_OPTIMIZATION_SOLVER_H
#define LIBMESH_TAO_OPTIMIZATION_SOLVER_H

#include "libmesh/libmesh_config.h"

// Petsc include files.
#if defined(LIBMESH_HAVE_PETSC_TAO) && !defined(LIBMESH_USE_COMPLEX_NUMBERS)

// Local includes
#include "libmesh/petsc_macro.h"
#include "libmesh/optimization_solver.h"

// Include header for the Tao optimization library
#include <petsctao.h>

namespace libMesh
{

// Allow users access to these functions in case they want to reuse them.  Note that users shouldn't
// need access to these most of the time as they are used internally by this object.
extern "C"
{
  PetscErrorCode __libmesh_tao_objective (Tao tao, Vec x, PetscReal * objective, void * ctx);
  PetscErrorCode __libmesh_tao_gradient(Tao tao, Vec x, Vec g, void * ctx);
  PetscErrorCode __libmesh_tao_hessian(Tao tao, Vec x, Mat h, Mat pc, void * ctx);
  PetscErrorCode __libmesh_tao_equality_constraints(Tao tao, Vec x, Vec ce, void * ctx);
  PetscErrorCode __libmesh_tao_equality_constraints_jacobian(Tao tao, Vec x, Mat J, Mat Jpre, void * ctx);
  PetscErrorCode __libmesh_tao_inequality_constraints(Tao tao, Vec x, Vec cineq, void * ctx);
  PetscErrorCode __libmesh_tao_inequality_constraints_jacobian(Tao tao, Vec x, Mat J, Mat Jpre, void * ctx);
}

/**
 * This class provides an interface to the Tao optimization solvers.
 *
 * \author David Knezevic
 * \date 2015
 */
template <typename T>
class TaoOptimizationSolver : public OptimizationSolver<T>
{
public:

  /**
   * The type of system that we use in conjunction with this solver.
   */
  typedef OptimizationSystem sys_type;

  /**
   *  Constructor. Initializes Tao data structures.
   */
  explicit
  TaoOptimizationSolver (sys_type & system);

  /**
   * Destructor.
   */
  ~TaoOptimizationSolver ();

  /**
   * Release all memory and clear data structures.
   */
  virtual void clear () libmesh_override;

  /**
   * Initialize data structures if not done so already.
   */
  virtual void init () libmesh_override;

  /**
   * Returns the raw PETSc Tao context pointer.
   */
  Tao tao() { this->init(); return _tao; }

  /**
   * Call the Tao solver.
   */
  virtual void solve () libmesh_override;

  /**
   * Get the current values of dual variables associated with
   * inequality and equality constraints. The variables will
   * be stored in _system.lambda_eq and _system.lambda_ineq.
   */
  virtual void get_dual_variables() libmesh_override;

  /**
   * Prints a useful message about why the latest optimization solve
   * con(di)verged.
   */
  virtual void print_converged_reason() libmesh_override;

  /**
   * Returns the currently-available (or most recently obtained, if the Tao object has
   * been destroyed) convergence reason.  Refer to Tao docs for the meaning of different
   * TaoConvergedReason.
   */
  virtual int get_converged_reason() libmesh_override;

protected:

  /**
   * Optimization solver context
   */
  Tao _tao;

  /**
   * Store the reason for Tao convergence/divergence for use even after _tao
   * has been cleared.  Note that print_converged_reason() will always *try* to
   * get the current reason with TaoGetConvergedReason(), but if the Tao object
   * has already been cleared, it will fall back on this stored value.  Note that
   * this value is therefore necessarily *not* cleared by the clear() function.
   */
  TaoConvergedReason _reason;

private:

  friend PetscErrorCode __libmesh_tao_objective (Tao tao, Vec x, PetscReal * objective, void * ctx);
  friend PetscErrorCode __libmesh_tao_gradient(Tao tao, Vec x, Vec g, void * ctx);
  friend PetscErrorCode __libmesh_tao_hessian(Tao tao, Vec x, Mat h, Mat pc, void * ctx);
  friend PetscErrorCode __libmesh_tao_equality_constraints(Tao tao, Vec x, Vec ce, void * ctx);
  friend PetscErrorCode __libmesh_tao_equality_constraints_jacobian(Tao tao, Vec x, Mat J, Mat Jpre, void * ctx);
  friend PetscErrorCode __libmesh_tao_inequality_constraints(Tao tao, Vec x, Vec cineq, void * ctx);
  friend PetscErrorCode __libmesh_tao_inequality_constraints_jacobian(Tao tao, Vec x, Mat J, Mat Jpre, void * ctx);
};



} // namespace libMesh


#endif // #if defined(LIBMESH_HAVE_PETSC_TAO) && !defined(LIBMESH_USE_COMPLEX_NUMBERS)
#endif // LIBMESH_TAO_OPTIMIZATION_SOLVER_H
