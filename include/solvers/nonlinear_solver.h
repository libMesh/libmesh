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



#ifndef LIBMESH_NONLINEAR_SOLVER_H
#define LIBMESH_NONLINEAR_SOLVER_H

// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/enum_solver_package.h"
#include "libmesh/reference_counted_object.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/libmesh.h"
#include "libmesh/parallel_object.h"
#include "libmesh/auto_ptr.h"

// C++ includes
#include <cstddef>

namespace libMesh
{

// forward declarations
template <typename T> class SparseMatrix;
template <typename T> class NumericVector;
template <typename T> class Preconditioner;
class SolverConfiguration;

/**
 * This class provides a uniform interface for nonlinear solvers.  This base
 * class is overloaded to provide nonlinear solvers from different packages
 * like PETSC.
 *
 * \author Benjamin Kirk
 * \date 2005
 */
template <typename T>
class NonlinearSolver : public ReferenceCountedObject<NonlinearSolver<T> >,
                        public ParallelObject
{
public:
  /**
   * The type of system
   */
  typedef NonlinearImplicitSystem sys_type;

  /**
   *  Constructor. Initializes Solver data structures
   */
  explicit
  NonlinearSolver (sys_type & s);

  /**
   * Destructor.
   */
  virtual ~NonlinearSolver ();

  /**
   * Builds a \p NonlinearSolver using the nonlinear solver package specified by
   * \p solver_package
   */
  static UniquePtr<NonlinearSolver<T> > build(sys_type & s,
                                              const SolverPackage solver_package = libMesh::default_solver_package());

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
   * May assign a name to the solver in some implementations
   */
  virtual void init (const char * name = libmesh_nullptr) = 0;

  /**
   * Solves the nonlinear system.
   */
  virtual std::pair<unsigned int, Real> solve (SparseMatrix<T> &,  // System Jacobian Matrix
                                               NumericVector<T> &, // Solution vector
                                               NumericVector<T> &, // Residual vector
                                               const double,      // Stopping tolerance
                                               const unsigned int) = 0; // N. Iterations

  /**
   * Prints a useful message about why the latest nonlinear solve
   * con(di)verged.
   */
  virtual void print_converged_reason() { libmesh_not_implemented(); }

  /**
   * Get the total number of linear iterations done in the last solve
   */
  virtual int get_total_linear_iterations() = 0;

  /**
   * If called *during* the solve(), for example by the user-specified
   * residual or Jacobian function, returns the current nonlinear iteration
   * number.  Must be redefined in derived classes.
   */
  virtual unsigned get_current_nonlinear_iteration_number() const = 0;

  /**
   * Function that computes the residual \p R(X) of the nonlinear system
   * at the input iterate \p X.
   */
  void (* residual) (const NumericVector<Number> & X,
                     NumericVector<Number> & R,
                     sys_type & S);

  /**
   * Object that computes the residual \p R(X) of the nonlinear system
   * at the input iterate \p X.
   */
  NonlinearImplicitSystem::ComputeResidual * residual_object;

  /**
   * Function that computes the Jacobian \p J(X) of the nonlinear system
   * at the input iterate \p X.
   */
  void (* jacobian) (const NumericVector<Number> & X,
                     SparseMatrix<Number> & J,
                     sys_type & S);

  /**
   * Object that computes the Jacobian \p J(X) of the nonlinear system
   * at the input iterate \p X.
   */
  NonlinearImplicitSystem::ComputeJacobian * jacobian_object;

  /**
   * Function that computes either the residual \f$ R(X) \f$ or the
   * Jacobian \f$ J(X) \f$ of the nonlinear system at the input
   * iterate \f$ X \f$.  Note that either \p R or \p J could be
   * \p NULL.
   */
  void (* matvec) (const NumericVector<Number> & X,
                   NumericVector<Number> * R,
                   SparseMatrix<Number> * J,
                   sys_type & S);

  /**
   * Object that computes either the residual \f$ R(X) \f$ or the
   * Jacobian \f$ J(X) \f$ of the nonlinear system at the input
   * iterate \f$ X \f$.  Note that either \p R or \p J could be
   * \p NULL.
   */
  NonlinearImplicitSystem::ComputeResidualandJacobian * residual_and_jacobian_object;

  /**
   * Function that computes the lower and upper bounds \p XL and \p XU on the solution of the nonlinear system.
   */
  void (* bounds) (NumericVector<Number> & XL,
                   NumericVector<Number> & XU,
                   sys_type & S);
  /**
   * Object that computes the bounds vectors  \f$ XL \f$ and \f$ XU \f$.
   */
  NonlinearImplicitSystem::ComputeBounds * bounds_object;

  /**
   * Function that computes a basis for the Jacobian's nullspace --
   * the kernel or the "zero energy modes" -- that can be used in
   * solving a degenerate problem iteratively, if the solver supports it
   * (e.g., PETSc's KSP).
   */
  void (* nullspace) (std::vector<NumericVector<Number> *> & sp, sys_type & S);

  /**
   * A callable object that computes a basis for the Jacobian's nullspace --
   * the kernel or the "zero energy modes" -- that can be used in
   * solving a degenerate problem iteratively, if the solver supports it
   * (e.g., PETSc's KSP).
   */
  NonlinearImplicitSystem::ComputeVectorSubspace * nullspace_object;

  /**
   * Function that computes a basis for the transpose Jacobian's nullspace --
   * when solving a degenerate problem iteratively, if the solver supports it
   * (e.g., PETSc's KSP), it is used to remove contributions outside of R(jac)
   */
  void (* transpose_nullspace) (std::vector<NumericVector<Number> *> & sp, sys_type & S);

  /**
   * A callable object that computes a basis for the transpose Jacobian's nullspace --
   * when solving a degenerate problem iteratively, if the solver supports it
   * (e.g., PETSc's KSP), it is used to remove contributions outside of R(jac)
   */
  NonlinearImplicitSystem::ComputeVectorSubspace * transpose_nullspace_object;

  /**
   * Function that computes a basis for the Jacobian's near nullspace --
   * the set of "low energy modes" -- that can be used for AMG coarsening,
   * if the solver supports it (e.g., ML, PETSc's GAMG).
   */
  void (* nearnullspace) (std::vector<NumericVector<Number> *> & sp, sys_type & S);

  /**
   * A callable object that computes a basis for the Jacobian's near nullspace --
   * the set of "low energy modes" -- that can be used for AMG coarsening,
   * if the solver supports it (e.g., ML, PETSc's GAMG).
   */
  NonlinearImplicitSystem::ComputeVectorSubspace * nearnullspace_object;

  /**
   * Customizable function pointer which users can attach to the
   * solver.  Gets called prior to every call to solve().
   */
  void (* user_presolve)(sys_type & S);

  /**
   * Function that performs a "check" on the Newton search direction
   * and solution after each nonlinear step. See documentation for the
   * NonlinearImplicitSystem::ComputePostCheck object for more
   * information about the calling sequence.
   */
  void (* postcheck) (const NumericVector<Number> & old_soln,
                      NumericVector<Number> & search_direction,
                      NumericVector<Number> & new_soln,
                      bool & changed_search_direction,
                      bool & changed_new_soln,
                      sys_type & S);

  /**
   * A callable object that is executed after each nonlinear
   * iteration. Allows the user to modify both the search direction
   * and the solution vector in an application-specific way.
   */
  NonlinearImplicitSystem::ComputePostCheck * postcheck_object;

  /**
   * @returns a constant reference to the system we are solving.
   */
  const sys_type & system () const { return _system; }

  /**
   * @returns a writeable reference to the system we are solving.
   */
  sys_type & system () { return _system; }

  /**
   * Attaches a Preconditioner object to be used during the linear solves.
   */
  void attach_preconditioner(Preconditioner<T> * preconditioner);

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

  /**
   * After a call to solve this will reflect whether or not the nonlinear
   * solve was successful.
   */
  bool converged;

  /**
   * Set the solver configuration object.
   */
  void set_solver_configuration(SolverConfiguration & solver_configuration);

protected:
  /**
   * A reference to the system we are solving.
   */
  sys_type & _system;

  /**
   * Flag indicating if the data structures have been initialized.
   */
  bool _is_initialized;

  /**
   * Holds the Preconditioner object to be used for the linear solves.
   */
  Preconditioner<T> * _preconditioner;

  /**
   * Optionally store a SolverOptions object that can be used
   * to set parameters like solver type, tolerances and iteration limits.
   */
  SolverConfiguration * _solver_configuration;
};




/*----------------------- inline functions ----------------------------------*/
template <typename T>
inline
NonlinearSolver<T>::NonlinearSolver (sys_type & s) :
  ParallelObject               (s),
  residual                     (libmesh_nullptr),
  residual_object              (libmesh_nullptr),
  jacobian                     (libmesh_nullptr),
  jacobian_object              (libmesh_nullptr),
  matvec                       (libmesh_nullptr),
  residual_and_jacobian_object (libmesh_nullptr),
  bounds                       (libmesh_nullptr),
  bounds_object                (libmesh_nullptr),
  nullspace                    (libmesh_nullptr),
  nullspace_object             (libmesh_nullptr),
  transpose_nullspace          (libmesh_nullptr),
  transpose_nullspace_object   (libmesh_nullptr),
  nearnullspace                (libmesh_nullptr),
  nearnullspace_object         (libmesh_nullptr),
  user_presolve                (libmesh_nullptr),
  postcheck                    (libmesh_nullptr),
  postcheck_object             (libmesh_nullptr),
  max_nonlinear_iterations(0),
  max_function_evaluations(0),
  absolute_residual_tolerance(0),
  relative_residual_tolerance(0),
  absolute_step_tolerance(0),
  relative_step_tolerance(0),
  max_linear_iterations(0),
  initial_linear_tolerance(0),
  minimum_linear_tolerance(0),
  converged(false),
  _system(s),
  _is_initialized (false),
  _preconditioner (libmesh_nullptr),
  _solver_configuration(libmesh_nullptr)
{
}



template <typename T>
inline
NonlinearSolver<T>::~NonlinearSolver ()
{
  this->clear ();
}


} // namespace libMesh


#endif // LIBMESH_NONLINEAR_SOLVER_H
