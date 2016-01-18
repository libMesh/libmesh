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



#ifndef LIBMESH_OPTIMIZATION_SOLVER_H
#define LIBMESH_OPTIMIZATION_SOLVER_H

// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/enum_solver_package.h"
#include "libmesh/reference_counted_object.h"
#include "libmesh/libmesh.h"
#include "libmesh/parallel_object.h"
#include "libmesh/auto_ptr.h"
#include "libmesh/optimization_system.h"

// C++ includes
#include <cstddef>

namespace libMesh
{

// forward declarations
template <typename T> class SparseMatrix;
template <typename T> class NumericVector;
template <typename T> class Preconditioner;

/**
 * This class provides a uniform interface for optimization solvers.  This base
 * class is overloaded to provide optimization solvers from different packages.
 *
 * \author David Knezevic
 * \date 2015
 */
template <typename T>
class OptimizationSolver : public ReferenceCountedObject<OptimizationSolver<T> >,
                           public ParallelObject
{
public:
  /**
   * The type of system
   */
  typedef OptimizationSystem sys_type;

  /**
   *  Constructor. Initializes Solver data structures
   */
  explicit
  OptimizationSolver (sys_type & s);

  /**
   * Destructor.
   */
  virtual ~OptimizationSolver ();

  /**
   * Builds an \p OptimizationSolver using the package specified by
   * \p solver_package
   */
  static UniquePtr<OptimizationSolver<T> > build(sys_type & s,
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
   */
  virtual void init () = 0;

  /**
   * Solves the optimization problem.
   */
  virtual void solve () = 0;

  /**
   * Get the current values of dual variables associated with
   * inequality and equality constraints. The variables will
   * be stored in _system.lambda_eq and _system.lambda_ineq.
   */
  virtual void get_dual_variables()
  { libmesh_not_implemented(); }

  /**
   * Prints a useful message about why the latest optimization solve
   * con(di)verged.
   */
  virtual void print_converged_reason() { libmesh_not_implemented(); }

  /**
   * Most optimization solver packages return an integral status
   * result of some kind.  This interface assumes they can be coerced
   * into an "int" type, which is usually safe since they are based on
   * enumerations.  Simply returns 0 if not implemented in derived
   * classes.
   */
  virtual int get_converged_reason() { return 0; }

  /**
   * Object that computes the objective function f(X)
   * at the input iterate X.
   */
  OptimizationSystem::ComputeObjective * objective_object;

  /**
   * Object that computes the gradient grad_f(X) of the objective function
   * at the input iterate X.
   */
  OptimizationSystem::ComputeGradient * gradient_object;

  /**
   * Object that computes the Hessian H_f(X) of the objective function
   * at the input iterate X.
   */
  OptimizationSystem::ComputeHessian * hessian_object;

  /**
   * Object that computes the equality constraints vector C_eq(X).
   * This will lead to the constraints C_eq(X) = 0 being imposed.
   */
  OptimizationSystem::ComputeEqualityConstraints * equality_constraints_object;

  /**
   * Object that computes the Jacobian of C_eq(X).
   */
  OptimizationSystem::ComputeEqualityConstraintsJacobian * equality_constraints_jacobian_object;

  /**
   * Object that computes the inequality constraints vector C_ineq(X).
   * This will lead to the constraints C_ineq(X) >= 0 being imposed.
   */
  OptimizationSystem::ComputeInequalityConstraints * inequality_constraints_object;

  /**
   * Object that computes the Jacobian of C_ineq(X).
   */
  OptimizationSystem::ComputeInequalityConstraintsJacobian * inequality_constraints_jacobian_object;

  /**
   * Object that computes the lower and upper bounds vectors.
   */
  OptimizationSystem::ComputeLowerAndUpperBounds * lower_and_upper_bounds_object;

  /**
   * @returns a constant reference to the system we are using to
   * define the optimization problem.
   */
  const sys_type & system () const { return _system; }

  /**
   * @returns a writeable reference to the system we are using to
   * define the optimization problem.
   */
  sys_type & system () { return _system; }

  /**
   * Maximum number of objective function evaluations allowed.
   */
  unsigned int max_objective_function_evaluations;

  /**
   * Required change in objective function which signals convergence.
   */
  Real objective_function_relative_tolerance;

  /**
   * Control how much is output from the OptimizationSolver as it's running.
   */
  bool verbose;

protected:

  /**
   * A reference to the system we are solving.
   */
  sys_type & _system;

  /**
   * Flag indicating if the data structures have been initialized.
   */
  bool _is_initialized;

};

} // namespace libMesh


#endif // LIBMESH_OPTIMIZATION_SOLVER_H
