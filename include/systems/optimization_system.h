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



#ifndef LIBMESH_OPTIMIZATION_SYSTEM_H
#define LIBMESH_OPTIMIZATION_SYSTEM_H

// Local Includes
#include "libmesh/implicit_system.h"

// C++ includes

namespace libMesh
{


// Forward declarations
template<typename T> class OptimizationSolver;


/**
 * This System subclass enables us to assemble an objective function,
 * gradient, Hessian and bounds for optimization problems.
 *
 * \author David Knezevic
 * \date 2015
 */
class OptimizationSystem : public ImplicitSystem
{
public:

  /**
   * Constructor.  Optionally initializes required
   * data structures.
   */
  OptimizationSystem (EquationSystems & es,
                      const std::string & name,
                      const unsigned int number);

  /**
   * Destructor.
   */
  virtual ~OptimizationSystem ();

  /**
   * The type of system.
   */
  typedef OptimizationSystem sys_type;

  /**
   * The type of the parent.
   */
  typedef ImplicitSystem Parent;

  /**
   * Abstract base class to be used to calculate the objective
   * function for optimization.
   */
  class ComputeObjective
  {
  public:
    virtual ~ComputeObjective () {}

    /**
     * This function will be called to compute the objective function
     * to be minimized, and must be implemented by the user in a
     * derived class. @return the value of the objective function at
     * the iterate \p X.
     */
    virtual Number objective (const NumericVector<Number> & X,
                              sys_type & S) = 0;
  };


  /**
   * Abstract base class to be used to calculate the gradient of
   * an objective function.
   */
  class ComputeGradient
  {
  public:
    virtual ~ComputeGradient () {}

    /**
     * This function will be called to compute the gradient of the
     * objective function, and must be implemented by the user in
     * a derived class. Set \p grad_f to be the gradient at the
     * iterate \p X.
     */
    virtual void gradient (const NumericVector<Number> & X,
                           NumericVector<Number> & grad_f,
                           sys_type & S) = 0;
  };


  /**
   * Abstract base class to be used to calculate the Hessian
   * of an objective function.
   */
  class ComputeHessian
  {
  public:
    virtual ~ComputeHessian () {}

    /**
     * This function will be called to compute the Hessian of
     * the objective function, and must be implemented by the
     * user in a derived class. Set \p H_f to be the gradient
     * at the iterate \p X.
     */
    virtual void hessian (const NumericVector<Number> & X,
                          SparseMatrix<Number> & H_f,
                          sys_type & S) = 0;
  };

  /**
   * Abstract base class to be used to calculate the equality constraints.
   */
  class ComputeEqualityConstraints
  {
  public:
    virtual ~ComputeEqualityConstraints () {}

    /**
     * This function will be called to evaluate the equality constraints
     * vector C_eq(X). This will impose the constraints C_eq(X) = 0.
     */
    virtual void equality_constraints (const NumericVector<Number> & X,
                                       NumericVector<Number> & C_eq,
                                       sys_type & S) = 0;
  };

  /**
   * Abstract base class to be used to calculate the Jacobian of the
   * equality constraints.
   */
  class ComputeEqualityConstraintsJacobian
  {
  public:
    virtual ~ComputeEqualityConstraintsJacobian () {}

    /**
     * This function will be called to evaluate the Jacobian of C_eq(X).
     */
    virtual void equality_constraints_jacobian (const NumericVector<Number> & X,
                                                SparseMatrix<Number> & C_eq_jac,
                                                sys_type & S) = 0;
  };

  /**
   * Abstract base class to be used to calculate the inequality constraints.
   */
  class ComputeInequalityConstraints
  {
  public:
    virtual ~ComputeInequalityConstraints () {}

    /**
     * This function will be called to evaluate the equality constraints
     * vector C_ineq(X). This will impose the constraints C_ineq(X) >= 0.
     */
    virtual void inequality_constraints (const NumericVector<Number> & X,
                                         NumericVector<Number> & C_ineq,
                                         sys_type & S) = 0;
  };

  /**
   * Abstract base class to be used to calculate the Jacobian of the
   * inequality constraints.
   */
  class ComputeInequalityConstraintsJacobian
  {
  public:
    virtual ~ComputeInequalityConstraintsJacobian () {}

    /**
     * This function will be called to evaluate the Jacobian of C_ineq(X).
     */
    virtual void inequality_constraints_jacobian (const NumericVector<Number> & X,
                                                  SparseMatrix<Number> & C_ineq_jac,
                                                  sys_type & S) = 0;
  };

  /**
   * Abstract base class to be used to calculate the lower and upper
   * bounds for all dofs in the system.
   */
  class ComputeLowerAndUpperBounds
  {
  public:
    virtual ~ComputeLowerAndUpperBounds () {}

    /**
     * This function should update the following two vectors:
     *   this->get_vector("lower_bounds"),
     *   this->get_vector("upper_bounds").
     */
    virtual void lower_and_upper_bounds (sys_type & S) = 0;
  };

  /**
   * @returns a clever pointer to the system.
   */
  sys_type & system () { return *this; }

  /**
   * Clear all the data structures associated with
   * the system.
   */
  virtual void clear () libmesh_override;

  /**
   * Initializes new data members of the system.
   */
  virtual void init_data () libmesh_override;

  /**
   * Reinitializes the member data fields associated with
   * the system, so that, e.g., \p assemble() may be used.
   */
  virtual void reinit () libmesh_override;

  /**
   * Solves the optimization problem.
   */
  virtual void solve () libmesh_override;

  /**
   * Initialize storage for the equality constraints, and the
   * corresponding Jacobian. The length of \p constraint_jac_sparsity
   * indicates the number of constraints that will be imposed,
   * and n_dofs_per_constraint[i] gives the indices that are non-zero
   * in row i of the Jacobian.
   */
  void initialize_equality_constraints_storage(const std::vector< std::set<numeric_index_type> > & constraint_jac_sparsity);

  /**
   * Initialize storage for the inequality constraints, as per
   * initialize_equality_constraints_storage.
   */
  void initialize_inequality_constraints_storage(const std::vector< std::set<numeric_index_type> > & constraint_jac_sparsity);

  /**
   * @returns \p "Optimization".  Helps in identifying
   * the system type in an equation system file.
   */
  virtual std::string system_type () const libmesh_override { return "Optimization"; }

  /**
   * The \p OptimizationSolver that is used for performing the optimization.
   */
  UniquePtr<OptimizationSolver<Number> > optimization_solver;

  /**
   * The vector that stores equality constraints.
   */
  UniquePtr<NumericVector<Number> > C_eq;

  /**
   * The sparse matrix that stores the Jacobian of C_eq.
   */
  UniquePtr<SparseMatrix<Number> > C_eq_jac;

  /**
   * The vector that stores inequality constraints.
   */
  UniquePtr<NumericVector<Number> > C_ineq;

  /**
   * The sparse matrix that stores the Jacobian of C_ineq.
   */
  UniquePtr<SparseMatrix<Number> > C_ineq_jac;

  /**
   * Vectors to store the dual variables associated with equality
   * and inequality constraints.
   */
  UniquePtr<NumericVector<Number> > lambda_eq;
  UniquePtr<NumericVector<Number> > lambda_ineq;

  /**
   * A copy of the equality and inequality constraint Jacobian sparsity
   * patterns.
   */
  std::vector< std::set<numeric_index_type> > eq_constraint_jac_sparsity;
  std::vector< std::set<numeric_index_type> > ineq_constraint_jac_sparsity;
};

} // namespace libMesh

#endif // LIBMESH_OPTIMIZATION_SYSTEM_H
