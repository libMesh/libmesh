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



#ifndef LIBMESH_NONLINEAR_IMPLICIT_SYSTEM_H
#define LIBMESH_NONLINEAR_IMPLICIT_SYSTEM_H

// Local Includes
#include "libmesh/implicit_system.h"

// C++ includes

namespace libMesh
{


// Forward declarations
class DiffSolver;
template<typename T> class NonlinearSolver;


/**
 * This class provides a specific system class.  It aims
 * at implicit systems, offering nothing more than just
 * the essentials needed to solve a system.  Note
 * that still additional vectors/matrices may be added,
 * as offered in the parent class \p ExplicitSystem.
 */
class NonlinearImplicitSystem : public ImplicitSystem
{
public:

  /**
   * Constructor.  Optionally initializes required
   * data structures.
   */
  NonlinearImplicitSystem (EquationSystems & es,
                           const std::string & name,
                           const unsigned int number);

  /**
   * Destructor.
   */
  virtual ~NonlinearImplicitSystem ();

  /**
   * The type of system.
   */
  typedef NonlinearImplicitSystem sys_type;

  /**
   * The type of the parent.
   */
  typedef ImplicitSystem Parent;

  /**
   * Abstract base class to be used to calculate the residual
   * of a nonlinear system.
   */
  class ComputeResidual
  {
  public:
    virtual ~ComputeResidual () {}
    /**
     * Residual function.  This function will be called to compute the
     * residual and must be implemented by the user in a derived class.
     */
    virtual void residual (const NumericVector<Number> & X,
                           NumericVector<Number> & R,
                           sys_type & S) = 0;
  };


  /**
   * Abstract base class to be used to calculate the Jacobian
   * of a nonlinear system.
   */
  class ComputeJacobian
  {
  public:
    virtual ~ComputeJacobian () {}

    /**
     * Jacobian function.  This function will be called to compute the
     * jacobian and must be implemented by the user in a derived class.
     */
    virtual void jacobian (const NumericVector<Number> & X,
                           SparseMatrix<Number> & J,
                           sys_type & S) = 0;
  };


  /**
   * Abstract base class to be used to calculate the bounds
   * on the degrees of freedom of a nonlinear system.
   */
  class ComputeBounds
  {
  public:
    virtual ~ComputeBounds () {}

    /**
     * This function will be called to compute the bounds vector and
     * must be implemented by the user in a derived class.
     */
    virtual void bounds (NumericVector<Number> & XL,
                         NumericVector<Number> & XU,
                         sys_type & S) = 0;
  };

  /**
   * Callable abstract base class to be used as a callback to provide
   * the solver with a basis for the system's Jacobian's nullspace
   * (the kernel or the "zero energy modes") or near-nullspace
   * (the "low energy modes").
   *
   * A nullspace can be used to solve a degenerate problem iteratively
   * (e.g., with a Krylov subspace method). A near nullspace can be used
   * by an Algebraic Multigrid (AMG) preconditioner to construct smoothed-
   * aggregation-like coarse spaces.
   */
  class ComputeVectorSubspace
  {
  public:
    virtual ~ComputeVectorSubspace () {}

    /**
     * This function will be called to compute the subspace basis
     * (e.g., nullspace or nearnullspace).
     * It must be implemented by the user in a derived class.
     */
    virtual void operator()(std::vector<NumericVector<Number> *> & sp,
                            sys_type & s);
  };

  /**
   * Abstract base class to be used to calculate the residual and Jacobian
   * simultaneously of a nonlinear system.
   */
  class ComputeResidualandJacobian
  {
  public:
    virtual ~ComputeResidualandJacobian () {}

    /**
     * Residual & Jacobian function, calculated simultaneously.
     * This function will be called to compute the residual and jacobian
     * simultaneously and must be implemented by the user in a derived class.
     */
    virtual void residual_and_jacobian (const NumericVector<Number> & X,
                                        NumericVector<Number> * R,
                                        SparseMatrix<Number> * J,
                                        sys_type & S) = 0;
  };

  /**
   * Abstract base class to be used for applying user modifications to
   * the solution vector and/or Newton update step after each
   * nonlinear step. A user should inherit from this class and define
   * a custom version of the postcheck() function.
   */
  class ComputePostCheck
  {
  public:
    virtual ~ComputePostCheck () {}

    /**
     * This interface, which is inspired by PETSc's, passes the user:
     * .) A constant reference to the "old" solution vector (previous Newton iterate),
     * .) A pointer to the search direction (update) vector, and
     * .) The "proposed" new solution vector.
     * The user can then modify either the search direction or the new
     * solution vector according to their own application-specific
     * criteria.  If they do, they must set the associated boolean
     * references appropriately.
     */
    virtual void postcheck (const NumericVector<Number> & old_soln,
                            NumericVector<Number> & search_direction,
                            NumericVector<Number> & new_soln,
                            bool & changed_search_direction,
                            bool & changed_new_soln,
                            sys_type & S) = 0;
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
   * Reinitializes the member data fields associated with
   * the system, so that, e.g., \p assemble() may be used.
   */
  virtual void reinit () libmesh_override;

  /**
   * Assembles & solves the nonlinear system R(x) = 0.
   */
  virtual void solve () libmesh_override;

  /**
   * Returns an integer corresponding to the upper iteration count
   * limit and a Real corresponding to the convergence tolerance to
   * be used in linear adjoint and/or sensitivity solves
   */
  virtual std::pair<unsigned int, Real>
  get_linear_solve_parameters() const libmesh_override;

  /**
   * Assembles a residual in \p rhs and/or a jacobian in \p matrix,
   * as requested.
   */
  virtual void assembly(bool get_residual,
                        bool get_jacobian,
                        bool apply_heterogeneous_constraints = false) libmesh_override;

  /**
   * @returns \p "NonlinearImplicit".  Helps in identifying
   * the system type in an equation system file.
   */
  virtual std::string system_type () const libmesh_override { return "NonlinearImplicit"; }

  /**
   * The \p NonlinearSolver defines the default interface used to
   * solve the nonlinear_implicit system.  This class handles all the
   * details of interfacing with various nonlinear algebra packages
   * like PETSc or LASPACK.
   */
  UniquePtr<NonlinearSolver<Number> > nonlinear_solver;

  /**
   * The \p DiffSolver defines an optional interface used to
   * solve the nonlinear_implicit system.
   */
  UniquePtr<DiffSolver> diff_solver;

  /**
   * Returns  the number of iterations
   * taken for the most recent nonlinear solve.
   */
  unsigned int n_nonlinear_iterations() const { return _n_nonlinear_iterations; }

  /**
   * Returns the final residual for the nonlinear system solve.
   */
  Real final_nonlinear_residual() const { return _final_nonlinear_residual; }

  /**
   * If called *during* the solve(), for example by the user-specified
   * residual or Jacobian function, return the current nonlinear iteration
   * number.
   */
  unsigned get_current_nonlinear_iteration_number() const;


protected:

  /**
   * Copies system parameters into nonlinear solver parameters
   */
  void set_solver_parameters();

  /**
   * The number of nonlinear iterations required to solve the nonlinear
   * system R(x)=0.
   */
  unsigned int _n_nonlinear_iterations;

  /**
   * The final residual for the nonlinear system R(x)
   */
  Real _final_nonlinear_residual;
};

} // namespace libMesh

#endif // LIBMESH_NONLINEAR_IMPLICIT_SYSTEM_H
