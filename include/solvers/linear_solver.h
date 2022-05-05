// The libMesh Finite Element Library.
// Copyright (C) 2002-2022 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_LINEAR_SOLVER_H
#define LIBMESH_LINEAR_SOLVER_H

// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/enum_subset_solve_mode.h" // SUBSET_ZERO
#include "libmesh/reference_counted_object.h"
#include "libmesh/libmesh.h"
#include "libmesh/parallel_object.h"

#ifdef LIBMESH_FORWARD_DECLARE_ENUMS
namespace libMesh
{
enum SolverPackage : int;
enum PreconditionerType : int;
enum SolverType : int;
enum LinearConvergenceReason : int;
}
#else
#include "libmesh/enum_solver_package.h"
#include "libmesh/enum_preconditioner_type.h"
#include "libmesh/enum_solver_type.h"
#include "libmesh/enum_convergence_flags.h"
#endif

// C++ includes
#include <cstddef>
#include <vector>
#include <memory>

namespace libMesh
{

// forward declarations
template <typename T> class SparseMatrix;
template <typename T> class NumericVector;
template <typename T> class ShellMatrix;
template <typename T> class Preconditioner;
class System;
class SolverConfiguration;

/**
 * This base class can be inherited from to provide interfaces to
 * linear solvers from different packages like PETSc and LASPACK.
 *
 * \author Benjamin Kirk
 * \date 2003
 */
template <typename T>
class LinearSolver : public ReferenceCountedObject<LinearSolver<T>>,
                     public ParallelObject
{
public:

  /**
   *  Constructor. Initializes Solver data structures
   */
  LinearSolver (const libMesh::Parallel::Communicator & comm_in);

  /**
   * Destructor.
   */
  virtual ~LinearSolver ();

  /**
   * Builds a \p LinearSolver using the linear solver package specified by
   * \p solver_package
   */
  static std::unique_ptr<LinearSolver<T>> build(const libMesh::Parallel::Communicator & comm_in,
                                                const SolverPackage solver_package = libMesh::default_solver_package());

  /**
   * \returns \p true if the data structures are
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
  virtual void init (const char * name = nullptr) = 0;

  /**
   * Apply names to the system to be solved.  For most packages this
   * is a no-op; for PETSc this sets an option prefix from the system
   * name and sets field names from the system's variable names.
   *
   * Since field names are applied to DoF numberings, this method must
   * be called again after any System reinit.
   */
  virtual void init_names (const System &) {}

  /**
   * \returns The type of solver to use.
   */
  SolverType solver_type () const { return _solver_type; }

  /**
   * Sets the type of solver to use.
   */
  void set_solver_type (const SolverType st)
  { _solver_type = st; }

  /**
   * \returns The type of preconditioner to use.
   */
  PreconditionerType preconditioner_type () const;

  /**
   * Sets the type of preconditioner to use.
   */
  void set_preconditioner_type (const PreconditionerType pct);

  /**
   * Attaches a Preconditioner object to be used
   */
  void attach_preconditioner(Preconditioner<T> * preconditioner);

  /**
   * Set the same_preconditioner flag, which indicates if we reuse the
   * same preconditioner for subsequent solves.
   */
  virtual void reuse_preconditioner(bool );

  /**
   * \returns \p same_preconditioner, which indicates if we reuse the
   * same preconditioner for subsequent solves.
   */
  bool get_same_preconditioner();

  /**
   * After calling this method, all successive solves will be
   * restricted to the given set of dofs, which must contain local
   * dofs on each processor only and not contain any duplicates.  This
   * mode can be disabled by calling this method with \p dofs being a
   * \p nullptr.
   */
  virtual void restrict_solve_to (const std::vector<unsigned int> * const dofs,
                                  const SubsetSolveMode subset_solve_mode=SUBSET_ZERO);

  /**
   * This function calls the solver \p _solver_type preconditioned
   * with the \p _preconditioner_type preconditioner.
   *
   * \note This method will compute the preconditioner from the system
   * matrix.
   */
  virtual std::pair<unsigned int, Real> solve (SparseMatrix<T> &,  // System Matrix
                                               NumericVector<T> &, // Solution vector
                                               NumericVector<T> &, // RHS vector
                                               const double,      // Stopping tolerance
                                               const unsigned int) = 0; // N. Iterations

  /**
   * Function to solve the adjoint system.
   *
   * \note This method will compute the preconditioner from the system
   * matrix. This is not a pure virtual function and is defined
   * linear_solver.C
   */
  virtual std::pair<unsigned int, Real> adjoint_solve (SparseMatrix<T> &,  // System Matrix
                                                       NumericVector<T> &, // Solution vector
                                                       NumericVector<T> &, // RHS vector
                                                       const double,      // Stopping tolerance
                                                       const unsigned int); // N. Iterations

  /**
   * This function calls the solver
   * "_solver_type" preconditioned with the
   * "_preconditioner_type" preconditioner.
   */
  virtual std::pair<unsigned int, Real> solve (SparseMatrix<T> &,  // System Matrix
                                               SparseMatrix<T> &,  // Preconditioning Matrix
                                               NumericVector<T> &, // Solution vector
                                               NumericVector<T> &, // RHS vector
                                               const double,      // Stopping tolerance
                                               const unsigned int) = 0; // N. Iterations

  /**
   * This function calls the solver "_solver_type" preconditioned with
   * the "_preconditioner_type" preconditioner.  The preconditioning
   * matrix is used if it is provided, or the system matrix is used if
   * \p precond_matrix is null
   */
  std::pair<unsigned int, Real> solve (SparseMatrix<T> & matrix,
                                       SparseMatrix<T> * precond_matrix,
                                       NumericVector<T> &, // Solution vector
                                       NumericVector<T> &, // RHS vector
                                       const double,      // Stopping tolerance
                                       const unsigned int); // N. Iterations



  /**
   * This function solves a system whose matrix is a shell matrix.
   */
  virtual std::pair<unsigned int, Real> solve (const ShellMatrix<T> & shell_matrix,
                                               NumericVector<T> &, // Solution vector
                                               NumericVector<T> &, // RHS vector
                                               const double,      // Stopping tolerance
                                               const unsigned int) = 0; // N. Iterations


  /**
   * This function solves a system whose matrix is a shell matrix, but
   * a sparse matrix is used as preconditioning matrix, this allowing
   * other preconditioners than JACOBI.
   */
  virtual std::pair<unsigned int, Real> solve (const ShellMatrix<T> & shell_matrix,
                                               const SparseMatrix<T> & precond_matrix,
                                               NumericVector<T> &, // Solution vector
                                               NumericVector<T> &, // RHS vector
                                               const double,      // Stopping tolerance
                                               const unsigned int) = 0; // N. Iterations


  /**
   * This function solves a system whose matrix is a shell matrix, but
   * an optional sparse matrix may be used as preconditioning matrix.
   */
  std::pair<unsigned int, Real> solve (const ShellMatrix<T> & matrix,
                                       const SparseMatrix<T> * precond_matrix,
                                       NumericVector<T> &, // Solution vector
                                       NumericVector<T> &, // RHS vector
                                       const double,      // Stopping tolerance
                                       const unsigned int); // N. Iterations


  /**
   * Prints a useful message about why the latest linear solve
   * con(di)verged.
   */
  virtual void print_converged_reason() const;

  /**
   * \returns The solver's convergence flag
   */
  virtual LinearConvergenceReason get_converged_reason() const = 0;

  /**
   * Set the solver configuration object.
   */
  void set_solver_configuration(SolverConfiguration & solver_configuration);

protected:


  /**
   * Enum stating which type of iterative solver to use.
   */
  SolverType _solver_type;

  /**
   * Enum stating with type of preconditioner to use.
   */
  PreconditionerType _preconditioner_type;

  /**
   * Flag indicating if the data structures have been initialized.
   */
  bool _is_initialized;

  /**
   * Holds the Preconditioner object to be used for the linear solves.
   */
  Preconditioner<T> * _preconditioner;

  /**
   * Boolean flag to indicate whether we want to use an identical
   * preconditioner to the previous solve. This can save
   * substantial work in the cases where the system matrix is
   * the same for successive solves.
   */
  bool same_preconditioner;

  /**
   * Optionally store a SolverOptions object that can be used
   * to set parameters like solver type, tolerances and iteration limits.
   */
  SolverConfiguration * _solver_configuration;
};




/*----------------------- inline functions ----------------------------------*/
template <typename T>
inline
LinearSolver<T>::~LinearSolver ()
{
  this->LinearSolver::clear ();
}

template <typename T>
inline
bool LinearSolver<T>::get_same_preconditioner()
{
  return same_preconditioner;
}

template <typename T>
inline
std::pair<unsigned int, Real>
LinearSolver<T>::solve (SparseMatrix<T> & mat,
                        SparseMatrix<T> * pc_mat,
                        NumericVector<T> & sol,
                        NumericVector<T> & rhs,
                        const double tol,
                        const unsigned int n_iter)
{
  if (pc_mat)
    return this->solve(mat, *pc_mat, sol, rhs, tol, n_iter);
  else
    return this->solve(mat, sol, rhs, tol, n_iter);
}


template <typename T>
inline
std::pair<unsigned int, Real>
LinearSolver<T>::solve (const ShellMatrix<T> & mat,
                        const SparseMatrix<T> * pc_mat,
                        NumericVector<T> & sol,
                        NumericVector<T> & rhs,
                        const double tol,
                        const unsigned int n_iter)
{
  if (pc_mat)
    return this->solve(mat, *pc_mat, sol, rhs, tol, n_iter);
  else
    return this->solve(mat, sol, rhs, tol, n_iter);
}

} // namespace libMesh


#endif // LIBMESH_LINEAR_SOLVER_H
