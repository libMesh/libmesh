// The libMesh Finite Element Library.
// Copyright (C) 2002-2026 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#ifndef LIBMESH_PETSC_SOLVE_REDUCTION_H
#define LIBMESH_PETSC_SOLVE_REDUCTION_H

#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_HAVE_PETSC

#include "libmesh/id_types.h"
#include "libmesh/petsc_macro.h"
#include "libmesh/wrapped_petsc.h"

// PETSc includes
#ifdef I
# define LIBMESH_SAW_I
#endif
#include <petscvec.h>
#ifndef LIBMESH_SAW_I
# undef I // Avoid complex.h contamination
#endif

#include <functional>
#include <memory>

namespace libMesh
{
class NonlinearImplicitSystem;
template <typename T>
class NumericVector;
template <typename T>
class PetscVector;
template <typename T>
class SparseMatrix;

/**
 * The map between the vector \f$ y \f$ that SNES solves for and the full solution \f$ x \f$ of a
 * \p NonlinearImplicitSystem.
 *
 * The full iterate is the prolongation \f$ x = P(y) \f$. Its derivative \f$ P' = dP/dy \f$ is the
 * step prolongation, which maps a step in \f$ y \f$ to a step in \f$ x \f$. SNES solves the
 * reduced residual \f$ P'^T F(P(y)) = 0 \f$, whose Jacobian is \f$ P'^T J P' \f$. Injection \f$ R
 * \f$ maps a full vector to \f$ y \f$ by taking its retained entries, so that \f$ R(P(y)) = y \f$.
 */
template <typename T>
class PetscSolveReduction
{
public:
  /// Assembles the full residual into its argument at the system's current_local_solution
  using AssembleResidual = std::function<void(NumericVector<T> &)>;
  /// Assembles the full Jacobian into its argument at the system's current_local_solution
  using AssembleJacobian = std::function<void(SparseMatrix<T> &)>;

  PetscSolveReduction(NonlinearImplicitSystem & sys);
  virtual ~PetscSolveReduction();

  /**
   * The number of entries of \f$ y \f$ on this processor
   */
  virtual numeric_index_type n_local_dofs() const = 0;

  /**
   * The global number of entries of \f$ y \f$
   */
  virtual numeric_index_type n_dofs() const = 0;

  /**
   * \returns The index in \f$ y \f$ of the local full degree of freedom \p full_dof, or \p
   * DofObject::invalid_id if \f$ y \f$ does not hold it
   */
  virtual dof_id_type reduced_dof(dof_id_type full_dof) const = 0;

  /**
   * \returns A new vector with the layout of \f$ y \f$. The caller owns it
   */
  virtual Vec create_vector() const = 0;

  /**
   * Prepares a solve whose initial guess is the full solution \p x and whose full residual is \p
   * r. Sets \p snes_solution() to \f$ R x \f$ and sizes \p snes_residual() like it
   */
  virtual void begin_solve(PetscVector<T> & x, PetscVector<T> & r) = 0;

  /**
   * Sets the full solution passed to \p begin_solve() to \f$ P(y) \f$, where \f$ y \f$ is
   * \p snes_solution()
   */
  virtual void end_solve() = 0;

  /**
   * The vector \f$ y \f$ that SNES solves for
   */
  Vec snes_solution() const { return _snes_solution; }

  /**
   * The vector that SNES evaluates the reduced residual into
   */
  Vec snes_residual() const { return _snes_residual; }

  /**
   * Sets the system's \p current_local_solution to the full iterate \f$ P(y) \f$, the point at
   * which the residual and Jacobian are assembled. The result goes into the existing storage of
   * the system's vectors
   */
  virtual void prolong(Vec y) = 0;

  /**
   * \returns A full vector that holds the iterate \f$ P(y) \f$, for a user callback such as a
   * line search postcheck. The
   * system's vectors are unaffected. The result may wrap \p y itself, for instance for the
   * identity reduction without constraint enforcement. In that case, users should be thoughtful
   * about how they use the returned vector for a couple reasons: 1) writes to the libMesh wrapped
   * vector will write to \p y itself 2) the libMesh wrapped vector should not outlive y itself
   * or else the users risks reading/writing into invalid memory.
   */
  virtual std::unique_ptr<NumericVector<T>> prolong_iterate(Vec y) const = 0;

  /**
   * \returns The full step \f$ P' dy \f$ for the step \p dy, for a user callback such as a line
   * search precheck or postcheck. \f$ P' \f$ is the derivative of \f$ P \f$, so it omits any
   * affine offset of \f$ P \f$, such as constraint values. When \f$ P' \f$ is the identity, the
   * result wraps \p dy itself and the user should be cautious with the returned libMesh wrapped
   * vector for the same reasons given in \p prolong_iterate().
   */
  virtual std::unique_ptr<NumericVector<T>> prolong_step(Vec dy) const = 0;

  /**
   * Sets \p y to the injection \f$ R x \f$ of the full vector \p x
   */
  virtual void inject(const NumericVector<T> & x, Vec y) const = 0;

  /**
   * Sets \p r to the reduced residual \f$ P'^T F(P(y)) \f$. \p assemble forms \f$ F \f$
   */
  virtual void residual(Vec y, Vec r, const AssembleResidual & assemble) = 0;

  /**
   * Forms the preconditioning matrix of the reduced system \f$ P'^T J P' \f$
   * where \f$ J \f$ is the Jacobian in the full degree of freedom space
   * evaluated by \p assemble using the prolongation of \p y. If \p jac does not
   * wrap the same PETSc matrix as \p pc, then \p jac should be a shell or
   * matrix-free finite difference matrix. Both \p pc and \p jac are closed
   * after assembly.
   */
  void jacobian(Vec y, SparseMatrix<T> & jac, SparseMatrix<T> & pc, const AssembleJacobian & assemble);

  virtual void finish_jacobian(SparseMatrix<T> & /*matrix*/) {}

protected:
  /**
   * Used for validity checking \p jac when it's different from the preconditioning matrix
   */
  virtual void check_jacobian_operator(SparseMatrix<T> & /*jac*/) const {}

  /// The nonlinear system in the full degree of freedom space
  NonlinearImplicitSystem & _sys;

  /// The full solution passed to \p begin_solve()
  PetscVector<T> * _solution;

  /// See \p snes_solution() and \p snes_residual()
  Vec _snes_solution;
  Vec _snes_residual;
};

} // namespace libMesh

#endif // LIBMESH_HAVE_PETSC
#endif // LIBMESH_PETSC_SOLVE_REDUCTION_H
