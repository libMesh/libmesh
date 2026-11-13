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
class ConstraintEliminationDofMap;
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
   * Sets the system's \p current_local_solution to the full iterate \f$ P(y) \f$
   */
  virtual void prolong(Vec y) = 0;

  /**
   * \returns A full vector that holds the iterate \f$ P(y) \f$. The result may share storage with
   * \p y
   */
  virtual std::unique_ptr<NumericVector<T>> prolong_iterate(Vec y) const = 0;

  /**
   * \returns The full step \f$ P' y \f$ for the step \p y. The result may share storage with \p y
   */
  virtual std::unique_ptr<NumericVector<T>> prolong_step(Vec y) const = 0;

  /**
   * Sets \p y to the injection \f$ R x \f$ of the full vector \p x
   */
  virtual void inject(const NumericVector<T> & x, Vec y) const = 0;

  /**
   * Sets \p r to the reduced residual \f$ P'^T F(P(y)) \f$. \p assemble forms \f$ F \f$
   */
  virtual void residual(Vec y, Vec r, const AssembleResidual & assemble) = 0;

  /**
   * Makes \p pc the reduced Jacobian \f$ P'^T J P' \f$ at \p y: prolongs \p y, has \p assemble
   * form \f$ J \f$ into \p pc, closes it, and calls \p finish_jacobian() on it. Closes \p jac if
   * it differs from \p pc, which requires that \p jac be matrix-free
   */
  void jacobian(Vec y, SparseMatrix<T> & jac, SparseMatrix<T> & pc, const AssembleJacobian & assemble);

  /**
   * Completes an assembled and closed Jacobian \p matrix as the reduced Jacobian. Does nothing by
   * default. The residual callback calls it on a Jacobian that it assembles together with the
   * residual
   */
  virtual void finish_jacobian(SparseMatrix<T> & /*matrix*/) {}

protected:
  /**
   * Checks \p jac, an operator that differs from the assembled preconditioning matrix. Does
   * nothing by default
   */
  virtual void check_jacobian_operator(SparseMatrix<T> & /*jac*/) const {}

  NonlinearImplicitSystem & _sys;

  /// The full solution passed to \p begin_solve()
  PetscVector<T> * _solution;

  /// See \p snes_solution() and \p snes_residual()
  Vec _snes_solution;
  Vec _snes_residual;
};

/**
 * The reduction \f$ P(y) = y \f$, which SNES uses when the system eliminates nothing. SNES then
 * iterates in the full solution vector itself. Constraint rows remain in the residual and the
 * Jacobian, where they are enforced if \p exact_constraint_enforcement is set
 */
template <typename T>
class PetscIdentityReduction : public PetscSolveReduction<T>
{
public:
  /**
   * \p exact_constraint_enforcement is read at every evaluation, so it may change after
   * construction
   */
  PetscIdentityReduction(NonlinearImplicitSystem & sys, const bool & exact_constraint_enforcement);

  virtual numeric_index_type n_local_dofs() const override;
  virtual numeric_index_type n_dofs() const override;
  virtual dof_id_type reduced_dof(dof_id_type full_dof) const override { return full_dof; }
  virtual Vec create_vector() const override;
  virtual void begin_solve(PetscVector<T> & x, PetscVector<T> & r) override;
  virtual void end_solve() override {}

  /**
   * Localizes \p y into the system's \p current_local_solution and enforces constraints on it if
   * requested. The system's \p solution keeps its value
   */
  virtual void prolong(Vec y) override;

  /**
   * \returns A vector that wraps \p y
   */
  virtual std::unique_ptr<NumericVector<T>> prolong_iterate(Vec y) const override;

  /**
   * \returns A vector that wraps \p y
   */
  virtual std::unique_ptr<NumericVector<T>> prolong_step(Vec y) const override;

  virtual void inject(const NumericVector<T> & x, Vec y) const override;

  /**
   * Assembles \f$ F \f$ into \p r and enforces its constraint rows if requested
   */
  virtual void residual(Vec y, Vec r, const typename PetscSolveReduction<T>::AssembleResidual & assemble) override;

  /**
   * Enforces the constraint rows of \p matrix if requested
   */
  virtual void finish_jacobian(SparseMatrix<T> & matrix) override;

protected:
  virtual void check_jacobian_operator(SparseMatrix<T> & jac) const override;

private:
  const bool & _exact_constraint_enforcement;
};

/**
 * The reduction \f$ P(y) = C y + h \f$ onto the unconstrained degrees of freedom of a \p
 * ConstraintEliminationDofMap, with \f$ P' = C \f$. Element assembly already applies \f$ C^T \f$,
 * so restricting a residual or a Jacobian drops its constrained rows and columns
 */
template <typename T>
class PetscConstraintEliminationReduction : public PetscSolveReduction<T>
{
public:
  PetscConstraintEliminationReduction(NonlinearImplicitSystem & sys);

  virtual numeric_index_type n_local_dofs() const override;
  virtual numeric_index_type n_dofs() const override;
  virtual dof_id_type reduced_dof(dof_id_type full_dof) const override;
  virtual Vec create_vector() const override;
  virtual void begin_solve(PetscVector<T> & x, PetscVector<T> & r) override;
  virtual void end_solve() override;

  /**
   * Sets the system's \p solution to \f$ C y + h \f$ and localizes it
   */
  virtual void prolong(Vec y) override;

  virtual std::unique_ptr<NumericVector<T>> prolong_iterate(Vec y) const override;
  virtual std::unique_ptr<NumericVector<T>> prolong_step(Vec y) const override;
  virtual void inject(const NumericVector<T> & x, Vec y) const override;

  /**
   * Assembles \f$ C^T F \f$ into the system's \p rhs and injects it into \p r
   */
  virtual void residual(Vec y, Vec r, const typename PetscSolveReduction<T>::AssembleResidual & assemble) override;

private:
  /**
   * Sets the parallel full vector \p x to \f$ C y \f$ if \p step, else to \f$ C y + h \f$
   */
  void prolong_into(Vec y, NumericVector<T> & x, bool step) const;

  const ConstraintEliminationDofMap & _dof_map;

  /// Entry \p i is the full index of the \p i-th local entry of \f$ y \f$
  WrappedPetsc<IS> _retained_is;

  /// The storage of \p _snes_solution and \p _snes_residual
  WrappedPetsc<Vec> _y, _r;
};

} // namespace libMesh

#endif // LIBMESH_HAVE_PETSC
#endif // LIBMESH_PETSC_SOLVE_REDUCTION_H
