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
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/wrapped_petsc.h"

// PETSc includes
#ifdef I
# define LIBMESH_SAW_I
#endif
#include <petscmat.h>
#ifndef LIBMESH_SAW_I
# undef I // Avoid complex.h contamination
#endif

#include <functional>
#include <memory>
#include <utility>

namespace libMesh
{
class ConstraintEliminationDofMap;
class StaticCondensation;
class StaticCondensationDofMap;
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
   * \returns Whether \f$ P(y) \f$ exists. It may not exist when forming it requires a nonlinear
   * solve, and then \p current_local_solution is unspecified
   */
  virtual bool prolong(Vec y) = 0;

  /**
   * \returns A full vector that holds the iterate \f$ P(y) \f$, or null if it does not exist. The
   * result may share storage with \p y. The system's solution vectors may change
   */
  virtual std::unique_ptr<NumericVector<T>> prolong_iterate(Vec y) = 0;

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
   * \returns Whether \f$ P(y) \f$ exists. If not, \p r is unspecified and \p assemble is not called
   */
  virtual bool residual(Vec y, Vec r, const AssembleResidual & assemble) = 0;

  /**
   * Makes the reduced Jacobian \f$ P'^T J P' \f$ at \p y: prolongs \p y, has \p assemble form
   * \f$ J \f$ into \p assembly_matrix(pc), closes it, and calls \p finish_jacobian() on it. Closes
   * \p jac if it differs from \p pc, which requires that \p jac be matrix-free
   * \returns Whether \f$ P(y) \f$ exists
   */
  bool jacobian(Vec y, SparseMatrix<T> & jac, SparseMatrix<T> & pc, const AssembleJacobian & assemble);

  /**
   * \returns The operator and preconditioning matrices that SNES uses, given the system's
   * matrices \p jac and \p pc. By default these are \p jac and \p pc
   */
  virtual std::pair<Mat, Mat> operators(Mat jac, Mat pc) const { return {jac, pc}; }

  /**
   * \returns The matrix that the user Jacobian is assembled into, given the preconditioning
   * matrix \p pc that SNES passes. By default this is \p pc
   */
  virtual SparseMatrix<T> & assembly_matrix(SparseMatrix<T> & pc) { return pc; }

  /**
   * Whether SNES's linear solves use the solver's attached \p Preconditioner
   */
  virtual bool uses_linear_preconditioner() const { return true; }

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
  virtual bool prolong(Vec y) override;

  /**
   * \returns A vector that wraps \p y
   */
  virtual std::unique_ptr<NumericVector<T>> prolong_iterate(Vec y) override;

  /**
   * \returns A vector that wraps \p y
   */
  virtual std::unique_ptr<NumericVector<T>> prolong_step(Vec y) const override;

  virtual void inject(const NumericVector<T> & x, Vec y) const override;

  /**
   * Assembles \f$ F \f$ into \p r and enforces its constraint rows if requested
   */
  virtual bool residual(Vec y, Vec r, const typename PetscSolveReduction<T>::AssembleResidual & assemble) override;

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
  virtual bool prolong(Vec y) override;

  virtual std::unique_ptr<NumericVector<T>> prolong_iterate(Vec y) override;
  virtual std::unique_ptr<NumericVector<T>> prolong_step(Vec y) const override;
  virtual void inject(const NumericVector<T> & x, Vec y) const override;

  /**
   * Assembles \f$ C^T F \f$ into the system's \p rhs and injects it into \p r
   */
  virtual bool residual(Vec y, Vec r, const typename PetscSolveReduction<T>::AssembleResidual & assemble) override;

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

#ifdef LIBMESH_HAVE_EIGEN
/**
 * Nonlinear elimination of the statically condensed degrees of freedom. \f$ y \f$ holds the
 * degrees of freedom that the system's \p StaticCondensationDofMap retains, and \f$ z = (z_1, y)
 * \f$ is an iterate of a parent reduction, which maps \f$ z \f$ to the full space. The
 * prolongation is \f$ P(y) = P_p(G(y), y) \f$, where \f$ G(y) \f$ solves the condensed equations
 * \f$ F_1(G(y), y) = 0 \f$ and \f$ P_p \f$ is the parent prolongation. A
 * \p ComputeEliminatedSolution object computes \f$ G \f$.
 *
 * The step prolongation is \f$ P' = P_p' [-A_{11}^{-1} A_{12}; I] \f$. Since \f$ F_1 = 0 \f$ at \f$
 * P(y) \f$, the reduced residual is the uncondensed part \f$ F_2 \f$ of the parent residual, and
 * the reduced Jacobian is the Schur complement \f$ S = A_{22} - A_{21} A_{11}^{-1} A_{12} \f$ of the
 * parent Jacobian, which the system's \p StaticCondensation matrix forms
 */
template <typename T>
class PetscNonlinearStaticCondensationReduction : public PetscSolveReduction<T>
{
public:
  /**
   * \p parent maps the parent space of the system's \p StaticCondensationDofMap to the full space
   */
  PetscNonlinearStaticCondensationReduction(
      NonlinearImplicitSystem & sys,
      NonlinearImplicitSystem::ComputeEliminatedSolution & eliminator,
      std::unique_ptr<PetscSolveReduction<T>> parent);

  virtual numeric_index_type n_local_dofs() const override;
  virtual numeric_index_type n_dofs() const override;
  virtual dof_id_type reduced_dof(dof_id_type full_dof) const override;
  virtual Vec create_vector() const override;
  virtual void begin_solve(PetscVector<T> & x, PetscVector<T> & r) override;

  /**
   * Sets the full solution to \f$ P(y) \f$ at the final \f$ y \f$, which need not be the last
   * iterate evaluated
   */
  virtual void end_solve() override;

  /**
   * Prolongs \f$ (z_1, y) \f$ through the parent, with the last converged condensed values \f$ z_1
   * \f$ as the initial guess, and then solves the condensed equations. On success the condensed
   * values become the initial guess of the next prolongation. After the first Jacobian, the initial
   * guess adds the tangent predictor \f$ -A_{11}^{-1} A_{12} (y - y_{old}) \f$, unless the
   * eliminator disables it
   */
  virtual bool prolong(Vec y) override;

  /**
   * Prolongs \p y and copies the system's solution
   */
  virtual std::unique_ptr<NumericVector<T>> prolong_iterate(Vec y) override;

  /**
   * Forms \f$ P_p' [-A_{11}^{-1} A_{12} y; y] \f$ with the condensed blocks of the last assembled
   * Jacobian
   */
  virtual std::unique_ptr<NumericVector<T>> prolong_step(Vec y) const override;

  virtual void inject(const NumericVector<T> & x, Vec y) const override;

  /**
   * Prolongs \p y and injects the uncondensed part of the parent residual into \p r
   */
  virtual bool residual(Vec y, Vec r, const typename PetscSolveReduction<T>::AssembleResidual & assemble) override;

  /**
   * The Schur complement for both
   */
  virtual std::pair<Mat, Mat> operators(Mat jac, Mat pc) const override;

  /**
   * The system's \p StaticCondensation matrix, whose closing forms the Schur complement
   */
  virtual SparseMatrix<T> & assembly_matrix(SparseMatrix<T> & pc) override;

  /**
   * Has the parent complete \p matrix
   */
  virtual void finish_jacobian(SparseMatrix<T> & matrix) override;

  /**
   * The linear solves are on the Schur complement directly, so the static condensation
   * preconditioner does not apply
   */
  virtual bool uses_linear_preconditioner() const override { return false; }

private:
  NonlinearImplicitSystem::ComputeEliminatedSolution & _eliminator;

  /// The reduction from the parent space to the full space
  std::unique_ptr<PetscSolveReduction<T>> _parent;

  const StaticCondensationDofMap & _dof_map;
  StaticCondensation & _sc;

  /// Entry \p i is the full index of the \p i-th local entry of \f$ y \f$
  WrappedPetsc<IS> _retained_is;

  /// Entry \p i is the parent index of the \p i-th local entry of \f$ y \f$
  WrappedPetsc<IS> _retained_parent_is;

  /// The storage of \p _snes_solution and \p _snes_residual
  WrappedPetsc<Vec> _y, _r;

  /// The parent iterate \f$ z \f$ and the parent residual
  WrappedPetsc<Vec> _z, _parent_r;

  /// Whether \p _sc holds the factored condensed blocks of an assembled Jacobian
  bool _have_factors = false;
};
#endif // LIBMESH_HAVE_EIGEN

} // namespace libMesh

#endif // LIBMESH_HAVE_PETSC
#endif // LIBMESH_PETSC_SOLVE_REDUCTION_H
