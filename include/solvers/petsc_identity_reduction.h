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

#ifndef LIBMESH_PETSC_IDENTITY_REDUCTION_H
#define LIBMESH_PETSC_IDENTITY_REDUCTION_H

#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_HAVE_PETSC

#include "libmesh/petsc_solve_reduction.h"

namespace libMesh
{
/**
 * The reduction whose prolongation is the identity, \f$ P(y) = y \f$. SNES uses it when the
 * system eliminates nothing. SNES then iterates in the full solution vector itself. Constraint
 * rows remain in the residual and the Jacobian, where they are enforced if \p
 * exact_constraint_enforcement is set
 */
template <typename T>
class PetscIdentityReduction : public PetscSolveReduction<T>
{
public:
  /**
   * \p exact_constraint_enforcement is read at every evaluation; it may change after
   * construction so we take it as a reference
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
   * requested. The system's \p solution is unaffected
   */
  virtual void prolong(Vec y) override;

  /**
   * \returns A copy of \p y with constraints enforced on it if requested, else a vector that wraps
   * \p y
   */
  virtual std::unique_ptr<NumericVector<T>> prolong_iterate(Vec y) const override;

  /**
   * \returns A vector that wraps \p dy
   */
  virtual std::unique_ptr<NumericVector<T>> prolong_step(Vec dy) const override;

  virtual void inject(const NumericVector<T> & x, Vec y) const override;

  /**
   * Assembles \f$ F \f$ into \p r and enforces its constraint rows if requested through
   * \p _exact_constraint_enforcement
   */
  virtual void residual(Vec y, Vec r, const typename PetscSolveReduction<T>::AssembleResidual & assemble) override;

  /**
   * Enforces the constraint rows of \p matrix if requested through \p _exact_constraint_enforcement
   */
  virtual void finish_jacobian(SparseMatrix<T> & matrix) override;

protected:
  virtual void check_jacobian_operator(SparseMatrix<T> & jac) const override;

private:
  /// Whether to perform exact constraint enforcement on the global residual and Jacobian
  /// structures after their assembly
  const bool & _exact_constraint_enforcement;
};

} // namespace libMesh

#endif // LIBMESH_HAVE_PETSC
#endif // LIBMESH_PETSC_IDENTITY_REDUCTION_H
