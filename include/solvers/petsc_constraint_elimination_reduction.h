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

#ifndef LIBMESH_PETSC_CONSTRAINT_ELIMINATION_REDUCTION_H
#define LIBMESH_PETSC_CONSTRAINT_ELIMINATION_REDUCTION_H

#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_HAVE_PETSC

#include "libmesh/petsc_solve_reduction.h"

namespace libMesh
{
class ConstraintEliminationDofMap;

/**
 * The reduction onto the unconstrained degrees of freedom of a \p ConstraintEliminationDofMap.
 * Its prolongation is \f$ P(y) = C y + h \f$, with \f$ P' = C \f$. Element assembly already
 * applies \f$ C^T \f$, so restricting a residual or a Jacobian simply drops its constrained rows
 * and columns
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
  virtual std::unique_ptr<NumericVector<T>> prolong_step(Vec dy) const override;
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
#endif // LIBMESH_PETSC_CONSTRAINT_ELIMINATION_REDUCTION_H
