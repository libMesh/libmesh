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

#include "libmesh/libmesh_common.h"

#ifdef LIBMESH_HAVE_PETSC

#include "libmesh/petsc_solve_reduction.h"
#include "libmesh/constraint_elimination_dof_map.h"
#include "libmesh/dof_map.h"
#include "libmesh/dof_object.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/petsc_matrix_base.h"
#include "libmesh/petsc_vector.h"

namespace libMesh
{

//------------------------------------------------------------------
// PetscSolveReduction

template <typename T>
PetscSolveReduction<T>::PetscSolveReduction(NonlinearImplicitSystem & sys)
  : _sys(sys), _solution(nullptr), _snes_solution(nullptr), _snes_residual(nullptr)
{
}

template <typename T>
PetscSolveReduction<T>::~PetscSolveReduction() = default;

template <typename T>
void PetscSolveReduction<T>::jacobian(Vec y,
                                      SparseMatrix<T> & jac,
                                      SparseMatrix<T> & pc,
                                      const AssembleJacobian & assemble)
{
  this->prolong(y);
  assemble(pc);
  pc.close();
  this->finish_jacobian(pc);

  if (&jac != &pc)
    {
      this->check_jacobian_operator(jac);
      jac.close();
    }
}

//------------------------------------------------------------------
// PetscIdentityReduction

template <typename T>
PetscIdentityReduction<T>::PetscIdentityReduction(NonlinearImplicitSystem & sys,
                                                  const bool & exact_constraint_enforcement)
  : PetscSolveReduction<T>(sys), _exact_constraint_enforcement(exact_constraint_enforcement)
{
}

template <typename T>
numeric_index_type PetscIdentityReduction<T>::n_local_dofs() const
{
  return this->_sys.get_dof_map().n_local_dofs();
}

template <typename T>
numeric_index_type PetscIdentityReduction<T>::n_dofs() const
{
  return this->_sys.get_dof_map().n_dofs();
}

template <typename T>
Vec PetscIdentityReduction<T>::create_vector() const
{
  Vec v;
  LibmeshPetscCall2(this->_sys.comm(),
                    VecDuplicate(cast_ref<PetscVector<T> &>(*this->_sys.solution).vec(), &v));
  return v;
}

template <typename T>
void PetscIdentityReduction<T>::begin_solve(PetscVector<T> & x, PetscVector<T> & r)
{
  // SNES iterates in x itself, so the system solution follows the solve
  this->_solution = &x;
  this->_snes_solution = x.vec();
  this->_snes_residual = r.vec();
}

template <typename T>
void PetscIdentityReduction<T>::prolong(Vec y)
{
  NonlinearImplicitSystem & sys = this->_sys;
  PetscVector<T> & X_sys = *cast_ptr<PetscVector<T> *>(sys.solution.get());
  PetscVector<T> X_global(y, sys.comm());

  // Use the system's update() to get a good local version of the
  // parallel solution.  This operation does not modify the incoming
  // "y" vector, it only localizes information from "y" into
  // sys.current_local_solution.
  X_global.swap(X_sys);
  sys.update();
  X_global.swap(X_sys);

  // Enforce constraints (if any) exactly on the
  // current_local_solution.  This is the solution vector that is
  // actually used in residual and Jacobian evaluations, and is
  // not locked by debug-enabled PETSc the way that "y" is.
  if (_exact_constraint_enforcement)
    sys.get_dof_map().enforce_constraints_exactly(sys, sys.current_local_solution.get());
}

template <typename T>
std::unique_ptr<NumericVector<T>> PetscIdentityReduction<T>::prolong_iterate(Vec y) const
{
  return this->prolong_step(y);
}

template <typename T>
std::unique_ptr<NumericVector<T>> PetscIdentityReduction<T>::prolong_step(Vec y) const
{
  return std::make_unique<PetscVector<T>>(y, this->_sys.comm());
}

template <typename T>
void PetscIdentityReduction<T>::inject(const NumericVector<T> & x, Vec y) const
{
  const auto & petsc_x = cast_ref<const PetscVector<T> &>(x);
  if (petsc_x.vec() != y)
    LibmeshPetscCall2(this->_sys.comm(), VecCopy(petsc_x.vec(), y));
}

template <typename T>
void PetscIdentityReduction<T>::residual(
    Vec y, Vec r, const typename PetscSolveReduction<T>::AssembleResidual & assemble)
{
  NonlinearImplicitSystem & sys = this->_sys;
  this->prolong(y);

  PetscVector<T> R(r, sys.comm());
  assemble(R);
  R.close();

  // Synchronize PETSc y to local solution since the local solution may be changed due to the constraints
  PetscVector<T> & X_sys = *cast_ptr<PetscVector<T> *>(sys.solution.get());
  PetscVector<T> X_global(y, sys.comm());

  X_global.swap(X_sys);
  sys.update();
  X_global.swap(X_sys);

  if (_exact_constraint_enforcement)
    {
      sys.get_dof_map().enforce_constraints_on_residual(
          sys, &R, sys.current_local_solution.get(), /*homogeneous=*/false);
      R.close();
    }
}

template <typename T>
void PetscIdentityReduction<T>::finish_jacobian(SparseMatrix<T> & matrix)
{
  if (_exact_constraint_enforcement)
    {
      this->_sys.get_dof_map().enforce_constraints_on_jacobian(this->_sys, &matrix);
      matrix.close();
    }
}

template <typename T>
void PetscIdentityReduction<T>::check_jacobian_operator(SparseMatrix<T> & jac) const
{
  // Constraint rows are enforced on the preconditioning matrix only. Assume that matrix-free
  // operators know what they're doing
#ifndef NDEBUG
  if (_exact_constraint_enforcement)
    {
      PetscBool j_is_mffd = PETSC_FALSE, j_is_shell = PETSC_FALSE;
      const Mat j = cast_ref<PetscMatrixBase<T> &>(jac).mat();
      LibmeshPetscCall2(this->_sys.comm(), PetscObjectTypeCompare((PetscObject)j, MATMFFD, &j_is_mffd));
      LibmeshPetscCall2(this->_sys.comm(), PetscObjectTypeCompare((PetscObject)j, MATSHELL, &j_is_shell));
      libmesh_assert(j_is_mffd == PETSC_TRUE || j_is_shell == PETSC_TRUE);
    }
#else
  libmesh_ignore(jac);
#endif
}

//------------------------------------------------------------------
// PetscConstraintEliminationReduction

template <typename T>
PetscConstraintEliminationReduction<T>::PetscConstraintEliminationReduction(
    NonlinearImplicitSystem & sys)
  : PetscSolveReduction<T>(sys), _dof_map(sys.get_dof_map().get_constraint_elimination())
{
  const auto & retained = _dof_map.local_retained_dofs();
  std::vector<PetscInt> indices(retained.begin(), retained.end());
  LibmeshPetscCall2(sys.comm(),
                    ISCreateGeneral(sys.comm().get(),
                                    cast_int<PetscInt>(indices.size()),
                                    indices.data(),
                                    PETSC_COPY_VALUES,
                                    _retained_is.get()));
}

template <typename T>
numeric_index_type PetscConstraintEliminationReduction<T>::n_local_dofs() const
{
  return _dof_map.n_local_dofs();
}

template <typename T>
numeric_index_type PetscConstraintEliminationReduction<T>::n_dofs() const
{
  return _dof_map.n_dofs();
}

template <typename T>
dof_id_type PetscConstraintEliminationReduction<T>::reduced_dof(const dof_id_type full_dof) const
{
  return _dof_map.find_reduced_dof(full_dof);
}

template <typename T>
Vec PetscConstraintEliminationReduction<T>::create_vector() const
{
  const auto & comm = this->_sys.comm();
  VecType vec_type;
  LibmeshPetscCall2(
      comm, VecGetType(cast_ref<PetscVector<T> &>(*this->_sys.solution).vec(), &vec_type));

  Vec v;
  LibmeshPetscCall2(comm, VecCreate(comm.get(), &v));
  LibmeshPetscCall2(comm,
                    VecSetSizes(v,
                                cast_int<PetscInt>(this->n_local_dofs()),
                                cast_int<PetscInt>(this->n_dofs())));
  LibmeshPetscCall2(comm, VecSetType(v, vec_type));
  LibmeshPetscCall2(comm, VecSetUp(v));
  return v;
}

template <typename T>
void PetscConstraintEliminationReduction<T>::begin_solve(PetscVector<T> & x, PetscVector<T> &)
{
  this->_solution = &x;
  _y.reset_to_zero();
  _r.reset_to_zero();
  *_y.get() = this->create_vector();
  *_r.get() = this->create_vector();
  this->_snes_solution = _y;
  this->_snes_residual = _r;

  // The initial guess is the unconstrained part of the current solution
  x.close();
  this->inject(x, _y);
}

template <typename T>
void PetscConstraintEliminationReduction<T>::end_solve()
{
  libmesh_assert(this->_solution);
  this->prolong_into(_y, *this->_solution, /*step=*/false);
}

template <typename T>
void PetscConstraintEliminationReduction<T>::prolong_into(Vec y,
                                                          NumericVector<T> & x,
                                                          const bool step) const
{
  auto & petsc_x = cast_ref<PetscVector<T> &>(x);
  libmesh_assert_equal_to(petsc_x.type(), PARALLEL);

  // The constrained entries of x are overwritten below, so their prior values do not matter
  LibmeshPetscCall2(this->_sys.comm(), VecISCopy(petsc_x.vec(), _retained_is, SCATTER_FORWARD, y));
  x.close();
  this->_sys.get_dof_map().enforce_constraints_exactly(this->_sys, &x, /*homogeneous=*/step);
}

template <typename T>
void PetscConstraintEliminationReduction<T>::prolong(Vec y)
{
  this->prolong_into(y, *this->_sys.solution, /*step=*/false);
  this->_sys.update();
}

template <typename T>
std::unique_ptr<NumericVector<T>> PetscConstraintEliminationReduction<T>::prolong_iterate(Vec y) const
{
  auto x = this->_sys.solution->zero_clone();
  this->prolong_into(y, *x, /*step=*/false);
  return x;
}

template <typename T>
std::unique_ptr<NumericVector<T>> PetscConstraintEliminationReduction<T>::prolong_step(Vec y) const
{
  auto x = this->_sys.solution->zero_clone();
  this->prolong_into(y, *x, /*step=*/true);
  return x;
}

template <typename T>
void PetscConstraintEliminationReduction<T>::inject(const NumericVector<T> & x, Vec y) const
{
  const auto & petsc_x = cast_ref<const PetscVector<T> &>(x);
  LibmeshPetscCall2(this->_sys.comm(),
                    VecISCopy(petsc_x.vec(), _retained_is, SCATTER_REVERSE, y));
}

template <typename T>
void PetscConstraintEliminationReduction<T>::residual(
    Vec y, Vec r, const typename PetscSolveReduction<T>::AssembleResidual & assemble)
{
  this->prolong(y);
  NumericVector<T> & F = *this->_sys.rhs;
  assemble(F);
  F.close();
  this->inject(F, r);
}

//------------------------------------------------------------------
// Explicit instantiations
template class LIBMESH_EXPORT PetscSolveReduction<Number>;
template class LIBMESH_EXPORT PetscIdentityReduction<Number>;
template class LIBMESH_EXPORT PetscConstraintEliminationReduction<Number>;

} // namespace libMesh

#endif // LIBMESH_HAVE_PETSC
