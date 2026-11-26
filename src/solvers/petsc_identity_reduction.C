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

#include "libmesh/petsc_identity_reduction.h"
#include "libmesh/dof_map.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/petsc_matrix_base.h"
#include "libmesh/petsc_vector.h"

namespace libMesh
{

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
  auto x = std::make_unique<PetscVector<T>>(y, this->_sys.comm());
  if (!_exact_constraint_enforcement)
    return x;

  // Match the iterate that prolong() gives to assembly
  auto enforced = x->clone();
  this->_sys.get_dof_map().enforce_constraints_exactly(this->_sys, enforced.get());
  return enforced;
}

template <typename T>
std::unique_ptr<NumericVector<T>> PetscIdentityReduction<T>::prolong_step(Vec dy) const
{
  return std::make_unique<PetscVector<T>>(dy, this->_sys.comm());
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
// Explicit instantiations
template class LIBMESH_EXPORT PetscIdentityReduction<Number>;

} // namespace libMesh

#endif // LIBMESH_HAVE_PETSC
