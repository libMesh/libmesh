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
#include "libmesh/libmesh_logging.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/petsc_matrix_base.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/static_condensation.h"
#include "libmesh/static_condensation_dof_map.h"

namespace
{
using namespace libMesh;

/**
 * \returns A new parallel vector of the type of \p sys's solution, with \p n_local entries on
 * this processor and \p n in total
 */
template <typename T>
Vec create_reduced_vector(const NonlinearImplicitSystem & sys,
                          const numeric_index_type n_local,
                          const numeric_index_type n)
{
  const auto & comm = sys.comm();
  VecType vec_type;
  LibmeshPetscCall2(comm, VecGetType(cast_ref<PetscVector<T> &>(*sys.solution).vec(), &vec_type));

  Vec v;
  LibmeshPetscCall2(comm, VecCreate(comm.get(), &v));
  LibmeshPetscCall2(comm, VecSetSizes(v, cast_int<PetscInt>(n_local), cast_int<PetscInt>(n)));
  LibmeshPetscCall2(comm, VecSetType(v, vec_type));
  LibmeshPetscCall2(comm, VecSetUp(v));
  return v;
}

/**
 * Sets \p is to the indices \p indices on the communicator of \p sys
 */
template <typename Indices>
void create_index_set(const NonlinearImplicitSystem & sys,
                      const Indices & indices,
                      WrappedPetsc<IS> & is)
{
  std::vector<PetscInt> petsc_indices(indices.begin(), indices.end());
  LibmeshPetscCall2(sys.comm(),
                    ISCreateGeneral(sys.comm().get(),
                                    cast_int<PetscInt>(petsc_indices.size()),
                                    petsc_indices.data(),
                                    PETSC_COPY_VALUES,
                                    is.get()));
}
}

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
bool PetscSolveReduction<T>::jacobian(Vec y,
                                      SparseMatrix<T> & jac,
                                      SparseMatrix<T> & pc,
                                      const AssembleJacobian & assemble)
{
  if (!this->prolong(y))
    return false;

  SparseMatrix<T> & matrix = this->assembly_matrix(pc);
  assemble(matrix);
  matrix.close();
  this->finish_jacobian(matrix);

  if (&jac != &pc)
    {
      this->check_jacobian_operator(jac);
      jac.close();
    }
  return true;
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
bool PetscIdentityReduction<T>::prolong(Vec y)
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

  return true;
}

template <typename T>
std::unique_ptr<NumericVector<T>> PetscIdentityReduction<T>::prolong_iterate(Vec y)
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
bool PetscIdentityReduction<T>::residual(
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
  return true;
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
  create_index_set(sys, _dof_map.local_retained_dofs(), _retained_is);
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
  return create_reduced_vector<T>(this->_sys, this->n_local_dofs(), this->n_dofs());
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
bool PetscConstraintEliminationReduction<T>::prolong(Vec y)
{
  this->prolong_into(y, *this->_sys.solution, /*step=*/false);
  this->_sys.update();
  return true;
}

template <typename T>
std::unique_ptr<NumericVector<T>> PetscConstraintEliminationReduction<T>::prolong_iterate(Vec y)
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
bool PetscConstraintEliminationReduction<T>::residual(
    Vec y, Vec r, const typename PetscSolveReduction<T>::AssembleResidual & assemble)
{
  this->prolong(y);
  NumericVector<T> & F = *this->_sys.rhs;
  assemble(F);
  F.close();
  this->inject(F, r);
  return true;
}

#ifdef LIBMESH_HAVE_EIGEN
//------------------------------------------------------------------
// PetscNonlinearStaticCondensationReduction

template <typename T>
PetscNonlinearStaticCondensationReduction<T>::PetscNonlinearStaticCondensationReduction(
    NonlinearImplicitSystem & sys,
    NonlinearImplicitSystem::ComputeEliminatedSolution & eliminator,
    std::unique_ptr<PetscSolveReduction<T>> parent)
  : PetscSolveReduction<T>(sys),
    _eliminator(eliminator),
    _parent(std::move(parent)),
    _dof_map(sys.get_dof_map().get_static_condensation()),
    _sc(sys.get_static_condensation())
{
  libmesh_assert(_parent);
  libmesh_assert_equal_to(_dof_map.parent_n_local_dofs(), _parent->n_local_dofs());
  create_index_set(sys, _dof_map.local_retained_dofs(), _retained_is);
  create_index_set(sys, _dof_map.local_retained_parent_dofs(), _retained_parent_is);
}

template <typename T>
numeric_index_type PetscNonlinearStaticCondensationReduction<T>::n_local_dofs() const
{
  return _dof_map.n_local_dofs();
}

template <typename T>
numeric_index_type PetscNonlinearStaticCondensationReduction<T>::n_dofs() const
{
  return _dof_map.n_dofs();
}

template <typename T>
dof_id_type PetscNonlinearStaticCondensationReduction<T>::reduced_dof(const dof_id_type full_dof) const
{
  return _dof_map.find_reduced_dof(full_dof);
}

template <typename T>
Vec PetscNonlinearStaticCondensationReduction<T>::create_vector() const
{
  return create_reduced_vector<T>(this->_sys, this->n_local_dofs(), this->n_dofs());
}

template <typename T>
void PetscNonlinearStaticCondensationReduction<T>::begin_solve(PetscVector<T> & x, PetscVector<T> &)
{
  this->_solution = &x;
  for (auto * v : {&_y, &_r, &_z, &_parent_r})
    v->reset_to_zero();
  *_y.get() = this->create_vector();
  *_r.get() = this->create_vector();
  *_z.get() = _parent->create_vector();
  *_parent_r.get() = _parent->create_vector();
  this->_snes_solution = _y;
  this->_snes_residual = _r;

  // The condensed entries of the current solution are the initial guess of the first local solves
  x.close();
  _parent->inject(x, _z);
  this->inject(x, _y);
}

template <typename T>
void PetscNonlinearStaticCondensationReduction<T>::end_solve()
{
  libmesh_assert(this->_solution);
  // The last evaluation may be a rejected trial iterate, so solve the condensed equations again at
  // the final one
  this->prolong(_y);
  *this->_solution = *this->_sys.solution;
}

template <typename T>
bool PetscNonlinearStaticCondensationReduction<T>::prolong(Vec y)
{
  NonlinearImplicitSystem & sys = this->_sys;

  // _z keeps the condensed values of the last successful prolongation. Once a Jacobian has been
  // assembled, the tangent predictor dz1 = -A11^{-1} A12 (y - y_old) moves them toward the new
  // solution of the condensed equations
  if (_have_factors && _eliminator.use_tangent_predictor())
    {
      // dz = [0; y - y_old], with y_old the uncondensed entries of _z
      WrappedPetsc<Vec> dz;
      LibmeshPetscCall2(sys.comm(), VecDuplicate(_z, dz.get()));
      LibmeshPetscCall2(sys.comm(), VecCopy(_z, dz));
      LibmeshPetscCall2(sys.comm(), VecISCopy(dz, _retained_parent_is, SCATTER_FORWARD, y));
      LibmeshPetscCall2(sys.comm(), VecAXPY(dz, -1, _z));
      {
        PetscVector<T> parent_step(dz, sys.comm());
        _sc.homogeneous_backwards_substitution(parent_step);
      }
      LibmeshPetscCall2(sys.comm(), VecAXPY(_z, 1, dz));
    }
  // Copy y exactly. y_old + (y - y_old) can differ from y by roundoff
  LibmeshPetscCall2(sys.comm(), VecISCopy(_z, _retained_parent_is, SCATTER_FORWARD, y));
  _parent->prolong(_z);

  // The eliminator reads current_local_solution and writes the condensed values into solution
  *sys.solution = *sys.current_local_solution;
  bool converged;
  {
    LOG_SCOPE("eliminate()", "PetscNonlinearStaticCondensationReduction");
    converged = _eliminator.eliminate(sys);
  }
  sys.solution->close();
  sys.comm().min(converged);
  if (!converged)
    return false;

  sys.update();

  // Keep the condensed values for the next prolongation. The uncondensed entries of _z stay y: the
  // parent may change them in the system solution, e.g. by enforcing constraints, but its residual
  // must see the iterate SNES holds
  _parent->inject(*sys.solution, _z);
  LibmeshPetscCall2(sys.comm(), VecISCopy(_z, _retained_parent_is, SCATTER_FORWARD, y));
  return true;
}

template <typename T>
std::unique_ptr<NumericVector<T>>
PetscNonlinearStaticCondensationReduction<T>::prolong_iterate(Vec y)
{
  if (!this->prolong(y))
    return nullptr;
  return this->_sys.solution->clone();
}

template <typename T>
std::unique_ptr<NumericVector<T>> PetscNonlinearStaticCondensationReduction<T>::prolong_step(Vec y) const
{
  // dz = [-A11^{-1} A12 y; y] in the parent space
  WrappedPetsc<Vec> dz;
  *dz.get() = _parent->create_vector();
  LibmeshPetscCall2(this->_sys.comm(), VecZeroEntries(dz));
  LibmeshPetscCall2(this->_sys.comm(), VecISCopy(dz, _retained_parent_is, SCATTER_FORWARD, y));
  PetscVector<T> parent_step(dz, this->_sys.comm());
  _sc.homogeneous_backwards_substitution(parent_step);

  // The parent's result may share storage with dz, which goes out of scope
  return _parent->prolong_step(dz)->clone();
}

template <typename T>
void PetscNonlinearStaticCondensationReduction<T>::inject(const NumericVector<T> & x, Vec y) const
{
  const auto & petsc_x = cast_ref<const PetscVector<T> &>(x);
  LibmeshPetscCall2(this->_sys.comm(), VecISCopy(petsc_x.vec(), _retained_is, SCATTER_REVERSE, y));
}

template <typename T>
bool PetscNonlinearStaticCondensationReduction<T>::residual(
    Vec y, Vec r, const typename PetscSolveReduction<T>::AssembleResidual & assemble)
{
  if (!this->prolong(y))
    return false;

  // _z now holds P(y) in the parent space. Its condensed rows of the parent residual vanish
  if (!_parent->residual(_z, _parent_r, assemble))
    return false;
  LibmeshPetscCall2(this->_sys.comm(),
                    VecISCopy(_parent_r, _retained_parent_is, SCATTER_REVERSE, r));
  return true;
}

template <typename T>
std::pair<Mat, Mat> PetscNonlinearStaticCondensationReduction<T>::operators(Mat, Mat) const
{
  const Mat schur = cast_ref<const PetscMatrixBase<T> &>(_sc.get_condensed_mat()).mat();
  return {schur, schur};
}

template <typename T>
SparseMatrix<T> & PetscNonlinearStaticCondensationReduction<T>::assembly_matrix(SparseMatrix<T> &)
{
  return _sc;
}

template <typename T>
void PetscNonlinearStaticCondensationReduction<T>::finish_jacobian(SparseMatrix<T> & matrix)
{
  _parent->finish_jacobian(matrix);
  // Closing the static condensation matrix factored the condensed blocks
  _have_factors = true;
}
#endif // LIBMESH_HAVE_EIGEN

//------------------------------------------------------------------
// Explicit instantiations
template class LIBMESH_EXPORT PetscSolveReduction<Number>;
template class LIBMESH_EXPORT PetscIdentityReduction<Number>;
template class LIBMESH_EXPORT PetscConstraintEliminationReduction<Number>;
#ifdef LIBMESH_HAVE_EIGEN
template class LIBMESH_EXPORT PetscNonlinearStaticCondensationReduction<Number>;
#endif

} // namespace libMesh

#endif // LIBMESH_HAVE_PETSC
