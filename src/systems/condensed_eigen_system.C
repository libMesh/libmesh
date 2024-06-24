// The libMesh Finite Element Library.
// Copyright (C) 2002-2024 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#include "libmesh/libmesh_config.h"

// Currently, the EigenSystem should only be available
// if SLEPc support is enabled.
#if defined(LIBMESH_HAVE_SLEPC)

#include "libmesh/condensed_eigen_system.h"

#include "libmesh/dof_map.h"
#include "libmesh/equation_systems.h"
#include "libmesh/int_range.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/parallel.h"

namespace libMesh
{

CondensedEigenSystem::CondensedEigenSystem (EquationSystems & es,
                                            const std::string & name_in,
                                            const unsigned int number_in)
  : Parent(es, name_in, number_in),
    condensed_dofs_initialized(false),
    _create_submatrices_in_solve(true)
{
}

CondensedEigenSystem::~CondensedEigenSystem() = default;

void
CondensedEigenSystem::initialize_condensed_dofs(const std::set<dof_id_type> & global_dirichlet_dofs_set)
{
  const DofMap & dof_map = this->get_dof_map();
  if (global_dirichlet_dofs_set.empty() && !dof_map.n_constrained_dofs())
    return;

  // First, put all unconstrained local dofs into non_dirichlet_dofs_set
  std::set<dof_id_type> local_non_condensed_dofs_set;
  for (auto i : make_range(dof_map.first_dof(), dof_map.end_dof()))
#if LIBMESH_ENABLE_CONSTRAINTS
    if (!dof_map.is_constrained_dof(i))
#endif
      local_non_condensed_dofs_set.insert(i);

  // Now erase the condensed dofs
  for (const auto & dof : global_dirichlet_dofs_set)
    if ((dof_map.first_dof() <= dof) && (dof < dof_map.end_dof()))
      local_non_condensed_dofs_set.erase(dof);

  // Finally, move local_non_condensed_dofs_set over to a vector for convenience in solve()
  this->local_non_condensed_dofs_vector.clear();

  for (const auto & dof : local_non_condensed_dofs_set)
    this->local_non_condensed_dofs_vector.push_back(dof);

  condensed_dofs_initialized = true;
}

dof_id_type CondensedEigenSystem::n_global_non_condensed_dofs() const
{
  if (!condensed_dofs_initialized)
    {
      return this->n_dofs();
    }
  else
    {
      dof_id_type n_global_non_condensed_dofs =
        cast_int<dof_id_type>(local_non_condensed_dofs_vector.size());
      this->comm().sum(n_global_non_condensed_dofs);

      return n_global_non_condensed_dofs;
    }
}

#ifdef LIBMESH_ENABLE_DEPRECATED
void
CondensedEigenSystem::set_raw_pointers()
{
  condensed_matrix_A = _condensed_matrix_A.get();
  condensed_matrix_B = _condensed_matrix_B.get();
}
#endif

void
CondensedEigenSystem::clear()
{
  Parent::clear();
  _condensed_matrix_A = nullptr;
  _condensed_matrix_B = nullptr;
  _condensed_precond_matrix = nullptr;
#ifdef LIBMESH_ENABLE_DEPRECATED
  set_raw_pointers();
#endif
}

void
CondensedEigenSystem::add_matrices()
{
  Parent::add_matrices();

  if (!this->use_shell_matrices())
  {
    if (!_condensed_matrix_A)
      _condensed_matrix_A = SparseMatrix<Number>::build(this->comm());
    if (!_condensed_matrix_B)
      _condensed_matrix_B = SparseMatrix<Number>::build(this->comm());

    // When not using shell matrices we use A for P as well so we don't need to add P
  }
  // we *are* using shell matrices but we might not be using a shell matrix for P
  else if (!this->use_shell_precond_matrix())
  {
    if (!_condensed_precond_matrix)
      _condensed_precond_matrix = SparseMatrix<Number>::build(this->comm());
  }

#ifdef LIBMESH_ENABLE_DEPRECATED
  set_raw_pointers();
#endif
}

void CondensedEigenSystem::solve()
{
  LOG_SCOPE("solve()", "CondensedEigenSystem");

  // If we haven't initialized any condensed dofs,
  // just use the default eigen_system
  if (!condensed_dofs_initialized)
    {
      Parent::solve();
      return;
    }

  // check that necessary parameters have been set
  libmesh_assert(
      this->get_equation_systems().parameters.have_parameter<unsigned int>("eigenpairs"));
  libmesh_assert(
      this->get_equation_systems().parameters.have_parameter<unsigned int>("basis vectors"));

  if (this->assemble_before_solve)
    {
      // Assemble the linear system
      this->assemble();

      // And close the assembled matrices; using a non-closed matrix
      // with create_submatrix() is deprecated.
      if (matrix_A)
        matrix_A->close();
      if (generalized() && matrix_B)
        matrix_B->close();
      if (precond_matrix)
        precond_matrix->close();
    }

  // If we reach here, then there should be some non-condensed dofs
  libmesh_assert(!local_non_condensed_dofs_vector.empty());

  if (_create_submatrices_in_solve)
    {
      if (matrix_A)
        matrix_A->create_submatrix(
            *_condensed_matrix_A, local_non_condensed_dofs_vector, local_non_condensed_dofs_vector);
      if (generalized() && matrix_B)
        matrix_B->create_submatrix(
            *_condensed_matrix_B, local_non_condensed_dofs_vector, local_non_condensed_dofs_vector);
      if (precond_matrix)
        precond_matrix->create_submatrix(*_condensed_precond_matrix,
                                         local_non_condensed_dofs_vector,
                                         local_non_condensed_dofs_vector);
    }
  else if (_condensed_precond_matrix.get() && !_condensed_precond_matrix->initialized())
    {
      const auto m = local_non_condensed_dofs_vector.size();
      auto M = m;
      _communicator.sum(M);
      _condensed_precond_matrix->init(M, M, m, m);
    }

  solve_helper(
      _condensed_matrix_A.get(), _condensed_matrix_B.get(), _condensed_precond_matrix.get());
}

std::pair<Real, Real> CondensedEigenSystem::get_eigenpair(dof_id_type i)
{
  LOG_SCOPE("get_eigenpair()", "CondensedEigenSystem");

  // If we haven't initialized any condensed dofs,
  // just use the default eigen_system
  if (!condensed_dofs_initialized)
    return Parent::get_eigenpair(i);

  // If we reach here, then there should be some non-condensed dofs
  libmesh_assert(!local_non_condensed_dofs_vector.empty());

  // This function assumes that condensed_solve has just been called.
  // If this is not the case, then we will trip an asset in get_eigenpair
  std::unique_ptr<NumericVector<Number>> temp = NumericVector<Number>::build(this->comm());
  const dof_id_type n_local =
    cast_int<dof_id_type>(local_non_condensed_dofs_vector.size());
  dof_id_type n       = n_local;
  this->comm().sum(n);

  temp->init (n, n_local, false, PARALLEL);

  std::pair<Real, Real> eval = eigen_solver->get_eigenpair (i, *temp);

  // Now map temp to solution. Loop over local entries of local_non_condensed_dofs_vector
  this->solution->zero();
  for (auto j : make_range(n_local))
    {
      const dof_id_type index = local_non_condensed_dofs_vector[j];
      solution->set(index,(*temp)(temp->first_local_index()+j));
    }

  // Enforcing constraints requires creating a ghosted version of the solution, which requires the
  // solution be assembled
  solution->close();
  get_dof_map().enforce_constraints_exactly(*this);
  this->update();

  return eval;
}



SparseMatrix<Number> & CondensedEigenSystem::get_condensed_matrix_A()
{
  libmesh_assert(_condensed_matrix_A);
  return *_condensed_matrix_A;
}

SparseMatrix<Number> & CondensedEigenSystem::get_condensed_matrix_B()
{
  libmesh_assert(_condensed_matrix_B);
  return *_condensed_matrix_B;
}

SparseMatrix<Number> & CondensedEigenSystem::get_condensed_precond_matrix()
{
  libmesh_assert(_condensed_precond_matrix);
  return *_condensed_precond_matrix;
}

void
CondensedEigenSystem::copy_sub_to_super(const NumericVector<Number> & sub,
                                        NumericVector<Number> & super)
{
  libmesh_assert_equal_to(sub.local_size(), local_non_condensed_dofs_vector.size());
  libmesh_assert_equal_to(sub.local_size() + this->get_dof_map().n_local_constrained_dofs(),
                          super.local_size());
  auto super_sub_view = super.get_subvector(local_non_condensed_dofs_vector);
  *super_sub_view = sub;
  super.restore_subvector(std::move(super_sub_view), local_non_condensed_dofs_vector);
}

void
CondensedEigenSystem::copy_super_to_sub(NumericVector<Number> & super,
                                        NumericVector<Number> & sub)
{
  libmesh_assert_equal_to(sub.local_size(), local_non_condensed_dofs_vector.size());
  libmesh_assert_equal_to(sub.local_size() + this->get_dof_map().n_local_constrained_dofs(),
                          super.local_size());
  auto super_sub_view = super.get_subvector(local_non_condensed_dofs_vector);
  sub = *super_sub_view;
  super.restore_subvector(std::move(super_sub_view), local_non_condensed_dofs_vector);
}

void
CondensedEigenSystem::copy_super_to_sub(const SparseMatrix<Number> & super,
                                        SparseMatrix<Number> & sub)
{
  libmesh_assert_equal_to(sub.local_m(), local_non_condensed_dofs_vector.size());
  libmesh_assert_equal_to(sub.local_m() + this->get_dof_map().n_local_constrained_dofs(),
                          super.local_m());
  auto super_sub_view = SparseMatrix<Number>::build(super.comm());
  super.create_submatrix(
      *super_sub_view, local_non_condensed_dofs_vector, local_non_condensed_dofs_vector);
  sub = *super_sub_view;
}

} // namespace libMesh


#endif // LIBMESH_HAVE_SLEPC
