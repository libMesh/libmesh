// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#include "libmesh/enum_parallel_type.h"
#include "libmesh/static_condensation.h"

#if defined(LIBMESH_HAVE_EIGEN) && defined(LIBMESH_HAVE_PETSC)

#include "libmesh/mesh_base.h"
#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/int_range.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/linear_solver.h"
#include "libmesh/static_condensation_preconditioner.h"
#include "libmesh/system.h"
#include "libmesh/petsc_matrix.h"
#include "libmesh/equation_systems.h"
#include "libmesh/static_condensation_dof_map.h"
#include "timpi/parallel_sync.h"
#include <unordered_set>

namespace libMesh
{
StaticCondensation::StaticCondensation(const MeshBase & mesh,
                                       System & system,
                                       const DofMap & full_dof_map,
                                       const StaticCondensationDofMap & reduced_dof_map)
  : PetscMatrixShellMatrix<Number>(full_dof_map.comm()),
    _mesh(mesh),
    _system(system),
    _full_dof_map(full_dof_map),
    _reduced_dof_map(reduced_dof_map),
    _current_elem_id(DofObject::invalid_id),
    _sc_is_initialized(false),
    _have_cached_values(false),
    _parallel_type(INVALID_PARALLELIZATION)
{
  _size_one_mat.resize(1, 1);
  _scp = std::make_unique<StaticCondensationPreconditioner>(*this);
}

StaticCondensation::~StaticCondensation() = default;

SparseMatrix<Number> & StaticCondensation::operator=(const SparseMatrix<Number> &)
{
  libmesh_not_implemented();
}

std::unique_ptr<SparseMatrix<Number>> StaticCondensation::zero_clone() const
{
  libmesh_not_implemented();
}

std::unique_ptr<SparseMatrix<Number>> StaticCondensation::clone() const
{
  libmesh_not_implemented();
}

void StaticCondensation::clear() noexcept
{
  PetscMatrixShellMatrix<Number>::clear();

  _elem_to_matrix_data.clear();
  _reduced_sys_mat.reset();
  _reduced_sol.reset();
  _reduced_rhs.reset();
  _reduced_solver.reset();
  _current_elem_id = DofObject::invalid_id;
  _have_cached_values = false;
  _sc_is_initialized = false;
}

void StaticCondensation::init(const numeric_index_type m,
                              const numeric_index_type n,
                              const numeric_index_type m_l,
                              const numeric_index_type n_l,
                              const numeric_index_type nnz,
                              const numeric_index_type noz,
                              const numeric_index_type blocksize)
{
  if (!this->initialized())
    {
      PetscMatrixShellMatrix<Number>::init(m, n, m_l, n_l, nnz, noz, blocksize);
      _parallel_type = ((m == m_l) && (this->n_processors() > 1)) ? SERIAL : PARALLEL;
      this->init();
    }
}

void StaticCondensation::init(const ParallelType type)
{
  if (!this->initialized())
    {
      PetscMatrixShellMatrix<Number>::init(type);
      _parallel_type = type;
      this->init();
    }
}

bool StaticCondensation::initialized() const
{
  return PetscMatrixShellMatrix<Number>::initialized() && _sc_is_initialized;
}

void StaticCondensation::init()
{
  if (this->initialized())
    return;

  libmesh_assert(_reduced_dof_map.initialized());

  // This API is public, so it can be called without going through the other init overloads which
  // would give us an indication of what kind of parallel type the user wants. If we've gotten no
  // indication so far, we will default to parallel
  if (_parallel_type == INVALID_PARALLELIZATION)
    _parallel_type = PARALLEL;
  libmesh_assert(_parallel_type != SERIAL);

  for (const auto & [elem_id, dof_data] : _reduced_dof_map._elem_to_dof_data)
    {
      auto & matrix_data = _elem_to_matrix_data[elem_id];

      const auto condensed_dof_size = dof_data.condensed_global_to_local_map.size();
      const auto uncondensed_dof_size = dof_data.uncondensed_global_to_local_map.size();

      matrix_data.Acc.setZero(condensed_dof_size, condensed_dof_size);
      matrix_data.Acu.setZero(condensed_dof_size, uncondensed_dof_size);
      matrix_data.Auc.setZero(uncondensed_dof_size, condensed_dof_size);
      matrix_data.Auu.setZero(uncondensed_dof_size, uncondensed_dof_size);
    }

  //
  // Build the reduced system data
  //
  const auto n = _reduced_dof_map.n_dofs();
  const auto n_local =
      (_parallel_type == SERIAL) ? _reduced_dof_map.n_dofs() : _reduced_dof_map.n_local_dofs();
  _reduced_solver = LinearSolver<Number>::build(this->comm());
  _reduced_solver->init((_system.name() + "_condensed_").c_str());
  _reduced_rhs = NumericVector<Number>::build(this->comm());
  // Init the RHS vector so we can conveniently get processor row offsets
  _reduced_rhs->init(n, n_local);

  // Initialize the reduced system matrix
  _sp = _reduced_dof_map._reduced_sp.get();
  _reduced_sys_mat = SparseMatrix<Number>::build(this->comm());
  if (auto * const petsc_mat = dynamic_cast<PetscMatrix<Number> *>(_reduced_sys_mat.get()))
    {
      // Optimization for PETSc. This is critical for problems in which there are SCALAR dofs that
      // introduce dense rows to avoid allocating a dense matrix
      petsc_mat->init(
          n, n, n_local, n_local, _reduced_dof_map._reduced_nnz, _reduced_dof_map._reduced_noz);
    }
  else
    {
      const auto & nnz = _sp->get_n_nz();
      const auto & noz = _sp->get_n_oz();
      const auto nz = nnz.empty() ? dof_id_type(0) : *std::max_element(nnz.begin(), nnz.end());
      const auto oz = noz.empty() ? dof_id_type(0) : *std::max_element(noz.begin(), noz.end());
      _reduced_sys_mat->init(n, n, n_local, n_local, nz, oz);
    }

  // Build ghosted full solution vector. Note that this is, in general, *not equal* to the system
  // solution, e.g. this may correspond to the solution for the Newton *update*
  _ghosted_full_sol = _system.current_local_solution->clone();

  _sc_is_initialized = true;
}

void StaticCondensation::close()
{
  _communicator.max(_have_cached_values);
  if (!_have_cached_values)
    {
      const bool closed = _reduced_sys_mat->closed();
      libmesh_assert(_communicator.verify(closed));
      if (!closed)
        _reduced_sys_mat->close();
      return;
    }

  DenseMatrix<Number> shim;
  std::vector<dof_id_type> reduced_space_indices;
  for (auto & [elem_id, matrix_data] : _elem_to_matrix_data)
    {
      const auto & dof_data = libmesh_map_find(_reduced_dof_map._elem_to_dof_data, elem_id);
      reduced_space_indices.clear();
      matrix_data.AccFactor = matrix_data.Acc.partialPivLu();
      const EigenMatrix S =
          matrix_data.Auu - matrix_data.Auc * matrix_data.AccFactor.solve(matrix_data.Acu);
      shim.resize(S.rows(), S.cols());
      for (const auto i : make_range(S.rows()))
        for (const auto j : make_range(S.cols()))
          shim(i, j) = S(i, j);
      for (const auto & var_reduced_space_indices : dof_data.reduced_space_indices)
        reduced_space_indices.insert(reduced_space_indices.end(),
                                     var_reduced_space_indices.begin(),
                                     var_reduced_space_indices.end());
      _reduced_sys_mat->add_matrix(shim, reduced_space_indices);
    }

  _reduced_sys_mat->close();

  _have_cached_values = false;
}

bool StaticCondensation::closed() const
{
  return _reduced_sys_mat->closed() && !_have_cached_values;
}

void StaticCondensation::zero()
{
  _reduced_sys_mat->zero();
  for (auto & [elem_id, matrix_data] : _elem_to_matrix_data)
    {
      libmesh_ignore(elem_id);
      matrix_data.Acc.setZero();
      matrix_data.Acu.setZero();
      matrix_data.Auc.setZero();
      matrix_data.Auu.setZero();
    }
}

void StaticCondensation::setup() { libmesh_assert(this->closed()); }

numeric_index_type StaticCondensation::m() const { return _full_dof_map.n_dofs(); }

numeric_index_type StaticCondensation::row_start() const { return _full_dof_map.first_dof(); }

numeric_index_type StaticCondensation::row_stop() const { return _full_dof_map.end_dof(); }

void StaticCondensation::set(const numeric_index_type, const numeric_index_type, const Number)
{
  libmesh_not_implemented();
}

void StaticCondensation::set_current_elem(const Elem & elem)
{
  libmesh_assert(!Threads::in_threads || libMesh::n_threads() == 1);
  _current_elem_id = elem.id();
}

void StaticCondensation::add(const numeric_index_type i,
                             const numeric_index_type j,
                             const Number value)
{
  _size_one_mat(0, 0) = value;
  this->add_matrix(_size_one_mat, {i}, {j});
}

void StaticCondensation::add_matrix(const DenseMatrix<Number> & dm,
                                    const std::vector<numeric_index_type> & rows,
                                    const std::vector<numeric_index_type> & cols)
{
  if (rows.empty() || cols.empty())
    return;

  libmesh_assert(_current_elem_id != DofObject::invalid_id);
  auto & matrix_data = libmesh_map_find(_elem_to_matrix_data, _current_elem_id);
  const auto & dof_data = libmesh_map_find(_reduced_dof_map._elem_to_dof_data, _current_elem_id);
  EigenMatrix * mat;

  auto info_from_index = [&dof_data](const auto global_index) {
    auto index_it = dof_data.condensed_global_to_local_map.find(global_index);
    const bool index_is_condensed = index_it != dof_data.condensed_global_to_local_map.end();
    if (!index_is_condensed)
      {
        index_it = dof_data.uncondensed_global_to_local_map.find(global_index);
        libmesh_assert(index_it != dof_data.uncondensed_global_to_local_map.end());
      }
    return std::make_pair(index_is_condensed, index_it->second);
  };

  for (const auto i : make_range(dm.m()))
    for (const auto j : make_range(dm.n()))
      {
        const auto global_i = rows[i];
        const auto global_j = cols[j];
        const auto [i_is_condensed, local_i] = info_from_index(global_i);
        const auto [j_is_condensed, local_j] = info_from_index(global_j);
        if (i_is_condensed)
          {
            if (j_is_condensed)
              mat = &matrix_data.Acc;
            else
              mat = &matrix_data.Acu;
          }
        else
          {
            if (j_is_condensed)
              mat = &matrix_data.Auc;
            else
              mat = &matrix_data.Auu;
          }
        (*mat)(local_i, local_j) += dm(i, j);
      }

  _have_cached_values = true;
}

void StaticCondensation::add_matrix(const DenseMatrix<Number> & dm,
                                    const std::vector<numeric_index_type> & dof_indices)
{
  this->add_matrix(dm, dof_indices, dof_indices);
}

void StaticCondensation::add(const Number, const SparseMatrix<Number> &)
{
  libmesh_not_implemented();
}

Number StaticCondensation::operator()(const numeric_index_type, const numeric_index_type) const
{
  libmesh_not_implemented();
}

Real StaticCondensation::l1_norm() const { libmesh_not_implemented(); }

Real StaticCondensation::linfty_norm() const { libmesh_not_implemented(); }

void StaticCondensation::print_personal(std::ostream &) const { libmesh_not_implemented(); }

void StaticCondensation::get_diagonal(NumericVector<Number> &) const { libmesh_not_implemented(); }

void StaticCondensation::get_transpose(SparseMatrix<Number> &) const { libmesh_not_implemented(); }

void StaticCondensation::get_row(numeric_index_type,
                                 std::vector<numeric_index_type> &,
                                 std::vector<Number> &) const
{
  libmesh_not_implemented();
}

void StaticCondensation::set_local_vectors(const NumericVector<Number> & global_vector,
                                           const std::vector<dof_id_type> & elem_dof_indices,
                                           std::vector<Number> & elem_dof_values_vec,
                                           EigenVector & elem_dof_values)
{
  global_vector.get(elem_dof_indices, elem_dof_values_vec);
  elem_dof_values.resize(elem_dof_indices.size());
  for (const auto i : index_range(elem_dof_indices))
    elem_dof_values(i) = elem_dof_values_vec[i];
}

void StaticCondensation::forward_elimination(const NumericVector<Number> & full_rhs)
{
  std::vector<dof_id_type> elem_condensed_dofs;
  std::vector<Number> elem_condensed_rhs_vec;
  EigenVector elem_condensed_rhs, elem_uncondensed_rhs;

  full_rhs.create_subvector(
      *_reduced_rhs, _reduced_dof_map._local_uncondensed_dofs, /*all_global_entries=*/false);

  std::vector<dof_id_type> reduced_space_indices;
  for (auto elem : _mesh.active_local_element_ptr_range())
    {
      auto & matrix_data = libmesh_map_find(_elem_to_matrix_data, elem->id());
      reduced_space_indices.clear();
      const auto & dof_data = libmesh_map_find(_reduced_dof_map._elem_to_dof_data, elem->id());
      for (const auto & var_reduced_space_indices : dof_data.reduced_space_indices)
        reduced_space_indices.insert(reduced_space_indices.end(),
                                     var_reduced_space_indices.begin(),
                                     var_reduced_space_indices.end());
      elem_condensed_dofs.resize(dof_data.condensed_global_to_local_map.size());
      for (const auto & [global_dof, local_dof] : dof_data.condensed_global_to_local_map)
        {
          libmesh_assert(local_dof < elem_condensed_dofs.size());
          elem_condensed_dofs[local_dof] = global_dof;
        }

      set_local_vectors(full_rhs, elem_condensed_dofs, elem_condensed_rhs_vec, elem_condensed_rhs);
      elem_uncondensed_rhs = -matrix_data.Auc * matrix_data.AccFactor.solve(elem_condensed_rhs);

      libmesh_assert(cast_int<std::size_t>(elem_uncondensed_rhs.size()) ==
                     reduced_space_indices.size());
      _reduced_rhs->add_vector(elem_uncondensed_rhs.data(), reduced_space_indices);
    }
  _reduced_rhs->close();
}

void StaticCondensation::backwards_substitution(const NumericVector<Number> & full_rhs,
                                                NumericVector<Number> & full_sol)
{
  std::vector<dof_id_type> elem_condensed_dofs, elem_uncondensed_dofs;
  std::vector<Number> elem_condensed_rhs_vec, elem_uncondensed_sol_vec;
  EigenVector elem_condensed_rhs, elem_uncondensed_sol, elem_condensed_sol;

  for (auto elem : _mesh.active_local_element_ptr_range())
    {
      auto & matrix_data = libmesh_map_find(_elem_to_matrix_data, elem->id());
      const auto & dof_data = libmesh_map_find(_reduced_dof_map._elem_to_dof_data, elem->id());
      elem_condensed_dofs.resize(dof_data.condensed_global_to_local_map.size());
      elem_uncondensed_dofs.resize(dof_data.uncondensed_global_to_local_map.size());
      for (const auto & [global_dof, local_dof] : dof_data.condensed_global_to_local_map)
        {
          libmesh_assert(local_dof < elem_condensed_dofs.size());
          elem_condensed_dofs[local_dof] = global_dof;
        }
      for (const auto & [global_dof, local_dof] : dof_data.uncondensed_global_to_local_map)
        {
          libmesh_assert(local_dof < elem_uncondensed_dofs.size());
          elem_uncondensed_dofs[local_dof] = global_dof;
        }

      set_local_vectors(full_rhs, elem_condensed_dofs, elem_condensed_rhs_vec, elem_condensed_rhs);
      set_local_vectors(*_ghosted_full_sol,
                        elem_uncondensed_dofs,
                        elem_uncondensed_sol_vec,
                        elem_uncondensed_sol);

      elem_condensed_sol =
          matrix_data.AccFactor.solve(elem_condensed_rhs - matrix_data.Acu * elem_uncondensed_sol);
      full_sol.insert(elem_condensed_sol.data(), elem_condensed_dofs);
    }

  full_sol.close();
}

void StaticCondensation::apply(const NumericVector<Number> & full_rhs,
                               NumericVector<Number> & full_parallel_sol)
{
  forward_elimination(full_rhs);
  // Apparently PETSc will send us the yvec without zeroing it ahead of time. This can be a poor
  // initial guess for the Krylov solve as well as lead to bewildered users who expect their initial
  // residual norm to equal the norm of the RHS
  full_parallel_sol.zero();
  _reduced_sol = full_parallel_sol.get_subvector(_reduced_dof_map._local_uncondensed_dofs);
  _reduced_solver->solve(*_reduced_sys_mat, *_reduced_sol, *_reduced_rhs, 1e-5, 300);
  // Must restore to the full solution because during backwards substitution we will need to be able
  // to read ghosted dofs and we don't support ghosting of subvectors
  full_parallel_sol.restore_subvector(std::move(_reduced_sol),
                                      _reduced_dof_map._local_uncondensed_dofs);
  *_ghosted_full_sol = full_parallel_sol;
  backwards_substitution(full_rhs, full_parallel_sol);
}

SolverPackage StaticCondensation::solver_package() { return libMesh::default_solver_package(); }

}
#else

#include "libmesh/dof_map.h"

namespace libMesh
{
StaticCondensation::StaticCondensation(const MeshBase &,
                                       const System &,
                                       const DofMap & full_dof_map)
  : SparseMatrix<Number>(full_dof_map.comm())
{
  libmesh_error_msg(
      "Static condensation requires configuring libMesh with PETSc and Eigen support");
}
}

#endif // LIBMESH_HAVE_EIGEN && LIBMESH_HAVE_PETSC
