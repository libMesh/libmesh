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

#include "libmesh/static_condensation.h"

#if defined(LIBMESH_HAVE_EIGEN) && defined(LIBMESH_HAVE_PETSC)

#include "libmesh/mesh_base.h"
#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/int_range.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/petsc_matrix_base.h"
#include "libmesh/linear_solver.h"
#include "libmesh/static_condensation_preconditioner.h"
#include "libmesh/system.h"
#include "libmesh/petsc_matrix.h"
#include "timpi/parallel_sync.h"
#include <unordered_set>

namespace libMesh
{
StaticCondensation::StaticCondensation(const MeshBase & mesh,
                                       const System & system,
                                       const DofMap & dof_map)
  : PetscMatrixShellMatrix<Number>(dof_map.comm()),
    _mesh(mesh),
    _system(system),
    _dof_map(dof_map),
    _current_elem_id(DofObject::invalid_id),
    _sc_is_initialized(false),
    _have_cached_values(false)
{
  _size_one_mat.resize(1, 1);
  _scp = std::make_unique<StaticCondensationPreconditioner>(*this);
}

StaticCondensation::~StaticCondensation() = default;

SparseMatrix<Number> &
StaticCondensation::operator=(const SparseMatrix<Number> &)
{
  libmesh_not_implemented();
}

std::unique_ptr<SparseMatrix<Number>>
StaticCondensation::zero_clone() const
{
  libmesh_not_implemented();
}

std::unique_ptr<SparseMatrix<Number>>
StaticCondensation::clone() const
{
  libmesh_not_implemented();
}

void
StaticCondensation::clear() noexcept
{
  PetscMatrixShellMatrix<Number>::clear();

  _elem_to_local_data.clear();
  _local_uncondensed_dofs.clear();
  _reduced_sys_mat.reset();
  _reduced_sol.reset();
  _reduced_rhs.reset();
  _reduced_solver.reset();
  _current_elem_id = DofObject::invalid_id;
  _have_cached_values = false;
  _sc_is_initialized = false;
}

void
StaticCondensation::init(const numeric_index_type m,
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
      this->init();
    }
}

void
StaticCondensation::init(const ParallelType type)
{
  if (!this->initialized())
    {
      PetscMatrixShellMatrix<Number>::init(type);
      this->init();
    }
}

bool
StaticCondensation::initialized() const
{
  return PetscMatrixShellMatrix<Number>::initialized() && _sc_is_initialized;
}

void
StaticCondensation::init()
{
  if (_sc_is_initialized)
    return;

  std::vector<dof_id_type> elem_dofs; // only used to satisfy API
  std::vector<dof_id_type> elem_uncondensed_dofs;
  std::unordered_set<dof_id_type> local_uncondensed_dofs;
  std::unordered_map<processor_id_type, std::unordered_set<dof_id_type>> nonlocal_uncondensed_dofs;
  dof_id_type condensed_local_dof_number = 0, uncondensed_local_dof_number = 0;
  std::unordered_map<dof_id_type, dof_id_type> *condensed_global_to_local_map = nullptr,
                                               *uncondensed_global_to_local_map = nullptr;

  // Handle SCALAR dofs
  for (const auto vg : make_range(_dof_map.n_variable_groups()))
    if (const auto & vg_description = _dof_map.variable_group(vg);
        vg_description.type().family == SCALAR)
      {
        std::vector<dof_id_type> scalar_dof_indices;
        const processor_id_type last_pid = this->comm().size() - 1;
        for (const auto vg_vn : make_range(vg_description.n_variables()))
          {
            const auto vn = vg_description.number(vg_vn);
            _dof_map.SCALAR_dof_indices(scalar_dof_indices, vn);
            if (this->comm().rank() == last_pid)
              local_uncondensed_dofs.insert(scalar_dof_indices.begin(), scalar_dof_indices.end());
            else
              nonlocal_uncondensed_dofs[last_pid].insert(scalar_dof_indices.begin(),
                                                         scalar_dof_indices.end());
          }
      }

  auto scalar_dofs_functor =
      [&elem_uncondensed_dofs, &uncondensed_local_dof_number, &uncondensed_global_to_local_map](
          const Elem & /*elem*/,
          std::vector<dof_id_type> & dof_indices,
          const std::vector<dof_id_type> & scalar_dof_indices) {
        dof_indices.insert(dof_indices.end(), scalar_dof_indices.begin(), scalar_dof_indices.end());
        elem_uncondensed_dofs.insert(
            elem_uncondensed_dofs.end(), scalar_dof_indices.begin(), scalar_dof_indices.end());
        for (const auto global_dof : scalar_dof_indices)
          (*uncondensed_global_to_local_map)[global_dof] = uncondensed_local_dof_number++;
      };

  auto field_dofs_functor = [this,
                             &local_uncondensed_dofs,
                             &nonlocal_uncondensed_dofs,
                             &elem_uncondensed_dofs,
                             &uncondensed_local_dof_number,
                             &condensed_local_dof_number,
                             &uncondensed_global_to_local_map,
                             &condensed_global_to_local_map](const Elem & elem,
                                                             const unsigned int node_num,
                                                             const unsigned int var_num,
                                                             std::vector<dof_id_type> & dof_indices,
                                                             const dof_id_type field_dof) {
    dof_indices.push_back(field_dof);

    bool uncondensed_dof = false;
    if (_uncondensed_vars.count(var_num))
      {
        libmesh_assert_msg(
            node_num == invalid_uint,
            "Users should not be providing continuous FEM variables to the uncondensed vars API");
        uncondensed_dof = true;
        elem_uncondensed_dofs.push_back(field_dof);
        if (elem.processor_id() == this->processor_id())
          local_uncondensed_dofs.insert(field_dof);
        else
          nonlocal_uncondensed_dofs[elem.processor_id()].insert(field_dof);
      }

    if (node_num != invalid_uint && !elem.is_internal(node_num))
      {
        uncondensed_dof = true;
        elem_uncondensed_dofs.push_back(field_dof);
        const auto & nd_ref = elem.node_ref(node_num);
        if (nd_ref.processor_id() == this->processor_id())
          local_uncondensed_dofs.insert(field_dof);
        else
          nonlocal_uncondensed_dofs[nd_ref.processor_id()].insert(field_dof);
      }

    if (uncondensed_dof)
      (*uncondensed_global_to_local_map)[field_dof] = uncondensed_local_dof_number++;
    else
      (*condensed_global_to_local_map)[field_dof] = condensed_local_dof_number++;
  };

  for (auto elem : _mesh.active_local_element_ptr_range())
    {
      auto & local_data = _elem_to_local_data[elem->id()];
      condensed_local_dof_number = 0;
      uncondensed_local_dof_number = 0;
      condensed_global_to_local_map = &local_data.condensed_global_to_local_map;
      uncondensed_global_to_local_map = &local_data.uncondensed_global_to_local_map;

      const auto sub_id = elem->subdomain_id();
      for (const auto vg : make_range(_dof_map.n_variable_groups()))
        {
          const auto & var_group = _dof_map.variable_group(vg);
          if (!var_group.active_on_subdomain(sub_id))
            continue;

          for (const auto v : make_range(var_group.n_variables()))
            {
              const auto var_num = var_group.number(v);
              elem_uncondensed_dofs.clear();
              _dof_map.dof_indices(elem,
                                   elem_dofs,
                                   var_num,
                                   scalar_dofs_functor,
                                   field_dofs_functor,
                                   elem->p_level());
              local_data.reduced_space_indices.insert(local_data.reduced_space_indices.end(),
                                                      elem_uncondensed_dofs.begin(),
                                                      elem_uncondensed_dofs.end());
            }
        }

      const auto condensed_dof_size = condensed_global_to_local_map->size();
      const auto uncondensed_dof_size = uncondensed_global_to_local_map->size();

      local_data.Acc.setZero(condensed_dof_size, condensed_dof_size);
      local_data.Acu.setZero(condensed_dof_size, uncondensed_dof_size);
      local_data.Auc.setZero(uncondensed_dof_size, condensed_dof_size);
      local_data.Auu.setZero(uncondensed_dof_size, uncondensed_dof_size);
    }

  _local_uncondensed_dofs.assign(local_uncondensed_dofs.begin(), local_uncondensed_dofs.end());
  local_uncondensed_dofs.clear();

  //
  // Build the reduced system data
  //

  const dof_id_type n_local = _local_uncondensed_dofs.size();
  dof_id_type n = n_local;
  this->comm().sum(n);
  _reduced_solver = LinearSolver<Number>::build(this->comm());
  _reduced_solver->init((_system.name() + "_condensed_").c_str());
  _reduced_rhs = NumericVector<Number>::build(this->comm());
  // Init the RHS vector so we can conveniently get processor row offsets
  _reduced_rhs->init(n, n_local);

  // Build a map from the full size problem uncondensed dof indices to the reduced problem
  // (uncondensed) dof indices
  std::unordered_map<dof_id_type, dof_id_type> full_dof_to_reduced_dof;
  const auto local_start = _reduced_rhs->first_local_index();
  for (const auto i : index_range(_local_uncondensed_dofs))
    full_dof_to_reduced_dof[_local_uncondensed_dofs[i]] = i + local_start;

  _reduced_sys_mat = SparseMatrix<Number>::build(this->comm());
  const auto & nnz = _sp->get_n_nz();
  const auto & noz = _sp->get_n_oz();
  libmesh_assert(nnz.size() == noz.size());
  if (auto * const petsc_mat = dynamic_cast<PetscMatrix<Number> *>(_reduced_sys_mat.get()))
    {
      // Optimization for PETSc. This is critical for problems in which there are SCALAR dofs that
      // introduce dense rows to avoid allocating a dense matrix
      std::vector<dof_id_type> reduced_nnz, reduced_noz;
      reduced_nnz.resize(_local_uncondensed_dofs.size());
      reduced_noz.resize(_local_uncondensed_dofs.size());
      for (const dof_id_type local_reduced_i : index_range(_local_uncondensed_dofs))
        {
          const dof_id_type full_i = _local_uncondensed_dofs[local_reduced_i];
          const dof_id_type local_full_i = full_i - _dof_map.first_dof();
          libmesh_assert(local_full_i < nnz.size());
          reduced_nnz[local_reduced_i] = nnz[local_full_i];
          reduced_noz[local_reduced_i] = noz[local_full_i];
        }
      petsc_mat->init(n, n, n_local, n_local, reduced_nnz, reduced_noz);
    }
  else
    {
      const auto nz = nnz.empty() ? dof_id_type(0) : *std::max_element(nnz.begin(), nnz.end());
      const auto oz = noz.empty() ? dof_id_type(0) : *std::max_element(noz.begin(), noz.end());
      _reduced_sys_mat->init(n, n, n_local, n_local, nz, oz);
    }

  //
  // Now we need to pull our nonlocal data
  //

  // build our queries
  std::unordered_map<processor_id_type, std::vector<dof_id_type>> nonlocal_uncondensed_dofs_mapvec;
  for (const auto & [pid, set] : nonlocal_uncondensed_dofs)
    {
      auto & vec = nonlocal_uncondensed_dofs_mapvec[pid];
      vec.assign(set.begin(), set.end());
    }
  // clear no longer needed memory
  nonlocal_uncondensed_dofs.clear();

  auto gather_functor = [&full_dof_to_reduced_dof](processor_id_type,
                                                   const std::vector<dof_id_type> & full_dof_ids,
                                                   std::vector<dof_id_type> & reduced_dof_ids) {
    reduced_dof_ids.resize(full_dof_ids.size());
    for (const auto i : index_range(full_dof_ids))
      reduced_dof_ids[i] = libmesh_map_find(full_dof_to_reduced_dof, full_dof_ids[i]);
  };

  auto action_functor =
      [&full_dof_to_reduced_dof](processor_id_type,
                                 const std::vector<dof_id_type> & full_dof_ids,
                                 const std::vector<dof_id_type> & reduced_dof_ids) {
        for (const auto i : index_range(full_dof_ids))
          {
            libmesh_assert(!full_dof_to_reduced_dof.count(full_dof_ids[i]));
            full_dof_to_reduced_dof[full_dof_ids[i]] = reduced_dof_ids[i];
          }
      };

  TIMPI::pull_parallel_vector_data(this->comm(),
                                   nonlocal_uncondensed_dofs_mapvec,
                                   gather_functor,
                                   action_functor,
                                   &DofObject::invalid_id);
  nonlocal_uncondensed_dofs_mapvec.clear();

  // Now we can finally set our element reduced dof indices
  std::vector<dof_id_type> full_dof_indices;
  for (auto & [elem, local_data] : _elem_to_local_data)
    {
      libmesh_ignore(elem);
      full_dof_indices = local_data.reduced_space_indices;
      local_data.reduced_space_indices.clear();
      for (const auto full_dof : full_dof_indices)
        local_data.reduced_space_indices.push_back(
            libmesh_map_find(full_dof_to_reduced_dof, full_dof));
    }

  // Build ghosted full solution vector. Note that this is, in general, *not equal* to the system
  // solution, e.g. this may correspond to the solution for the Newton *update*
  _ghosted_full_sol = _system.current_local_solution->clone();
  _sc_is_initialized = true;
}

void
StaticCondensation::close()
{
  _communicator.max(_have_cached_values);
  if (!_have_cached_values)
    {
      const bool closed = _reduced_sys_mat->closed();
#ifndef NDEBUG
      _communicator.verify(closed);
#endif
      if (!closed)
        _reduced_sys_mat->close();
      return;
    }

  DenseMatrix<Number> shim;
  for (auto & [elem_id, local_data] : _elem_to_local_data)
    {
      libmesh_ignore(elem_id);
      local_data.AccFactor = local_data.Acc.partialPivLu();
      const EigenMatrix S =
          local_data.Auu - local_data.Auc * local_data.AccFactor.solve(local_data.Acu);
      shim.resize(S.rows(), S.cols());
      for (const auto i : make_range(S.rows()))
        for (const auto j : make_range(S.cols()))
          shim(i, j) = S(i, j);
      _reduced_sys_mat->add_matrix(shim, local_data.reduced_space_indices);
    }

  _reduced_sys_mat->close();

  _have_cached_values = false;
}

bool
StaticCondensation::closed() const
{
  return _reduced_sys_mat->closed() && !_have_cached_values;
}

void
StaticCondensation::zero()
{
  _reduced_sys_mat->zero();
  for (auto & [elem_id, local_data] : _elem_to_local_data)
    {
      libmesh_ignore(elem_id);
      local_data.Acc.setZero();
      local_data.Acu.setZero();
      local_data.Auc.setZero();
      local_data.Auu.setZero();
    }
}

void
StaticCondensation::setup()
{
  libmesh_assert(this->closed());
}

numeric_index_type
StaticCondensation::m() const
{
  return _dof_map.n_dofs();
}

numeric_index_type
StaticCondensation::row_start() const
{
  return _dof_map.first_dof();
}

numeric_index_type
StaticCondensation::row_stop() const
{
  return _dof_map.end_dof();
}

void
StaticCondensation::set(const numeric_index_type, const numeric_index_type, const Number)
{
  libmesh_not_implemented();
}

void
StaticCondensation::set_current_elem(const Elem & elem)
{
  libmesh_assert(!Threads::in_threads || libMesh::n_threads() == 1);
  _current_elem_id = elem.id();
}

void
StaticCondensation::add(const numeric_index_type i, const numeric_index_type j, const Number value)
{
  _size_one_mat(0, 0) = value;
  this->add_matrix(_size_one_mat, {i}, {j});
}

void
StaticCondensation::add_matrix(const DenseMatrix<Number> & dm,
                               const std::vector<numeric_index_type> & rows,
                               const std::vector<numeric_index_type> & cols)
{
  libmesh_assert(_current_elem_id != DofObject::invalid_id);
  auto & local_data = libmesh_map_find(_elem_to_local_data, _current_elem_id);
  EigenMatrix * mat;

  auto info_from_index = [&local_data](const auto global_index) {
    auto index_it = local_data.condensed_global_to_local_map.find(global_index);
    const bool index_is_condensed = index_it != local_data.condensed_global_to_local_map.end();
    if (!index_is_condensed)
      {
        index_it = local_data.uncondensed_global_to_local_map.find(global_index);
        libmesh_assert(index_it != local_data.uncondensed_global_to_local_map.end());
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
              mat = &local_data.Acc;
            else
              mat = &local_data.Acu;
          }
        else
          {
            if (j_is_condensed)
              mat = &local_data.Auc;
            else
              mat = &local_data.Auu;
          }
        (*mat)(local_i, local_j) += dm(i, j);
      }

  _have_cached_values = true;
}

void
StaticCondensation::add_matrix(const DenseMatrix<Number> & dm,
                               const std::vector<numeric_index_type> & dof_indices)
{
  this->add_matrix(dm, dof_indices, dof_indices);
}

void
StaticCondensation::add(const Number, const SparseMatrix<Number> &)
{
  libmesh_not_implemented();
}

Number
StaticCondensation::operator()(const numeric_index_type, const numeric_index_type) const
{
  libmesh_not_implemented();
}

Real
StaticCondensation::l1_norm() const
{
  libmesh_not_implemented();
}

Real
StaticCondensation::linfty_norm() const
{
  libmesh_not_implemented();
}

void
StaticCondensation::print_personal(std::ostream &) const
{
  libmesh_not_implemented();
}

void
StaticCondensation::get_diagonal(NumericVector<Number> &) const
{
  libmesh_not_implemented();
}

void
StaticCondensation::get_transpose(SparseMatrix<Number> &) const
{
  libmesh_not_implemented();
}

void
StaticCondensation::get_row(numeric_index_type,
                            std::vector<numeric_index_type> &,
                            std::vector<Number> &) const
{
  libmesh_not_implemented();
}

void
StaticCondensation::set_local_vectors(const NumericVector<Number> & global_vector,
                                      const std::vector<dof_id_type> & elem_dof_indices,
                                      std::vector<Number> & elem_dof_values_vec,
                                      EigenVector & elem_dof_values)
{
  global_vector.get(elem_dof_indices, elem_dof_values_vec);
  elem_dof_values.resize(elem_dof_indices.size());
  for (const auto i : index_range(elem_dof_indices))
    elem_dof_values(i) = elem_dof_values_vec[i];
}

void
StaticCondensation::forward_elimination(const NumericVector<Number> & full_rhs)
{
  std::vector<dof_id_type> elem_condensed_dofs;
  std::vector<Number> elem_condensed_rhs_vec;
  EigenVector elem_condensed_rhs, elem_uncondensed_rhs;

  full_rhs.create_subvector(*_reduced_rhs, _local_uncondensed_dofs, /*all_global_entries=*/false);

  for (auto elem : _mesh.active_local_element_ptr_range())
    {
      auto & local_data = _elem_to_local_data[elem->id()];
      elem_condensed_dofs.resize(local_data.condensed_global_to_local_map.size());
      for (const auto & [global_dof, local_dof] : local_data.condensed_global_to_local_map)
        {
          libmesh_assert(local_dof < elem_condensed_dofs.size());
          elem_condensed_dofs[local_dof] = global_dof;
        }

      set_local_vectors(full_rhs, elem_condensed_dofs, elem_condensed_rhs_vec, elem_condensed_rhs);
      elem_uncondensed_rhs = -local_data.Auc * local_data.AccFactor.solve(elem_condensed_rhs);

      libmesh_assert(cast_int<std::size_t>(elem_uncondensed_rhs.size()) ==
                     local_data.reduced_space_indices.size());
      _reduced_rhs->add_vector(elem_uncondensed_rhs.data(), local_data.reduced_space_indices);
    }
  _reduced_rhs->close();
}

void
StaticCondensation::backwards_substitution(const NumericVector<Number> & full_rhs,
                                           NumericVector<Number> & full_sol)
{
  std::vector<dof_id_type> elem_condensed_dofs, elem_uncondensed_dofs;
  std::vector<Number> elem_condensed_rhs_vec, elem_uncondensed_sol_vec;
  EigenVector elem_condensed_rhs, elem_uncondensed_sol, elem_condensed_sol;

  for (auto elem : _mesh.active_local_element_ptr_range())
    {
      auto & local_data = _elem_to_local_data[elem->id()];
      elem_condensed_dofs.resize(local_data.condensed_global_to_local_map.size());
      elem_uncondensed_dofs.resize(local_data.uncondensed_global_to_local_map.size());
      for (const auto & [global_dof, local_dof] : local_data.condensed_global_to_local_map)
        {
          libmesh_assert(local_dof < elem_condensed_dofs.size());
          elem_condensed_dofs[local_dof] = global_dof;
        }
      for (const auto & [global_dof, local_dof] : local_data.uncondensed_global_to_local_map)
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
          local_data.AccFactor.solve(elem_condensed_rhs - local_data.Acu * elem_uncondensed_sol);
      full_sol.insert(elem_condensed_sol.data(), elem_condensed_dofs);
    }

  full_sol.close();
}

void
StaticCondensation::apply(const NumericVector<Number> & full_rhs,
                          NumericVector<Number> & full_parallel_sol)
{
  forward_elimination(full_rhs);
  // Apparently PETSc will send us the yvec without zeroing it ahead of time. This can be a poor
  // initial guess for the Krylov solve as well as lead to bewildered users who expect their initial
  // residual norm to equal the norm of the RHS
  full_parallel_sol.zero();
  _reduced_sol = full_parallel_sol.get_subvector(_local_uncondensed_dofs);
  _reduced_solver->solve(*_reduced_sys_mat, *_reduced_sol, *_reduced_rhs, 1e-5, 300);
  // Must restore to the full solution because during backwards substitution we will need to be able
  // to read ghosted dofs and we don't support ghosting of subvectors
  full_parallel_sol.restore_subvector(std::move(_reduced_sol), _local_uncondensed_dofs);
  *_ghosted_full_sol = full_parallel_sol;
  backwards_substitution(full_rhs, full_parallel_sol);
}

SolverPackage
StaticCondensation::solver_package()
{
  return libMesh::default_solver_package();
}

}

#else

#include "libmesh/dof_map.h"

namespace libMesh
{
StaticCondensation::StaticCondensation(const MeshBase &, const System &, const DofMap & dof_map)
  : SparseMatrix<Number>(dof_map.comm())
{
  libmesh_error_msg(
      "Static condensation requires configuring libMesh with PETSc and Eigen support");
}
}

#endif // LIBMESH_HAVE_EIGEN && LIBMESH_HAVE_PETSC
