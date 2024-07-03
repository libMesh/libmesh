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

#include "libmesh/static_condensation.h"

#ifdef LIBMESH_HAVE_PETSC

#include "libmesh/mesh_base.h"
#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/int_range.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/linear_solver.h"
#include "timpi/parallel_sync.h"
#include <unordered_set>

namespace libMesh
{
StaticCondensation::StaticCondensation(const MeshBase & mesh, const DofMap & dof_map)
  : Preconditioner<Number>(dof_map.comm()), _mesh(mesh), _dof_map(dof_map)
{
}

StaticCondensation::~StaticCondensation() = default;

auto
StaticCondensation::total_and_condensed_from_scalar_dofs_functor(
    std::vector<dof_id_type> & /*elem_condensed_dofs*/)
{
  return [](const Elem &,
            std::vector<dof_id_type> & elem_dof_indices,
            const std::vector<dof_id_type> & scalar_dof_indices)
  {
    elem_dof_indices.insert(
        elem_dof_indices.end(), scalar_dof_indices.begin(), scalar_dof_indices.end());
  };
}

auto
StaticCondensation::total_and_condensed_from_field_dofs_functor(
    std::vector<dof_id_type> & elem_condensed_dofs)
{
  return [this, &elem_condensed_dofs](const Elem & elem,
                                      const unsigned int node_num,
                                      const unsigned int var_num,
                                      std::vector<dof_id_type> & elem_dof_indices,
                                      const dof_id_type field_dof)
  {
    elem_dof_indices.push_back(field_dof);
    if (_uncondensed_vars.count(var_num))
      // If we're going to keep all the dofs for this var (e.g. we don't want to condense any out)
      // then we can return now
      return;

    if (node_num != invalid_uint)
    {
      // This is a nodal dof
      if (elem.is_internal(node_num))
        elem_condensed_dofs.push_back(field_dof);
    }
    else
      elem_condensed_dofs.push_back(field_dof);
  };
}

void
StaticCondensation::clear()
{
  _elem_to_local_data.clear();
  _local_uncondensed_dofs.clear();
  _reduced_sys_mat.reset();
  _reduced_sol.reset();
  _reduced_rhs.reset();
  _reduced_solver.reset();
}

void
StaticCondensation::init()
{
  if (_is_initialized)
    return;

  // Some APIs that mark the preconditioner as uninitialized may not clear data, so we do it to be
  // safe
  clear();

  std::vector<dof_id_type> elem_dofs; // only used to satisfy API
  std::vector<dof_id_type> elem_condensed_dofs, elem_uncondensed_dofs;
  std::unordered_set<dof_id_type> local_uncondensed_dofs;
  std::unordered_map<processor_id_type, std::unordered_set<dof_id_type>> nonlocal_uncondensed_dofs;

  const Elem * first_elem = nullptr;
  const auto first_elem_it = _mesh.active_local_elements_begin();
  if (first_elem_it != _mesh.active_local_elements_end())
    first_elem = *first_elem_it;

  auto scalar_dofs_functor = [this,
                              &local_uncondensed_dofs,
                              &nonlocal_uncondensed_dofs,
                              &elem_condensed_dofs,
                              &elem_uncondensed_dofs,
                              first_elem](const Elem & elem,
                                          std::vector<dof_id_type> & dof_indices,
                                          const std::vector<dof_id_type> & scalar_dof_indices)
  {
    total_and_condensed_from_scalar_dofs_functor(elem_condensed_dofs)(
        elem, dof_indices, scalar_dof_indices);
    elem_uncondensed_dofs.insert(
        elem_uncondensed_dofs.end(), scalar_dof_indices.begin(), scalar_dof_indices.end());

    // Only need to do this for the first element we encounter
    if (&elem == first_elem)
    {
      const processor_id_type last_pid = _dof_map.comm().size() - 1;
      if (_dof_map.comm().rank() == last_pid)
        local_uncondensed_dofs.insert(scalar_dof_indices.begin(), scalar_dof_indices.end());
      else
        nonlocal_uncondensed_dofs[last_pid].insert(scalar_dof_indices.begin(),
                                                   scalar_dof_indices.end());
    }
  };

  auto field_dofs_functor = [this,
                             &local_uncondensed_dofs,
                             &nonlocal_uncondensed_dofs,
                             &elem_condensed_dofs,
                             &elem_uncondensed_dofs](const Elem & elem,
                                                     const unsigned int node_num,
                                                     const unsigned int var_num,
                                                     std::vector<dof_id_type> & dof_indices,
                                                     const dof_id_type field_dof)
  {
    total_and_condensed_from_field_dofs_functor(elem_condensed_dofs)(
        elem, node_num, var_num, dof_indices, field_dof);

    if (node_num != invalid_uint && !elem.is_internal(node_num))
    {
      elem_uncondensed_dofs.push_back(field_dof);
      const auto & nd_ref = elem.node_ref(node_num);
      if (nd_ref.processor_id() == _dof_map.processor_id())
        local_uncondensed_dofs.insert(field_dof);
      else
        nonlocal_uncondensed_dofs[nd_ref.processor_id()].insert(field_dof);
    }
  };

  for (auto elem : _mesh.active_local_element_ptr_range())
  {
    unsigned int uncondensed_dof_size = 0;
    unsigned int condensed_dof_size = 0;
    auto & local_data = _elem_to_local_data[elem->id()];

    const auto sub_id = elem->subdomain_id();
    for (const auto vg : make_range(_dof_map.n_variable_groups()))
    {
      const auto & var_group = _dof_map.variable_group(vg);
      if (!var_group.active_on_subdomain(sub_id))
        continue;

      for (const auto v : make_range(var_group.n_variables()))
      {
        const auto var_num = var_group.number(v);
        auto & var_data = local_data.var_to_data[var_num];
        elem_condensed_dofs.clear();
        elem_uncondensed_dofs.clear();
        _dof_map.dof_indices(
            elem, elem_dofs, var_num, scalar_dofs_functor, field_dofs_functor, elem->p_level());
        var_data.uncondensed_dofs_offset = uncondensed_dof_size;
        var_data.condensed_dofs_offset = condensed_dof_size;
        var_data.num_uncondensed_dofs = elem_uncondensed_dofs.size();
        var_data.num_condensed_dofs = elem_condensed_dofs.size();
        uncondensed_dof_size += var_data.num_uncondensed_dofs;
        condensed_dof_size += var_data.num_condensed_dofs;
        local_data.reduced_space_indices.insert(local_data.reduced_space_indices.end(),
                                                elem_uncondensed_dofs.begin(),
                                                elem_uncondensed_dofs.end());
      }
    }

    local_data.Acc.resize(condensed_dof_size, condensed_dof_size);
    local_data.Acu.resize(condensed_dof_size, uncondensed_dof_size);
    local_data.Auc.resize(uncondensed_dof_size, condensed_dof_size);
    local_data.Auu.resize(uncondensed_dof_size, uncondensed_dof_size);
  }
  _local_uncondensed_dofs.assign(local_uncondensed_dofs.begin(), local_uncondensed_dofs.end());
  local_uncondensed_dofs.clear();

  //
  // Build the reduced system data
  //

  const dof_id_type n_local = _local_uncondensed_dofs.size();
  dof_id_type n = n_local;
  _dof_map.comm().sum(n);
  _reduced_sys_mat = SparseMatrix<Number>::build(_dof_map.comm());
  auto sp = _dof_map.build_sparsity(
      _mesh, /*calculate_constrained=*/false, /*uncondensed_dofs_only=*/true);
  const auto & nnz = sp->get_n_nz();
  const auto & noz = sp->get_n_oz();
  const auto nz = nnz.empty() ? dof_id_type(0) : *std::max_element(nnz.begin(), nnz.end());
  const auto oz = noz.empty() ? dof_id_type(0) : *std::max_element(noz.begin(), noz.end());
  _reduced_sys_mat->init(n, n, n_local, n_local, nz, oz);
  _reduced_solver = LinearSolver<Number>::build(_dof_map.comm());
  _reduced_solver->init("condensed_");
  _reduced_rhs = NumericVector<Number>::build(_dof_map.comm());

  // Build a map from the full size problem uncondensed dof indices to the reduced problem
  // (uncondensed) dof indices
  std::unordered_map<dof_id_type, dof_id_type> full_dof_to_reduced_dof;
  const auto local_start = _reduced_sys_mat->row_start();
  for (const auto i : index_range(_local_uncondensed_dofs))
    full_dof_to_reduced_dof[_local_uncondensed_dofs[i]] = i + local_start;

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
                                                   std::vector<dof_id_type> & reduced_dof_ids)
  {
    reduced_dof_ids.resize(full_dof_ids.size());
    for (const auto i : index_range(full_dof_ids))
      reduced_dof_ids[i] = libmesh_map_find(full_dof_to_reduced_dof, full_dof_ids[i]);
  };

  auto action_functor = [&full_dof_to_reduced_dof](processor_id_type,
                                                   const std::vector<dof_id_type> & full_dof_ids,
                                                   const std::vector<dof_id_type> & reduced_dof_ids)
  {
    for (const auto i : index_range(full_dof_ids))
    {
      libmesh_assert(!full_dof_to_reduced_dof.count(full_dof_ids[i]));
      full_dof_to_reduced_dof[full_dof_ids[i]] = reduced_dof_ids[i];
    }
  };

  TIMPI::pull_parallel_vector_data(_dof_map.comm(),
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

  _is_initialized = true;
}

void
StaticCondensation::add_matrix(const Elem & elem,
                               const unsigned int i_var,
                               const unsigned int j_var,
                               const DenseMatrix<Number> & k_var_ij)
{
  auto & local_data = libmesh_map_find(_elem_to_local_data, elem.id());
  auto & ivar_data = libmesh_map_find(local_data.var_to_data, i_var);
  auto & jvar_data = libmesh_map_find(local_data.var_to_data, j_var);
  auto copy_mat = [&k_var_ij](EigenMatrix & k_sc,
                              const unsigned int var_i_num_split_dofs,
                              const unsigned int k_sc_i_offset,
                              const unsigned int k_var_ij_i_offset,
                              const unsigned int var_j_num_split_dofs,
                              const unsigned int k_sc_j_offset,
                              const unsigned int k_var_ij_j_offset)
  {
    for (const auto i : make_range(var_i_num_split_dofs))
      for (const auto j : make_range(var_j_num_split_dofs))
        k_sc(k_sc_i_offset + i, k_sc_j_offset + j) =
            k_var_ij(k_var_ij_i_offset + i, k_var_ij_j_offset + j);
  };

  copy_mat(local_data.Auu,
           ivar_data.num_uncondensed_dofs,
           ivar_data.uncondensed_dofs_offset,
           0,
           jvar_data.num_uncondensed_dofs,
           jvar_data.uncondensed_dofs_offset,
           0);
  copy_mat(local_data.Auc,
           ivar_data.num_uncondensed_dofs,
           ivar_data.uncondensed_dofs_offset,
           0,
           jvar_data.num_condensed_dofs,
           jvar_data.condensed_dofs_offset,
           jvar_data.num_uncondensed_dofs);
  copy_mat(local_data.Acu,
           ivar_data.num_condensed_dofs,
           ivar_data.condensed_dofs_offset,
           ivar_data.num_uncondensed_dofs,
           jvar_data.num_uncondensed_dofs,
           jvar_data.uncondensed_dofs_offset,
           0);
  copy_mat(local_data.Acc,
           ivar_data.num_condensed_dofs,
           ivar_data.condensed_dofs_offset,
           ivar_data.num_uncondensed_dofs,
           jvar_data.num_condensed_dofs,
           jvar_data.condensed_dofs_offset,
           jvar_data.num_uncondensed_dofs);
}

void
StaticCondensation::setup()
{
  _reduced_sys_mat->zero();

  const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision, 0, ", ", "\n");
  auto original_flags = libMesh::out.flags();
  libMesh::out << std::fixed << std::setprecision(2);

  for (auto & [elem_id, local_data] : _elem_to_local_data)
  {
    libmesh_ignore(elem_id);
    libMesh::out << "Matrix for elem ID " << elem_id << ":\n"
                 << local_data.Acc.format(CSVFormat) << std::endl;
    libMesh::out.flags(original_flags);
    libMesh::out << "For elem id " << elem_id << " the A determinant is "
                 << local_data.Acc.determinant() << std::endl;
    local_data.AccFactor = local_data.Acc.partialPivLu();
    const EigenMatrix S =
        local_data.Auu - local_data.Auc * local_data.AccFactor.solve(local_data.Acu);
    DenseMatrix<Number> shim(S.rows(), S.cols());
    for (const auto i : make_range(S.rows()))
      for (const auto j : make_range(S.cols()))
        shim(i, j) = S(i, j);
    _reduced_sys_mat->add_matrix(shim, local_data.reduced_space_indices);
  }

  _reduced_sys_mat->close();
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
  std::vector<dof_id_type> elem_dofs; // only used to satisfy API
  std::vector<dof_id_type> elem_condensed_dofs;
  std::vector<Number> elem_condensed_rhs_vec;
  EigenVector elem_condensed_rhs, elem_uncondensed_rhs;

  full_rhs.create_subvector(*_reduced_rhs, _local_uncondensed_dofs, /*all_global_entries=*/false);

  //
  // Forward elimination step
  //

  for (auto elem : _mesh.active_local_element_ptr_range())
  {
    elem_condensed_dofs.clear();
    auto & local_data = _elem_to_local_data[elem->id()];

    const auto sub_id = elem->subdomain_id();
    for (const auto vg : make_range(_dof_map.n_variable_groups()))
    {
      const auto & var_group = _dof_map.variable_group(vg);
      if (!var_group.active_on_subdomain(sub_id))
        continue;

      for (const auto v : make_range(var_group.n_variables()))
        _dof_map.dof_indices(elem,
                             elem_dofs,
                             var_group.number(v),
                             total_and_condensed_from_scalar_dofs_functor(elem_condensed_dofs),
                             total_and_condensed_from_field_dofs_functor(elem_condensed_dofs),
                             elem->p_level());
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
  std::vector<dof_id_type> elem_dofs; // only used to satisfy API
  std::vector<dof_id_type> elem_condensed_dofs, elem_uncondensed_dofs;
  std::vector<Number> elem_condensed_rhs_vec, elem_uncondensed_sol_vec;
  EigenVector elem_condensed_rhs, elem_uncondensed_sol, elem_condensed_sol;

  auto scalar_dofs_functor = [&elem_condensed_dofs, &elem_uncondensed_dofs](
                                 const Elem & elem,
                                 std::vector<dof_id_type> & dof_indices,
                                 const std::vector<dof_id_type> & scalar_dof_indices)
  {
    total_and_condensed_from_scalar_dofs_functor(elem_condensed_dofs)(
        elem, dof_indices, scalar_dof_indices);
    elem_uncondensed_dofs.insert(
        elem_uncondensed_dofs.end(), scalar_dof_indices.begin(), scalar_dof_indices.end());
  };

  auto field_dofs_functor =
      [this, &elem_condensed_dofs, &elem_uncondensed_dofs](const Elem & elem,
                                                           const unsigned int node_num,
                                                           const unsigned int var_num,
                                                           std::vector<dof_id_type> & dof_indices,
                                                           const dof_id_type field_dof)
  {
    total_and_condensed_from_field_dofs_functor(elem_condensed_dofs)(
        elem, node_num, var_num, dof_indices, field_dof);

    if (node_num != invalid_uint && !elem.is_internal(node_num))
      elem_uncondensed_dofs.push_back(field_dof);
  };

  for (auto elem : _mesh.active_local_element_ptr_range())
  {
    elem_condensed_dofs.clear();
    elem_uncondensed_dofs.clear();
    auto & local_data = _elem_to_local_data[elem->id()];

    const auto sub_id = elem->subdomain_id();
    for (const auto vg : make_range(_dof_map.n_variable_groups()))
    {
      const auto & var_group = _dof_map.variable_group(vg);
      if (!var_group.active_on_subdomain(sub_id))
        continue;

      for (const auto v : make_range(var_group.n_variables()))
        _dof_map.dof_indices(elem,
                             elem_dofs,
                             var_group.number(v),
                             scalar_dofs_functor,
                             field_dofs_functor,
                             elem->p_level());
    }

    set_local_vectors(full_rhs, elem_condensed_dofs, elem_condensed_rhs_vec, elem_condensed_rhs);
    set_local_vectors(
        full_sol, elem_uncondensed_dofs, elem_uncondensed_sol_vec, elem_uncondensed_sol);

    elem_condensed_sol =
        local_data.AccFactor.solve(elem_condensed_rhs - local_data.Acu * elem_uncondensed_sol);
    full_sol.insert(elem_condensed_sol.data(), elem_condensed_dofs);
  }

  full_sol.close();
}

void
StaticCondensation::apply(const NumericVector<Number> & full_rhs, NumericVector<Number> & full_sol)
{
  forward_elimination(full_rhs);
  _reduced_sol = full_sol.get_subvector(_local_uncondensed_dofs);
  _reduced_solver->solve(*_reduced_sys_mat, *_reduced_sol, *_reduced_rhs, 1e-5, 300);
  // Must restore to the full solution because during backwards substitution we will need to be able
  // to read ghosted dofs and we don't support ghosting of subvectors
  full_sol.restore_subvector(std::move(_reduced_sol), _local_uncondensed_dofs);
  backwards_substitution(full_rhs, full_sol);
}
}

#endif // LIBMESH_HAVE_PETSC
