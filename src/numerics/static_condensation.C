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
  : _mesh(mesh), _dof_map(dof_map)
{
  init();
}

auto
StaticCondensation::computeElemDofsScalar(std::vector<dof_id_type> & /*elem_interior_dofs*/,
                                          std::vector<dof_id_type> & elem_trace_dofs)
    -> decltype(auto)
{
  return [&elem_trace_dofs](const Elem &,
                            std::vector<dof_id_type> & elem_dof_indices,
                            const std::vector<dof_id_type> & scalar_dof_indices)
  {
    elem_dof_indices.insert(
        elem_trace_dofs.end(), scalar_dof_indices.begin(), scalar_dof_indices.end());
    elem_trace_dofs.insert(
        elem_trace_dofs.end(), scalar_dof_indices.begin(), scalar_dof_indices.end());
  };
}

auto
StaticCondensation::computeElemDofsField(std::vector<dof_id_type> & elem_interior_dofs,
                                         std::vector<dof_id_type> & elem_trace_dofs)
    -> decltype(auto)
{
  return [&elem_interior_dofs, &elem_trace_dofs](const Elem & elem,
                                                 const unsigned int node_num,
                                                 std::vector<dof_id_type> & elem_dof_indices,
                                                 const dof_id_type field_dof)
  {
    elem_dof_indices.push_back(field_dof);
    if (node_num != invalid_uint)
    {
      // This is a nodal dof
      if (elem.is_internal(node_num))
        elem_interior_dofs.push_back(field_dof);
      else
        elem_trace_dofs.push_back(field_dof);
    }
    else
      elem_interior_dofs.push_back(field_dof);
  };
}

void
StaticCondensation::init()
{
  std::vector<dof_id_type> elem_dofs; // only used to satisfy API
  std::vector<dof_id_type> elem_interior_dofs, elem_trace_dofs;
  std::unordered_set<dof_id_type> local_trace_dofs;
  std::unordered_map<processor_id_type, std::unordered_set<dof_id_type>> nonlocal_trace_dofs;

  const Elem * first_elem = nullptr;
  const auto first_elem_it = _mesh.active_local_elements_begin();
  if (first_elem_it != _mesh.active_local_elements_end())
    first_elem = *first_elem_it;

  auto scalar_dofs_functor = [this,
                              &local_trace_dofs,
                              &nonlocal_trace_dofs,
                              &elem_interior_dofs,
                              &elem_trace_dofs,
                              first_elem](const Elem & elem,
                                          std::vector<dof_id_type> & dof_indices,
                                          const std::vector<dof_id_type> & scalar_dof_indices)
  {
    computeElemDofsScalar(elem_interior_dofs,
                          elem_trace_dofs)(elem, dof_indices, scalar_dof_indices);

    // Only need to do this for the first element we encounter
    if (&elem == first_elem)
    {
      const processor_id_type last_pid = _dof_map.comm().size() - 1;
      if (_dof_map.comm().rank() == last_pid)
        local_trace_dofs.insert(scalar_dof_indices.begin(), scalar_dof_indices.end());
      else
        nonlocal_trace_dofs[last_pid].insert(scalar_dof_indices.begin(), scalar_dof_indices.end());
    }
  };

  auto field_dofs_functor =
      [this, &local_trace_dofs, &nonlocal_trace_dofs, &elem_interior_dofs, &elem_trace_dofs](
          const Elem & elem,
          const unsigned int node_num,
          std::vector<dof_id_type> & dof_indices,
          const dof_id_type field_dof)
  {
    computeElemDofsField(elem_interior_dofs,
                         elem_trace_dofs)(elem, node_num, dof_indices, field_dof);

    if (node_num != invalid_uint && !elem.is_internal(node_num))
    {
      const auto & nd_ref = elem.node_ref(node_num);
      if (nd_ref.processor_id() == _dof_map.processor_id())
        local_trace_dofs.insert(field_dof);
      else
        nonlocal_trace_dofs[nd_ref.processor_id()].insert(field_dof);
    }
  };

  for (auto elem : _mesh.active_local_element_ptr_range())
  {
    elem_interior_dofs.clear();
    elem_trace_dofs.clear();
    unsigned int boundary_dof_size = 0;
    unsigned int interior_dof_size = 0;
    auto & local_data = _elem_to_local_data[elem->id()];

    const auto sub_id = elem->subdomain_id();
    for (const auto vg : make_range(_dof_map.n_variable_groups()))
    {
      const auto & var_group = _dof_map.variable_group(vg);
      if (!var_group.active_on_subdomain(sub_id))
        continue;

      for (const auto v : make_range(var_group.n_variables()))
      {
        auto & var_data = local_data.var_to_data[v];
        _dof_map.dof_indices(
            elem, elem_dofs, v, scalar_dofs_functor, field_dofs_functor, elem->p_level());
        var_data.boundary_dofs_offset = boundary_dof_size;
        var_data.interior_dofs_offset = interior_dof_size;
        var_data.num_boundary_dofs = elem_trace_dofs.size();
        var_data.num_interior_dofs = elem_interior_dofs.size();
        boundary_dof_size += var_data.num_boundary_dofs;
        interior_dof_size += var_data.num_interior_dofs;
        local_data.reduced_space_indices = elem_trace_dofs;
      }
    }

    local_data.Aii.resize(interior_dof_size, interior_dof_size);
    local_data.Aib.resize(interior_dof_size, boundary_dof_size);
    local_data.Abi.resize(boundary_dof_size, interior_dof_size);
    local_data.Abb.resize(boundary_dof_size, boundary_dof_size);
  }

  //
  // Build the reduced soln and rhs vectors and system mat
  //

  _reduced_sol = NumericVector<Number>::build(_dof_map.comm());
  const dof_id_type n_local = local_trace_dofs.size();
  dof_id_type n = n_local;
  _dof_map.comm().sum(n);
  _reduced_sol->init(n, n_local, /*fast=*/false, PARALLEL);
  _reduced_rhs = _reduced_sol->clone();
  _reduced_sys_mat = SparseMatrix<Number>::build(_dof_map.comm());
  // // Make some assumptions about the sparsity pattern
  // dof_id_type num_couplings = 0;
  // if (_dof.comm().rank() == 0 && !_elem_to_local_data.empty())
  //   num_couplings = _elem_to_local_data.begin()->second.Abb.n();
  // _dof.comm().broadcast(num_couplings);
  // // We multiply by 2 for process on-diagonal since the facet dofs are on faces between 2 elements
  // // But actually for general static condensation we could be talking about continuous Galerkin.
  // // For hex elements this would mean an 8x multiplier! So for now we just accept the libmesh
  // // defaults and prepare ourselves for possible nonzero allocations to get started
  _reduced_sys_mat->init(n, n, n_local, n_local /*, 2 * num_couplings, num_couplings*/);
  _reduced_solver = LinearSolver<Number>::build(_dof_map.comm());
  _reduced_solver->init("condensed_");

  // Build a map from the full size problem trace dof indices to the reduced problem (trace) dof
  // indices
  std::unordered_map<dof_id_type, dof_id_type> full_dof_to_reduced_dof;
  auto set_it = local_trace_dofs.begin();
  for (const auto i :
       make_range(_reduced_sol->first_local_index(), _reduced_sol->last_local_index()))
  {
    full_dof_to_reduced_dof[*set_it] = i;
    ++set_it;
  }

  //
  // Now we need to pull our nonlocal data
  //

  // build our queries
  std::unordered_map<processor_id_type, std::vector<dof_id_type>> nonlocal_trace_dofs_mapvec;
  for (const auto & [pid, set] : nonlocal_trace_dofs)
  {
    auto & vec = nonlocal_trace_dofs_mapvec[pid];
    vec.assign(set.begin(), set.end());
  }
  // clear no longer needed memory
  nonlocal_trace_dofs.clear();

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
                                   nonlocal_trace_dofs_mapvec,
                                   gather_functor,
                                   action_functor,
                                   &DofObject::invalid_id);
  nonlocal_trace_dofs_mapvec.clear();

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

  copy_mat(local_data.Abb,
           ivar_data.num_boundary_dofs,
           ivar_data.boundary_dofs_offset,
           0,
           jvar_data.num_boundary_dofs,
           jvar_data.boundary_dofs_offset,
           0);
  copy_mat(local_data.Abi,
           ivar_data.num_boundary_dofs,
           ivar_data.boundary_dofs_offset,
           0,
           jvar_data.num_interior_dofs,
           jvar_data.interior_dofs_offset,
           jvar_data.num_boundary_dofs);
  copy_mat(local_data.Aib,
           ivar_data.num_interior_dofs,
           ivar_data.interior_dofs_offset,
           ivar_data.num_boundary_dofs,
           jvar_data.num_boundary_dofs,
           jvar_data.boundary_dofs_offset,
           0);
  copy_mat(local_data.Aii,
           ivar_data.num_interior_dofs,
           ivar_data.interior_dofs_offset,
           ivar_data.num_boundary_dofs,
           jvar_data.num_interior_dofs,
           jvar_data.interior_dofs_offset,
           jvar_data.num_boundary_dofs);
}

void
StaticCondensation::assemble_reduced_mat()
{
  for (auto & [elem_id, local_data] : _elem_to_local_data)
  {
    libmesh_ignore(elem_id);
    local_data.AiiFactor = local_data.Aii.partialPivLu();
    const auto S = local_data.Abb - local_data.Abi * local_data.AiiFactor.solve(local_data.Aib);
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
  std::vector<dof_id_type> elem_interior_dofs, elem_trace_dofs_full;
  std::vector<Number> elem_interior_rhs_vec, elem_trace_rhs_vec;
  EigenVector elem_interior_rhs, elem_trace_rhs;

  _reduced_rhs->zero();

  //
  // Forward elimination step
  //

  for (auto elem : _mesh.active_local_element_ptr_range())
  {
    elem_interior_dofs.clear();
    elem_trace_dofs_full.clear();
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
                             v,
                             computeElemDofsScalar(elem_interior_dofs, elem_trace_dofs_full),
                             computeElemDofsField(elem_interior_dofs, elem_trace_dofs_full),
                             elem->p_level());
    }

    set_local_vectors(full_rhs, elem_interior_dofs, elem_interior_rhs_vec, elem_interior_rhs);
    set_local_vectors(full_rhs, elem_trace_dofs_full, elem_trace_rhs_vec, elem_trace_rhs);

    elem_trace_rhs -= local_data.Abi * local_data.AiiFactor.solve(elem_interior_rhs);

    libmesh_assert(cast_int<std::size_t>(elem_trace_rhs.size()) ==
                   local_data.reduced_space_indices.size());
    _reduced_rhs->add_vector(elem_trace_rhs.data(), local_data.reduced_space_indices);
  }
  _reduced_rhs->close();
}

void
StaticCondensation::backwards_substitution(const NumericVector<Number> & full_rhs,
                                           NumericVector<Number> & full_sol)
{
  std::vector<dof_id_type> elem_dofs; // only used to satisfy API
  std::vector<dof_id_type> elem_interior_dofs, elem_trace_dofs_full;
  std::vector<Number> elem_interior_rhs_vec, elem_trace_sol_vec;
  EigenVector elem_interior_rhs, elem_trace_sol, elem_interior_sol;

  for (auto elem : _mesh.active_local_element_ptr_range())
  {
    elem_interior_dofs.clear();
    elem_trace_dofs_full.clear();
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
                             v,
                             computeElemDofsScalar(elem_interior_dofs, elem_trace_dofs_full),
                             computeElemDofsField(elem_interior_dofs, elem_trace_dofs_full),
                             elem->p_level());
    }

    set_local_vectors(full_rhs, elem_interior_dofs, elem_interior_rhs_vec, elem_interior_rhs);
    set_local_vectors(
        *_reduced_sol, local_data.reduced_space_indices, elem_trace_sol_vec, elem_trace_sol);

    elem_interior_sol =
        local_data.AiiFactor.solve(elem_interior_rhs - local_data.Aib * elem_trace_sol);
    full_sol.insert(elem_interior_sol.data(), elem_interior_dofs);
    // We are doing this redundantly since multiple elements share trace dofs. An alternative would
    // be to hold a global reduced to full index map
    full_sol.insert(elem_trace_sol.data(), elem_trace_dofs_full);
  }

  full_sol.close();
}

void
StaticCondensation::solve(const NumericVector<Number> & full_rhs, NumericVector<Number> & full_sol)
{
  assemble_reduced_mat();
  forward_elimination(full_rhs);
  _reduced_sol->zero();
  _reduced_solver->solve(*_reduced_sys_mat, *_reduced_sol, *_reduced_rhs, 1e-5, 300);
  backwards_substitution(full_rhs, full_sol);
}
}
