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
#include "libmesh/petsc_matrix.h"
#include "libmesh/linear_solver.h"
#include "libmesh/implicit_system.h"
#include "timpi/parallel_sync.h"
#include <unordered_set>

namespace libMesh
{
StaticCondensation::StaticCondensation(const MeshBase & mesh,
                                       ImplicitSystem & system,
                                       const DofMap & dof_map)
  : Preconditioner<Number>(dof_map.comm()), _mesh(mesh), _system(system), _dof_map(dof_map)
{
}

StaticCondensation::~StaticCondensation() = default;

void
StaticCondensation::total_dofs_from_scalar_dofs(std::vector<dof_id_type> & dofs,
                                                const std::vector<dof_id_type> & scalar_dofs)
{
  dofs.insert(dofs.end(), scalar_dofs.begin(), scalar_dofs.end());
}

void
StaticCondensation::condensed_dofs_from_scalar_dofs(
    std::vector<dof_id_type> & /*condensed_dofs*/, const std::vector<dof_id_type> & /*scalar_dofs*/)
{
}

void
StaticCondensation::uncondensed_dofs_from_scalar_dofs(std::vector<dof_id_type> & uncondensed_dofs,
                                                      const std::vector<dof_id_type> & scalar_dofs)
{
  uncondensed_dofs.insert(uncondensed_dofs.end(), scalar_dofs.begin(), scalar_dofs.end());
}

void
StaticCondensation::total_dofs_from_field_dof(std::vector<dof_id_type> & dofs,
                                              const Elem & /*elem*/,
                                              const unsigned int /*node_num*/,
                                              const unsigned int /*var_num*/,
                                              const dof_id_type field_dof)
{
  dofs.push_back(field_dof);
}

void
StaticCondensation::condensed_dofs_from_field_dof(std::vector<dof_id_type> & condensed_dofs,
                                                  const Elem & elem,
                                                  const unsigned int node_num,
                                                  const unsigned int var_num,
                                                  const dof_id_type field_dof) const
{
  if (_uncondensed_vars.count(var_num))
    // If we're going to keep all the dofs for this var (e.g. we don't want to condense any out)
    // then we can return now
    return;

  if (node_num != invalid_uint)
  {
    // This is a nodal dof
    if (elem.is_internal(node_num))
      condensed_dofs.push_back(field_dof);
  }
  else
    condensed_dofs.push_back(field_dof);
}

void
StaticCondensation::uncondensed_dofs_from_field_dof(std::vector<dof_id_type> & uncondensed_dofs,
                                                    const Elem & elem,
                                                    const unsigned int node_num,
                                                    const unsigned int var_num,
                                                    const dof_id_type field_dof) const
{
  if (_uncondensed_vars.count(var_num))
  {
    libmesh_assert_msg(
        node_num == invalid_uint,
        "Users should not be providing continuous FEM variables to the uncondensed vars API");
    uncondensed_dofs.push_back(field_dof);
  }

  if (node_num != invalid_uint && !elem.is_internal(node_num))
    uncondensed_dofs.push_back(field_dof);
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

  const Elem * first_elem = nullptr;
  const auto first_elem_it = _mesh.active_local_elements_begin();
  if (first_elem_it != _mesh.active_local_elements_end())
    first_elem = *first_elem_it;

  auto scalar_dofs_functor =
      [this, &local_uncondensed_dofs, &elem_condensed_dofs, &elem_uncondensed_dofs, first_elem](
          const Elem & elem,
          std::vector<dof_id_type> & dof_indices,
          const std::vector<dof_id_type> & scalar_dof_indices)
  {
    total_dofs_from_scalar_dofs(dof_indices, scalar_dof_indices);
    condensed_dofs_from_scalar_dofs(elem_condensed_dofs, scalar_dof_indices);
    uncondensed_dofs_from_scalar_dofs(elem_uncondensed_dofs, scalar_dof_indices);

    // Only need to do this for the first element we encounter
    if (&elem == first_elem)
    {
      const processor_id_type last_pid = this->comm().size() - 1;
      if (this->comm().rank() == last_pid)
        local_uncondensed_dofs.insert(scalar_dof_indices.begin(), scalar_dof_indices.end());
    }
  };

  auto field_dofs_functor =
      [this, &local_uncondensed_dofs, &elem_condensed_dofs, &elem_uncondensed_dofs](
          const Elem & elem,
          const unsigned int node_num,
          const unsigned int var_num,
          std::vector<dof_id_type> & dof_indices,
          const dof_id_type field_dof)
  {
    total_dofs_from_field_dof(dof_indices, elem, node_num, var_num, field_dof);
    condensed_dofs_from_field_dof(elem_condensed_dofs, elem, node_num, var_num, field_dof);
    uncondensed_dofs_from_field_dof(elem_uncondensed_dofs, elem, node_num, var_num, field_dof);

    if (_uncondensed_vars.count(var_num) && elem.processor_id() == this->processor_id())
      local_uncondensed_dofs.insert(field_dof);

    if (node_num != invalid_uint && !elem.is_internal(node_num))
    {
      const auto & nd_ref = elem.node_ref(node_num);
      if (nd_ref.processor_id() == this->processor_id())
        local_uncondensed_dofs.insert(field_dof);
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
      }
    }

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

  _reduced_solver = LinearSolver<Number>::build(this->comm());
  _reduced_solver->init("condensed_");
  _reduced_rhs = NumericVector<Number>::build(this->comm());
  _reduced_sys_mat = SparseMatrix<Number>::build(this->comm());

  // Build ghosted full solution vector. Note that this is, in general, *not equal* to the system
  // solution, e.g. this may correspond to the solution for the Newton *update*
  _ghosted_full_sol = _system.current_local_solution->clone();
  // Need a RHS for storing the eliminated version that still has all the dof indices
  _eliminated_rhs = _system.rhs->clone();

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
  auto & full_mat = _system.get_system_matrix();
  if (full_mat.closed())
    full_mat.zero();

  std::vector<dof_id_type> elem_dofs, elem_uncondensed_dofs;

  auto scalar_dofs_functor =
      [this, &elem_uncondensed_dofs](const Elem & /*elem*/,
                                     std::vector<dof_id_type> & dof_indices,
                                     const std::vector<dof_id_type> & scalar_dof_indices)
  {
    total_dofs_from_scalar_dofs(dof_indices, scalar_dof_indices);
    uncondensed_dofs_from_scalar_dofs(elem_uncondensed_dofs, scalar_dof_indices);
  };

  auto field_dofs_functor = [this, &elem_uncondensed_dofs](const Elem & elem,
                                                           const unsigned int node_num,
                                                           const unsigned int var_num,
                                                           std::vector<dof_id_type> & dof_indices,
                                                           const dof_id_type field_dof)
  {
    total_dofs_from_field_dof(dof_indices, elem, node_num, var_num, field_dof);
    uncondensed_dofs_from_field_dof(elem_uncondensed_dofs, elem, node_num, var_num, field_dof);
  };

  DenseMatrix<Number> shim;
  for (auto elem : _mesh.active_local_element_ptr_range())
  {
    elem_uncondensed_dofs.clear();
    auto & local_data = libmesh_map_find(_elem_to_local_data, elem->id());

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

    local_data.AccFactor = local_data.Acc.partialPivLu();
    const EigenMatrix S =
        local_data.Auu - local_data.Auc * local_data.AccFactor.solve(local_data.Acu);
    shim.resize(S.rows(), S.cols());
    for (const auto i : make_range(S.rows()))
      for (const auto j : make_range(S.cols()))
        shim(i, j) = S(i, j);
    full_mat.add_matrix(shim, elem_uncondensed_dofs);
  }

  full_mat.close();
  full_mat.create_submatrix(*_reduced_sys_mat, _local_uncondensed_dofs, _local_uncondensed_dofs);
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
StaticCondensation::forward_elimination(NumericVector<Number> & rhs)
{
  std::vector<dof_id_type> elem_dofs; // only used to satisfy API
  std::vector<dof_id_type> elem_condensed_dofs, elem_uncondensed_dofs;
  std::vector<Number> elem_condensed_rhs_vec;
  EigenVector elem_condensed_rhs, elem_uncondensed_rhs;

  auto scalar_dofs_functor = [this, &elem_condensed_dofs, &elem_uncondensed_dofs](
                                 const Elem & /*elem*/,
                                 std::vector<dof_id_type> & dof_indices,
                                 const std::vector<dof_id_type> & scalar_dof_indices)
  {
    total_dofs_from_scalar_dofs(dof_indices, scalar_dof_indices);
    condensed_dofs_from_scalar_dofs(elem_condensed_dofs, scalar_dof_indices);
    uncondensed_dofs_from_scalar_dofs(elem_uncondensed_dofs, scalar_dof_indices);
  };

  auto field_dofs_functor =
      [this, &elem_condensed_dofs, &elem_uncondensed_dofs](const Elem & elem,
                                                           const unsigned int node_num,
                                                           const unsigned int var_num,
                                                           std::vector<dof_id_type> & dof_indices,
                                                           const dof_id_type field_dof)
  {
    total_dofs_from_field_dof(dof_indices, elem, node_num, var_num, field_dof);
    condensed_dofs_from_field_dof(elem_condensed_dofs, elem, node_num, var_num, field_dof);
    uncondensed_dofs_from_field_dof(elem_uncondensed_dofs, elem, node_num, var_num, field_dof);
  };

  //
  // Forward elimination step
  //

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

    set_local_vectors(rhs, elem_condensed_dofs, elem_condensed_rhs_vec, elem_condensed_rhs);
    elem_uncondensed_rhs = -local_data.Auc * local_data.AccFactor.solve(elem_condensed_rhs);

    rhs.add_vector(elem_uncondensed_rhs.data(), elem_uncondensed_dofs);
  }
  rhs.close();
}

void
StaticCondensation::backwards_substitution(const NumericVector<Number> & rhs,
                                           NumericVector<Number> & sol)
{
  std::vector<dof_id_type> elem_dofs; // only used to satisfy API
  std::vector<dof_id_type> elem_condensed_dofs, elem_uncondensed_dofs;
  std::vector<Number> elem_condensed_rhs_vec, elem_uncondensed_sol_vec;
  EigenVector elem_condensed_rhs, elem_uncondensed_sol, elem_condensed_sol;

  auto scalar_dofs_functor = [&elem_condensed_dofs, &elem_uncondensed_dofs](
                                 const Elem & /*elem*/,
                                 std::vector<dof_id_type> & dof_indices,
                                 const std::vector<dof_id_type> & scalar_dof_indices)
  {
    total_dofs_from_scalar_dofs(dof_indices, scalar_dof_indices);
    condensed_dofs_from_scalar_dofs(elem_condensed_dofs, scalar_dof_indices);
    uncondensed_dofs_from_scalar_dofs(elem_uncondensed_dofs, scalar_dof_indices);
  };

  auto field_dofs_functor =
      [this, &elem_condensed_dofs, &elem_uncondensed_dofs](const Elem & elem,
                                                           const unsigned int node_num,
                                                           const unsigned int var_num,
                                                           std::vector<dof_id_type> & dof_indices,
                                                           const dof_id_type field_dof)
  {
    total_dofs_from_field_dof(dof_indices, elem, node_num, var_num, field_dof);
    condensed_dofs_from_field_dof(elem_condensed_dofs, elem, node_num, var_num, field_dof);
    uncondensed_dofs_from_field_dof(elem_uncondensed_dofs, elem, node_num, var_num, field_dof);
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

    set_local_vectors(rhs, elem_condensed_dofs, elem_condensed_rhs_vec, elem_condensed_rhs);
    set_local_vectors(
        *_ghosted_full_sol, elem_uncondensed_dofs, elem_uncondensed_sol_vec, elem_uncondensed_sol);

    elem_condensed_sol =
        local_data.AccFactor.solve(elem_condensed_rhs - local_data.Acu * elem_uncondensed_sol);
    sol.insert(elem_condensed_sol.data(), elem_condensed_dofs);
  }

  sol.close();
}

void
StaticCondensation::apply(const NumericVector<Number> & full_rhs,
                          NumericVector<Number> & full_parallel_sol)
{
  *_eliminated_rhs = full_rhs;
  forward_elimination(*_eliminated_rhs);
  // Apparently PETSc will send us the yvec without zeroing it ahead of time. This can be a poor
  // initial guess for the Krylov solve as well as lead to bewildered users who expect their initial
  // residual norm to equal the norm of the RHS
  full_parallel_sol.zero();
  _reduced_sol = full_parallel_sol.get_subvector(_local_uncondensed_dofs);
  _reduced_rhs = _eliminated_rhs->get_subvector(_local_uncondensed_dofs);
  _reduced_solver->solve(*_reduced_sys_mat, *_reduced_sol, *_reduced_rhs, 1e-5, 300);
  // Must restore to the full solution because during backwards substitution we will need to be able
  // to read ghosted dofs and we don't support ghosting of subvectors
  full_parallel_sol.restore_subvector(std::move(_reduced_sol), _local_uncondensed_dofs);
  _eliminated_rhs->restore_subvector(std::move(_reduced_rhs), _local_uncondensed_dofs);
  *_ghosted_full_sol = full_parallel_sol;
  backwards_substitution(full_rhs, full_parallel_sol);
}
}

#endif // LIBMESH_HAVE_PETSC
