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

#include "libmesh/reduced_dof_map.h"

#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/equation_systems.h"
#include "libmesh/int_range.h"
#include "libmesh/mesh_base.h"
#include "libmesh/system.h"
#include "timpi/parallel_sync.h"

#include <algorithm>
#include <set>

namespace libMesh
{
ReducedDofMap::ReducedDofMap(MeshBase & mesh,
                             System & system,
                             const DofMap & dof_map,
                             std::string reduced_system_name)
  : DofMapBase(dof_map.comm()),
    _mesh(mesh),
    _system(system),
    _dof_map(dof_map),
    _reduced_system_name(std::move(reduced_system_name)),
    _initialized(false),
    _reduced_system(nullptr)
{
}

ReducedDofMap::~ReducedDofMap() = default;

void ReducedDofMap::retained_dof_added(const Elem &, dof_id_type) {}

void ReducedDofMap::unretained_dof_found(const Elem &, dof_id_type) {}

void ReducedDofMap::retain_dof(const dof_id_type full_dof)
{
  if (_dof_map.local_index(full_dof))
    _local_retained_set.insert(full_dof);
  else
    _nonlocal_retained[_dof_map.dof_owner(full_dof)].insert(full_dof);
}

void ReducedDofMap::retain_dof_and_constraining_dofs(const Elem * const elem,
                                                     const dof_id_type full_dof,
                                                     std::unordered_set<dof_id_type> & elem_seen,
                                                     std::vector<dof_id_type> & elem_var_dofs)
{
  if (!elem_seen.insert(full_dof).second)
    return;

  this->retain_dof(full_dof);
  elem_var_dofs.push_back(full_dof);
  this->retained_dof_added(*elem, full_dof);

#ifdef LIBMESH_ENABLE_CONSTRAINTS
  const auto & constraints = _dof_map.get_dof_constraints();
  if (const auto it = constraints.find(full_dof); it != constraints.end())
    for (const auto & [constraining_dof, weight] : it->second)
      {
        libmesh_ignore(weight);
        // Constraining dofs may themselves be constrained
        this->retain_dof_and_constraining_dofs(elem, constraining_dof, elem_seen, elem_var_dofs);
      }
#endif
}

void ReducedDofMap::reinit()
{
  if (this->initialized())
    this->clear();

  // SCALAR dofs are always retained. They are owned by the last rank
  for (const auto vg : make_range(_dof_map.n_variable_groups()))
    if (const auto & vg_description = _dof_map.variable_group(vg);
        vg_description.type().family == SCALAR)
      {
        std::vector<dof_id_type> scalar_dof_indices;
        for (const auto vg_vn : make_range(vg_description.n_variables()))
          {
            _dof_map.SCALAR_dof_indices(scalar_dof_indices, vg_description.number(vg_vn));
            for (const auto dof : scalar_dof_indices)
              this->retain_dof(dof);
          }
      }

  std::vector<dof_id_type> elem_dofs; // only used to satisfy the dof_indices API
  std::vector<dof_id_type> * elem_var_dofs = nullptr;
  std::unordered_set<dof_id_type> elem_seen;
  std::set<unsigned int> full_vars_present;

  auto scalar_dofs_functor = [this, &elem_seen, &elem_var_dofs](
                                 const Elem & elem,
                                 std::vector<dof_id_type> & dof_indices,
                                 const std::vector<dof_id_type> & scalar_dof_indices) {
    dof_indices.insert(dof_indices.end(), scalar_dof_indices.begin(), scalar_dof_indices.end());
    for (const auto dof : scalar_dof_indices)
      this->retain_dof_and_constraining_dofs(&elem, dof, elem_seen, *elem_var_dofs);
  };

  auto field_dofs_functor = [this, &elem_seen, &elem_var_dofs](const Elem & elem,
                                                              const unsigned int node_num,
                                                              const unsigned int var_num,
                                                              std::vector<dof_id_type> & dof_indices,
                                                              const dof_id_type field_dof) {
    dof_indices.push_back(field_dof);
    if (this->retain_field_dof(elem, node_num, var_num))
      this->retain_dof_and_constraining_dofs(&elem, field_dof, elem_seen, *elem_var_dofs);
    else
      this->unretained_dof_found(elem, field_dof);
  };

  for (const auto * const elem : _mesh.active_local_element_ptr_range())
    {
      auto & elem_dofs_by_var = _elem_reduced_dofs[elem->id()];
      elem_dofs_by_var.resize(_dof_map.n_variables());
      elem_seen.clear();

      const auto sub_id = elem->subdomain_id();
      for (const auto vg : make_range(_dof_map.n_variable_groups()))
        {
          const auto & var_group = _dof_map.variable_group(vg);
          if (!var_group.active_on_subdomain(sub_id))
            continue;
          for (const auto v : make_range(var_group.n_variables()))
            {
              const auto var_num = var_group.number(v);
              elem_var_dofs = &elem_dofs_by_var[var_num];
              _dof_map.dof_indices(
                  elem, elem_dofs, var_num, scalar_dofs_functor, field_dofs_functor, elem->p_level());
              if (!elem_var_dofs->empty())
                full_vars_present.insert(var_num);
            }
        }
    }

  //
  // A rank may own a retained dof that none of its local elements use, for instance a mid-edge
  // node of a coarse element that is a vertex of a finer element owned elsewhere. Every rank that
  // retains a nonlocal dof therefore tells the dof's owner
  //
  std::unordered_map<processor_id_type, std::vector<dof_id_type>> nonlocal_retained;
  for (const auto & [pid, set] : _nonlocal_retained)
    nonlocal_retained[pid].assign(set.begin(), set.end());
  _nonlocal_retained.clear();

  TIMPI::push_parallel_vector_data(
      this->comm(),
      nonlocal_retained,
      [this](processor_id_type, const std::vector<dof_id_type> & dofs) {
        _local_retained_set.insert(dofs.begin(), dofs.end());
      });

  _local_retained_dofs.assign(_local_retained_set.begin(), _local_retained_set.end());
  _local_retained_set.clear();
  std::sort(_local_retained_dofs.begin(), _local_retained_dofs.end());

  //
  // Number the reduced space
  //
  this->compute_dof_info(_local_retained_dofs.size());

  const auto full_first = _dof_map.first_dof();
  _local_full_to_reduced.assign(_dof_map.n_local_dofs(), DofObject::invalid_id);
  for (const auto i : index_range(_local_retained_dofs))
    _local_full_to_reduced[_local_retained_dofs[i] - full_first] = this->first_dof() + i;

  auto gather_reduced = [this](processor_id_type,
                               const std::vector<dof_id_type> & full_dofs,
                               std::vector<dof_id_type> & reduced_dofs) {
    reduced_dofs.resize(full_dofs.size());
    for (const auto i : index_range(full_dofs))
      reduced_dofs[i] = this->reduced_dof(full_dofs[i]);
  };
  auto store_reduced = [this](processor_id_type,
                              const std::vector<dof_id_type> & full_dofs,
                              const std::vector<dof_id_type> & reduced_dofs) {
    for (const auto i : index_range(full_dofs))
      _nonlocal_full_to_reduced[full_dofs[i]] = reduced_dofs[i];
  };
  TIMPI::pull_parallel_vector_data(this->comm(),
                                   nonlocal_retained,
                                   gather_reduced,
                                   store_reduced,
                                   &DofObject::invalid_id);

  //
  // Determine the variables with any dofs in the reduced space
  //
  this->comm().set_union(full_vars_present);
  _reduced_vars.reserve(full_vars_present.size());
  unsigned int first_local_number = 0;
  for (const auto full_var_num : full_vars_present)
    {
      const auto & full_var = _dof_map.variable(full_var_num);
      _reduced_vars.push_back(Variable{nullptr,
                                       full_var.name(),
                                       cast_int<unsigned int>(_reduced_vars.size()),
                                       first_local_number,
                                       full_var.type()});
      first_local_number += _reduced_vars.back().n_components(_mesh);
    }

  // Renumber the element dofs, keeping only the variables in the reduced space
  for (auto & [elem_id, elem_dofs_by_var] : _elem_reduced_dofs)
    {
      libmesh_ignore(elem_id);
      std::vector<std::vector<dof_id_type>> reduced_dofs_by_var;
      reduced_dofs_by_var.reserve(full_vars_present.size());
      for (const auto full_var_num : full_vars_present)
        {
          auto & var_dofs = reduced_dofs_by_var.emplace_back(
              std::move(elem_dofs_by_var[full_var_num]));
          for (auto & dof : var_dofs)
            dof = this->reduced_dof(dof);
        }
      elem_dofs_by_var = std::move(reduced_dofs_by_var);
    }

  // The sparsity pattern of the reduced space, with rows indexed by the local full dofs
  _reduced_sp = _dof_map.build_sparsity(_mesh, /*calculate_constrained=*/false, this);
  const auto & nnz = _reduced_sp->get_n_nz();
  const auto & noz = _reduced_sp->get_n_oz();
  libmesh_assert_equal_to(nnz.size(), noz.size());

  // Exact per-row preallocation counts. These are critical when SCALAR dofs introduce dense rows,
  // since the maximum row count would then allocate a dense matrix
  _reduced_nnz.resize(_local_retained_dofs.size());
  _reduced_noz.resize(_local_retained_dofs.size());
  for (const auto i : index_range(_local_retained_dofs))
    {
      const dof_id_type local_full_i = _local_retained_dofs[i] - full_first;
      libmesh_assert_less(local_full_i, nnz.size());
      _reduced_nnz[i] = nnz[local_full_i];
      _reduced_noz[i] = noz[local_full_i];
    }

  //
  // The stand-in system gives the reduced variables a system number. Nodes hold no reduced dofs,
  // since reduced dofs are queried by element only
  //
  const std::vector<unsigned int> nvpg(_reduced_vars.size(), 1);

  // add_system returns the existing system if there is one, so this is safe on every reinit
  _reduced_system = &_system.get_equation_systems().add_system("Basic", _reduced_system_name);
  if (!_reduced_system->is_initialized())
    {
      _reduced_system->init();
      for (auto * const nd : _mesh.active_node_ptr_range())
        {
          nd->set_n_vars_per_group(_reduced_system->number(), nvpg);
          for (const auto g : index_range(nvpg))
            nd->set_n_comp_group(_reduced_system->number(), g, 0);
        }
    }

  // The stand-in system holds no solution data worth writing
  _reduced_system->hide_output() = true;

  _initialized = true;
}

dof_id_type ReducedDofMap::reduced_dof(const dof_id_type full_dof) const
{
  if (_dof_map.local_index(full_dof))
    {
      const auto reduced = _local_full_to_reduced[full_dof - _dof_map.first_dof()];
      libmesh_assert_msg(reduced != DofObject::invalid_id,
                         "Full dof " << full_dof << " is not retained in the reduced space");
      return reduced;
    }
  return libmesh_map_find(_nonlocal_full_to_reduced, full_dof);
}

const std::vector<std::vector<dof_id_type>> &
ReducedDofMap::elem_dof_indices(const dof_id_type elem_id) const
{
  return libmesh_map_find(_elem_reduced_dofs, elem_id);
}

unsigned int ReducedDofMap::n_variables() const { return _reduced_vars.size(); }

const Variable & ReducedDofMap::variable(const unsigned int c) const
{
  libmesh_assert_less(c, _reduced_vars.size());
  return _reduced_vars[c];
}

void ReducedDofMap::dof_indices(const Elem * const elem,
                                std::vector<dof_id_type> & di,
                                const unsigned int vn,
                                int /*p_level*/) const
{
  di = libmesh_map_find(_elem_reduced_dofs, elem->id())[vn];
}

void ReducedDofMap::dof_indices(const Node * const, std::vector<dof_id_type> &, const unsigned int) const
{
  libmesh_error_msg("ReducedDofMap dof indices are only meant to be queried with elements, not "
                    "nodes");
}

void ReducedDofMap::clear()
{
  DofMapBase::clear();
  _local_retained_dofs.clear();
  _local_full_to_reduced.clear();
  _nonlocal_full_to_reduced.clear();
  _local_retained_set.clear();
  _nonlocal_retained.clear();
  _elem_reduced_dofs.clear();
  _reduced_vars.clear();
  _reduced_sp = nullptr;
  _reduced_nnz.clear();
  _reduced_noz.clear();
  _initialized = false;
}
}
