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

#include "libmesh/static_condensation_dof_map.h"

#include "libmesh/mesh_base.h"
#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/int_range.h"
#include "libmesh/system.h"
#include "libmesh/equation_systems.h"
#include "timpi/parallel_sync.h"
#include <unordered_set>

namespace libMesh
{
StaticCondensationDofMap::StaticCondensationDofMap(MeshBase & mesh,
                                                   System & system,
                                                   const DofMap & dof_map)
  : DofMapBase(dof_map.comm()),
    _mesh(mesh),
    _system(system),
    _dof_map(dof_map),
    _sc_is_initialized(false)
{
  libmesh_experimental();
}

StaticCondensationDofMap::~StaticCondensationDofMap() = default;

void StaticCondensationDofMap::add_uncondensed_dof(
    const dof_id_type full_dof_number,
    const bool involved_in_constraints,
    std::unordered_map<dof_id_type, dof_id_type> & uncondensed_global_to_local_map,
    std::unordered_set<dof_id_type> & local_uncondensed_dofs_set,
    std::unordered_map<processor_id_type, std::unordered_set<dof_id_type>> &
        nonlocal_uncondensed_dofs,
    std::vector<dof_id_type> & elem_uncondensed_dofs,
    dof_id_type & uncondensed_local_dof_number,
    std::unordered_set<dof_id_type> & constraint_dofs)
{
  if (uncondensed_global_to_local_map.count(full_dof_number))
    // We've already seen this dof on this element
    return;

  if (_dof_map.local_index(full_dof_number))
    local_uncondensed_dofs_set.insert(full_dof_number);
  else
    nonlocal_uncondensed_dofs[_dof_map.dof_owner(full_dof_number)].insert(full_dof_number);

  elem_uncondensed_dofs.push_back(full_dof_number);
  (uncondensed_global_to_local_map)[full_dof_number] = uncondensed_local_dof_number++;
  if (involved_in_constraints)
    constraint_dofs.insert(full_dof_number);
}

void StaticCondensationDofMap::add_uncondensed_dof_plus_constraint_dofs(
    const dof_id_type full_dof_number,
    bool involved_in_constraints,
    std::unordered_map<dof_id_type, dof_id_type> & uncondensed_global_to_local_map,
    std::unordered_set<dof_id_type> & local_uncondensed_dofs_set,
    std::unordered_map<processor_id_type, std::unordered_set<dof_id_type>> &
        nonlocal_uncondensed_dofs,
    std::vector<dof_id_type> & elem_uncondensed_dofs,
    dof_id_type & uncondensed_local_dof_number,
    std::unordered_set<dof_id_type> & constraint_dofs)
{
  const auto & full_dof_constraints = _dof_map.get_dof_constraints();
  auto it = full_dof_constraints.find(full_dof_number);
  const bool is_constrained = it != full_dof_constraints.end();
  involved_in_constraints = involved_in_constraints || is_constrained;

  this->add_uncondensed_dof(full_dof_number,
                            involved_in_constraints,
                            uncondensed_global_to_local_map,
                            local_uncondensed_dofs_set,
                            nonlocal_uncondensed_dofs,
                            elem_uncondensed_dofs,
                            uncondensed_local_dof_number,
                            constraint_dofs);
  if (is_constrained)
    for (const auto [full_constraining_dof, weight] : it->second)
      {
        libmesh_ignore(weight);
        // Our constraining dofs may themselves be constrained
        this->add_uncondensed_dof_plus_constraint_dofs(full_constraining_dof,
                                                       /*involved_in_constraints=*/true,
                                                       uncondensed_global_to_local_map,
                                                       local_uncondensed_dofs_set,
                                                       nonlocal_uncondensed_dofs,
                                                       elem_uncondensed_dofs,
                                                       uncondensed_local_dof_number,
                                                       constraint_dofs);
      }
}

void StaticCondensationDofMap::reinit()
{
  if (this->initialized())
    this->clear();

  std::vector<dof_id_type> elem_dofs, elem_uncondensed_dofs; // only used to satisfy API
  dof_id_type condensed_local_dof_number = 0, uncondensed_local_dof_number = 0;
  std::unordered_map<dof_id_type, dof_id_type> *condensed_global_to_local_map = nullptr,
                                               *uncondensed_global_to_local_map = nullptr;
  std::set<unsigned int> full_vars_present_in_reduced_sys;
  std::unordered_set<dof_id_type> local_uncondensed_dofs_set, constraint_dofs;
  std::unordered_map<processor_id_type, std::unordered_set<dof_id_type>> nonlocal_uncondensed_dofs;

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
              local_uncondensed_dofs_set.insert(scalar_dof_indices.begin(),
                                                scalar_dof_indices.end());
            else
              nonlocal_uncondensed_dofs[last_pid].insert(scalar_dof_indices.begin(),
                                                         scalar_dof_indices.end());
          }
      }

  auto scalar_dofs_functor =
      [this,
       &uncondensed_global_to_local_map,
       &local_uncondensed_dofs_set,
       &nonlocal_uncondensed_dofs,
       &elem_uncondensed_dofs,
       &uncondensed_local_dof_number,
       &constraint_dofs](const Elem & /*elem*/,
                         std::vector<dof_id_type> & dof_indices,
                         const std::vector<dof_id_type> & scalar_dof_indices) {
        dof_indices.insert(dof_indices.end(), scalar_dof_indices.begin(), scalar_dof_indices.end());
        for (const auto global_dof : scalar_dof_indices)
          this->add_uncondensed_dof_plus_constraint_dofs(global_dof,
                                                         false,
                                                         *uncondensed_global_to_local_map,
                                                         local_uncondensed_dofs_set,
                                                         nonlocal_uncondensed_dofs,
                                                         elem_uncondensed_dofs,
                                                         uncondensed_local_dof_number,
                                                         constraint_dofs);
      };

  auto field_dofs_functor = [this,
                             &condensed_local_dof_number,
                             &condensed_global_to_local_map,
                             &uncondensed_global_to_local_map,
                             &local_uncondensed_dofs_set,
                             &nonlocal_uncondensed_dofs,
                             &elem_uncondensed_dofs,
                             &uncondensed_local_dof_number,
                             &constraint_dofs](const Elem & elem,
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
      }

    if (node_num != invalid_uint && !elem.is_internal(node_num))
      uncondensed_dof = true;

    if (uncondensed_dof)
      this->add_uncondensed_dof_plus_constraint_dofs(field_dof,
                                                     false,
                                                     *uncondensed_global_to_local_map,
                                                     local_uncondensed_dofs_set,
                                                     nonlocal_uncondensed_dofs,
                                                     elem_uncondensed_dofs,
                                                     uncondensed_local_dof_number,
                                                     constraint_dofs);
    else
      (*condensed_global_to_local_map)[field_dof] = condensed_local_dof_number++;
  };

  for (auto elem : _mesh.active_local_element_ptr_range())
    {
      auto & dof_data = _elem_to_dof_data[elem->id()];
      condensed_local_dof_number = 0;
      uncondensed_local_dof_number = 0;
      condensed_global_to_local_map = &dof_data.condensed_global_to_local_map;
      uncondensed_global_to_local_map = &dof_data.uncondensed_global_to_local_map;

      const auto sub_id = elem->subdomain_id();
      for (const auto vg : make_range(_dof_map.n_variable_groups()))
        {
          const auto & var_group = _dof_map.variable_group(vg);
          for (const auto v : make_range(var_group.n_variables()))
            {
              const auto var_num = var_group.number(v);
              dof_data.reduced_space_indices.resize(var_num + 1);
              if (!var_group.active_on_subdomain(sub_id))
                continue;
              elem_uncondensed_dofs.clear();
              _dof_map.dof_indices(elem,
                                   elem_dofs,
                                   var_num,
                                   scalar_dofs_functor,
                                   field_dofs_functor,
                                   elem->p_level());
              if (!elem_uncondensed_dofs.empty())
                {
                  auto & var_reduced_space_indices = dof_data.reduced_space_indices[var_num];
                  var_reduced_space_indices.insert(var_reduced_space_indices.end(),
                                                   elem_uncondensed_dofs.begin(),
                                                   elem_uncondensed_dofs.end());
                  full_vars_present_in_reduced_sys.insert(var_num);
                }
            }
        }
    }

  //
  // We've built our local uncondensed dofs container ... but only using local element dof_indices
  // calls. It can be the case that we own a degree of freedom that is not actually needed by our
  // local element assembly but is needed by other processes element assembly. One example we've run
  // into of this is a mid-edge coarse element node holding side hierarchic dofs which is also a
  // fine element's vertex node. This node may be owned by the process holding the fine element
  // which doesn't need those side hierarchic dofs for its assembly
  //

  // Build supported query type. Has to be map to contiguous data for calls to MPI
  std::unordered_map<processor_id_type, std::vector<dof_id_type>> nonlocal_uncondensed_dofs_mapvec;
  for (const auto & [pid, set] : nonlocal_uncondensed_dofs)
    {
      auto & vec = nonlocal_uncondensed_dofs_mapvec[pid];
      vec.assign(set.begin(), set.end());
    }
  // clear no longer needed memory
  nonlocal_uncondensed_dofs.clear();

  auto receive_needed_local_dofs =
      [&local_uncondensed_dofs_set](processor_id_type,
                                    const std::vector<dof_id_type> & local_dofs_to_insert) {
        local_uncondensed_dofs_set.insert(local_dofs_to_insert.begin(), local_dofs_to_insert.end());
      };

  TIMPI::push_parallel_vector_data(
      _mesh.comm(), nonlocal_uncondensed_dofs_mapvec, receive_needed_local_dofs);

  _local_uncondensed_dofs.assign(local_uncondensed_dofs_set.begin(),
                                 local_uncondensed_dofs_set.end());
  local_uncondensed_dofs_set.clear();

  //
  // Build the reduced system data
  //

  const dof_id_type n_local = _local_uncondensed_dofs.size();
  dof_id_type n = n_local;
  this->comm().sum(n);

  // Get DOF counts on all processors
  this->compute_dof_info(n_local);

  // Build a map from the full size problem uncondensed dof indices to the reduced problem
  // (uncondensed) dof indices
  std::unordered_map<dof_id_type, dof_id_type> full_dof_to_reduced_dof;
  const auto local_start = _first_df[this->processor_id()];
  for (const auto i : index_range(_local_uncondensed_dofs))
    full_dof_to_reduced_dof[_local_uncondensed_dofs[i]] = i + local_start;

  // Build the condensed system sparsity pattern
  _reduced_sp = this->_dof_map.build_sparsity(
      this->_mesh, /*calculate_constrained=*/false, /*use_condensed_system=*/true);
  const auto & nnz = _reduced_sp->get_n_nz();
  const auto & noz = _reduced_sp->get_n_oz();
  libmesh_assert(nnz.size() == noz.size());

  // Optimization for PETSc. This is critical for problems in which there are SCALAR dofs that
  // introduce dense rows to avoid allocating a dense matrix
  _reduced_nnz.resize(_local_uncondensed_dofs.size());
  _reduced_noz.resize(_local_uncondensed_dofs.size());
  for (const dof_id_type local_reduced_i : index_range(_local_uncondensed_dofs))
    {
      const dof_id_type full_i = _local_uncondensed_dofs[local_reduced_i];
      const dof_id_type local_full_i = full_i - _dof_map.first_dof();
      libmesh_assert(local_full_i < nnz.size());
      _reduced_nnz[local_reduced_i] = nnz[local_full_i];
      _reduced_noz[local_reduced_i] = noz[local_full_i];
    }

  //
  // Now we need to pull our nonlocal data
  //

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

  // Determine the variables with any degrees of freedom present in the reduced system
  _communicator.set_union(full_vars_present_in_reduced_sys);
  _reduced_vars.reserve(full_vars_present_in_reduced_sys.size());
  unsigned int first_local_number = 0;
  for (const auto i : index_range(full_vars_present_in_reduced_sys))
    {
      const auto full_var_num = *std::next(full_vars_present_in_reduced_sys.begin(), i);
      const auto & full_var = _dof_map.variable(full_var_num);
      _reduced_vars.push_back(Variable{nullptr,
                                       full_var.name(),
                                       cast_int<unsigned int>(i),
                                       first_local_number,
                                       full_var.type()});
      first_local_number += _reduced_vars.back().n_components(_mesh);
    }

  // Now we can finally set our element reduced dof indices
  std::vector<dof_id_type> var_full_dof_indices;
  for (auto & [elem, dof_data] : _elem_to_dof_data)
    {
      libmesh_ignore(elem);
      auto & reduced_space_indices = dof_data.reduced_space_indices;
      // Keep around only those variables which are present in our reduced system
      {
        std::size_t i = 0;
        reduced_space_indices.erase(
            std::remove_if(
                reduced_space_indices.begin(),
                reduced_space_indices.end(),
                [&full_vars_present_in_reduced_sys, &i](const std::vector<dof_id_type> &) {
                  return !full_vars_present_in_reduced_sys.count(i++);
                }),
            reduced_space_indices.end());
      }
      libmesh_assert(reduced_space_indices.size() == full_vars_present_in_reduced_sys.size());

      for (auto & var_dof_indices : reduced_space_indices)
        {
          var_full_dof_indices = var_dof_indices;
          var_dof_indices.clear();
          for (const auto full_dof : var_full_dof_indices)
            var_dof_indices.push_back(libmesh_map_find(full_dof_to_reduced_dof, full_dof));
        }
    }

  // Build our dof constraints map
  for (const auto full_dof : constraint_dofs)
    _full_to_reduced_constraint_dofs[full_dof] =
        libmesh_map_find(full_dof_to_reduced_dof, full_dof);
  constraint_dofs.clear();

  // Prevent querying Nodes for dof indices
  std::vector<unsigned int> nvpg(_reduced_vars.size());
  for (auto & elem : nvpg)
    elem = 1;

  // add_system returns a system if it already exists instead of erroring so there's no harm if
  // we do this multiple times
  _reduced_system = &_system.get_equation_systems().add_system("Basic", "reduced");
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

  _sc_is_initialized = true;
}

unsigned int StaticCondensationDofMap::n_variables() const { return _reduced_vars.size(); }

const Variable & StaticCondensationDofMap::variable(const unsigned int c) const
{
  return _reduced_vars[c];
}

void StaticCondensationDofMap::dof_indices(const Elem * const elem,
                                           std::vector<dof_id_type> & di,
                                           const unsigned int vn,
                                           int /*p_level*/) const
{
  di.clear();
  di = libmesh_map_find(_elem_to_dof_data, elem->id()).reduced_space_indices[vn];
}

void StaticCondensationDofMap::dof_indices(const Node * const,
                                           std::vector<dof_id_type> &,
                                           const unsigned int) const
{
  libmesh_error_msg("StaticCondensationDofMap dof indices are only meant to be queried with "
                    "elements, not nodes");
}

void StaticCondensationDofMap::clear()
{
  DofMapBase::clear();
  _elem_to_dof_data.clear();
  _uncondensed_vars.clear();
  _reduced_vars.clear();
  _reduced_sp = nullptr;
  _reduced_nnz.clear();
  _reduced_noz.clear();
  _sc_is_initialized = false;
}
}
