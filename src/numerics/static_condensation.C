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

namespace libMesh
{
template <typename T>
StaticCondensation<T>::StaticCondensation(const MeshBase & mesh, const DofMap & dof_map)
  : _mesh(mesh), _dof_map(dof_map)
{
  build_idx_maps();
}

template <typename T>
void
StaticCondensation<T>::build_idx_maps()
{
  std::vector<dof_id_type> di, dii, dib;
  for (auto elem : _mesh.active_local_element_ptr_range())
  {
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
        _dof_map.dof_indices(elem, di, v, elem->p_level(), &dib, &dii);
#ifdef DEBUG
        auto merge_di = dib;
        merge_di.insert(merge_di.end(), dii.begin(), dii.end());
        libmesh_assert(merge_di == di);
#endif
        var_data.boundary_dofs_offset = boundary_dof_size;
        var_data.interior_dofs_offset = interior_dof_size;
        var_data.num_boundary_dofs = dib.size();
        var_data.num_interior_dofs = dii.size();
        boundary_dof_size += var_data.num_boundary_dofs;
        interior_dof_size += var_data.num_interior_dofs;
      }
    }

    local_data.Aii.resize(interior_dof_size, interior_dof_size);
    local_data.Aib.resize(interior_dof_size, boundary_dof_size);
    local_data.Abi.resize(boundary_dof_size, interior_dof_size);
    local_data.Abb.resize(boundary_dof_size, boundary_dof_size);
  }
}

template <typename T>
void
StaticCondensation<T>::add_matrix(const Elem & elem,
                                  const unsigned int i_var,
                                  const unsigned int j_var,
                                  const DenseMatrix<T> & k)
{
  auto & local_data = libmesh_map_find(_elem_to_local_data, elem.id());
  auto & ivar_data = libmesh_map_find(local_data.var_to_data, i_var);
  auto & jvar_data = libmesh_map_find(local_data.var_to_data, j_var);
  auto copy_mat = [&k](DenseMatrix<T> & dest_mat,
                       const unsigned int i_offset,
                       const unsigned int i_start,
                       const unsigned int i_stop,
                       const unsigned int j_offset,
                       const unsigned int j_start,
                       const unsigned int j_stop)
  {
    for (const auto i : make_range(i_start, i_stop))
      for (const auto j : make_range(j_start, j_stop))
        dest_mat(i_offset + i, j_offset + j) = k(i, j);
  };

  copy_mat(local_data.Abb,
           ivar_data.boundary_dofs_offset,
           0,
           ivar_data.num_boundary_dofs,
           jvar_data.boundary_dofs_offset,
           0,
           jvar_data.num_boundary_dofs);
  copy_mat(local_data.Abi,
           ivar_data.boundary_dofs_offset,
           0,
           ivar_data.num_boundary_dofs,
           jvar_data.interior_dofs_offset,
           jvar_data.num_boundary_dofs,
           jvar_data.num_boundary_dofs + jvar_data.num_interior_dofs);
  copy_mat(local_data.Aib,
           ivar_data.interior_dofs_offset,
           ivar_data.num_boundary_dofs,
           ivar_data.num_boundary_dofs + ivar_data.num_interior_dofs,
           jvar_data.boundary_dofs_offset,
           0,
           jvar_data.num_boundary_dofs);
  copy_mat(local_data.Aii,
           ivar_data.interior_dofs_offset,
           ivar_data.num_boundary_dofs,
           ivar_data.num_boundary_dofs + ivar_data.num_interior_dofs,
           jvar_data.interior_dofs_offset,
           jvar_data.num_boundary_dofs,
           jvar_data.num_boundary_dofs + jvar_data.num_interior_dofs);
}
}
