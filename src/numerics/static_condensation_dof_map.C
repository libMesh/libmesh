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

#include "libmesh/static_condensation_dof_map.h"

#include "libmesh/elem.h"
#include "libmesh/mesh_base.h"

namespace libMesh
{
StaticCondensationDofMap::StaticCondensationDofMap(MeshBase & mesh,
                                                   System & system,
                                                   const DofMap & dof_map)
  : ReducedDofMap(mesh, system, dof_map, "reduced")
{
  libmesh_experimental();
}

StaticCondensationDofMap::~StaticCondensationDofMap() = default;

bool StaticCondensationDofMap::retain_field_dof(const Elem & elem,
                                                const unsigned int node_num,
                                                const unsigned int var_num) const
{
  if (_uncondensed_vars.count(var_num))
    {
      libmesh_assert_msg(
          node_num == invalid_uint,
          "Users should not be providing continuous FEM variables to the uncondensed vars API");
      return true;
    }

  return node_num != invalid_uint && !elem.is_internal(node_num);
}

void StaticCondensationDofMap::retained_dof_added(const Elem & elem, const dof_id_type full_dof)
{
  auto & map = _elem_to_dof_data[elem.id()].uncondensed_global_to_local_map;
  const auto local_number = cast_int<dof_id_type>(map.size());
  map.emplace(full_dof, local_number);
}

void StaticCondensationDofMap::unretained_dof_found(const Elem & elem, const dof_id_type full_dof)
{
  auto & map = _elem_to_dof_data[elem.id()].condensed_global_to_local_map;
  const auto local_number = cast_int<dof_id_type>(map.size());
  map.emplace(full_dof, local_number);
}

void StaticCondensationDofMap::reinit()
{
  ReducedDofMap::reinit();

  // StaticCondensation expects an entry for every active local element, including those without
  // any degrees of freedom
  for (const auto * const elem : _mesh.active_local_element_ptr_range())
    _elem_to_dof_data[elem->id()];
}

void StaticCondensationDofMap::clear()
{
  ReducedDofMap::clear();
  _elem_to_dof_data.clear();
  _uncondensed_vars.clear();
}
}
