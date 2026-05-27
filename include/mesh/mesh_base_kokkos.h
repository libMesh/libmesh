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

#ifndef LIBMESH_MESH_BASE_KOKKOS_H
#define LIBMESH_MESH_BASE_KOKKOS_H

#include "libmesh/mesh_base.h"

#ifdef LIBMESH_HAVE_KOKKOS

#include "libmesh/kokkos_storage_policy.h"

namespace libMesh
{

struct MeshBase::KokkosGeometryCache
{
  using node_id_view = ::Kokkos::View<dof_id_type *>;
  using elem_id_view = ::Kokkos::View<dof_id_type *>;
  using node_coord_view = ::Kokkos::View<Real **>;
  using elem_node_id_view = ::Kokkos::View<unsigned int **>;
  using elem_type_view = ::Kokkos::View<ElemType *>;
  using elem_mapping_type_view = ::Kokkos::View<ElemMappingType *>;
  using elem_n_nodes_view = ::Kokkos::View<unsigned int *>;
  using elem_p_level_view = ::Kokkos::View<unsigned int *>;
  using elem_subdomain_view = ::Kokkos::View<subdomain_id_type *>;

  node_id_view node_ids;
  elem_id_view element_ids;
  node_coord_view node_coordinates;
  elem_node_id_view element_node_ids;
  elem_type_view element_types;
  elem_mapping_type_view element_mapping_types;
  elem_n_nodes_view element_n_nodes;
  elem_p_level_view element_p_levels;
  elem_subdomain_view element_subdomains;
  std::vector<dof_id_type> host_node_ids;
  std::vector<dof_id_type> host_element_ids;
  std::vector<ElemType> host_element_types;
  std::vector<ElemMappingType> host_element_mapping_types;
  std::vector<unsigned int> host_element_n_nodes;
  std::vector<unsigned int> host_element_p_levels;
  std::vector<subdomain_id_type> host_element_subdomains;
  std::unordered_map<dof_id_type, unsigned int> node_lookup;
  std::unordered_map<dof_id_type, unsigned int> element_lookup;
  unsigned int max_nodes = 0;
};

} // namespace libMesh

#endif

#endif
