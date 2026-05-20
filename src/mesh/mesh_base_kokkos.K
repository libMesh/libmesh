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

#include "libmesh/mesh_base_kokkos.h"

#include "libmesh/elem.h"
#include "libmesh/utility.h"

#include <algorithm>

namespace libMesh
{

#ifdef LIBMESH_HAVE_KOKKOS
void MeshBase::delete_kokkos_geometry_cache(KokkosGeometryCache * cache)
{
  delete cache;
}

const MeshBase::KokkosGeometryCache &
MeshBase::get_kokkos_geometry_cache() const
{
  if (_kokkos_geometry_cache)
    return *_kokkos_geometry_cache;

  auto cache = kokkos_geometry_cache_ptr(new KokkosGeometryCache(),
                                         &MeshBase::delete_kokkos_geometry_cache);
  cache->host_element_ids.reserve(this->n_active_local_elem());

  for (const auto & elem : this->active_local_element_ptr_range())
    {
      cache->element_lookup.emplace(elem->id(), cast_int<unsigned int>(cache->host_element_ids.size()));
      cache->host_element_ids.push_back(elem->id());
      cache->max_nodes = std::max(cache->max_nodes, elem->n_nodes());

      for (unsigned int n = 0; n != elem->n_nodes(); ++n)
        {
          const dof_id_type node_id = elem->node_id(n);
          if (!cache->node_lookup.count(node_id))
            {
              cache->node_lookup.emplace(node_id, cast_int<unsigned int>(cache->host_node_ids.size()));
              cache->host_node_ids.push_back(node_id);
            }
        }
    }

  cache->node_ids =
    KokkosGeometryCache::node_id_view("mesh_kokkos_node_ids", cache->host_node_ids.size());
  cache->element_ids =
    KokkosGeometryCache::elem_id_view("mesh_kokkos_element_ids", cache->host_element_ids.size());
  cache->node_coordinates =
    KokkosGeometryCache::node_coord_view("mesh_kokkos_node_coordinates",
                                         cache->host_node_ids.size(),
                                         LIBMESH_DIM);
  cache->element_node_ids =
    KokkosGeometryCache::elem_node_id_view("mesh_kokkos_element_node_ids",
                                           cache->host_element_ids.size(),
                                           cache->max_nodes);
  cache->element_types =
    KokkosGeometryCache::elem_type_view("mesh_kokkos_element_types", cache->host_element_ids.size());
  cache->element_mapping_types =
    KokkosGeometryCache::elem_mapping_type_view("mesh_kokkos_element_mapping_types",
                                                cache->host_element_ids.size());
  cache->element_n_nodes =
    KokkosGeometryCache::elem_n_nodes_view("mesh_kokkos_element_n_nodes", cache->host_element_ids.size());
  cache->element_p_levels =
    KokkosGeometryCache::elem_p_level_view("mesh_kokkos_element_p_levels", cache->host_element_ids.size());
  cache->element_subdomains =
    KokkosGeometryCache::elem_subdomain_view("mesh_kokkos_element_subdomains",
                                             cache->host_element_ids.size());

  auto h_node_ids = ::Kokkos::create_mirror_view(cache->node_ids);
  auto h_element_ids = ::Kokkos::create_mirror_view(cache->element_ids);
  auto h_node_coordinates = ::Kokkos::create_mirror_view(cache->node_coordinates);
  auto h_element_node_ids = ::Kokkos::create_mirror_view(cache->element_node_ids);
  auto h_element_types = ::Kokkos::create_mirror_view(cache->element_types);
  auto h_element_mapping_types = ::Kokkos::create_mirror_view(cache->element_mapping_types);
  auto h_element_n_nodes = ::Kokkos::create_mirror_view(cache->element_n_nodes);
  auto h_element_p_levels = ::Kokkos::create_mirror_view(cache->element_p_levels);
  auto h_element_subdomains = ::Kokkos::create_mirror_view(cache->element_subdomains);

  for (auto node_index : index_range(cache->host_node_ids))
    {
      const dof_id_type node_id = cache->host_node_ids[node_index];
      const Node & node = *this->query_node_ptr(node_id);
      h_node_ids(cast_int<unsigned int>(node_index)) = node_id;
      for (unsigned int component = 0; component != LIBMESH_DIM; ++component)
        h_node_coordinates(cast_int<unsigned int>(node_index), component) = node(component);
    }

  for (auto elem_index : index_range(cache->host_element_ids))
    {
      const dof_id_type elem_id = cache->host_element_ids[elem_index];
      const Elem & elem = *this->query_elem_ptr(elem_id);
      h_element_ids(cast_int<unsigned int>(elem_index)) = elem_id;
      h_element_types(cast_int<unsigned int>(elem_index)) = elem.type();
      h_element_mapping_types(cast_int<unsigned int>(elem_index)) = elem.mapping_type();
      h_element_n_nodes(cast_int<unsigned int>(elem_index)) = elem.n_nodes();
      h_element_p_levels(cast_int<unsigned int>(elem_index)) = elem.p_level();
      h_element_subdomains(cast_int<unsigned int>(elem_index)) = elem.subdomain_id();
      for (unsigned int n = 0; n != elem.n_nodes(); ++n)
        h_element_node_ids(cast_int<unsigned int>(elem_index), n) =
          libmesh_map_find(cache->node_lookup, elem.node_id(n));
    }

  ::Kokkos::deep_copy(cache->node_ids, h_node_ids);
  ::Kokkos::deep_copy(cache->element_ids, h_element_ids);
  ::Kokkos::deep_copy(cache->node_coordinates, h_node_coordinates);
  ::Kokkos::deep_copy(cache->element_node_ids, h_element_node_ids);
  ::Kokkos::deep_copy(cache->element_types, h_element_types);
  ::Kokkos::deep_copy(cache->element_mapping_types, h_element_mapping_types);
  ::Kokkos::deep_copy(cache->element_n_nodes, h_element_n_nodes);
  ::Kokkos::deep_copy(cache->element_p_levels, h_element_p_levels);
  ::Kokkos::deep_copy(cache->element_subdomains, h_element_subdomains);
  _kokkos_geometry_cache = std::move(cache);
  return *_kokkos_geometry_cache;
}

void MeshBase::prepare_kokkos_geometry_cache() const
{
  libmesh_ignore(this->get_kokkos_geometry_cache());
}

unsigned int MeshBase::get_kokkos_elem_index(const Elem & elem) const
{
  const auto & cache = this->get_kokkos_geometry_cache();
  if (auto it = cache.element_lookup.find(elem.id()); it != cache.element_lookup.end())
    return it->second;

  return libMesh::invalid_uint;
}
#endif

} // namespace libMesh
