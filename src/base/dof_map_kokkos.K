// The libMesh Finite Element Library.
// Copyright (C) 2002-2026 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public  License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

#include "libmesh/dof_map_kokkos.h"

#include "libmesh/elem.h"
#include "libmesh/mesh_base_kokkos.h"

#ifdef LIBMESH_HAVE_PETSC
#include "libmesh/petsc_vector.h"
#endif

#include <algorithm>
#include <memory>
#include <vector>

namespace libMesh
{

#ifdef LIBMESH_HAVE_KOKKOS
void DofMap::KokkosDofIndexCacheDeleter::operator()(KokkosDofIndexCache * cache) const
{
  delete cache;
}

void DofMap::KokkosLocalIndexCacheDeleter::operator()(KokkosLocalIndexCache * cache) const
{
  delete cache;
}

const DofMap::KokkosDofIndexCache *
DofMap::get_kokkos_dof_index_cache(const unsigned int vn) const
{
  if (auto it = _kokkos_dof_index_caches.find(vn); it != _kokkos_dof_index_caches.end())
    return it->second.get();

  const auto & geometry_cache = _mesh.get_kokkos_geometry_cache();
  auto cache = kokkos_dof_index_cache_ptr(new KokkosDofIndexCache());
  std::vector<dof_id_type> dof_indices;
  cache->host_element_ids = geometry_cache.host_element_ids;

  for (const auto elem_id : cache->host_element_ids)
    {
      const Elem * elem = _mesh.query_elem_ptr(elem_id);
      libmesh_assert(elem);
      this->dof_indices(elem, dof_indices, vn);
      cache->max_dofs = std::max(cache->max_dofs, cast_int<unsigned int>(dof_indices.size()));
    }

  cache->element_ids =
    KokkosDofIndexCache::elem_id_view("dof_map_kokkos_element_ids",
                                      cache->host_element_ids.size());
  cache->element_dof_indices =
    KokkosDofIndexCache::elem_dof_id_view("dof_map_kokkos_dof_indices",
                                          cache->host_element_ids.size(),
                                          cache->max_dofs);
  cache->element_n_dofs =
    KokkosDofIndexCache::elem_dof_count_view("dof_map_kokkos_n_dofs",
                                             cache->host_element_ids.size());
  cache->element_subdomains =
    KokkosDofIndexCache::elem_subdomain_view("dof_map_kokkos_elem_subdomains",
                                             cache->host_element_ids.size());

  auto h_element_ids = ::Kokkos::create_mirror_view(cache->element_ids);
  auto h_dof_indices = ::Kokkos::create_mirror_view(cache->element_dof_indices);
  auto h_n_dofs = ::Kokkos::create_mirror_view(cache->element_n_dofs);
  auto h_subdomains = ::Kokkos::create_mirror_view(cache->element_subdomains);
  cache->host_element_dof_indices.resize(cache->host_element_ids.size() * cache->max_dofs,
                                         DofObject::invalid_id);
  cache->host_element_n_dofs.resize(cache->host_element_ids.size(), 0);
  cache->host_element_subdomains.resize(cache->host_element_ids.size(),
                                        static_cast<subdomain_id_type>(Elem::invalid_subdomain_id));

  for (auto elem_index : index_range(cache->host_element_ids))
    {
      const dof_id_type elem_id = cache->host_element_ids[elem_index];
      const Elem * elem = _mesh.query_elem_ptr(elem_id);
      libmesh_assert(elem);
      h_element_ids(cast_int<unsigned int>(elem_index)) = elem_id;
      h_subdomains(cast_int<unsigned int>(elem_index)) = elem->subdomain_id();
      cache->host_element_subdomains[elem_index] = elem->subdomain_id();
      this->dof_indices(elem, dof_indices, vn);
      h_n_dofs(cast_int<unsigned int>(elem_index)) = cast_int<unsigned int>(dof_indices.size());
      cache->host_element_n_dofs[elem_index] = cast_int<unsigned int>(dof_indices.size());
      for (auto i : index_range(dof_indices))
        {
          h_dof_indices(cast_int<unsigned int>(elem_index), cast_int<unsigned int>(i)) =
            dof_indices[i];
          cache->host_element_dof_indices[elem_index * cache->max_dofs + i] = dof_indices[i];
        }
    }

  ::Kokkos::deep_copy(cache->element_ids, h_element_ids);
  ::Kokkos::deep_copy(cache->element_dof_indices, h_dof_indices);
  ::Kokkos::deep_copy(cache->element_n_dofs, h_n_dofs);
  ::Kokkos::deep_copy(cache->element_subdomains, h_subdomains);

  auto [it, inserted] = _kokkos_dof_index_caches.emplace(vn, std::move(cache));
  libmesh_ignore(inserted);
  return it->second.get();
}

const DofMap::KokkosLocalIndexCache *
DofMap::get_kokkos_local_index_cache(const NumericVector<Number> & local_vector,
                                     const unsigned int vn) const
{
#ifdef LIBMESH_HAVE_PETSC
  const auto key = std::make_pair(vn, &local_vector);
  if (auto it = _kokkos_local_index_caches.find(key); it != _kokkos_local_index_caches.end())
    return it->second.get();

  const auto * petsc_vector = dynamic_cast<const PetscVector<Number> *>(&local_vector);
  if (!petsc_vector)
    return nullptr;

  const auto * dof_index_cache = this->get_kokkos_dof_index_cache(vn);
  if (!dof_index_cache)
    return nullptr;

  const auto & geometry_cache = _mesh.get_kokkos_geometry_cache();
  auto cache = kokkos_local_index_cache_ptr(new KokkosLocalIndexCache());
  cache->max_dofs = dof_index_cache->max_dofs;

  cache->element_local_indices =
    KokkosLocalIndexCache::elem_local_index_view("dof_map_kokkos_local_indices",
                                                 geometry_cache.host_element_ids.size(),
                                                 cache->max_dofs);
  auto h_local_indices = ::Kokkos::create_mirror_view(cache->element_local_indices);

  for (auto elem_index : index_range(geometry_cache.host_element_ids))
    {
      const unsigned int n_dofs = dof_index_cache->host_element_n_dofs[elem_index];
      for (unsigned int i = 0; i != n_dofs; ++i)
        h_local_indices(cast_int<unsigned int>(elem_index), i) =
          cast_int<unsigned int>(
            petsc_vector->map_global_to_local_index(
              dof_index_cache->host_element_dof_indices[elem_index * cache->max_dofs + i]));
    }

  ::Kokkos::deep_copy(cache->element_local_indices, h_local_indices);
  auto [it, inserted] = _kokkos_local_index_caches.emplace(key, std::move(cache));
  libmesh_ignore(inserted);
  return it->second.get();
#else
  libmesh_ignore(local_vector, vn);
  return nullptr;
#endif
}

const DofMap::KokkosLocalIndexCache *
DofMap::require_kokkos_local_index_cache(const NumericVector<Number> & local_vector,
                                         const unsigned int vn) const
{
  this->prepare_kokkos_local_index_cache(local_vector, vn);
  return this->get_kokkos_local_index_cache(local_vector, vn);
}

void DofMap::prepare_kokkos_dof_index_caches() const
{
  for (auto vn : make_range(this->n_variables()))
    libmesh_ignore(this->get_kokkos_dof_index_cache(vn));
}

void
DofMap::prepare_kokkos_local_index_cache(const NumericVector<Number> & local_vector,
                                         const unsigned int vn) const
{
  libmesh_ignore(this->get_kokkos_local_index_cache(local_vector, vn));
}
#endif

} // namespace libMesh
