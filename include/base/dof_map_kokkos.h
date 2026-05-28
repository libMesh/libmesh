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

#ifndef LIBMESH_DOF_MAP_KOKKOS_H
#define LIBMESH_DOF_MAP_KOKKOS_H

#include "libmesh/dof_map.h"

#ifdef LIBMESH_HAVE_KOKKOS

#include "libmesh/kokkos_storage_policy.h"

namespace libMesh
{

struct DofMap::KokkosDofIndexCache
{
  using memory_space = typename ::Kokkos::DefaultExecutionSpace::memory_space;
  using elem_id_view = ::Kokkos::View<dof_id_type *, memory_space>;
  using elem_dof_id_view = ::Kokkos::View<dof_id_type **, memory_space>;
  using elem_dof_count_view = ::Kokkos::View<unsigned int *, memory_space>;
  using elem_subdomain_view = ::Kokkos::View<subdomain_id_type *, memory_space>;

  elem_id_view element_ids;
  elem_dof_id_view element_dof_indices;
  elem_dof_count_view element_n_dofs;
  elem_subdomain_view element_subdomains;
  std::vector<dof_id_type> host_element_ids;
  std::vector<dof_id_type> host_element_dof_indices;
  std::vector<unsigned int> host_element_n_dofs;
  std::vector<subdomain_id_type> host_element_subdomains;
  unsigned int max_dofs = 0;
};

struct DofMap::KokkosLocalIndexCache
{
  using memory_space = typename ::Kokkos::DefaultExecutionSpace::memory_space;
  using elem_local_index_view = ::Kokkos::View<unsigned int **, memory_space>;

  elem_local_index_view element_local_indices;
  unsigned int max_dofs = 0;
};

} // namespace libMesh

#endif

#endif
