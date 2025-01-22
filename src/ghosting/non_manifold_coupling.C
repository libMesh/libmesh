// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

// libMesh includes
#include "libmesh/elem.h"
#include "libmesh/non_manifold_coupling.h"
#include "libmesh/simple_range.h" // as_range()

namespace libMesh
{

NonManifoldGhostingFunctor::
NonManifoldGhostingFunctor(const MeshBase & mesh) :
  GhostingFunctor(mesh),
  _stem(MeshTools::SidesToElemMap::build(mesh))
{
  // We call MeshBase::reinit_ghosting_functors() in MeshBase::prepare_for_use(),
  // so it _should_ be safe to omit this call to SidesToElemMap::build() here,
  // but to be on the safe side, we'll construct it both here and during any call
  // to mesh_reinit().
}

std::unique_ptr<GhostingFunctor>
NonManifoldGhostingFunctor::clone () const
{
  return std::make_unique<NonManifoldGhostingFunctor>(*this);
}

void
NonManifoldGhostingFunctor::mesh_reinit ()
{
  _stem = MeshTools::SidesToElemMap::build(*_mesh);
}

void
NonManifoldGhostingFunctor::redistribute ()
{
  this->mesh_reinit();
}

void
NonManifoldGhostingFunctor::delete_remote_elements()
{
  this->mesh_reinit();
}

void
NonManifoldGhostingFunctor::operator() (
  const MeshBase::const_element_iterator & range_begin,
  const MeshBase::const_element_iterator & range_end,
  processor_id_type p,
  map_type & coupled_elements)
{
  // This code is just for geometric coupling, so we use a null
  // CouplingMatrix pointer.  We'll declare that here so as to avoid
  // confusing the emplace() calls later.
  CouplingMatrix * nullcm = nullptr;

  for (const auto & elem : as_range(range_begin, range_end))
    {
      // The idea here is that the input range [range_begin,
      // range_end] is a list of Elems whose DOFs we already know we
      // want to ghost (if they aren't already owned by processor p,
      // we don't need to worry about ghosting "local" Elem DOFs)
      // Therefore, we always put *all* the non-local Elems in the
      // input range into the output "coupled_elements" map.  The
      // purpose of the rest of this of this function is then to
      // determine which *additional* elements should be ghosted,
      // given that we care about the entire input range.
      if (elem->processor_id() != p)
        coupled_elements.emplace(elem, nullcm);

      // For each side of "elem", look up that (elem, side) pair in
      // the SidesToElemMap and then add all the Elem pointers from
      // that list to the coupled_elements map.
      for (auto s : elem->side_index_range())
        for (const auto & neigh : as_range(_stem.get_connected_elems(elem, s)))
          if (neigh->processor_id() != p)
            coupled_elements.emplace(neigh, nullcm);
    }
}

} // namespace libMesh
