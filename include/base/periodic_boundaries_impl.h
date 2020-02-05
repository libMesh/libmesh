// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#ifndef LIBMESH_PERIODIC_BOUNDARIES_IMPL_H
#define LIBMESH_PERIODIC_BOUNDARIES_IMPL_H

// Local Includes
#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_ENABLE_PERIODIC

#include "libmesh/periodic_boundaries.h"
#include "libmesh/point_locator_base.h"
#include "libmesh/elem.h"
#include "libmesh/periodic_boundary.h"
#include "libmesh/mesh_base.h"

namespace libMesh
{

template <typename RealType>
PeriodicBoundariesTempl<RealType>::~PeriodicBoundariesTempl()
{
  for (auto & pr : *this)
    delete pr.second;
}



template <typename RealType>
PeriodicBoundaryBaseTempl<RealType> * PeriodicBoundariesTempl<RealType>::boundary(boundary_id_type id)
{
  iterator i = this->find(id);
  if (i == this->end())
    return nullptr;
  return i->second;
}



template <typename RealType>
const PeriodicBoundaryBaseTempl<RealType> * PeriodicBoundariesTempl<RealType>::boundary(boundary_id_type id) const
{
  const_iterator i = this->find(id);
  if (i == this->end())
    return nullptr;
  return i->second;
}

template <typename RealType>
const ElemTempl<RealType> *
PeriodicBoundariesTempl<RealType>::neighbor(boundary_id_type boundary_id,
                                            const PointLocatorBase & point_locator,
                                            const Elem * e,
                                            unsigned int side) const
{
  // Find a point on that side (and only that side)

  Point p = e->build_side_ptr(side)->centroid();

  const PeriodicBoundaryBase * b = this->boundary(boundary_id);
  libmesh_assert (b);
  p = b->get_corresponding_pos(p);

  std::set<const Elem *> candidate_elements;
  point_locator.operator()(p, candidate_elements);

  // We might have found multiple elements, e.g. if two distinct periodic
  // boundaries are overlapping (see systems_of_equations_ex9, for example).
  // As a result, we need to search for the element that has boundary_id.
  const MeshBase & mesh = point_locator.get_mesh();
  for(const Elem * elem_it : candidate_elements)
    {
      unsigned int s_neigh =
        mesh.get_boundary_info().side_with_boundary_id(elem_it, b->pairedboundary);

      // If s_neigh is not invalid then we have found an element that contains
      // boundary_id, so return this element
      if(s_neigh != libMesh::invalid_uint)
        {
          return elem_it;
        }
    }

  libmesh_error_msg("Periodic boundary neighbor not found");
  return nullptr;
}

} // namespace libMesh





#endif // LIBMESH_ENABLE_PERIODIC

#endif // LIBMESH_PERIODIC_BOUNDARIES_IMPL_H
