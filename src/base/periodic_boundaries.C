// The libMesh Finite Element Library.
// Copyright (C) 2002-2020 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

// Local Includes
#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_ENABLE_PERIODIC

#include "libmesh/periodic_boundaries.h"

#include "libmesh/boundary_info.h"
#include "libmesh/elem.h"
#include "libmesh/mesh_base.h"
#include "libmesh/periodic_boundary.h"
#include "libmesh/point_locator_base.h"
#include "libmesh/remote_elem.h"

namespace libMesh
{

PeriodicBoundaries::~PeriodicBoundaries()
{
  for (auto & pr : *this)
    delete pr.second;
}



PeriodicBoundaryBase * PeriodicBoundaries::boundary(boundary_id_type id)
{
  iterator i = this->find(id);
  if (i == this->end())
    return nullptr;
  return i->second;
}



const PeriodicBoundaryBase * PeriodicBoundaries::boundary(boundary_id_type id) const
{
  const_iterator i = this->find(id);
  if (i == this->end())
    return nullptr;
  return i->second;
}




const Elem * PeriodicBoundaries::neighbor(boundary_id_type boundary_id,
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

  // If we should have found a periodic neighbor but didn't then
  // either we're on a ghosted element with a remote periodic neighbor
  // or we're on a mesh with an inconsistent periodic boundary.
  libmesh_assert_msg(!mesh.is_serial() &&
                     (e->processor_id() != mesh.processor_id()),
                     "Periodic boundary neighbor not found");
  return remote_elem;
}

} // namespace libMesh





#endif // LIBMESH_ENABLE_PERIODIC
