// The libMesh Finite Element Library.
// Copyright (C) 2002-2021 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

PeriodicBoundaries::~PeriodicBoundaries() = default;



PeriodicBoundaryBase * PeriodicBoundaries::boundary(boundary_id_type id)
{
  iterator i = this->find(id);
  if (i == this->end())
    return nullptr;
  return i->second.get();
}



const PeriodicBoundaryBase * PeriodicBoundaries::boundary(boundary_id_type id) const
{
  const_iterator i = this->find(id);
  if (i == this->end())
    return nullptr;
  return i->second.get();
}




const Elem * PeriodicBoundaries::neighbor(boundary_id_type boundary_id,
                                          const PointLocatorBase & point_locator,
                                          const Elem * e,
                                          unsigned int side,
                                          unsigned int * neigh_side) const
{
  std::unique_ptr<const Elem> neigh_side_proxy;

  // Find a point on that side (and only that side)
  Point p = e->build_side_ptr(side)->vertex_average();

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
      std::vector<unsigned int> neigh_sides =
        mesh.get_boundary_info().sides_with_boundary_id(elem_it, b->pairedboundary);

      for (auto ns : neigh_sides)
        {
          if (neigh_side)
            {
              elem_it->build_side_ptr(neigh_side_proxy, ns);
              if (neigh_side_proxy->contains_point(p))
                {
                  *neigh_side = ns;
                  return elem_it;
                }
            }
          else
            // checking contains_point is too expensive if we don't
            // definitely need it to find neigh_side
            return elem_it;
        }
    }

  // If we should have found a periodic neighbor but didn't then
  // either we're on a ghosted element with a remote periodic neighbor
  // or we're on a mesh with an inconsistent periodic boundary.
  libmesh_error_msg_if(mesh.is_serial() ||
                       (e->processor_id() == mesh.processor_id()),
                       "Periodic boundary neighbor not found");

  if (neigh_side)
    *neigh_side = libMesh::invalid_uint;
  return remote_elem;
}

} // namespace libMesh





#endif // LIBMESH_ENABLE_PERIODIC
