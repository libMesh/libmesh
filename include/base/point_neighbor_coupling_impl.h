// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#ifndef LIBMESH_POINT_NEIGHBOR_COUPLING_IMPL_H
#define LIBMESH_POINT_NEIGHBOR_COUPLING_IMPL_H


// Local Includes
#include "libmesh/point_neighbor_coupling.h"

#include "libmesh/elem.h"
#include "libmesh/periodic_boundaries.h"
#include "libmesh/remote_elem.h"
#include "libmesh/libmesh_logging.h"

// C++ Includes
#include <unordered_set>

namespace libMesh
{

template <typename RealType>
void PointNeighborCouplingTempl<RealType>::mesh_reinit()
{
  // Unless we have periodic boundary conditions, we don't need
  // anything precomputed.
#ifdef LIBMESH_ENABLE_PERIODIC
  if (_periodic_bcs && !_periodic_bcs->empty())
#endif
    {
      // If we do have periodic boundary conditions, we'll need a master
      // point locator, so we'd better have a mesh to build it on.
      libmesh_assert(_mesh);

      // Make sure an up-to-date master point locator has been
      // constructed; we'll need to grab sub-locators soon.
      _mesh->sub_point_locator();
    }
}



template <typename RealType>
void PointNeighborCouplingTempl<RealType>::operator()
  (const typename MeshBaseTempl<RealType>::const_element_iterator & range_begin,
   const typename MeshBaseTempl<RealType>::const_element_iterator & range_end,
   processor_id_type p,
   map_type & coupled_elements)
{
  LOG_SCOPE("operator()", "PointNeighborCoupling");

#ifdef LIBMESH_ENABLE_PERIODIC
  bool check_periodic_bcs =
    (_periodic_bcs && !_periodic_bcs->empty());
  std::unique_ptr<PointLocatorBase> point_locator;
  if (check_periodic_bcs)
    {
      libmesh_assert(_mesh);
      point_locator = _mesh->sub_point_locator();
    }
#endif

  if (!this->_n_levels)
    {
      for (const auto & elem : as_range(range_begin, range_end))
        if (elem->processor_id() != p)
          coupled_elements.insert (std::make_pair(elem,_dof_coupling));

      return;
    }

  typedef std::unordered_set<const Elem*> set_type;
  set_type next_elements_to_check(range_begin, range_end);
  set_type elements_to_check;
  set_type elements_checked;

  for (unsigned int i=0; i != this->_n_levels; ++i)
    {
      elements_to_check.swap(next_elements_to_check);
      next_elements_to_check.clear();
      elements_checked.insert(elements_to_check.begin(), elements_to_check.end());

      for (const auto & elem : elements_to_check)
        {
          std::set<const Elem *> point_neighbors;

          if (elem->processor_id() != p)
            coupled_elements.insert (std::make_pair(elem,_dof_coupling));

#ifdef LIBMESH_ENABLE_PERIODIC
          // We might have a periodic neighbor here
          if (check_periodic_bcs)
            {
              libmesh_not_implemented();
            }
          else
#endif
            {
              elem->find_point_neighbors(point_neighbors);
            }

          for (const auto & neighbor : point_neighbors)
            {
              if (!elements_checked.count(neighbor))
                next_elements_to_check.insert(neighbor);

              if (neighbor->processor_id() != p)
                coupled_elements.insert
                  (std::make_pair(neighbor, _dof_coupling));
            }
        }
    }
}


} // namespace libMesh

#endif // LIBMESH_POINT_NEIGHBOR_COUPLING_IMPL_H
