// The libMesh Finite Element Library.
// Copyright (C) 2002-2016 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



// Local Includes -----------------------------------
#include "libmesh/default_coupling.h"

#include "libmesh/elem.h"
#include "libmesh/remote_elem.h"

namespace libMesh
{

//-------------------------------------------------------
// we need to implement these constructors here so that
// a full DofMap definition is available.
void DefaultCoupling::operator()
  (const MeshBase::const_element_iterator & range_begin,
   const MeshBase::const_element_iterator & range_end,
   processor_id_type p,
   map_type & coupled_elements)
{
  for (MeshBase::const_element_iterator elem_it = range_begin;
       elem_it != range_end; ++elem_it)
    {
      std::vector<const Elem *> active_neighbors;

      const Elem * const elem = *elem_it;

      if (elem->processor_id() != p)
        coupled_elements.insert (std::make_pair(elem,_dof_coupling));

      if (_couple_neighbor_dofs)
        for (unsigned int s=0; s<elem->n_sides(); s++)
          if (elem->neighbor_ptr(s) != libmesh_nullptr)
            {
              const Elem * const neighbor_0 = elem->neighbor_ptr(s);

              // Mesh ghosting might ask us about what we want to
              // distribute along with non-local elements, and those
              // non-local elements might have remote neighbors, and
              // if they do then we can't say anything about them.
              if (neighbor_0 == remote_elem)
                continue;

#ifdef LIBMESH_ENABLE_AMR
              neighbor_0->active_family_tree_by_neighbor(active_neighbors,elem);
#else
              active_neighbors.clear();
              active_neighbors.push_back(neighbor_0);
#endif

              for (std::size_t a=0; a != active_neighbors.size(); ++a)
                {
                  const Elem * neighbor = active_neighbors[a];

                  if (neighbor->processor_id() != p)
                    coupled_elements.insert
                      (std::make_pair(neighbor, _dof_coupling));
                }
            }
    }
}


} // namespace libMesh
