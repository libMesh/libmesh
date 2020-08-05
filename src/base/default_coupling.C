// The libMesh Finite Element Library.
// Copyright (C) 2002-2020 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



// Local Includes
#include "libmesh/default_coupling.h"
#include "libmesh/coupling_matrix.h"
#include "libmesh/elem.h"
#include "libmesh/periodic_boundaries.h"
#include "libmesh/remote_elem.h"
#include "libmesh/int_range.h"
#include "libmesh/libmesh_logging.h"

// C++ Includes
#include <unordered_set>

namespace libMesh
{

void DefaultCoupling::set_dof_coupling(const CouplingMatrix * dof_coupling)
{
  // We used to treat an empty 0x0 _dof_coupling matrix as if it
  // were an NxN all-ones matrix.  We'd like to stop supporting this
  // behavior, but for now we'll just warn about it, while supporting
  // it via the preferred mechanism: a nullptr _dof_coupling
  // matrix pointer is interpreted as a full coupling matrix.
  if (dof_coupling && dof_coupling->empty())
    {
      libmesh_deprecated();
      _dof_coupling = nullptr;
    }
  else
    _dof_coupling = dof_coupling;
}



void DefaultCoupling::mesh_reinit()
{
  // Unless we have periodic boundary conditions, we don't need
  // anything precomputed.
#ifdef LIBMESH_ENABLE_PERIODIC
  if (!_periodic_bcs || _periodic_bcs->empty())
    return;
#endif

  // If we do have periodic boundary conditions, we'll need a master
  // point locator, so we'd better have a mesh to build it on.
  libmesh_assert(_mesh);

  // Make sure an up-to-date master point locator has been
  // constructed; we'll need to grab sub-locators soon.
  _mesh->sub_point_locator();
}



void DefaultCoupling::operator()
  (const MeshBase::const_element_iterator & range_begin,
   const MeshBase::const_element_iterator & range_end,
   processor_id_type p,
   map_type & coupled_elements)
{
  LOG_SCOPE("operator()", "DefaultCoupling");

  libmesh_assert(_mesh);

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
      {
        libmesh_assert(_mesh->query_elem_ptr(elem->id()) ==elem);
        if (elem->processor_id() != p)
          coupled_elements.emplace(elem, _dof_coupling);
      }
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
          std::vector<const Elem *> active_neighbors;

          libmesh_assert(_mesh->query_elem_ptr(elem->id()) ==elem);

          if (elem->processor_id() != p)
            coupled_elements.emplace(elem, _dof_coupling);

          for (auto s : elem->side_index_range())
            {
              const Elem * neigh = elem->neighbor_ptr(s);

              // If we have a neighbor here
              if (neigh)
                {
                  // Mesh ghosting might ask us about what we want to
                  // distribute along with non-local elements, and those
                  // non-local elements might have remote neighbors, and
                  // if they do then we can't say anything about them.
                  if (neigh == remote_elem)
                    continue;
                }
#ifdef LIBMESH_ENABLE_PERIODIC
              // We might still have a periodic neighbor here
              else if (check_periodic_bcs)
                {
                  libmesh_assert(_mesh);

                  neigh = elem->topological_neighbor
                    (s, *_mesh, *point_locator, _periodic_bcs);
                }
#endif

              // With no regular *or* periodic neighbors we have nothing
              // to do.
              if (!neigh)
                continue;

              // With any kind of neighbor, we need to couple to all the
              // active descendants on our side.
#ifdef LIBMESH_ENABLE_AMR
              if (neigh == elem->neighbor_ptr(s))
                neigh->active_family_tree_by_neighbor(active_neighbors,elem);
#  ifdef LIBMESH_ENABLE_PERIODIC
              else
                neigh->active_family_tree_by_topological_neighbor
                  (active_neighbors,elem,*_mesh,*point_locator,_periodic_bcs);
#  endif
#else
              active_neighbors.clear();
              active_neighbors.push_back(neigh);
#endif

              for (const auto & neighbor : active_neighbors)
                {
                  if (!elements_checked.count(neighbor))
                    next_elements_to_check.insert(neighbor);

                  if (neighbor->processor_id() != p)
                    coupled_elements.emplace(neighbor, _dof_coupling);
                }
            }
        }
    }
}


} // namespace libMesh
