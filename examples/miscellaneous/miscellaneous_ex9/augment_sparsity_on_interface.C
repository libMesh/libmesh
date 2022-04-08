// local includes
#include "augment_sparsity_on_interface.h"

// libMesh includes
#include "libmesh/elem.h"
#include "libmesh/boundary_info.h"
#include "libmesh/elem_side_builder.h"

using namespace libMesh;

AugmentSparsityOnInterface::AugmentSparsityOnInterface(MeshBase & mesh,
                                                       boundary_id_type crack_boundary_lower,
                                                       boundary_id_type crack_boundary_upper) :
  _mesh(mesh),
  _crack_boundary_lower(crack_boundary_lower),
  _crack_boundary_upper(crack_boundary_upper),
  _initialized(false)
{
  this->mesh_reinit();
}

const ElementSideMap & AugmentSparsityOnInterface::get_lower_to_upper() const
{
  libmesh_assert(this->_initialized);
  return _lower_to_upper;
}

void AugmentSparsityOnInterface::mesh_reinit ()
{
  this->_initialized = true;

  // Loop over all elements (not just local elements) to make sure we find
  // "neighbor" elements on opposite sides of the crack.

  // Map from (elem, side) to vertex average
  std::map<std::pair<const Elem *, unsigned char>, Point> lower;
  std::map<std::pair<const Elem *, unsigned char>, Point> upper;

  // To avoid extraneous allocation when getting element side vertex averages
  std::unique_ptr<const Elem> side_elem;

  for (const auto & elem : _mesh.active_element_ptr_range())
    for (auto side : elem->side_index_range())
      if (elem->neighbor_ptr(side) == nullptr)
        {
          if (_mesh.get_boundary_info().has_boundary_id(elem, side, _crack_boundary_lower))
            {
              elem->build_side_ptr(side_elem, side);

              lower[std::make_pair(elem, side)] = side_elem->vertex_average();
            }

          if (_mesh.get_boundary_info().has_boundary_id(elem, side, _crack_boundary_upper))
            {
              elem->build_side_ptr(side_elem, side);

              upper[std::make_pair(elem, side)] = side_elem->vertex_average();
            }
        }

  // If we're doing a reinit on a distributed mesh then we may not see
  // all the vertex averages, or even a matching number of vertex averages.
  // std::size_t n_lower = lower.size();
  // std::size_t n_upper = upper.size();
  // libmesh_assert(n_lower == n_upper);

  // Clear _lower_to_upper. This map will be used for matrix assembly later on.
  _lower_to_upper.clear();

  // Clear _upper_to_lower. This map will be used for element ghosting
  // on distributed meshes, communication send_list construction in
  // parallel, and sparsity calculations
  _upper_to_lower.clear();

  // We do an N^2 search to find elements with matching vertex averages. This could be optimized,
  // e.g. by first sorting the vertex averages based on their (x,y,z) location.
  {
    for (const auto & [lower_key, lower_val] : lower)
      {
        // find closest vertex average in upper
        Real min_distance = std::numeric_limits<Real>::max();

        for (const auto & [upper_key, upper_val] : upper)
          {
            Real distance = (upper_val - lower_val).norm();
            if (distance < min_distance)
              {
                min_distance = distance;
                _lower_to_upper[lower_key] = upper_key.first;
              }
          }

        // For pairs with local elements, we should have found a
        // matching pair by now.
        const Elem * elem     = lower_key.first;
        const Elem * neighbor = _lower_to_upper[lower_key];
        if (min_distance < TOLERANCE)
          {
            // fill up the inverse map
            _upper_to_lower[neighbor] = elem;
          }
        else
          {
            libmesh_assert_not_equal_to(elem->processor_id(), _mesh.processor_id());
            // This must have a false positive; a remote element would
            // have been closer.
            _lower_to_upper.erase(lower_key);
          }
      }

    // Let's make sure we didn't miss any upper elements either
#ifndef NDEBUG
    for (const auto & upper_pr : upper)
      {
        const Elem * neighbor = upper_pr.first.first;
        if (neighbor->processor_id() != _mesh.processor_id())
          continue;
        ElementMap::const_iterator utl_it =
          _upper_to_lower.find(neighbor);
        libmesh_assert(utl_it != _upper_to_lower.end());
      }
#endif
  }
}

void AugmentSparsityOnInterface::operator()
  (const MeshBase::const_element_iterator & range_begin,
   const MeshBase::const_element_iterator & range_end,
   processor_id_type p,
   map_type & coupled_elements)
{
  libmesh_assert(this->_initialized);

  const CouplingMatrix * const null_mat = nullptr;

  for (const auto & elem : as_range(range_begin, range_end))
    {
      if (elem->processor_id() != p)
        coupled_elements.emplace(elem, null_mat);

      for (auto side : elem->side_index_range())
        if (elem->neighbor_ptr(side) == nullptr)
          {
            ElementSideMap::const_iterator ltu_it =
              _lower_to_upper.find(std::make_pair(elem, side));
            if (ltu_it != _lower_to_upper.end())
              {
                const Elem * neighbor = ltu_it->second;
                if (neighbor->processor_id() != p)
                  coupled_elements.emplace(neighbor, null_mat);
              }
          }

      ElementMap::const_iterator utl_it =
        _upper_to_lower.find(elem);
      if (utl_it != _upper_to_lower.end())
        {
          const Elem * neighbor = utl_it->second;
          if (neighbor->processor_id() != p)
            coupled_elements.emplace(neighbor, null_mat);
        }

    }
}
