// local includes
#include "augment_sparsity_on_interface.h"

// libMesh includes
#include "libmesh/elem.h"

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

  // Map from (elem, side) to centroid
  std::map<std::pair<const Elem *, unsigned char>, Point> lower_centroids;
  std::map<std::pair<const Elem *, unsigned char>, Point> upper_centroids;

  for (const auto & elem : _mesh.active_element_ptr_range())
    for (auto side : elem->side_index_range())
      if (elem->neighbor_ptr(side) == libmesh_nullptr)
        {
          if (_mesh.get_boundary_info().has_boundary_id(elem, side, _crack_boundary_lower))
            {
              std::unique_ptr<const Elem> side_elem = elem->build_side_ptr(side);

              lower_centroids[std::make_pair(elem, side)] = side_elem->centroid();
            }

          if (_mesh.get_boundary_info().has_boundary_id(elem, side, _crack_boundary_upper))
            {
              std::unique_ptr<const Elem> side_elem = elem->build_side_ptr(side);

              upper_centroids[std::make_pair(elem, side)] = side_elem->centroid();
            }
        }

  // If we're doing a reinit on a distributed mesh then we may not see
  // all the centroids, or even a matching number of centroids.
  // std::size_t n_lower_centroids = lower_centroids.size();
  // std::size_t n_upper_centroids = upper_centroids.size();
  // libmesh_assert(n_lower_centroids == n_upper_centroids);

  // Clear _lower_to_upper. This map will be used for matrix assembly later on.
  _lower_to_upper.clear();

  // Clear _upper_to_lower. This map will be used for element ghosting
  // on distributed meshes, communication send_list construction in
  // parallel, and sparsity calculations
  _upper_to_lower.clear();

  // We do an N^2 search to find elements with matching centroids. This could be optimized,
  // e.g. by first sorting the centroids based on their (x,y,z) location.
  {
    std::map<std::pair<const Elem *, unsigned char>, Point>::iterator it     = lower_centroids.begin();
    std::map<std::pair<const Elem *, unsigned char>, Point>::iterator it_end = lower_centroids.end();
    for ( ; it != it_end; ++it)
      {
        Point lower_centroid = it->second;

        // find closest centroid in upper_centroids
        Real min_distance = std::numeric_limits<Real>::max();

        std::map<std::pair<const Elem *, unsigned char>, Point>::iterator inner_it     = upper_centroids.begin();
        std::map<std::pair<const Elem *, unsigned char>, Point>::iterator inner_it_end = upper_centroids.end();

        for ( ; inner_it != inner_it_end; ++inner_it)
          {
            Point upper_centroid = inner_it->second;

            Real distance = (upper_centroid - lower_centroid).norm();
            if (distance < min_distance)
              {
                min_distance = distance;
                _lower_to_upper[it->first] = inner_it->first.first;
              }
          }

        // For pairs with local elements, we should have found a
        // matching pair by now.
        const Elem * elem     = it->first.first;
        const Elem * neighbor = _lower_to_upper[it->first];
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
            _lower_to_upper.erase(it->first);
          }
      }

    // Let's make sure we didn't miss any upper elements either
#ifndef NDEBUG
    std::map<std::pair<const Elem *, unsigned char>, Point>::iterator inner_it     = upper_centroids.begin();
    std::map<std::pair<const Elem *, unsigned char>, Point>::iterator inner_it_end = upper_centroids.end();

    for ( ; inner_it != inner_it_end; ++inner_it)
      {
        const Elem * neighbor = inner_it->first.first;
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

  const CouplingMatrix * const null_mat = libmesh_nullptr;

  for (MeshBase::const_element_iterator elem_it = range_begin;
       elem_it != range_end; ++elem_it)
    {
      const Elem * const elem = *elem_it;

      if (elem->processor_id() != p)
        coupled_elements.insert (std::make_pair(elem, null_mat));

      for (auto side : elem->side_index_range())
        if (elem->neighbor_ptr(side) == libmesh_nullptr)
          {
            ElementSideMap::const_iterator ltu_it =
              _lower_to_upper.find(std::make_pair(elem, side));
            if (ltu_it != _lower_to_upper.end())
              {
                const Elem * neighbor = ltu_it->second;
                if (neighbor->processor_id() != p)
                  coupled_elements.insert (std::make_pair(neighbor, null_mat));
              }
          }

      ElementMap::const_iterator utl_it =
        _upper_to_lower.find(elem);
      if (utl_it != _upper_to_lower.end())
        {
          const Elem * neighbor = utl_it->second;
          if (neighbor->processor_id() != p)
            coupled_elements.insert (std::make_pair(neighbor, null_mat));
        }

    }
}
