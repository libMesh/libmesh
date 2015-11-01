// local includes
#include "augment_sparsity_on_interface.h"

// libMesh includes
#include "libmesh/linear_implicit_system.h"
#include "libmesh/elem.h"

using namespace libMesh;

AugmentSparsityOnInterface::AugmentSparsityOnInterface(EquationSystems& es,
                                                       boundary_id_type crack_boundary_lower,
                                                       boundary_id_type crack_boundary_upper)
  :
  _es(es),
  _crack_boundary_lower(crack_boundary_lower),
  _crack_boundary_upper(crack_boundary_upper)
{}

const ElementIdMap& AugmentSparsityOnInterface::get_lower_to_upper() const
{
  return _lower_to_upper;
}

void AugmentSparsityOnInterface::augment_sparsity_pattern (SparsityPattern::Graph & ,
                                                           std::vector<dof_id_type> & n_nz,
                                                           std::vector<dof_id_type> & n_oz)
{
  // get a constant reference to the mesh object
  const MeshBase& mesh = _es.get_mesh();

  // Loop over all elements (not just local elements) to make sure we find
  // "neighbor" elements on opposite sides of the crack.

  // TODO: This won't work with ParallelMesh.
  MeshBase::const_element_iterator       el     = mesh.active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_elements_end();

  // Map from (elem,side) to centroid
  std::map< std::pair<dof_id_type,unsigned char>, Point> lower_centroids;
  std::map< std::pair<dof_id_type,unsigned char>, Point> upper_centroids;

  for ( ; el != end_el; ++el)
    {
      const Elem* elem = *el;

      {
        for (unsigned char side=0; side<elem->n_sides(); side++)
          if (elem->neighbor(side) == NULL)
            {
              if( mesh.get_boundary_info().has_boundary_id
                    (elem,side, _crack_boundary_lower) )
                {
                  UniquePtr<Elem> side_elem = elem->build_side(side);

                  lower_centroids[std::make_pair(elem->id(), side)] = side_elem->centroid();
                }
            }
      }

      {
        for (unsigned char side=0; side<elem->n_sides(); side++)
          if (elem->neighbor(side) == NULL)
            {
              if( mesh.get_boundary_info().has_boundary_id
                    (elem,side, _crack_boundary_upper) )
                {
                  UniquePtr<Elem> side_elem = elem->build_side(side);

                  upper_centroids[std::make_pair(elem->id(), side)] = side_elem->centroid();
                }
            }
      }
    }
  std::size_t n_lower_centroids = lower_centroids.size();
  std::size_t n_upper_centroids = upper_centroids.size();
  libmesh_assert(n_lower_centroids == n_upper_centroids);

  // Clear _lower_to_upper. This map will be used for matrix assembly later on.
  _lower_to_upper.clear();

  // Initialize the two-way element ID map
  std::map<dof_id_type, dof_id_type> two_way_element_ID_map;

  // We do an N^2 search to find elements with matching centroids. This could be optimized,
  // e.g. by first sorting the centroids based on their (x,y,z) location.
  {
    std::map< std::pair<dof_id_type,unsigned char>, Point>::iterator it     = lower_centroids.begin();
    std::map< std::pair<dof_id_type,unsigned char>, Point>::iterator it_end = lower_centroids.end();
    for( ; it != it_end; ++it)
      {
        Point lower_centroid = it->second;

        // find centroid in upper_centroids
        std::map< std::pair<dof_id_type,unsigned char>, Point>::iterator inner_it     = upper_centroids.begin();
        std::map< std::pair<dof_id_type,unsigned char>, Point>::iterator inner_it_end = upper_centroids.end();

        Real min_distance = std::numeric_limits<Real>::max();
        for( ; inner_it != inner_it_end; ++inner_it)
          {
            Point upper_centroid = inner_it->second;

            Real distance = (upper_centroid - lower_centroid).size();
            if( distance < min_distance )
              {
                min_distance = distance;
                _lower_to_upper[it->first] = inner_it->first.first;
              }
          }

        // We should've found matching elements on either side of the crack by now
        libmesh_assert_less(min_distance, TOLERANCE);

        // fill up a "two way element ID map" which we iterate over below
        dof_id_type this_id     = it->first.first;
        dof_id_type neighbor_id = _lower_to_upper[it->first];
        two_way_element_ID_map[this_id]     = neighbor_id;
        two_way_element_ID_map[neighbor_id] = this_id;
      }
  }

  {
    const LinearImplicitSystem &system =
      _es.get_system<LinearImplicitSystem> ("Poisson");
    const DofMap &dof_map = system.get_dof_map();
    const unsigned int sys_num = system.number();
    const unsigned int var_num = system.variable_number("u");

    // Create a set that will have the DoF numbers we need to augment the sparsity pattern
    std::set<dof_id_type> local_coupled_dofs, remote_coupled_dofs;

    std::map<dof_id_type, dof_id_type>::iterator it     = two_way_element_ID_map.begin();
    std::map<dof_id_type, dof_id_type>::iterator it_end = two_way_element_ID_map.end();

    for( ; it != it_end; ++it)
      {
        local_coupled_dofs.clear();
        remote_coupled_dofs.clear();

        // get a pointer to the one of the elements on the crack
        const Elem* this_elem = mesh.elem(it->first);

        // get a pointer to the neighbor element on the other side
        // of the "crack" in the mesh
        const Elem* neighbor_elem = mesh.elem(it->second);

        // Now loop over the nodes and get the DoFs on neighbor_elem.
        // These will be used to augment the sparsity patttern.
        for (unsigned int n=0; n<neighbor_elem->n_nodes(); n++)
          {
            const dof_id_type global_dof_number = neighbor_elem->get_node(n)->dof_number(sys_num,var_num,0);

            // and finally insert it into one of the sets
            if ((global_dof_number <  dof_map.first_dof()) ||
                (global_dof_number >= dof_map.end_dof()))
              {
                remote_coupled_dofs.insert(global_dof_number);
              }
            else
              {
                local_coupled_dofs.insert(global_dof_number);
              }

          } // end loop over nodes on neighbor element

        // conservatively increase the preallocation for the implicit matrix contribution
        // to account for these dofs
        for (unsigned int nl=0; nl<this_elem->n_nodes(); nl++)
          {
            const dof_id_type
              global_dof_number = this_elem->get_node(nl)->dof_number(sys_num,var_num,0),
              n_local_dofs      = dof_map.n_local_dofs(),
              n_remote_dofs     = dof_map.n_dofs() - n_local_dofs;

            // only monkey with the sparsity pattern for local dofs!
            if ((global_dof_number >= dof_map.first_dof()) &&
                (global_dof_number  < dof_map.end_dof()))
              {
                const dof_id_type
                  dof_offset = global_dof_number - dof_map.first_dof();

                libmesh_assert_less (dof_offset, n_nz.size());
                libmesh_assert_less (dof_offset, n_oz.size());

                n_nz[dof_offset] += local_coupled_dofs.size();
                n_oz[dof_offset] += remote_coupled_dofs.size();

                // for crazy coarse problems on many processors we need to impose sane limits
                // since we shouldn't allow n_nz > n_local_dofs, for example
                n_nz[dof_offset] = std::min(n_nz[dof_offset], n_local_dofs);
                n_oz[dof_offset] = std::min(n_oz[dof_offset], n_remote_dofs);
              }
          } // end loop over nodes on target element
      }
  }

}
