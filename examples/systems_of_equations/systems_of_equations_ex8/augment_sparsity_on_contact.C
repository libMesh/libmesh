// local includes
#include "augment_sparsity_on_contact.h"

// libMesh includes
#include "libmesh/linear_implicit_system.h"
#include LIBMESH_INCLUDE_UNORDERED_SET

using namespace libMesh;

AugmentSparsityOnContact::AugmentSparsityOnContact(
  System& sys)
  :
  _sys(sys)
{}

void AugmentSparsityOnContact::clear_contact_element_map()
{
  _contact_element_map.clear();
}

void AugmentSparsityOnContact::add_contact_element(
  dof_id_type element_id,
  dof_id_type other_element_id)
{
  LIBMESH_BEST_UNORDERED_MAP<dof_id_type, std::set<dof_id_type> >::iterator it =
    _contact_element_map.find(element_id);

  if(it == _contact_element_map.end())
  {
    std::set<dof_id_type> new_set;
    new_set.insert(other_element_id);

    _contact_element_map[element_id] = new_set;
  }
  else
  {
    std::set<dof_id_type> existing_set = it->second;
    existing_set.insert(other_element_id);

    _contact_element_map[element_id] = existing_set;
  }
}

void AugmentSparsityOnContact::augment_sparsity_pattern(
  SparsityPattern::Graph & ,
  std::vector<dof_id_type> & n_nz,
  std::vector<dof_id_type> & n_oz)
{
  // get a constant reference to the mesh object
  const MeshBase& mesh = _sys.get_mesh();

  const DofMap& dof_map = _sys.get_dof_map();

  // Create a set that will have the DoF numbers we need to augment the sparsity pattern
  LIBMESH_BEST_UNORDERED_SET<dof_id_type> local_coupled_dofs, remote_coupled_dofs;

  LIBMESH_BEST_UNORDERED_MAP<dof_id_type, std::set<dof_id_type> >::iterator it     = _contact_element_map.begin();
  LIBMESH_BEST_UNORDERED_MAP<dof_id_type, std::set<dof_id_type> >::iterator it_end = _contact_element_map.end();

  for( ; it != it_end; ++it)
  {
    const Elem* this_elem = mesh.elem(it->first);
    std::set<dof_id_type> neighbor_elem_ids = it->second;

    local_coupled_dofs.clear();
    remote_coupled_dofs.clear();

    std::set<dof_id_type>::iterator set_it = neighbor_elem_ids.begin();
    std::set<dof_id_type>::iterator set_it_end = neighbor_elem_ids.end();
    for( ; set_it != set_it_end; ++set_it)
    {
      const Elem* neighbor_elem = mesh.elem( *set_it );

      // Now loop over the nodes and get the DoFs on neighbor_elem.
      // These will be used to augment the sparsity pattern.
      // Note, we should also account for element-based dofs,
      // but for now we're just being lazy and considering only
      // node-based dofs.
      for(unsigned int var_num=0; var_num<_sys.n_vars(); var_num++)
      {
        for (unsigned int n=0; n<neighbor_elem->n_nodes(); n++)
        {
          unsigned int var_dofs_on_node =
            neighbor_elem->get_node(n)->n_dofs(_sys.number(), var_num);

          if(var_dofs_on_node > 0)
          {
            const dof_id_type global_dof_number =
              neighbor_elem->get_node(n)->dof_number(_sys.number(),var_num,0);

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
          }
        } // end loop over nodes on neighbor element
      }
    }

    // Conservatively increase the preallocation for the implicit matrix contribution
    // to account for remote_coupled_dofs and local_coupled_dofs.
    for(unsigned int var_num=0; var_num<_sys.n_vars(); var_num++)
    {
      for (unsigned int nl=0; nl<this_elem->n_nodes(); nl++)
      {
        const dof_id_type
          global_dof_number = this_elem->get_node(nl)->dof_number(_sys.number(),var_num,0),
          n_local_dofs      = dof_map.n_local_dofs(),
          n_remote_dofs     = dof_map.n_dofs() - n_local_dofs;

        // Only monkey with the sparsity pattern for local dofs!
        // The other dofs will be handled on other processors.
        if ((global_dof_number >= dof_map.first_dof()) &&
            (global_dof_number  < dof_map.end_dof()))
        {
          const dof_id_type
            dof_offset = global_dof_number - dof_map.first_dof();

          libmesh_assert_less (dof_offset, n_nz.size());
          libmesh_assert_less (dof_offset, n_oz.size());

          // TODO: This factor of 5 should not be needed. We must be making
          // a mistake somewhere with keeping track of the dof coupling.
          n_nz[dof_offset] += 5*local_coupled_dofs.size();
          n_oz[dof_offset] += 5*remote_coupled_dofs.size();

          // for crazy coarse problems on many processors we need to impose sane limits
          // since we shouldn't allow n_nz > n_local_dofs, for example
          n_nz[dof_offset] = std::min(n_nz[dof_offset], n_local_dofs);
          n_oz[dof_offset] = std::min(n_oz[dof_offset], n_remote_dofs);
        }
      }
    }

  }

}
