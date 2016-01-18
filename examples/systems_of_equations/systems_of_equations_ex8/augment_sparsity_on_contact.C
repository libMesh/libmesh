// local includes
#include "augment_sparsity_on_contact.h"

// libMesh includes
#include "libmesh/linear_implicit_system.h"
#include LIBMESH_INCLUDE_UNORDERED_SET
#include "libmesh/elem.h"

using namespace libMesh;

AugmentSparsityOnContact::AugmentSparsityOnContact(System & sys) :
  _sys(sys)
{}

void AugmentSparsityOnContact::clear_contact_node_map()
{
  _contact_node_map.clear();
}

void AugmentSparsityOnContact::add_contact_node(dof_id_type lower_node_id,
                                                dof_id_type upper_node_id)
{
  _contact_node_map[lower_node_id] = upper_node_id;
}

void AugmentSparsityOnContact::augment_sparsity_pattern(SparsityPattern::Graph &,
                                                        std::vector<dof_id_type> & n_nz,
                                                        std::vector<dof_id_type> & n_oz)
{
  // get a constant reference to the mesh object
  const MeshBase & mesh = _sys.get_mesh();

  LIBMESH_BEST_UNORDERED_MAP<dof_id_type, dof_id_type>::iterator it     = _contact_node_map.begin();
  LIBMESH_BEST_UNORDERED_MAP<dof_id_type, dof_id_type>::iterator it_end = _contact_node_map.end();

  for ( ; it != it_end; ++it)
    {
      const Node & this_node = mesh.node(it->first);
      const Node & other_node = mesh.node(it->second);

      set_sparsity_values(this_node, other_node, n_nz, n_oz);
      set_sparsity_values(other_node, this_node, n_nz, n_oz);
    }
}

void AugmentSparsityOnContact::set_sparsity_values(const Node & this_node,
                                                   const Node & other_node,
                                                   std::vector<dof_id_type> & n_nz,
                                                   std::vector<dof_id_type> & n_oz)
{
  const DofMap & dof_map = _sys.get_dof_map();

  for (unsigned int var_i=0; var_i<3; var_i++)
    {
      dof_id_type dof_index_on_this_node =
        this_node.dof_number(_sys.number(), var_i, 0);

      unsigned int n_local_coupled_dofs = 0;
      unsigned int n_remote_coupled_dofs = 0;

      for (unsigned int var_j=0; var_j<3; var_j++)
        {
          dof_id_type dof_index_on_other_node =
            other_node.dof_number(_sys.number(), var_j, 0);

          if ((dof_index_on_other_node >= dof_map.first_dof()) &&
              (dof_index_on_other_node  < dof_map.end_dof()))
            n_local_coupled_dofs++;
          else
            n_remote_coupled_dofs++;
        }

      // only monkey with the sparsity pattern for local dofs!
      if ((dof_index_on_this_node >= dof_map.first_dof()) &&
          (dof_index_on_this_node  < dof_map.end_dof()))
        {
          const unsigned int
            dof_offset = dof_index_on_this_node - dof_map.first_dof();

          libmesh_assert_less (dof_offset, n_nz.size());
          libmesh_assert_less (dof_offset, n_oz.size());

          n_nz[dof_offset] += n_local_coupled_dofs;
          n_oz[dof_offset] += n_remote_coupled_dofs;
        }
    }
}
