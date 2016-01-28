#include "libmesh/parallel_ghost_sync.h"

namespace libMesh
{

SyncNodalPositions::SyncNodalPositions(MeshBase & m)
  : mesh(m)
{}



void SyncNodalPositions::gather_data (const std::vector<dof_id_type> & ids,
                                      std::vector<datum> & data)
{
  data.resize(ids.size());

  // Gather (x,y,z) data for all node IDs in the ids vector
  for (std::size_t i=0; i<ids.size(); ++i)
    {
      // Look for this node in the mesh
      Node * node = mesh.node_ptr(ids[i]);

      if (node == libmesh_nullptr)
        libmesh_error_msg("Error! Mesh returned a NULL node pointer in SyncNodalPosition::gather_data().");

      // Store this node's position in the data array.
      // This should call Point::op=
      data[i] = *node;
    } // end for
} // gather_data()



void SyncNodalPositions::act_on_data (const std::vector<dof_id_type> & ids,
                                      std::vector<datum> & data)
{
  for (std::size_t i=0; i<ids.size(); ++i)
    {

      // Get a pointer to the node whose position is to be updated.
      Node * node = mesh.node_ptr(ids[i]);

      if (node == libmesh_nullptr)
        libmesh_error_msg("Error! Mesh returned a NULL node pointer in SyncNodalPosition::act_on_data().");

      // Update this node's position.  Should call Point::op=
      *node = data[i];
    } // end for
} // act_on_data()

} // namespace
