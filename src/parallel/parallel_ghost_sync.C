#include "libmesh/parallel_ghost_sync.h"

namespace libMesh
{

SyncNodalPositions::SyncNodalPositions(MeshBase & m)
  : mesh(m)
{}



void SyncNodalPositions::gather_data (const std::vector<dof_id_type> & ids,
                                      std::vector<datum> & data) const
{
  data.resize(ids.size());

  // Gather (x,y,z) data for all node IDs in the ids vector
  for (std::size_t i=0; i<ids.size(); ++i)
    {
      // Look for this node in the mesh
      const Point & pt = mesh.point(ids[i]);

      // Store this node's position in the data array.
      // This should call Point::op=
      data[i] = pt;
    } // end for
} // gather_data()



void SyncNodalPositions::act_on_data (const std::vector<dof_id_type> & ids,
                                      std::vector<datum> & data) const
{
  for (std::size_t i=0; i<ids.size(); ++i)
    {

      // Get a pointer to the node whose position is to be updated.
      Node & node = mesh.node_ref(ids[i]);

      // Update this node's position.  Should call Point::op=
      node = data[i];
    } // end for
} // act_on_data()

} // namespace
