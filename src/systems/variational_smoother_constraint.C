// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

// Local Includes
#include "libmesh/variational_smoother_constraint.h"
#include "libmesh/mesh_tools.h"

namespace libMesh
{

VariationalSmootherConstraint::VariationalSmootherConstraint(System & sys, const bool & preserve_subdomain_boundaries)
  :
    Constraint(),
    _sys(sys),
    _preserve_subdomain_boundaries(preserve_subdomain_boundaries)
  {}

VariationalSmootherConstraint::~VariationalSmootherConstraint() = default;

void VariationalSmootherConstraint::constrain()
{
  auto & mesh = _sys.get_mesh();

  // Only compute the node to elem map once
  std::unordered_map<dof_id_type, std::vector<const Elem *>> nodes_to_elem_map;
  MeshTools::build_nodes_to_elem_map(mesh, nodes_to_elem_map);

  const auto boundary_node_ids = MeshTools::find_boundary_nodes (mesh);
  for (const auto & bid : boundary_node_ids)
  {
    //if (bid == 2)
    //  continue;
    const auto & node = mesh.node_ref(bid);
    // Find all the nodal neighbors... that is the nodes directly connected
    // to this node through one edge
    std::vector<const Node *> neighbors;
    MeshTools::find_nodal_neighbors(mesh, node, nodes_to_elem_map, neighbors);

    // Remove any neighbors that are not boundary nodes
    neighbors.erase(
      std::remove_if(neighbors.begin(), neighbors.end(),
        [&boundary_node_ids](const Node * node) {
          return boundary_node_ids.find(node->id()) == boundary_node_ids.end();
        }
      ),
      neighbors.end()
    );

    // if 2D, determine if node and neighbors lie on the same line
    // if 3D, determine if node and neighbors lie on the same plane
    // if not same line/plane, then node is either part of a curved surface or
    // it is the vertex where two boundary surfaces meet. In the first case,
    // we should just fix the node. For the latter case:
    //   - 2D: just fix the node
    //   - 3D: If the node is at the intersection of 3 surfaces (i.e., the
    //   vertex of a cube), fix the node. If the node is at the intersection
    //   of 2 surfaces (i.e., the edge of a cube), constrain it to slide along
    //   this edge.

    //   But for now, just fix all the boundary nodes to not move
    this->fix_node(node);
  }

  // Constrain subdomain boundary nodes, if requested
  if (_preserve_subdomain_boundaries)
  {
    auto already_constrained_node_ids = boundary_node_ids;
    for (const auto * elem : mesh.active_element_ptr_range())
    {
      const auto & subdomain_id = elem->subdomain_id();
      for (const auto side : elem->side_index_range())
      {
        const auto * neighbor = elem->neighbor_ptr(side);
        if (neighbor == nullptr)
          continue;

        const auto & neighbor_subdomain_id = neighbor->subdomain_id();
        if (subdomain_id != neighbor_subdomain_id)
        {
          for (const auto local_node_id : elem->nodes_on_side(side))
          {
            const auto & node = mesh.node_ref(elem->node_id(local_node_id));
            // Make sure we haven't already constrained this node
            if (
                std::find(already_constrained_node_ids.begin(),
                          already_constrained_node_ids.end(),
                          node.id()) == already_constrained_node_ids.end()
            )
            {
              this->fix_node(node);
              already_constrained_node_ids.insert(node.id());
            }
          }//for local_node_id
        }
      }// for side
    }// for elem
  }

}

void VariationalSmootherConstraint::fix_node(const Node & node)
{
  for (const auto d : make_range(_sys.get_mesh().mesh_dimension()))
  {
    const auto constrained_dof_index = node.dof_number(_sys.number(), d, 0);
    DofConstraintRow constraint_row;
    // Leave the constraint row as all zeros so this dof is independent from other dofs
    const auto constrained_value = node(d);
    // Simply constrain this dof to retain it's current value
    _sys.get_dof_map().add_constraint_row( constrained_dof_index, constraint_row, constrained_value, true);
  }
}

} // namespace libMesh
