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
  const auto & mesh = _sys.get_mesh();
  const auto dim = mesh.mesh_dimension();

  // Only compute the node to elem map once
  std::unordered_map<dof_id_type, std::vector<const Elem *>> nodes_to_elem_map;
  MeshTools::build_nodes_to_elem_map(mesh, nodes_to_elem_map);

  const auto boundary_node_ids = MeshTools::find_boundary_nodes (mesh);
  for (const auto & bid : boundary_node_ids)
  {
    const auto & node = mesh.node_ref(bid);
    // Find all the nodal neighbors... that is the nodes directly connected
    // to this node through one edge
    std::vector<const Node *> neighbors;
    MeshTools::find_nodal_neighbors(mesh, node, nodes_to_elem_map, neighbors);

    // Remove any neighbors that are not boundary nodes
    neighbors.erase(
      std::remove_if(neighbors.begin(), neighbors.end(),
        [&boundary_node_ids](const Node * neigh) {
          return boundary_node_ids.find(neigh->id()) == boundary_node_ids.end();
        }
      ),
      neighbors.end()
    );

    // Determine whether the current node is colinear (coplanar) with its boundary
    // neighbors Start by computing the vectors from the current node to each boundary
    // neighbor node
    //std::cout << "Node " << node.id() << ":" << std::endl;
    std::vector<Point> dist_vecs;
    for (const auto & neighbor : neighbors)
    {
      dist_vecs.push_back((*neighbor) - node);
      //std::cout << "  Neighbor " << neighbor->id() << ", dist_vec = " << dist_vecs.back() << std::endl;
    }

    // 2D: If the current node and all (two) boundary neighbor nodes lie on the same line,
    // the magnitude of the dot product of the distance vectors will be equal
    // to the product of the magnitudes of the vectors. This is because the distance
    // vectors lie on the same line, so the cos(theta) term in the dot product
    // evaluates to -1 or 1.
    if (dim == 2)
    {
      // Physically, the boundary of a 2D mesh is a 1D curve. By the
      // definition of a "neighbor", it is only possible for a node
      // to have 2 neighbors on the boundary.
      libmesh_assert_equal_to(dist_vecs.size(), dim);
      const Real dot_product = dist_vecs[0] * dist_vecs[1];
      const Real norm_product = dist_vecs[0].norm() * dist_vecs[1].norm();

      // node is not colinear with its boundary neighbors and is thus immovable
      if (!relative_fuzzy_equals(std::abs(dot_product), norm_product))
      {
        this->fix_node(node);
        continue;
      }

      // else, determine equation of line: c_x * x + c_y * y + c = 0
      Real c_x, c_y, c;
      const auto & vec = dist_vecs[0];
      if (vec(0) == 0) // vertical line
      {
        c_y = 0.;
        c_x = 1.;
        c = -node(0);
      }
      else // not a vertical line
      {
        c_y = 1.;
        c_x = -vec(1) / vec(0); // c_x = -m from y = mx + b
        c = -c_x * node(0) - node(1);
      }

      const std::vector<Real> xy_coefs{c_x, c_y};

      // Constrain the dimension with the largest coefficient
      const unsigned int constrained_dim = (std::abs(c_x) > std::abs(c_y)) ? 0 : 1;
      const unsigned int free_dim = (std::abs(c_x) > std::abs(c_y)) ? 1 : 0;

      const auto constrained_dof_index = node.dof_number(_sys.number(), constrained_dim, 0);
      const auto free_dof_index = node.dof_number(_sys.number(), free_dim, 0);
      DofConstraintRow constraint_row;
      constraint_row[free_dof_index] =  -xy_coefs[free_dim] / xy_coefs[constrained_dim];
      const auto constrained_value = -c / xy_coefs[constrained_dim];
      _sys.get_dof_map().add_constraint_row( constrained_dof_index, constraint_row, constrained_value, true);
    }

    // 3D: If the current node and all boundary neighbor nodes lie on the same plane,
    // all the distance vectors from the current node to the boundary nodes will be
    // orthogonal to the plane normal. We can obtain a reference normal by normalizing
    // the cross product between two of the distance vectors. If the normalized cross
    // products of all other combinations (excluding self combinations) match this
    // reference normal, then the current node is coplanar with all of its boundary nodes.
    else if (dim == 3)
    {
      // We should have at least 2 distance vectors to compute a normal with in 3D
      libmesh_assert_greater_equal(dist_vecs.size(), 2);

      // Compute the reference normal by taking the cross product of two vectors in
      // dist_vecs. We will use dist_vecs[0] as the first vector and the next available
      // vector in dist_vecs that is not (anti)parallel to dist_vecs[0]. Without this
      // check we may end up with a zero vector for the reference vector.
      unsigned int vec_index;
      const Point vec_0_normalized = dist_vecs[0] / dist_vecs[0].norm();
      for (const auto ii : make_range(size_t(1), dist_vecs.size()))
      {
        // (anti)parallel check
        const bool is_parallel =
            vec_0_normalized.relative_fuzzy_equals(dist_vecs[ii] / dist_vecs[ii].norm());
        const bool is_antiparallel =
            vec_0_normalized.relative_fuzzy_equals(-dist_vecs[ii] / dist_vecs[ii].norm());
        if (!(is_parallel || is_antiparallel))
        {
          vec_index = ii;
          break;
        }
      }

      const Point reference_cross_prod = dist_vecs[0].cross(dist_vecs[vec_index]);
      const Point reference_normal = reference_cross_prod / reference_cross_prod.norm();
      
      bool node_is_coplanar = true;
      for (const auto ii : index_range(dist_vecs))
      {
        const Point vec_ii_normalized = dist_vecs[ii] / dist_vecs[ii].norm();
        for (const auto jj : make_range(ii + 1, dist_vecs.size()))
        {
          // No need to compute the cross product for this case, as it is by
          // definition equal to the reference normal computed above.
          // Also check for dist_vecs that are (anti)parallel, and skip,
          // as their cross product will be zero.
          const Point vec_jj_normalized = dist_vecs[jj] / dist_vecs[jj].norm();
          const bool is_parallel =
              vec_ii_normalized.relative_fuzzy_equals(vec_jj_normalized);
          const bool is_antiparallel = vec_ii_normalized.relative_fuzzy_equals(
              -vec_jj_normalized);

          if ((ii == 0 and jj == vec_index) || (is_parallel || is_antiparallel))
            continue;

          const Point cross_prod = dist_vecs[ii].cross(dist_vecs[jj]);
          const Point normal = cross_prod / cross_prod.norm();

          // node is not coplanar with its boundary neighbors and is thus immovable
          if (!(reference_normal.relative_fuzzy_equals(normal) ||
                reference_normal.relative_fuzzy_equals(-normal)))
          {
            node_is_coplanar = false;
            break;
          }
        }
      }
      
      // TODO: Need to add check for sliding edge node

      if (!node_is_coplanar)
      {
        this->fix_node(node);
        continue;
      }

      // else, determine equation of plane: c_x * x + c_y * y + c_z * z + c = 0
      const Real c_x = reference_normal(0);
      const Real c_y = reference_normal(1);
      const Real c_z = reference_normal(2);
      const Real c = -(c_x * node(0) + c_y * node(1) + c_z * node(2));

      const std::vector<Real> xyz_coefs{c_x, c_y, c_z};

      // Find the dimension with the largest nonzero magnitude coefficient
      auto it = std::max_element(xyz_coefs.begin(), xyz_coefs.end(),
          [](double a, double b) {
              return std::abs(a) < std::abs(b);
          });
      const unsigned int constrained_dim = std::distance(xyz_coefs.begin(), it);

      //std::cout << "Node " << node.id() << " (" << node(0) << ", " << node(1) << ", " << node(2)
      //          << ") lies on the plane " << std::endl
      //          << c_x << " * x + " << c_y << " * y + " << c_z << " * z + " << c << " = 0" << std::endl;


      //const std::vector<std::string> dim_names{"x", "y", "z"};
      //std::cout << "Constraining " << dim_names[constrained_dim] << " = ";

      DofConstraintRow constraint_row;
      auto constrained_value = -c / xyz_coefs[constrained_dim];
      for (const auto free_dim : index_range(xyz_coefs))
      {
        if (free_dim == constrained_dim)
          continue;
        const auto free_dof_index = node.dof_number(_sys.number(), free_dim, 0);
        constraint_row[free_dof_index] =  -xyz_coefs[free_dim] / xyz_coefs[constrained_dim];
        //std::cout << constraint_row[free_dof_index] << " * " << dim_names[free_dim] << " + ";
      }

      //std::cout << constrained_value << std::endl;
      const auto constrained_dof_index = node.dof_number(_sys.number(), constrained_dim, 0);
      _sys.get_dof_map().add_constraint_row( constrained_dof_index, constraint_row, constrained_value, true);

    }

    // if not same line/plane, then node is either part of a curved surface or
    // it is the vertex where two boundary surfaces meet. In the first case,
    // we should just fix the node. For the latter case:
    //   - 2D: just fix the node
    //   - 3D: If the node is at the intersection of 3 surfaces (i.e., the
    //   vertex of a cube), fix the node. If the node is at the intersection
    //   of 2 surfaces (i.e., the edge of a cube), constrain it to slide along
    //   this edge.

    //  1D
    else
      this->fix_node(node);

  }// end bid

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
