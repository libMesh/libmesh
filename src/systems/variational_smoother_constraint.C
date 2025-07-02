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
#include "libmesh/boundary_info.h"

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
  const auto &mesh = _sys.get_mesh();

  // Only compute the node to elem map once
  std::unordered_map<dof_id_type, std::vector<const Elem *>> nodes_to_elem_map;
  MeshTools::build_nodes_to_elem_map(mesh, nodes_to_elem_map);

  const auto & boundary_info = mesh.get_boundary_info();

  const auto boundary_node_ids = MeshTools::find_boundary_nodes (mesh);
  for (const auto & bid : boundary_node_ids)
  {
    const auto & node = mesh.node_ref(bid);

    // Find all the nodal neighbors... that is the nodes directly connected
    // to this node through one edge
    std::vector<const Node *> neighbors;
    MeshTools::find_nodal_neighbors(mesh, node, nodes_to_elem_map, neighbors);

    // Remove any neighbors that are not boundary nodes OR boundary neighbor nodes
    // that don't share a boundary id with node
    auto remove_neighbor = [&boundary_node_ids, &node, &nodes_to_elem_map, &boundary_info]
      (const Node * neigh) -> bool
    {
      const bool is_neighbor_boundary_node = boundary_node_ids.find(neigh->id()) != boundary_node_ids.end();

      // Determine whether nodes share a common boundary id
      // First, find the common element that both node and neigh belong to
      const auto & elems_containing_node = nodes_to_elem_map[node.id()];
      const auto & elems_containing_neigh = nodes_to_elem_map[neigh->id()];
      const Elem * common_elem = nullptr;
      bool nodes_have_common_bid = false;
      for (const auto * neigh_elem : elems_containing_neigh)
        if (std::find(elems_containing_node.begin(), elems_containing_node.end(), neigh_elem) != elems_containing_node.end())
        {
          common_elem = neigh_elem;
          // Keep this in the loop because there can be multiple common elements
          // Now, determine whether node and neigh share a common boundary id
          nodes_have_common_bid = nodes_share_boundary_id(node, *neigh, *common_elem, boundary_info) || nodes_have_common_bid;
        }

      // remove if neighbor is not boundary node or nodes don't share a common bid
      return (is_neighbor_boundary_node && nodes_have_common_bid) ? false : true;
    };

    neighbors.erase(
      std::remove_if(neighbors.begin(), neighbors.end(), remove_neighbor),
      neighbors.end()
    );

    this->impose_constraints(node, neighbors);

  }// end bid

  // Constrain subdomain boundary nodes, if requested
  if (_preserve_subdomain_boundaries)
  {
    auto already_constrained_node_ids = boundary_node_ids;
    for (const auto * elem : mesh.active_element_ptr_range())
    {
      const auto & sub_id1 = elem->subdomain_id();
      for (const auto side : elem->side_index_range())
      {
        const auto * neighbor = elem->neighbor_ptr(side);
        if (neighbor == nullptr)
          continue;

        const auto & sub_id2 = neighbor->subdomain_id();
        if (sub_id1 == sub_id2)
          continue;

        // elem and neighbor are in difference subdomains, and share nodes
        // that need to be constrained
        for (const auto local_node_id : elem->nodes_on_side(side))
        {
          const auto & node = mesh.node_ref(elem->node_id(local_node_id));
          // Make sure we haven't already constrained this node
          if (
              std::find(already_constrained_node_ids.begin(),
                        already_constrained_node_ids.end(),
                        node.id()) != already_constrained_node_ids.end()
          )
            continue;

          // Find all the nodal neighbors... that is the nodes directly connected
          // to this node through one edge
          std::vector<const Node *> neighbors;
          MeshTools::find_nodal_neighbors(mesh, node, nodes_to_elem_map, neighbors);

          // Remove any neighbors that are not on the subdomain boundary
          auto remove_neighbor = [&node, &nodes_to_elem_map, &sub_id1, &sub_id2]
            (const Node * neigh) -> bool
          {
            // Determine whether the neighbor is on the subdomain boundary
            // First, find the common element that both node and neigh belong to
            const auto & elems_containing_node = nodes_to_elem_map[node.id()];
            const auto & elems_containing_neigh = nodes_to_elem_map[neigh->id()];
            const Elem * common_elem = nullptr;
            for (const auto * neigh_elem : elems_containing_neigh)
            {
              if (std::find(elems_containing_node.begin(), elems_containing_node.end(), neigh_elem) != elems_containing_node.end())
              {
                common_elem = neigh_elem;
                break;
              }
            }

            libmesh_assert(common_elem != nullptr);
            const auto common_sub_id = common_elem->subdomain_id();
            libmesh_assert(common_sub_id == sub_id1 || common_sub_id == sub_id2);

            // Define this allias for convenience
            const auto & other_sub_id = (common_sub_id == sub_id1) ? sub_id2: sub_id1;

            // Now, determine whether node and neigh are on a side coincident
            // with the interval boundary
            for (const auto common_side : common_elem->side_index_range())
            {
              bool node_found_on_side = false;
              bool neigh_found_on_side = false;
              for (const auto local_node_id : common_elem->nodes_on_side(common_side))
              {
                if (common_elem->node_id(local_node_id) == node.id())
                  node_found_on_side = true;
                else if (common_elem->node_id(local_node_id) == neigh->id())
                  neigh_found_on_side = true;
              }

              if (node_found_on_side && neigh_found_on_side)
              {
                const auto matched_side = common_side;
                // There could be multiple matched sides, so keep this next part
                // inside the loop
                //
                // Does matched_side, containing both node and neigh, lie on the
                // subdomain boundary between sub_id1 (= common_sub_id or other_sub_id)
                // and sub_id2 (= other sub_id or common_sub_id)?
                const auto matched_neighbor_sub_id = common_elem->neighbor_ptr(matched_side)->subdomain_id();
                const bool is_matched_side_on_subdomain_boundary = matched_neighbor_sub_id == other_sub_id;
                if (is_matched_side_on_subdomain_boundary)
                  return false; // Don't remove the neighbor node
              }
            }

            return true; // Remove the neighbor node
          };

          neighbors.erase(
            std::remove_if(neighbors.begin(), neighbors.end(), remove_neighbor),
            neighbors.end()
          );

          this->impose_constraints(node, neighbors);
          already_constrained_node_ids.insert(node.id());

        }//for local_node_id

      }// for side
    }// for elem
  }
}

void VariationalSmootherConstraint::impose_constraints(
    const Node &node, const std::vector<const Node *> neighbors) {
  const auto &mesh = _sys.get_mesh();
  const auto dim = mesh.mesh_dimension();

  // Determine whether the current node is colinear (2D) or coplanar 3D with
  // its boundary neighbors. Start by computing the vectors from the current
  // node to each boundary neighbor node
  std::vector<Point> dist_vecs;
  for (const auto &neighbor : neighbors)
    dist_vecs.push_back((*neighbor) - node);

  // 2D: If the current node and all (two) boundary neighbor nodes lie on the
  // same line, the magnitude of the dot product of the distance vectors will be
  // equal to the product of the magnitudes of the vectors. This is because the
  // distance vectors lie on the same line, so the cos(theta) term in the dot
  // product evaluates to -1 or 1.
  if (dim == 2) {
    // Physically, the boundary of a 2D mesh is a 1D curve. By the
    // definition of a "neighbor", it is only possible for a node
    // to have 2 neighbors on the boundary.
    libmesh_assert_equal_to(dist_vecs.size(), dim);
    const Real dot_product = dist_vecs[0] * dist_vecs[1];
    const Real norm_product = dist_vecs[0].norm() * dist_vecs[1].norm();

    // node is not colinear with its boundary neighbors and is thus immovable
    if (!relative_fuzzy_equals(std::abs(dot_product), norm_product))
      this->fix_node(node);

    else {
      // Yes, yes, we are using a function called "constrain_node_to_plane" to
      // constrain a node to a line in a 2D mesh. However, the line
      // c_x * x + c_y * y + c = 0 is equivalent to the plane
      // c_x * x + c_y * y + 0 * z + c = 0, so the same logic applies here.
      // Since all dist_vecs reside in the xy plane, and are parallel to the
      // line we are constraining to, crossing one of the dist_vecs with the
      // unit vector in the z direction should give us a vector normal to the
      // constraining line. This reference normal vector also resides in the xy
      // plane.
      //
      // TODO: what if z is not the inactive dimension in 2D?
      // Would this even happen!?!?
      const auto reference_normal = dist_vecs[0].cross(Point(0., 0., 1.));
      this->constrain_node_to_plane(node, reference_normal);
    }
  }

  // 3D: If the current node and all boundary neighbor nodes lie on the same
  // plane, all the distance vectors from the current node to the boundary nodes
  // will be orthogonal to the plane normal. We can obtain a reference normal by
  // normalizing the cross product between two of the distance vectors. If the
  // normalized cross products of all other combinations (excluding self
  // combinations) match this reference normal, then the current node is
  // coplanar with all of its boundary nodes and should be constrained to the
  // plane. If not same line/plane, then node is either part of a curved surface
  // or it is the vertex where two boundary surfaces meet. In the first case, we
  // should just fix the node. For the latter case, if the node is at the
  // intersection of 3 surfaces (i.e., the vertex of a cube), fix the node. If
  // the node is at the intersection of 2 surfaces (i.e., the edge of a cube),
  // constrain it to slide along this edge.
  else if (dim == 3) {
    // We should have at least 2 distance vectors to compute a normal with in 3D
    libmesh_assert_greater_equal(dist_vecs.size(), 2);

    // Compute the reference normal by taking the cross product of two vectors
    // in dist_vecs. We will use dist_vecs[0] as the first vector and the next
    // available vector in dist_vecs that is not (anti)parallel to dist_vecs[0].
    // Without this check we may end up with a zero vector for the reference
    // vector.
    unsigned int vec_index;
    const Point vec_0_normalized = dist_vecs[0] / dist_vecs[0].norm();
    for (const auto ii : make_range(size_t(1), dist_vecs.size())) {
      // (anti)parallel check
      const bool is_parallel = vec_0_normalized.relative_fuzzy_equals(
          dist_vecs[ii] / dist_vecs[ii].norm());
      const bool is_antiparallel = vec_0_normalized.relative_fuzzy_equals(
          -dist_vecs[ii] / dist_vecs[ii].norm());
      if (!(is_parallel || is_antiparallel)) {
        vec_index = ii;
        break;
      }
    }

    const Point reference_cross_prod = dist_vecs[0].cross(dist_vecs[vec_index]);
    const Point reference_normal =
        reference_cross_prod / reference_cross_prod.norm();

    // Does the node lie within a 2D surface? (Not on the edge)
    bool node_is_coplanar = true;

    // Does the node lie on the intersection of two 2D surfaces (i.e., a
    // line)? If so, and the node is not located at the intersection of three
    // 2D surfaces (i.e., a single point or vertex of the mesh), then some of
    // the dist_vecs will be parallel.

    // Each entry will be a vector from dist_vecs that has a corresponding
    // (anti)parallel vector, also from dist_vecs
    std::vector<Point> parallel_pairs;

    for (const auto ii : index_range(dist_vecs)) {
      const Point vec_ii_normalized = dist_vecs[ii] / dist_vecs[ii].norm();
      for (const auto jj : make_range(ii + 1, dist_vecs.size())) {
        const Point vec_jj_normalized = dist_vecs[jj] / dist_vecs[jj].norm();
        const bool is_parallel =
            vec_ii_normalized.relative_fuzzy_equals(vec_jj_normalized);
        const bool is_antiparallel =
            vec_ii_normalized.relative_fuzzy_equals(-vec_jj_normalized);

        if (is_parallel || is_antiparallel) {
          parallel_pairs.push_back(vec_ii_normalized);
          // Don't bother computing the cross product of parallel vector below,
          // it will be zero and cannot be used to define a normal vector.
          continue;
        }

        const Point cross_prod = dist_vecs[ii].cross(dist_vecs[jj]);
        const Point normal = cross_prod / cross_prod.norm();

        // node is not coplanar with its boundary neighbors
        if (!(reference_normal.relative_fuzzy_equals(normal) ||
              reference_normal.relative_fuzzy_equals(-normal)))
          node_is_coplanar = false;
      }
    }

    if (node_is_coplanar)
      this->constrain_node_to_plane(node, reference_normal);
    else if (parallel_pairs.size())
      this->constrain_node_to_line(node, parallel_pairs[0]);
    else
      this->fix_node(node);
  }

  //  1D
  else
    this->fix_node(node);
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

void VariationalSmootherConstraint::constrain_node_to_plane(const Node & node, const Point & ref_normal_vec)
{
  const auto dim = _sys.get_mesh().mesh_dimension();
  // determine equation of plane: c_x * x + c_y * y + c_z * z + c = 0
  std::vector<Real> xyz_coefs; // vector to hold c_x, c_y, c_z
  Real c = 0.;

  // We choose to constrain the dimension with the largest magnitude coefficient
  // This approach ensures the coefficients added to the constraint_row
  // (i.e., -c_xyz / c_max) have as small magnitude as possible
  unsigned int constrained_dim;
  Real max_abs_coef = 0.;
  for (const auto d : make_range(dim))
  {
    const auto coef = ref_normal_vec(d);
    xyz_coefs.push_back(coef);
    c -= coef * node(d);

    const auto coef_abs = std::abs(coef);
    if (coef_abs > max_abs_coef)
    {
      max_abs_coef = coef_abs;
      constrained_dim = d;
    }
  }

  DofConstraintRow constraint_row;
  for (const auto free_dim : make_range(dim))
  {
    if (free_dim == constrained_dim)
      continue;

    const auto free_dof_index = node.dof_number(_sys.number(), free_dim, 0);
    constraint_row[free_dof_index] =  -xyz_coefs[free_dim] / xyz_coefs[constrained_dim];
  }

  const auto inhomogeneous_part = -c / xyz_coefs[constrained_dim];
  const auto constrained_dof_index = node.dof_number(_sys.number(), constrained_dim, 0);
  _sys.get_dof_map().add_constraint_row( constrained_dof_index, constraint_row, inhomogeneous_part, true);
}

void VariationalSmootherConstraint::constrain_node_to_line(const Node & node, const Point & line_vec)
{
  const auto dim = _sys.get_mesh().mesh_dimension();

  // We will free the dimension most paralle to line_vec to keep the
  // constraint coefficients small
  const std::vector<Real> line_vec_coefs{line_vec(0), line_vec(1), line_vec(2)};
  auto it = std::max_element(line_vec_coefs.begin(), line_vec_coefs.end(),
      [](double a, double b) {
          return std::abs(a) < std::abs(b);
      });
  const unsigned int free_dim = std::distance(line_vec_coefs.begin(), it);
  const auto free_dof_index = node.dof_number(_sys.number(), free_dim, 0);

  // A line is parameterized as r(t) = node + t * line_vec, so
  // x(t) = node(x) + t * line_vec(x)
  // y(t) = node(y) + t * line_vec(y)
  // z(t) = node(z) + t * line_vec(z)
  // Let's say we leave x free. Then t = (x(t) - node(x)) / line_vec(x)
  // Then y and z can be constrained as
  // y = node(y) + line_vec_y * (x(t) - node(x)) / line_vec(x)
  //   = x(t) * line_vec(y) / line_vec(x) + (node(y) - node(x) * line_vec(y) / line_vec(x))
  // z = x(t) * line_vec(z) / line_vec(x) + (node(z) - node(x) * line_vec(z) / line_vec(x))

  libmesh_assert(!relative_fuzzy_equals(line_vec(free_dim), 0.));
  for (const auto constrained_dim : make_range(dim))
  {
    if (constrained_dim == free_dim)
      continue;

    DofConstraintRow constraint_row;
    constraint_row[free_dof_index] = line_vec(constrained_dim) / line_vec(free_dim);
    const auto inhomogeneous_part = node(constrained_dim) - node(free_dim) * line_vec(constrained_dim) / line_vec(free_dim);
    const auto constrained_dof_index = node.dof_number(_sys.number(), constrained_dim, 0);
    _sys.get_dof_map().add_constraint_row( constrained_dof_index, constraint_row, inhomogeneous_part, true);
  }
}

// Utility function to determine whether two nodes share a boundary ID.
// The motivation for this is that a sliding boundary node on a triangular
// element can have a neighbor boundary node in the same element that is not
// part of the same boundary
// Consider the below example with nodes A, C, D, F, G that comprise elements
// E1, E2, E3, with boundaries B1, B2, B3, B4. To determine the constraint
// equations for the sliding node C, the neighboring nodes A and D need to
// be identified to construct the line that C is allowed to slide along.
// Note that neighbors A and D both share the boundary B1 with C.
// Without ensuring that neighbors share the same boundary as the current
// node, a neighboring node that does not lie on the same boundary
// (i.e. F and G) might be selected to define the constraining line,
// resulting in an incorrect constraint.
// Note that if, for example, boundaries B1 and B2 were to be combined
// into boundary B12, the node F would share a boundary id with node C
// and result in an incorrect constraint. It would be useful to design
// additional checks to detect cases like this.

//         B3
//    G-----------F
//    | \       / |
// B4 |  \  E2 /  | B2
//    |   \   /   |
//    | E1 \ / E3 |
//    A-----C-----D
//         B1

bool VariationalSmootherConstraint::nodes_share_boundary_id(
  const Node & boundary_node,
  const Node & neighbor_node,
  const Elem & containing_elem,
  const BoundaryInfo & boundary_info)
{
  bool nodes_share_bid = false;

  // Node ids local to containing_elem
  const auto node_id = containing_elem.get_node_index(&boundary_node);
  const auto neighbor_id = containing_elem.get_node_index(&neighbor_node);

  for (const auto side_id : containing_elem.side_index_range())
  {
    // We don't care about this side if it doesn't contain our boundary and neighbor nodes
    if (!(containing_elem.is_node_on_side(node_id, side_id) && containing_elem.is_node_on_side(neighbor_id, side_id)))
      continue;

    // If the current side, containing boundary_node and neighbor_node, lies on a boundary,
    // we can say that boundary_node and neighbor_node have a common boundary id.
    std::vector<boundary_id_type> boundary_ids;
    boundary_info.boundary_ids(&containing_elem, side_id, boundary_ids);
    if (boundary_ids.size())
    {
      nodes_share_bid = true;
      break;
    }
  }
  return nodes_share_bid;
}

} // namespace libMesh
