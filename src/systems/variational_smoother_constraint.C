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

PointConstraint::PointConstraint(const Point &p) : location(p) {}

bool PointConstraint::operator<(const PointConstraint &other) const {
  return location < other.location;
}

bool PointConstraint::operator==(const PointConstraint &other) const {
  return location == other.location;
}

ConstraintVariant
PointConstraint::intersect(const ConstraintVariant &other) const {
  return std::visit(
      [&](auto &&o) -> ConstraintVariant {
        if constexpr (std::is_same_v<std::decay_t<decltype(o)>,
                                     PointConstraint>) {
          libmesh_error_msg_if(!(this->location == o.location),
                               "Points do not match.");
          return *this;
        } else {
          libmesh_error_msg_if(!o.contains_point(*this),
                               "Point is not on the constraint.");
          return *this;
        }
      },
      other);
}

LineConstraint::LineConstraint(const Point &p, const Point &d) {
  r0 = p;
  libmesh_error_msg_if(
      d.norm() < TOLERANCE,
      "Can't define a line with zero magnitude direction vector.");
  // Flip direction vector if necessary so it points in the positive x/y/z
  // direction This helps to eliminate duplicate lines
  Point canonical{0, 0, 0};
  // Choose the canonical dimension to ensure the dot product below is nonzero
  for (const auto dim_id : make_range(3))
    if (!absolute_fuzzy_equals(d(dim_id), 0.)) {
      canonical(dim_id) = 1.;
      break;
    }

  const auto dot_prod = d * canonical;
  libmesh_assert(!absolute_fuzzy_equals(dot_prod, 0.));
  dir = (dot_prod > 0) ? d.unit() : -d.unit();
}

bool LineConstraint::operator<(const LineConstraint &other) const {
  if (!(dir.absolute_fuzzy_equals(other.dir, TOLERANCE)))
    return dir < other.dir;
  return (dir * r0) < (other.dir * other.r0) - TOLERANCE;
}

bool LineConstraint::operator==(const LineConstraint &other) const {
  if (!(dir.absolute_fuzzy_equals(other.dir, TOLERANCE)))
    return false;
  return this->contains_point(other.r0);
}

bool LineConstraint::contains_point(const PointConstraint &p) const {
  // If the point lies on the line, then the vector p - r0 is parallel to the
  // line In that case, the cross product of p - r0 with the line's direction
  // will be zero.
  return dir.cross(p.location - r0).norm() < TOLERANCE;
}

bool LineConstraint::is_parallel(const LineConstraint &l) const {
  return dir.absolute_fuzzy_equals(l.dir, TOLERANCE);
}

bool LineConstraint::is_parallel(const PlaneConstraint &p) const {
  return dir * p.normal < TOLERANCE;
}

ConstraintVariant
LineConstraint::intersect(const ConstraintVariant &other) const {
  return std::visit(
      [&](auto &&o) -> ConstraintVariant {
        using T = std::decay_t<decltype(o)>;
        if constexpr (std::is_same_v<T, LineConstraint>) {
          if (*this == o)
            return *this;
          libmesh_error_msg_if(this->is_parallel(o),
                               "Lines are parallel and do not intersect.");

          // Solve for t in the equation p1 + t·d1 = p2 + s·d2
          // The shortest vector between skew lines lies along the normal vector
          // (d1 × d2). Projecting the vector (p2 - p1) onto this normal gives a
          // scalar proportional to the distance. This is equivalent to solving:
          //   ((p2 - p1) × d2) · (d1 × d2) = t · |d1 × d2|²
          //   ⇒ t = ((delta × d2) · (d1 × d2)) / |d1 × d2|²

          const Point delta = o.r0 - r0;
          const Point cross_d1_d2 = dir.cross(o.dir);
          const Real cross_dot = (delta.cross(o.dir)) * cross_d1_d2;
          const Real denom = cross_d1_d2.norm_sq();

          const Real t = cross_dot / denom;
          const Point intersection = r0 + t * dir;

          // Verify that intersection lies on both lines
          libmesh_error_msg_if(o.dir.cross(intersection - o.r0).norm() >
                                   TOLERANCE,
                               "Lines do not intersect at a single point.");

          return PointConstraint{intersection};
        } else if constexpr (std::is_same_v<T, PlaneConstraint>) {
          return o.intersect(*this);
        } else if constexpr (std::is_same_v<T, PointConstraint>) {
          libmesh_error_msg_if(!this->contains_point(o),
                               "Point is not on the line.");
          return o;
        } else
          libmesh_error_msg("Unsupported constraint type in Line::intersect.");
      },
      other);
}

PlaneConstraint::PlaneConstraint(const Point &p, const Point &n) {
  point = p;
  libmesh_error_msg_if(
      n.norm() < TOLERANCE,
      "Can't define a plane with zero magnitude direction vector.");
  // Flip normal vector if necessary so it points in the positive x/y/z
  // direction This helps to eliminate duplicate points
  Point canonical{0, 0, 0};
  // Choose the canonical dimension to ensure the dot product below is nonzero
  for (const auto dim_id : make_range(3))
    if (!absolute_fuzzy_equals(n(dim_id), 0.)) {
      canonical(dim_id) = 1.;
      break;
    }

  const auto dot_prod = n * canonical;
  libmesh_assert(!absolute_fuzzy_equals(dot_prod, 0.));
  normal = (dot_prod > 0) ? n.unit() : -n.unit();
}

bool PlaneConstraint::operator<(const PlaneConstraint &other) const {
  if (!(normal.absolute_fuzzy_equals(other.normal, TOLERANCE)))
    return normal < other.normal;
  return (normal * point) < (other.normal * other.point) - TOLERANCE;
}

bool PlaneConstraint::operator==(const PlaneConstraint &other) const {
  if (!(normal.absolute_fuzzy_equals(other.normal, TOLERANCE)))
    return false;
  return this->contains_point(other.point);
}

bool PlaneConstraint::is_parallel(const PlaneConstraint &p) const {
  return normal.absolute_fuzzy_equals(p.normal, TOLERANCE);
}

bool PlaneConstraint::is_parallel(const LineConstraint &l) const {
  return l.is_parallel(*this);
}

bool PlaneConstraint::contains_point(const PointConstraint &p) const {
  // distance between the point and the plane
  const Real dist = (p.location - point) * normal;
  return std::abs(dist) < TOLERANCE;
}

bool PlaneConstraint::contains_line(const LineConstraint &l) const {
  const bool base_on_plane = this->contains_point(PointConstraint(l.r0));
  const bool dir_orthogonal = std::abs(normal * l.dir) < TOLERANCE;
  return base_on_plane && dir_orthogonal;
}

ConstraintVariant
PlaneConstraint::intersect(const ConstraintVariant &other) const {
  return std::visit(
      [&](auto &&o) -> ConstraintVariant {
        using T = std::decay_t<decltype(o)>;
        if constexpr (std::is_same_v<T, PlaneConstraint>) {
          // If planes are identical, return one of them
          if (*this == o)
            return *this;
          libmesh_error_msg_if(this->is_parallel(o),
                               "Planes are parallel and do not intersect.");

          // Solve for a point on the intersection line of two planes.
          // Given planes:
          //   Plane 1: n1 · (x - p1) = 0
          //   Plane 2: n2 · (x - p2) = 0
          // The line of intersection has direction dir = n1 × n2.
          // To find a point on this line, we assume:
          //   x = p1 + s·n1 = p2 + t·n2
          //   ⇒ p1 - p2 = t·n2 - s·n1
          // Taking dot products with n1 and n2 leads to:
          //   [-n1·n1   n1·n2] [s] = [n1 · (p1 - p2)]
          //   [-n1·n2   n2·n2] [t]   [n2 · (p1 - p2)]

          const Point dir =
              this->normal.cross(o.normal); // direction of line of intersection
          libmesh_assert(dir.norm() > TOLERANCE);
          const Point w = this->point - o.point;

          // Dot product terms used in 2x2 system
          const Real n1_dot_n1 = normal * normal;
          const Real n1_dot_n2 = normal * o.normal;
          const Real n2_dot_n2 = o.normal * o.normal;
          const Real n1_dot_w = normal * w;
          const Real n2_dot_w = o.normal * w;

          const Real denom = -(n1_dot_n1 * n2_dot_n2 - n1_dot_n2 * n1_dot_n2);
          libmesh_assert(std::abs(denom) > TOLERANCE);

          const Real s = -(n1_dot_n2 * n2_dot_w - n2_dot_n2 * n1_dot_w) / denom;
          const Point p0 = point + s * normal;

          return LineConstraint{p0, dir};
        } else if constexpr (std::is_same_v<T, LineConstraint>) {
          if (this->contains_line(o))
            return o;
          libmesh_error_msg_if(
              this->is_parallel(o),
              "Line is parallel and does not intersect the plane.");

          // Solve for t in the parametric equation:
          //   p(t) = r0 + t·d
          // such that this point also satisfies the plane equation:
          //   n · (p(t) - p0) = 0
          // which leads to:
          //   t = (n · (p0 - r0)) / (n · d)

          const Real denom = normal * o.dir;
          libmesh_assert(std::abs(denom) > TOLERANCE);
          const Real t = (normal * (point - o.r0)) / denom;
          return PointConstraint{o.r0 + t * o.dir};
        } else if constexpr (std::is_same_v<T, PointConstraint>) {
          libmesh_error_msg_if(!this->contains_point(o),
                               "Point is not on the plane.");
          return o;
        } else
          libmesh_error_msg("Unsupported constraint type in Plane::intersect.");
      },
      other);
}

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
  const auto dim = mesh.mesh_dimension();

  // Only compute the node to elem map once
  std::unordered_map<dof_id_type, std::vector<const Elem *>> nodes_to_elem_map;
  MeshTools::build_nodes_to_elem_map(mesh, nodes_to_elem_map);

  const auto & boundary_info = mesh.get_boundary_info();

  const auto boundary_node_ids = MeshTools::find_boundary_nodes(mesh);

  // Identify/constrain subdomain boundary nodes, if requested
  std::unordered_map<dof_id_type, ConstraintVariant> subdomain_boundary_map;
  if (_preserve_subdomain_boundaries)
  {
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
          // Make sure we haven't already processed this node
          if (subdomain_boundary_map.count(node.id()))
            continue;

          // Get the relevant nodal neighbors for the subdomain constraint
          const auto side_grouped_boundary_neighbors =
              get_neighbors_for_subdomain_constraint(
                  mesh, node, sub_id1, nodes_to_elem_map);

          // Determine which constraint should be imposed
          const auto subdomain_constraint =
              determine_constraint(node, dim, side_grouped_boundary_neighbors);

          // This subdomain boundary node does not lie on an external boundary,
          // go ahead and impose constraint
          if (boundary_node_ids.find(node.id()) == boundary_node_ids.end())
            this->impose_constraint(node, subdomain_constraint);

          // This subdomain boundary node could lie on an external boundary, save it
          // for later to combine with the external boundary constraint.
          // We also save constraints for non-boundary nodes so we don't try to
          // re-constrain the node when accessed from the neighboring elem.
          // See subdomain_boundary_map.count call above.
          subdomain_boundary_map[node.id()] = subdomain_constraint;

        }//for local_node_id

      }// for side
    }// for elem
  }


  // Loop through boundary nodes and impose constraints
  for (const auto & bid : boundary_node_ids)
  {
    const auto & node = mesh.node_ref(bid);

    // Get the relevant nodal neighbors for the boundary constraint
    const auto side_grouped_boundary_neighbors =
        get_neighbors_for_boundary_constraint(mesh, node, boundary_node_ids,
                                              boundary_info, nodes_to_elem_map);

    // Determine which constraint should be imposed
    const auto boundary_constraint =
        determine_constraint(node, dim, side_grouped_boundary_neighbors);

    // Check for the case where this boundary node is also part of a subdomain id boundary
    const auto it = subdomain_boundary_map.find(bid);
    if (it != subdomain_boundary_map.end())
    {
      const auto &subdomain_constraint = it->second;
      // Combine current boundary constraint with previously determined
      // subdomain_constraint
      try
      {
        const auto combined_constraint =
            intersect_constraints(subdomain_constraint, boundary_constraint);
        this->impose_constraint(node, combined_constraint);
      }
      catch (const std::exception & e)
      {
        // This will catch cases where constraints have no intersection
        // Fall back to fixed node constraint
        this->impose_constraint(node, PointConstraint(node));
      }

    } else
      this->impose_constraint(node, boundary_constraint);

  } // end bid
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

  // We will free the dimension most parallel to line_vec to keep the
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

std::set<std::set<const Node *>>
VariationalSmootherConstraint::get_neighbors_for_subdomain_constraint(
    const MeshBase &mesh, const Node &node, const subdomain_id_type sub_id,
    const std::unordered_map<dof_id_type, std::vector<const Elem *>>
        &nodes_to_elem_map) {

  // Find all the nodal neighbors... that is the nodes directly connected
  // to this node through one edge
  std::vector<const Node *> neighbors;
  MeshTools::find_nodal_neighbors(mesh, node, nodes_to_elem_map, neighbors);

  // Each constituent set corresponds to neighbors sharing a face on the
  // subdomain boundary
  std::set<std::set<const Node *>> side_grouped_boundary_neighbors;

  for (const auto *neigh : neighbors) {
    // Determine whether the neighbor is on the subdomain boundary
    // First, find the common elements that both node and neigh belong to
    const auto &elems_containing_node = nodes_to_elem_map.at(node.id());
    const auto &elems_containing_neigh = nodes_to_elem_map.at(neigh->id());
    const Elem * common_elem = nullptr;
    for (const auto * neigh_elem : elems_containing_neigh)
    {
      if (
          (std::find(elems_containing_node.begin(), elems_containing_node.end(), neigh_elem)
          != elems_containing_node.end())
          // We should be able to find a common element on the sub_id boundary
          && (neigh_elem->subdomain_id() == sub_id)
      )
        common_elem = neigh_elem;
      else
        continue;

      // Now, determine whether node and neigh are on a side coincident
      // with the subdomain boundary
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

        if (!(node_found_on_side && neigh_found_on_side &&
              common_elem->neighbor_ptr(common_side)))
          continue;

        const auto matched_side = common_side;
        // There could be multiple matched sides, so keep this next part
        // inside the common_side loop
        //
        // Does matched_side, containing both node and neigh, lie on the
        // sub_id subdomain boundary?
        const auto matched_neighbor_sub_id =
            common_elem->neighbor_ptr(matched_side)->subdomain_id();
        const bool is_matched_side_on_subdomain_boundary =
            matched_neighbor_sub_id != sub_id;

        if (is_matched_side_on_subdomain_boundary) {
          // Store all nodes that live on this side
          const auto nodes_on_side = common_elem->nodes_on_side(common_side);
          std::set<const Node *> node_ptrs_on_side;
          for (const auto local_node_id : nodes_on_side)
            node_ptrs_on_side.insert(common_elem->node_ptr(local_node_id));
          node_ptrs_on_side.erase(node_ptrs_on_side.find(&node));
          side_grouped_boundary_neighbors.insert(node_ptrs_on_side);

          continue;
        }

      }// for common_side

    }// for neigh_elem
  }

  libmesh_assert_msg(!side_grouped_boundary_neighbors.empty(),
      "No boundary neighbors found for node " << node << " on the subdomain "
      << "boundary for subdomain " << sub_id);

  return side_grouped_boundary_neighbors;
}

std::set<std::set<const Node *>>
VariationalSmootherConstraint::get_neighbors_for_boundary_constraint(
    const MeshBase &mesh, const Node &node,
    const std::unordered_set<dof_id_type> &boundary_node_ids,
    const BoundaryInfo &boundary_info,
    const std::unordered_map<dof_id_type, std::vector<const Elem *>>
        &nodes_to_elem_map) {

  // Find all the nodal neighbors... that is the nodes directly connected
  // to this node through one edge
  std::vector<const Node *> neighbors;
  MeshTools::find_nodal_neighbors(mesh, node, nodes_to_elem_map, neighbors);

  // Each constituent set corresponds to neighbors sharing a face on the
  // boundary
  std::set<std::set<const Node *>> side_grouped_boundary_neighbors;

  for (const auto *neigh : neighbors) {
    const bool is_neighbor_boundary_node = boundary_node_ids.find(neigh->id()) != boundary_node_ids.end();
    if (!is_neighbor_boundary_node)
      continue;

    // Determine whether nodes share a common boundary id
    // First, find the common element that both node and neigh belong to
    const auto &elems_containing_node = nodes_to_elem_map.at(node.id());
    const auto &elems_containing_neigh = nodes_to_elem_map.at(neigh->id());
    const Elem *common_elem = nullptr;
    for (const auto *neigh_elem : elems_containing_neigh) {
      const bool is_neigh_common =
          std::find(elems_containing_node.begin(), elems_containing_node.end(),
                    neigh_elem) != elems_containing_node.end();
      if (!is_neigh_common)
        continue;
      common_elem = neigh_elem;
      // Keep this in the neigh_elem loop because there can be multiple common
      // elements Now, determine whether node and neigh share a common boundary
      // id
      const bool nodes_have_common_bid =
          VariationalSmootherConstraint::nodes_share_boundary_id(
              node, *neigh, *common_elem, boundary_info);
      if (nodes_have_common_bid) {
        // Find the side coinciding with the shared boundary
        for (const auto side : common_elem->side_index_range()) {
          // We only care about external boundaries here, make sure side doesn't
          // have a neighbor
          if (common_elem->neighbor_ptr(side))
            continue;

          bool node_found_on_side = false;
          bool neigh_found_on_side = false;
          const auto nodes_on_side = common_elem->nodes_on_side(side);
          for (const auto local_node_id : nodes_on_side) {
            if (common_elem->node_id(local_node_id) == node.id())
              node_found_on_side = true;
            else if (common_elem->node_id(local_node_id) == neigh->id())
              neigh_found_on_side = true;
          }
          if (!(node_found_on_side && neigh_found_on_side))
            continue;

          std::set<const Node *> node_ptrs_on_side;
          for (const auto local_node_id : nodes_on_side)
            node_ptrs_on_side.insert(common_elem->node_ptr(local_node_id));
          node_ptrs_on_side.erase(node_ptrs_on_side.find(&node));
          side_grouped_boundary_neighbors.insert(node_ptrs_on_side);
        }
        continue;
      }
    }
  }

  libmesh_assert_msg(!side_grouped_boundary_neighbors.empty(),
      "No boundary neighbors found for node " << node << " on the external boundary");

  return side_grouped_boundary_neighbors;
}

ConstraintVariant VariationalSmootherConstraint::determine_constraint(
    const Node &node, const unsigned int dim,
    const std::set<std::set<const Node *>> &side_grouped_boundary_neighbors) {
  // Determines the appropriate geometric constraint for a node based on its
  // neighbors.

  // Extract neighbors in flat vector
  std::vector<const Node *> neighbors;
  for (const auto &side : side_grouped_boundary_neighbors)
    neighbors.insert(neighbors.end(), side.begin(), side.end());
  libmesh_assert_greater_equal(neighbors.size(), 1);

  // Constrain the node to it's current location
  if (dim == 1 || neighbors.size() == 1)
    return PointConstraint{node};

  if (dim == 2 || neighbors.size() == 2) {
    // Determine whether node + all neighbors are colinear
    bool all_colinear = true;
    const Point ref_dir = (*neighbors[0] - node).unit();
    for (const auto i : make_range(size_t(1), neighbors.size())) {
      const Point delta = *(neighbors[i]) - node;
      libmesh_assert(delta.norm() >= TOLERANCE);
      const Point dir = delta.unit();
      if (!dir.relative_fuzzy_equals(ref_dir) &&
          !dir.relative_fuzzy_equals(-ref_dir)) {
        all_colinear = false;
        break;
      }
    }
    if (all_colinear)
      return LineConstraint{node, ref_dir};

    return PointConstraint{node};
  }

  // dim == 3, neighbors.size() >= 3
  std::set<PlaneConstraint> valid_planes;
  for (const auto &side_nodes : side_grouped_boundary_neighbors) {
    std::vector<const Node *> side_nodes_vec(side_nodes.begin(),
                                             side_nodes.end());
    for (const auto i : index_range(side_nodes_vec)) {
      const Point vec_i = (*side_nodes_vec[i] - node);
      for (const auto j : make_range(i)) {
        const Point vec_j = (*side_nodes_vec[j] - node);
        Point candidate_normal = vec_i.cross(vec_j);
        if (candidate_normal.norm() <= TOLERANCE)
          continue;

        const PlaneConstraint candidate_plane{node, candidate_normal};
        valid_planes.emplace(candidate_plane);
      }
    }
  }

  // Fall back to point constraint
  if (valid_planes.empty())
    return PointConstraint(node);

  // Combine all the planes together to get a common constraint
  auto it = valid_planes.begin();
  ConstraintVariant current = *it++;
  for (; it != valid_planes.end(); ++it)
  {
    try
    {
      current = intersect_constraints(current, *it);
    }
    catch (const std::exception & e)
    {
      // This will catch cases where constraints have no intersection
      // Fall back to fixed node constraint
      current = PointConstraint(node);
      break;
    }
  }

  return current;
}

// Applies the computed constraint (PointConstraint, LineConstraint, or
// PlaneConstraint) to a node.
void VariationalSmootherConstraint::impose_constraint(
    const Node &node, const ConstraintVariant &constraint) {
  if (std::holds_alternative<PointConstraint>(constraint))
    fix_node(node);
  else if (std::holds_alternative<LineConstraint>(constraint))
    constrain_node_to_line(node, std::get<LineConstraint>(constraint).dir);
  else if (std::holds_alternative<PlaneConstraint>(constraint))
    constrain_node_to_plane(node, std::get<PlaneConstraint>(constraint).normal);
  else
    libmesh_assert_msg(false, "Unknown constraint type.");
}

} // namespace libMesh
