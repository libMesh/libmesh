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

#ifndef LIBMESH_VARIATIONAL_SMOOTHER_CONSTRAINT_H
#define LIBMESH_VARIATIONAL_SMOOTHER_CONSTRAINT_H

// Local Includes
#include "libmesh/system.h"
#include "libmesh/dof_map.h"

// C++ includes
#include <variant>

namespace libMesh
{

// Forward declarations
class PointConstraint;
class LineConstraint;
class PlaneConstraint;
class InvalidConstraint;

/**
 * Type used to store a constraint that may be a PlaneConstraint,
 * LineConstraint, or PointConstraint. std::variant is an alternative to using
 * the classic polymorphic approach where these constraints inherit from an
 * common base class.
 */
using ConstraintVariant = std::variant<PointConstraint, LineConstraint,
                                       PlaneConstraint, InvalidConstraint>;

/**
 * Represents a fixed point constraint.
 */
class PointConstraint
{

public:
  PointConstraint() = default;

  /**
   * Constructor
   * @param point The point defining the constraint.
   * @param tol The tolerance to use for numerical comparisons.
   */
  PointConstraint(const Point &point, const Real &tol = TOLERANCE * TOLERANCE);

  /**
   * Comparison operator for ordering PointConstraint objects.
   * A PointConstraint is considered less than another if its location
   * is lexicographically less than the other's location.
   * @param other The PointConstraint to compare with.
   * @param tol The tolerance to use for numerical comparisons.
   * @return True if this PointConstraint is less than the other.
   */
  bool operator<(const PointConstraint &other) const;

  /**
   * Equality operator.
   * @param other The PointConstraint to compare with.
   * @return True if both PointConstraints have the same location.
   */
  bool operator==(const PointConstraint &other) const;

  /**
   * Query whether a point lies on another point.
   * @param p The point in question
   * @return bool indicating whether p lies on this point.
   */
  bool contains_point(const PointConstraint &p) const { return *this == p; }

  /**
   * Computes the intersection of this point with another constraint.
   * Handles intersection with PointConstraint, LineConstraint, or
   * PlaneConstraint.
   * @param other The constraint to intersect with.
   * @return The most specific ConstraintVariant that satisfies both
   * constraints. constraints. If no intersection exists, return an
   * InvalidConstraint.
   */
  ConstraintVariant intersect(const ConstraintVariant &other) const;

  /**
   * Const getter for the _point attribute
   */
  const Point &point() const { return _point; }

  /**
   * Const getter for the _tol attribute
   */
  const Real &tol() const { return _tol; }

private:
  // Life is easier if we don't make this const
  /**
   * Location of constraint
   */
  Point _point;

  /**
   * Tolerance to use for numerical comparisons
   */
  Real _tol;
};

/**
 * Represents a line constraint defined by a base point and direction vector.
 */
class LineConstraint
{
public:
  LineConstraint() = default;

  /**
   * Constructor
   * @param point A point on the constraining line.
   * @param direction the direction of the constraining line.
   * @param tol The tolerance to use for numerical comparisons.
   */
  LineConstraint(const Point &point, const Point &direction,
                 const Real &tol = TOLERANCE * TOLERANCE);

  /**
   * Comparison operator for ordering LineConstraint objects.
   * The comparison is primarily based on the direction vector. If the direction
   * vectors are equal (within tolerance), the tie is broken using the dot
   * product of the direction with the base point.
   * @param other The LineConstraint to compare with.
   * @return True if this LineConstraint is less than the other.
   */
  bool operator<(const LineConstraint &other) const;

  /**
   * Equality operator.
   * @param other The LineConstraint to compare with.
   * @return True if both LineConstraints represent the same line.
   */
  bool operator==(const LineConstraint &other) const;

  /**
   * Query whether a point lies on the line.
   * @param p The point in question
   * @return bool indicating whether p lies on the line.
   */
  bool contains_point(const PointConstraint &p) const;

  /**
   * Query whether a line is parallel to this line
   * @param l The line in question
   * @return bool indicating whether l is parallel to this line.
   */
  bool is_parallel(const LineConstraint &l) const;

  /**
   * Query whether a plane is parallel to this line
   * @param p The plane in question
   * @return bool indicating whether p is parallel to this line.
   */
  bool is_parallel(const PlaneConstraint &p) const;

  /**
   * Computes the intersection of this line with another constraint.
   * Handles intersection with LineConstraint, PlaneConstraint, or
   * PointConstraint.
   * @param other The constraint to intersect with.
   * @return The most specific ConstraintVariant that satisfies both
   * constraints. constraints. If no intersection exists, return an
   * InvalidConstraint.
   */
  ConstraintVariant intersect(const ConstraintVariant &other) const;

  /**
   * Const getter for the _point attribute
   */
  const Point &point() const { return _point; }

  /**
   * Const getter for the _direction attribute
   */
  const Point &direction() const { return _direction; }

  /**
   * Const getter for the _tol attribute
   */
  const Real &tol() const { return _tol; }

private:
  // Life is easier if we don't make these const
  /**
   * A point on the constraining line
   */
  Point _point;

  /**
   * Direction of the constraining line
   */
  Point _direction;

  /**
   * Tolerance to use for numerical comparisons
   */
  Real _tol;
};

/**
 * Represents a plane constraint defined by a point and normal vector.
 */
class PlaneConstraint
{

public:
  PlaneConstraint() = default;

  /**
   * Constructor
   * @param point A point on the constraining plane.
   * @param normal the direction normal to the constraining plane.
   * @param tol The tolerance to use for numerical comparisons.
   */
  PlaneConstraint(const Point &point, const Point &normal,
                  const Real &tol = TOLERANCE * TOLERANCE);

  /**
   * Comparison operator for ordering PlaneConstraint objects.
   * The comparison is primarily based on the normal vector. If the normal
   * vectors are equal (within tolerance), the tie is broken using the dot
   * product of the normal with the point on the plane.
   * @param other The PlaneConstraint to compare with.
   * @return True if this PlaneConstraint is less than the other.
   */
  bool operator<(const PlaneConstraint &other) const;

  /**
   * Equality operator.
   * @param other The PlaneConstraint to compare with.
   * @return True if both PlaneConstraints represent the same plane.
   */
  bool operator==(const PlaneConstraint &other) const;

  /**
   * Query whether a point lies on the plane.
   * @param p The point in question
   * @return bool indicating whether p lies on the plane.
   */
  bool contains_point(const PointConstraint &p) const;

  /**
   * Query whether a line lies on the plane.
   * @param l The line in question
   * @return bool indicating whether l lies on the plane.
   */
  bool contains_line(const LineConstraint &l) const;

  /**
   * Query whether a plane is parallel to this plane
   * @param p The plane in question
   * @return bool indicating whether p is parallel to this plane.
   */
  bool is_parallel(const PlaneConstraint &p) const;

  /**
   * Query whether a line is parallel to this plane
   * @param l The line in question
   * @return bool indicating whether l is parallel to this plane.
   */
  bool is_parallel(const LineConstraint &l) const;

  /**
   * Computes the intersection of this plane with another constraint.
   * Handles intersection with PlaneConstraint, LineConstraint, or
   * PointConstraint.
   * @param other The constraint to intersect with.
   * @return The most specific ConstraintVariant that satisfies both
   * constraints. constraints. If no intersection exists, return an
   * InvalidConstraint.
   */
  ConstraintVariant intersect(const ConstraintVariant &other) const;

  /**
   * Const getter for the _point attribute
   */
  const Point &point() const { return _point; }

  /**
   * Const getter for the _normal attribute
   */
  const Point &normal() const { return _normal; }

  /**
   * Const getter for the _tol attribute
   */
  const Real &tol() const { return _tol; }

private:
  // Life is easier if we don't make these const
  /**
   * A point on the constraining plane
   */
  Point _point;

  /**
   * The direction normal to the constraining plane
   */
  Point _normal;

  /**
   * Tolerance to use for numerical comparisons
   */
  Real _tol;
};

/**
 * Represents an invalid constraint (i.e., when the two constraints don't
 * intersect)
 */
class InvalidConstraint
{

public:
  InvalidConstraint()
    : _err_msg("We should never get here! The InvalidConstraint object should be "
               "detected and replaced with a valid ConstraintVariant prior to calling "
               "any class methods.")
  {
  }

  /**
   * Dummy intersect method that should never be called.
   */
  ConstraintVariant intersect(const ConstraintVariant &) const {
    libmesh_assert_msg(false, _err_msg);
    return *this;
  }

  /**
   * Dummy contains_point method that should never be called.
   */
  bool contains_point(const PointConstraint &) const {
    libmesh_assert_msg(false, _err_msg);
    return false;
  }

private:
  std::string _err_msg;
};

/**
 * Dispatch intersection between two constraint variants.
 * Resolves to the appropriate method based on the type of the first operand.
 * @param a First constraint.
 * @param b Constraint to combine with a.
 * @return Combination (intersection) of constraint a and b.
 */
inline ConstraintVariant intersect_constraints(const ConstraintVariant &a,
                                               const ConstraintVariant &b) {
  // std::visit applies the visitor v (a Callable that can be called with any
  // combination of types from Variants) to the active value inside a
  // std::Variant. This circumvents the issue that the literal ConstraintVariant
  // type does not have a method called 'intersect' (but the types defining
  // ConstraintVariant do)
  return std::visit(
      [](const auto &lhs, const auto &rhs) -> ConstraintVariant {
        return lhs.intersect(rhs);
      },
      a, b);
}

/**
 * Constraint class for the VariationalMeshSmoother.
 *
 * Currently, all mesh boundary nodes are constrained to not move during
 * smoothing. If requested (preserve_subdomain_boundaries = true), nodes on
 * subdomain boundaries are also constrained to not move.
 */
class VariationalSmootherConstraint : public System::Constraint
{
private:

  System & _sys;

  /**
   * Whether nodes on subdomain boundaries are subject to change via smoothing
   */
  const bool _preserve_subdomain_boundaries;

  /**
   * Constrain (i.e., fix) a node to not move during mesh smoothing.
   * @param node Node to fix.
   */
  void fix_node(const Node & node);

  /**
   * Constrain a node to remain in the given plane during mesh smoothing.
   * @param node Node to constrain
   * @param ref_normal_vec Reference normal vector to the constraining plane.
   * This, along with the coordinates of node, are used to define the
   * constraining plane.
   */
  void constrain_node_to_plane(const Node & node, const Point & ref_normal_vec);

  /**
   * Constrain a node to remain on the given line during mesh smoothing.
   * @param node Node to constrain
   * @param line_vec vector parallel to the constraining line.
   * This, along with the coordinates of node, are used to define the
   * constraining line.
   */
  void constrain_node_to_line(const Node & node, const Point & line_vec);

  /**
   * Given a mesh and a node in the mesh, the vector will be filled with
   * every node directly attached to the given one. IF NO neighbors are found,
   * all nodes on the containing side are treated as neighbors. This is useful
   * when the node does not lie on an edge, such as the central face node in
   * HEX27 elements.
   */
  static void find_nodal_or_face_neighbors(
      const MeshBase & mesh,
      const Node & node,
      const std::unordered_map<dof_id_type, std::vector<const Elem *>> & nodes_to_elem_map,
      std::vector<const Node *> & neighbors);

  /**
   * Determines whether two neighboring nodes share a common boundary id.
   * @param boundary_node The first of the two prospective nodes.
   * @param neighbor_node The second of the two prospective nodes.
   * @param containing_elem The element containing node1 and node2.
   * @param boundary_info The mesh's BoundaryInfo.
   * @return nodes_share_bid Whether node1 and node2 share a common boundary id.
   */
  static bool nodes_share_boundary_id(
      const Node & boundary_node,
      const Node & neighbor_node,
      const Elem & containing_elem,
      const BoundaryInfo & boundary_info);

  /**
   * Get the relevant nodal neighbors for a subdomain constraint.
   * @param mesh The mesh being smoothed.
   * @param node The node (on the subdomain boundary) being constrained.
   * @param sub_id The subdomain id of the block on one side of the subdomain
   * boundary.
   * @param nodes_to_elem_map A mapping from node id to containing element ids.
   * @return A set of node pointer sets containing nodal neighbors to 'node' on
   * the sub_id1-sub_id2 boundary. The subsets are grouped by element faces
   * that form the subdomain boundary. Note that 'node' itself does not appear
   * in this set.
   */
  static std::set<std::set<const Node *>>
  get_neighbors_for_subdomain_constraint(
      const MeshBase &mesh, const Node &node, const subdomain_id_type sub_id,
      const std::unordered_map<dof_id_type, std::vector<const Elem *>>
          &nodes_to_elem_map);

  /**
   * Get the relevant nodal neighbors for an external boundary constraint.
   * @param mesh The mesh being smoothed.
   * @param node The node (on the external boundary) being constrained.
   * @param boundary_node_ids The set of mesh's external boundary node ids.
   * @param boundary_info The mesh's BoundaryInfo.
   * @param nodes_to_elem_map A mapping from node id to containing element ids.
   * @return A set of node pointer sets containing nodal neighbors to 'node' on
   * the external boundary. The subsets are grouped by element faces that form
   * the external boundary. Note that 'node' itself does not appear in this
   * set.
   */
  static std::set<std::set<const Node *>> get_neighbors_for_boundary_constraint(
      const MeshBase &mesh, const Node &node,
      const std::unordered_set<dof_id_type> &boundary_node_ids,
      const BoundaryInfo &boundary_info,
      const std::unordered_map<dof_id_type, std::vector<const Elem *>>
          &nodes_to_elem_map);

  /**
   * Determines the appropriate constraint (PointConstraint, LineConstraint, or
   * PlaneConstraint) for a node based on its neighbors.
   * @param node The node to constrain.
   * @param dim The mesh dimension.
   * @return The best-fit constraint for the given geometry.
   */
  static ConstraintVariant determine_constraint(
      const Node &node, const unsigned int dim,
      const std::set<std::set<const Node *>> &side_grouped_boundary_neighbors);

  /**
   * Applies a given constraint to a node (e.g., fixing it, restricting it to a
   * line or plane).
   * @param node The node to constrain.
   * @param constraint The geometric constraint variant to apply.
   * @throw libMesh::logicError If the constraint cannot be imposed.
   */
  void impose_constraint(const Node &node, const ConstraintVariant &constraint);

public:
  /**
   * Constructor
   * @param sys System to constrain.
   * @param preserve_subdomain_boundaries Whether to constrain nodes on
   * subdomain boundaries to not move.
   */
  VariationalSmootherConstraint(System & sys, const bool & preserve_subdomain_boundaries);

  virtual ~VariationalSmootherConstraint() override;

  virtual void constrain() override;
};


} // namespace libMesh

#endif // LIBMESH_VARIATIONAL_SMOOTHER_CONSTRAINT_H
