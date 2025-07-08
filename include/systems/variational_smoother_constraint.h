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

namespace libMesh
{
/*
 * Constraint class for the VariationalMeshSmoother.
 *
 * Currently, all mesh boundary nodes are constrained to not move during smoothing.
 * If requested (preserve_subdomain_boundaries = true), nodes on subdomain boundaries
 * are also constrained to not move.
 */
class VariationalSmootherConstraint : public System::Constraint
{
private:

  System & _sys;

  /// Whether nodes on subdomain boundaries are subject to change via smoothing
  const bool _preserve_subdomain_boundaries;

  /*
   * Identifies and imposes the appropriate constraints on a node.
   * @param node The node to constrain.
   * @param neighbors Vector of neighbors to use to identify the constraint to
   * impose.
   */
  void impose_constraints(const Node &node,
                          const std::vector<const Node *> neighbors);

  /*
   * Constrain (i.e., fix) a node to not move during mesh smoothing.
   * @param node Node to fix.
   */
  void fix_node(const Node & node);

  /*
   * Constrain a node to remain in the given plane during mesh smoothing.
   * @param node Node to constrain
   * @param ref_normal_vec Reference normal vector to the constraining plane.
   * This, along with the coordinates of node, are used to define the
   * constraining plane.
   */
  void constrain_node_to_plane(const Node & node, const Point & ref_normal_vec);

  /*
   * Constrain a node to remain on the given line during mesh smoothing.
   * @param node Node to constrain
   * @param line_vec vector parallel to the constraining line.
   * This, along with the coordinates of node, are used to define the
   * constraining line.
   */
  void constrain_node_to_line(const Node & node, const Point & line_vec);

  /*
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

  static void filter_neighbors_for_subdomain_constraint(
      const Node & node,
      std::vector<const Node *> & neighbors,
      const subdomain_id_type sub_id1,
      const subdomain_id_type sub_id2,
      std::unordered_map<dof_id_type, std::vector<const Elem *>> & nodes_to_elem_map
    );

  static void filter_neighbors_for_boundary_constraint(
      const Node & node,
      std::vector<const Node *> & neighbors,
      std::unordered_map<dof_id_type, std::vector<const Elem *>> & nodes_to_elem_map,
      const std::unordered_set<dof_id_type> & boundary_node_ids,
      const BoundaryInfo & boundary_info
    );

public:

  /*
   * Constructor
   * @param sys System to constrain.
   * @param preserve_subdomain_boundaries Whether to constrain nodes on subdomain boundaries to not move.
   */
  VariationalSmootherConstraint(System & sys, const bool & preserve_subdomain_boundaries);

  virtual ~VariationalSmootherConstraint() override;

  virtual void constrain() override;
};


} // namespace libMesh

#endif // LIBMESH_VARIATIONAL_SMOOTHER_CONSTRAINT_H
