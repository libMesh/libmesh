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
