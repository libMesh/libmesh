// The libMesh Finite Element Library.
// Copyright (C) 2002-2024 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


#include "libmesh/libmesh_config.h"

// libmesh includes
#include "libmesh/mesh_triangle_interface.h"
#include "libmesh/unstructured_mesh.h"
#include "libmesh/face_tri3.h"
#include "libmesh/face_tri6.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/mesh_smoother_laplace.h"
#include "libmesh/boundary_info.h"
#include "libmesh/mesh_triangle_holes.h"
#include "libmesh/mesh_triangle_wrapper.h"
#include "libmesh/enum_elem_type.h"
#include "libmesh/enum_to_string.h"
#include "libmesh/utility.h"

// C/C++ includes
#include <sstream>


namespace libMesh
{
//
// Function definitions for the TriangulatorInterface class
//

// Constructor
TriangulatorInterface::TriangulatorInterface(UnstructuredMesh & mesh)
  : _mesh(mesh),
    _holes(nullptr),
    _markers(nullptr),
    _regions(nullptr),
    _elem_type(TRI3),
    _desired_area(0.1),
    _minimum_angle(20.0),
    _triangulation_type(GENERATE_CONVEX_HULL),
    _insert_extra_points(false),
    _smooth_after_generating(true),
    _quiet(true)
{}


void TriangulatorInterface::set_interpolate_boundary_points (int n_points)
{
  // Maybe we'll reserve a meaning for negatives later?
  libmesh_assert(n_points >= 0);

  _interpolate_boundary_points = n_points;

  // backwards compatibility - someone (including us) might want to
  // query this via the old API.
  _insert_extra_points = n_points;
}



int TriangulatorInterface::get_interpolate_boundary_points () const
{
  // backwards compatibility - someone might have turned this off via
  // the old API
  if (!_insert_extra_points)
    return 0;

  return _interpolate_boundary_points;
}



void TriangulatorInterface::elems_to_segments()
{
  // Don't try to override manually specified segments
  if (!this->segments.empty())
    return;

  // If we have edges, they should form the polyline with the ordering
  // we want.  Let's turn them into segments for later use, because
  // we're going to delete the original elements to replace with our
  // triangulation.
  if (_mesh.n_elem())
    {
      // Mapping from points to node ids, to back those out from
      // MeshedHole results later
      std::map<Point, dof_id_type> point_id_map;

      for (Node * node : _mesh.node_ptr_range())
        {
          // We're not going to support overlapping nodes on the boundary
          libmesh_error_msg_if
            (point_id_map.count(*node),
             "TriangulatorInterface does not support overlapping nodes found at "
             << static_cast<Point&>(*node));

          point_id_map.emplace(*node, node->id());
        }

      // We don't support directly generating Tri6, so for
      // compatibility with future stitching we need to be working
      // with first-order elements.  Let's get rid of any non-vertex
      // nodes we just added.
      for (Elem * elem : _mesh.element_ptr_range())
        for (auto n : make_range(elem->n_vertices(), elem->n_nodes()))
          point_id_map.erase(elem->point(n));

      // We'll steal the ordering calculation from
      // the MeshedHole code
      const TriangulatorInterface::MeshedHole mh { _mesh, this->_bdy_ids };

      // If we've specified only a subset of the mesh as our outer
      // boundary, then we may have nodes that don't actually fall
      // inside that boundary.  Triangulator code doesn't like Steiner
      // points that aren't inside the triangulation domain, so we
      // need to get rid of them.
      //
      // Also, if we're using Edge3 elements to define our outer
      // boundary, we're only dealing with their 2 end nodes and we'll
      // need to get rid of their central nodes.
      std::unordered_set<Node *> nodes_to_delete;

      for (Elem * elem : _mesh.element_ptr_range())
        for (auto n : make_range(elem->n_vertices(), elem->n_nodes()))
          nodes_to_delete.insert(elem->node_ptr(n));

      if (!this->_bdy_ids.empty())
        {
          for (auto & node : _mesh.node_ptr_range())
            if (!mh.contains(*node))
              nodes_to_delete.insert(node);
        }

      // And now we're done with elements.  Delete them lest they have
      // dangling pointers to nodes we'll be deleting.
      _mesh.clear_elems();

      // Make segments from boundary nodes; also make sure we don't
      // delete them.
      const std::size_t np = mh.n_points();
      for (auto i : make_range(np))
        {
          const Point pt = mh.point(i);
          const dof_id_type id0 = libmesh_map_find(point_id_map, pt);
          nodes_to_delete.erase(_mesh.node_ptr(id0));
          const Point next_pt = mh.point((np+i+1)%np);
          const dof_id_type id1 = libmesh_map_find(point_id_map, next_pt);
          this->segments.emplace_back(id0, id1);
        }

      for (Node * node : nodes_to_delete)
        _mesh.delete_node(node);

      if (this->_verify_hole_boundaries && _holes)
        this->verify_holes(mh);
    }
}



void TriangulatorInterface::nodes_to_segments(dof_id_type max_node_id)
{
  // Don't try to override manually specified segments, or try to add
  // segments if we're doing a convex hull
  if (!this->segments.empty() || _triangulation_type != PSLG)
    return;

  for (auto node_it = _mesh.nodes_begin(),
       node_end = _mesh.nodes_end();
       node_it != node_end;)
    {
      Node * node = *node_it;

      // If we're out of boundary nodes, the rest are going to be
      // Steiner points or hole points
      if (node->id() >= max_node_id)
        break;

      ++node_it;

      Node * next_node = (node_it == node_end) ?
        *_mesh.nodes_begin() : *node_it;

      this->segments.emplace_back(node->id(), next_node->id());
    }

  if (this->_verify_hole_boundaries && _holes)
    {
      std::vector<Point> outer_pts;
      for (auto segment : this->segments)
        outer_pts.push_back(_mesh.point(segment.first));

      ArbitraryHole ah(outer_pts);
      this->verify_holes(ah);
    }
}



void TriangulatorInterface::insert_any_extra_boundary_points()
{
  // If the initial PSLG is really simple, e.g. an L-shaped domain or
  // a square/rectangle, the resulting triangulation may be very
  // "structured" looking.  Sometimes this is a problem if your
  // intention is to work with an "unstructured" looking grid.  We can
  // attempt to work around this limitation by inserting midpoints
  // into the original PSLG.  Inserting additional points into a
  // set of points meant to be a convex hull usually makes less sense.

  const int n_interpolated = this->get_interpolate_boundary_points();
  if ((_triangulation_type==PSLG) && n_interpolated)
    {
      // If we were lucky enough to start with contiguous node ids,
      // let's keep them that way.
      dof_id_type nn = _mesh.max_node_id();

      std::vector<std::pair<unsigned int, unsigned int>> old_segments =
        std::move(this->segments);

      // We expect to have converted any elems and/or nodes into
      // segments by now.
      libmesh_assert(!old_segments.empty());

      this->segments.clear();

      // Insert a new point on each segment at evenly spaced locations
      // between existing boundary points.
      // np=index into new points vector
      // n =index into original points vector
      for (auto old_segment : old_segments)
        {
          Node * begin_node = _mesh.node_ptr(old_segment.first);
          Node * end_node = _mesh.node_ptr(old_segment.second);
          dof_id_type current_id = begin_node->id();
          for (auto i : make_range(n_interpolated))
            {
              // new points are equispaced along the original segments
              const Point new_point =
                ((n_interpolated-i) * *(Point *)(begin_node) +
                 (i+1) * *(Point *)(end_node)) /
                (n_interpolated + 1);
              Node * next_node = _mesh.add_point(new_point, nn++);
              this->segments.emplace_back(current_id,
                                          next_node->id());
              current_id = next_node->id();
            }
          this->segments.emplace_back(current_id,
                                      end_node->id());
        }
    }
}


void TriangulatorInterface::verify_holes(const Hole & outer_bdy)
{
  for (const Hole * hole : *_holes)
    {
      for (const Hole * hole2 : *_holes)
        {
          if (hole == hole2)
            continue;

          for (auto i : make_range(hole2->n_points()))
            if (hole->contains(hole2->point(i)))
              libmesh_error_msg
                ("Found point " << hole2->point(i) <<
                 " on one hole boundary and another's interior");
        }

      for (auto i : make_range(hole->n_points()))
        if (!outer_bdy.contains(hole->point(i)))
          libmesh_error_msg
            ("Found point " << hole->point(i) <<
             " on hole boundary but outside outer boundary");
    }
}


unsigned int TriangulatorInterface::total_hole_points()
{
  // If the holes vector is non-nullptr (and non-empty) we need to determine
  // the number of additional points which the holes will add to the
  // triangulation.
  // Note that the number of points is always equal to the number of segments
  // that form the holes.
  unsigned int n_hole_points = 0;

  if (_holes)
    for (const auto & hole : *_holes)
    {
      n_hole_points += hole->n_points();
      // A hole at least has one enclosure.
      // Points on enclosures are ordered so that we can add segments implicitly.
      // Elements in segment_indices() indicates the starting points of all enclosures.
      // The last element in segment_indices() is the number of total points.
      libmesh_assert_greater(hole->segment_indices().size(), 1);
      libmesh_assert_equal_to(hole->segment_indices().back(), hole->n_points());
    }

  return n_hole_points;
}

} // namespace libMesh

