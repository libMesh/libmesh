// The libMesh Finite Element Library.
// Copyright (C) 2002-2021 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


void TriangulatorInterface::insert_any_extra_boundary_points()
{
  // If the initial PSLG is really simple, e.g. an L-shaped domain or
  // a square/rectangle, the resulting triangulation may be very
  // "structured" looking.  Sometimes this is a problem if your
  // intention is to work with an "unstructured" looking grid.  We can
  // attempt to work around this limitation by inserting midpoints
  // into the original PSLG.  Inserting additional points into a
  // set of points meant to be a convex hull usually makes less sense.

  if ((_triangulation_type==PSLG) && (_insert_extra_points))
    {
      // Make a copy of the original points from the Mesh
      std::vector<Point> original_points;
      original_points.reserve (_mesh.n_nodes());
      for (auto & node : _mesh.node_ptr_range())
        original_points.push_back(*node);

      // Clear out the mesh
      _mesh.clear();

      // Make sure the new Mesh will be 2D
      _mesh.set_mesh_dimension(2);

      // Insert a new point on each PSLG at some random location
      // np=index into new points vector
      // n =index into original points vector
      for (std::size_t np=0, n=0, tops=2*original_points.size(); np<tops; ++np)
        {
          // the even entries are the original points
          if (np%2==0)
            _mesh.add_point(original_points[n++]);

          else // the odd entries are the midpoints of the original PSLG segments
            _mesh.add_point ((original_points[n] + original_points[n-1])/2);
        }
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

