// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#ifdef LIBMESH_HAVE_TRIANGLE

// C/C++ includes
#include <sstream>

#include "libmesh/mesh_triangle_interface.h"
#include "libmesh/unstructured_mesh.h"
#include "libmesh/face_tri3.h"
#include "libmesh/face_tri6.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/mesh_smoother_laplace.h"
#include "libmesh/boundary_info.h"
#include "libmesh/mesh_triangle_holes.h"
#include "libmesh/mesh_triangle_wrapper.h"

namespace libMesh
{
//
// Function definitions for the TriangleInterface class
//

// Constructor
TriangleInterface::TriangleInterface(UnstructuredMesh & mesh)
  : _mesh(mesh),
    _holes(libmesh_nullptr),
    _elem_type(TRI3),
    _desired_area(0.1),
    _minimum_angle(20.0),
    _extra_flags(""),
    _triangulation_type(GENERATE_CONVEX_HULL),
    _insert_extra_points(false),
    _smooth_after_generating(true),
    _serializer(_mesh)
{}



// Primary function responsible for performing the triangulation
void TriangleInterface::triangulate()
{
  // Will the triangulation have holes?
  const bool have_holes = ((_holes != libmesh_nullptr) && (!_holes->empty()));

  // If the initial PSLG is really simple, e.g. an L-shaped domain or
  // a square/rectangle, the resulting triangulation may be very
  // "structured" looking.  Sometimes this is a problem if your
  // intention is to work with an "unstructured" looking grid.  We can
  // attempt to work around this limitation by inserting midpoints
  // into the original PSLG.  Inserting additional points into a
  // set of points meant to be a convex hull usually makes less sense.

  // May or may not need to insert new points ...
  if ((_triangulation_type==PSLG) && (_insert_extra_points))
    {
      // Make a copy of the original points from the Mesh
      std::vector<Point> original_points (_mesh.n_nodes());

      MeshBase::node_iterator       node_it  = _mesh.nodes_begin();
      const MeshBase::node_iterator node_end = _mesh.nodes_end();

      for (unsigned int ctr=0; node_it != node_end; ++node_it)
        original_points[ctr++] = **node_it;

      // Clear out the mesh
      _mesh.clear();

      // Make sure the new Mesh will be 2D
      _mesh.set_mesh_dimension(2);

      // Insert a new point on each PSLG at some random location
      // np=index into new points vector
      // n =index into original points vector
      for (std::size_t np=0, n=0; np<2*original_points.size(); ++np)
        {
          // the even entries are the original points
          if (np%2==0)
            _mesh.add_point(original_points[n++]);

          else // the odd entries are the midpoints of the original PSLG segments
            _mesh.add_point ((original_points[n] + original_points[n-1])/2);
        }
    }

  // Regardless of whether we added additional points, the set of points to
  // triangulate is now sitting in the mesh.

  // If the holes vector is non-NULL (and non-empty) we need to determine
  // the number of additional points which the holes will add to the
  // triangulation.
  unsigned int n_hole_points = 0;

  if (have_holes)
    {
      for (std::size_t i=0; i<_holes->size(); ++i)
        n_hole_points += (*_holes)[i]->n_points();
    }

  // Triangle data structure for the mesh
  TriangleWrapper::triangulateio initial;
  TriangleWrapper::triangulateio final;

  // Pseudo-Constructor for the triangle io structs
  TriangleWrapper::init(initial);
  TriangleWrapper::init(final);

  initial.numberofpoints = _mesh.n_nodes() + n_hole_points;
  initial.pointlist      = static_cast<REAL*>(std::malloc(initial.numberofpoints * 2 * sizeof(REAL)));

  if (_triangulation_type==PSLG)
    {
      // Implicit segment ordering: One segment per point, including hole points
      if (this->segments.empty())
        initial.numberofsegments = initial.numberofpoints;

      // User-defined segment ordering: One segment per entry in the segments vector
      else
        initial.numberofsegments = this->segments.size();
    }

  else if (_triangulation_type==GENERATE_CONVEX_HULL)
    initial.numberofsegments = n_hole_points; // One segment for each hole point

  // Debugging
  // libMesh::out << "Number of segments set to: " << initial.numberofsegments << std::endl;

  // Allocate space for the segments (2 int per segment)
  if (initial.numberofsegments > 0)
    {
      initial.segmentlist = static_cast<int *> (std::malloc(initial.numberofsegments * 2 * sizeof(int)));
    }


  // Copy all the holes' points and segments into the triangle struct.

  // The hole_offset is a constant offset into the points vector which points
  // past the end of the last hole point added.
  unsigned int hole_offset=0;

  if (have_holes)
    for (std::size_t i=0; i<_holes->size(); ++i)
      {
        for (unsigned int ctr=0, h=0; h<(*_holes)[i]->n_points(); ctr+=2, ++h)
          {
            Point p = (*_holes)[i]->point(h);

            const unsigned int index0 = 2*hole_offset+ctr;
            const unsigned int index1 = 2*hole_offset+ctr+1;

            // Save the x,y locations in the triangle struct.
            initial.pointlist[index0] = p(0);
            initial.pointlist[index1] = p(1);

            // Set the points which define the segments
            initial.segmentlist[index0] = hole_offset+h;
            initial.segmentlist[index1] = (h==(*_holes)[i]->n_points()-1) ? hole_offset : hole_offset+h+1; // wrap around
          }

        // Update the hole_offset for the next hole
        hole_offset += (*_holes)[i]->n_points();
      }


  // Copy all the non-hole points and segments into the triangle struct.
  {
    MeshBase::node_iterator it = _mesh.nodes_begin();
    const MeshBase::node_iterator end = _mesh.nodes_end();

    for (dof_id_type ctr=0; it != end; ctr+=2, ++it)
      {
        dof_id_type index = 2*hole_offset + ctr;

        // Get pointer to the current node
        Node * node = *it;

        // Set x,y values in pointlist
        initial.pointlist[index] = (*node)(0);
        initial.pointlist[index+1] = (*node)(1);

        // If the user requested a PSLG, the non-hole points are also segments
        if (_triangulation_type==PSLG)
          {
            // Use implicit ordering to define segments
            if (this->segments.empty())
              {
                dof_id_type n = ctr/2; // ctr is always even
                initial.segmentlist[index] = hole_offset+n;
                initial.segmentlist[index+1] = (n==_mesh.n_nodes()-1) ? hole_offset : hole_offset+n+1; // wrap around
              }
          }
      }
  }


  // If the user provided it, use his ordering to define the segments
  for (std::size_t ctr=0, s=0; s<this->segments.size(); ctr+=2, ++s)
    {
      const unsigned int index0 = 2*hole_offset+ctr;
      const unsigned int index1 = 2*hole_offset+ctr+1;

      initial.segmentlist[index0] = hole_offset + this->segments[s].first;
      initial.segmentlist[index1] = hole_offset + this->segments[s].second;
    }



  // Tell the input struct about the holes
  if (have_holes)
    {
      initial.numberofholes = _holes->size();
      initial.holelist      = static_cast<REAL*>(std::malloc(initial.numberofholes * 2 * sizeof(REAL)));
      for (std::size_t i=0, ctr=0; i<_holes->size(); ++i, ctr+=2)
        {
          Point inside_point = (*_holes)[i]->inside();
          initial.holelist[ctr]   = inside_point(0);
          initial.holelist[ctr+1] = inside_point(1);
        }
    }

  // Set the triangulation flags.
  // c ~ enclose convex hull with segments
  // z ~ use zero indexing
  // B ~ Suppresses boundary markers in the output
  // Q ~ run in "quiet" mode
  // p ~ Triangulates a Planar Straight Line Graph
  //     If the `p' switch is used, `segmentlist' must point to a list of
  //     segments, `numberofsegments' must be properly set, and
  //     `segmentmarkerlist' must either be set to NULL (in which case all
  //     markers default to zero), or must point to a list of markers.
  // D ~ Conforming Delaunay: use this switch if you want all triangles
  //     in the mesh to be Delaunay, and not just constrained Delaunay
  // q ~  Quality mesh generation with no angles smaller than 20 degrees.
  //      An alternate minimum angle may be specified after the q
  // a ~ Imposes a maximum triangle area constraint.
  // -P  Suppresses the output .poly file. Saves disk space, but you lose the ability to maintain
  //     constraining segments on later refinements of the mesh.
  // Create the flag strings, depends on element type
  std::ostringstream flags;

  // Default flags always used
  flags << "zBPQ";

  // Flags which are specific to the type of triangulation
  switch (_triangulation_type)
    {
    case GENERATE_CONVEX_HULL:
      {
        flags << "c";
        break;
      }

    case PSLG:
      {
        flags << "p";
        break;
      }

    case INVALID_TRIANGULATION_TYPE:
      libmesh_error_msg("ERROR: INVALID_TRIANGULATION_TYPE selected!");

    default:
      libmesh_error_msg("Unrecognized _triangulation_type");
    }


  // Flags specific to the type of element
  switch (_elem_type)
    {
    case TRI3:
      {
        // do nothing.
        break;
      }

    case TRI6:
      {
        flags << "o2";
        break;
      }

    default:
      libmesh_error_msg("ERROR: Unrecognized triangular element type.");
    }


  // If we do have holes and the user asked to GENERATE_CONVEX_HULL,
  // need to add the p flag so the triangulation respects those segments.
  if ((_triangulation_type==GENERATE_CONVEX_HULL) && (have_holes))
    flags << "p";

  // Finally, add the area constraint
  if (_desired_area > TOLERANCE)
    flags << "a" << std::fixed << _desired_area;

  // add minimum angle constraint
  if (_minimum_angle > TOLERANCE)
    flags << "q" << std::fixed << _minimum_angle;

  // add user provided extra flags
  if (_extra_flags.size() > 0)
    flags << _extra_flags;

  // Refine the initial output to conform to the area constraint
  TriangleWrapper::triangulate(const_cast<char *>(flags.str().c_str()),
                               &initial,
                               &final,
                               libmesh_nullptr); // voronoi ouput -- not used


  // Send the information computed by Triangle to the Mesh.
  TriangleWrapper::copy_tri_to_mesh(final,
                                    _mesh,
                                    _elem_type);

  // To the naked eye, a few smoothing iterations usually looks better,
  // so we do this by default unless the user says not to.
  if (this->_smooth_after_generating)
    LaplaceMeshSmoother(_mesh).smooth(2);


  // Clean up.
  TriangleWrapper::destroy(initial,      TriangleWrapper::INPUT);
  TriangleWrapper::destroy(final,        TriangleWrapper::OUTPUT);

  // Prepare the mesh for use before returning.  This ensures (among
  // other things) that it is partitioned and therefore users can
  // iterate over local elements, etc.
  _mesh.prepare_for_use();
}

} // namespace libMesh







#endif // LIBMESH_HAVE_TRIANGLE
