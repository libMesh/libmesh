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

#ifdef LIBMESH_HAVE_TRIANGLE

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
// Function definitions for the TrianguleInterface class
//

// Constructor
TriangleInterface::TriangleInterface(UnstructuredMesh & mesh)
  : TriangulatorInterface(mesh),
    _extra_flags(""),
    _serializer(_mesh)
{}




// Primary function responsible for performing the triangulation
void TriangleInterface::triangulate()
{
  // Will the triangulation have holes?
  const bool have_holes = ((_holes != nullptr) && (!_holes->empty()));

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

  // Regardless of whether we added additional points, the set of points to
  // triangulate is now sitting in the mesh.

  // If the holes vector is non-nullptr (and non-empty) we need to determine
  // the number of additional points which the holes will add to the
  // triangulation.
  // Note that the number of points is always equal to the number of segments
  // that form the holes.
  unsigned int n_hole_points = 0;

  if (have_holes)
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

  // Triangle data structure for the mesh
  TriangleWrapper::triangulateio initial;
  TriangleWrapper::triangulateio final;
  TriangleWrapper::triangulateio voronoi;

  // Pseudo-Constructor for the triangle io structs
  TriangleWrapper::init(initial);
  TriangleWrapper::init(final);
  TriangleWrapper::init(voronoi);

  initial.numberofpoints = _mesh.n_nodes() + n_hole_points;
  initial.pointlist      = static_cast<REAL*>(std::malloc(initial.numberofpoints * 2 * sizeof(REAL)));

  if (_triangulation_type==PSLG)
    {
      // Implicit segment ordering: One segment per point, including hole points
      if (this->segments.empty())
        initial.numberofsegments = initial.numberofpoints;

      // User-defined segment ordering: One segment per entry in the segments vector
      else
        initial.numberofsegments = this->segments.size() + n_hole_points;
    }

  else if (_triangulation_type==GENERATE_CONVEX_HULL)
    initial.numberofsegments = n_hole_points; // One segment for each hole point

  // Allocate space for the segments (2 int per segment)
  if (initial.numberofsegments > 0)
    {
      initial.segmentlist = static_cast<int *> (std::malloc(initial.numberofsegments * 2 * sizeof(int)));
      if (_markers)
        initial.segmentmarkerlist = static_cast<int *> (std::malloc(initial.numberofsegments * sizeof(int)));
    }


  // Copy all the holes' points and segments into the triangle struct.

  // The hole_offset is a constant offset into the points vector which points
  // past the end of the last hole point added.
  unsigned int hole_offset=0;

  if (have_holes)
    for (const auto & hole : *_holes)
      {
        for (unsigned int ctr=0, h=0, i=0, hsism=hole->segment_indices().size()-1; i<hsism; ++i)
          {
            unsigned int begp = hole_offset + hole->segment_indices()[i];
            unsigned int endp = hole->segment_indices()[i+1];

            for (; h<endp; ctr+=2, ++h)
              {
                Point p = hole->point(h);

                const unsigned int index0 = 2*hole_offset+ctr;
                const unsigned int index1 = 2*hole_offset+ctr+1;

                // Save the x,y locations in the triangle struct.
                initial.pointlist[index0] = p(0);
                initial.pointlist[index1] = p(1);

                // Set the points which define the segments
                initial.segmentlist[index0] = hole_offset+h;
                initial.segmentlist[index1] = (h == endp - 1) ? begp : hole_offset + h + 1; // wrap around
                if (_markers)
                  // 1 is reserved for boundaries of holes
                  initial.segmentmarkerlist[hole_offset+h] = 1;
              }
          }

        // Update the hole_offset for the next hole
        hole_offset += hole->n_points();
      }


  // Copy all the non-hole points and segments into the triangle struct.
  {
    dof_id_type ctr=0;
    for (auto & node : _mesh.node_ptr_range())
      {
        dof_id_type index = 2*hole_offset + ctr;

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
                if (_markers)
                  initial.segmentmarkerlist[hole_offset + n] = (*_markers)[n];
              }
          }

        ctr +=2;
      }
  }


  // If the user provided it, use his ordering to define the segments
  for (std::size_t ctr=0, s=0, ss=this->segments.size(); s<ss; ctr+=2, ++s)
    {
      const unsigned int index0 = 2*hole_offset+ctr;
      const unsigned int index1 = 2*hole_offset+ctr+1;

      initial.segmentlist[index0] = hole_offset + this->segments[s].first;
      initial.segmentlist[index1] = hole_offset + this->segments[s].second;
      if (_markers)
        initial.segmentmarkerlist[hole_offset + s] = (*_markers)[s];
    }



  // Tell the input struct about the holes
  if (have_holes)
    {
      initial.numberofholes = _holes->size();
      initial.holelist      = static_cast<REAL*>(std::malloc(initial.numberofholes * 2 * sizeof(REAL)));
      for (std::size_t i=0, ctr=0, hs=_holes->size(); i<hs; ++i, ctr+=2)
        {
          Point inside_point = (*_holes)[i]->inside();
          initial.holelist[ctr]   = inside_point(0);
          initial.holelist[ctr+1] = inside_point(1);
        }
    }

  if (_regions)
    {
      initial.numberofregions = _regions->size();
      initial.regionlist      = static_cast<REAL*>(std::malloc(initial.numberofregions * 4 * sizeof(REAL)));
      for (std::size_t i=0, ctr=0, rs=_regions->size(); i<rs; ++i, ctr+=4)
        {
          Point inside_point = (*_regions)[i]->inside();
          initial.regionlist[ctr]   = inside_point(0);
          initial.regionlist[ctr+1] = inside_point(1);
          initial.regionlist[ctr+2] = (*_regions)[i]->attribute();
          initial.regionlist[ctr+3] = (*_regions)[i]->max_area();
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
  //     `segmentmarkerlist' must either be set to nullptr (in which case all
  //     markers default to zero), or must point to a list of markers.
  // D ~ Conforming Delaunay: use this switch if you want all triangles
  //     in the mesh to be Delaunay, and not just constrained Delaunay
  // q ~  Quality mesh generation with no angles smaller than 20 degrees.
  //      An alternate minimum angle may be specified after the q
  // a ~ Imposes a maximum triangle area constraint.
  // -P  Suppresses the output .poly file. Saves disk space, but you lose the ability to maintain
  //     constraining segments on later refinements of the mesh.
  // -e  Outputs (to an .edge file) a list of edges of the triangulation.
  // -v  Outputs the Voronoi diagram associated with the triangulation.
  // Create the flag strings, depends on element type
  std::ostringstream flags;

  // Default flags always used
  flags << "z";

  if (_quiet)
    flags << "QP";
  else
    flags << "V";

  if (_markers)
    flags << "ev";

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
      libmesh_error_msg("ERROR: Unrecognized triangular element type == " << Utility::enum_to_string(_elem_type));
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

  if (_regions)
    flags << "Aa";

  // add user provided extra flags
  if (_extra_flags.size() > 0)
    flags << _extra_flags;

  // Refine the initial output to conform to the area constraint
  if (_markers)
  {
    // need Voronoi to generate boundary information
    TriangleWrapper::triangulate(const_cast<char *>(flags.str().c_str()),
                                 &initial,
                                 &final,
                                 &voronoi);

    // Send the information computed by Triangle to the Mesh.
    TriangleWrapper::copy_tri_to_mesh(final,
                                      _mesh,
                                      _elem_type,
                                      &voronoi);
  }
  else
  {
    TriangleWrapper::triangulate(const_cast<char *>(flags.str().c_str()),
                                 &initial,
                                 &final,
                                 nullptr);
    // Send the information computed by Triangle to the Mesh.
    TriangleWrapper::copy_tri_to_mesh(final,
                                      _mesh,
                                      _elem_type);
  }



  // To the naked eye, a few smoothing iterations usually looks better,
  // so we do this by default unless the user says not to.
  if (this->_smooth_after_generating)
    LaplaceMeshSmoother(_mesh).smooth(2);


  // Clean up.
  TriangleWrapper::destroy(initial,      TriangleWrapper::INPUT);
  TriangleWrapper::destroy(final,        TriangleWrapper::OUTPUT);
  TriangleWrapper::destroy(voronoi,      TriangleWrapper::OUTPUT);

  // Prepare the mesh for use before returning.  This ensures (among
  // other things) that it is partitioned and therefore users can
  // iterate over local elements, etc.
  _mesh.prepare_for_use();
}

} // namespace libMesh







#endif // LIBMESH_HAVE_TRIANGLE
