// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


#ifndef LIBMESH_MESH_TRIANGLE_INTERFACE_H
#define LIBMESH_MESH_TRIANGLE_INTERFACE_H


#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_HAVE_TRIANGLE

// Local Includes
#include "libmesh/enum_elem_type.h"
#include "libmesh/libmesh.h"
#include "libmesh/mesh_serializer.h"

// C++ includes
#include <cstddef>
#include <vector>

namespace libMesh
{

  // Forward Declarations

  class UnstructuredMesh;

  /**
   * A C++ interface between LibMesh and the Triangle library written by
   * J.R. Shewchuk.
   *
   * @author John W. Peterson
   */
  class TriangleInterface
  {
  public:
    /**
     * The constructor.  A reference to the mesh containing the points
     * which are to be triangulated must be provided.  Unless otherwise
     * specified, a convex hull will be computed for the set of input points
     * and the convex hull will be meshed.
     */
    explicit
    TriangleInterface(UnstructuredMesh& mesh);

    /**
     * Empty destructor.
     */
    ~TriangleInterface() {}

    /**
     * The TriangulationType is used with the general triangulate function
     * defind below.
     */
    enum TriangulationType
      {
	/**
	 * Uses the triangle library to first generate a convex hull from the set
	 * of points passed in, and then triangulate this set of points.  This
	 * is probably the most common type of usage.
	 */
	GENERATE_CONVEX_HULL = 0,

	/**
	 * Use the triangle library to triangulate a Planar Straight Line
	 * Graph which is defined implicitly by the order of the "points" vector:
	 * a straight line is assumed to lie between each successive pair of
	 * points, with an additional line joining the final and first points.
	 * In case your triangulation is a little too "structured" looking
	 * (which can happen when the initial PSLG is really simple) you can try to
	 * make the resulting triangulation a little more "unstructured" looking
	 * by setting insert_points to true in the triangulate() function.
	 */
	PSLG = 1,

	/**
	 * Does nothing, used as a default value.
	 */
	INVALID_TRIANGULATION_TYPE
      };

    /**
     * The hole class and its several subclasses define the interface
     * and functionality of a "hole" which appears in a 2D mesh.
     * See mesh_triangle_holes.C/h for definitions.
     */
    class Hole;
    class PolygonHole;
    class ArbitraryHole;

    /**
     * This is the main public interface for this function.
     * Internally, it calls Triangle's triangulate routine.
     */
    void triangulate();

    /**
     * Sets and/or gets the desired element type.
     */
    ElemType& elem_type() {return _elem_type;}

    /**
     * Sets and/or gets the desired triangle area. Set to zero to disable
     * area constraint.
     */
    Real& desired_area() {return _desired_area;}

    /**
     * Sets and/or gets the minimum angle. Set to zero to disable area
     * constraint.
     */
    Real& minimum_angle() {return _minimum_angle;}

    /**
     * Sets and/or gets additional flags to be passed to triangle
     */
    std::string& extra_flags() {return _extra_flags;}

    /**
     * Sets and/or gets the desired triangulation type.
     */
    TriangulationType& triangulation_type() {return _triangulation_type;}

    /**
     * Sets and/or gets the flag for inserting add'l points.
     */
    bool& insert_extra_points() {return _insert_extra_points;}

    /**
     * Sets/gets flag which tells whether to do Delaunay mesh
     * smoothing after generating the grid.
     */
    bool& smooth_after_generating() {return _smooth_after_generating;}

    /**
     * Attaches a vector of Hole* pointers which will be
     * meshed around.
     */
    void attach_hole_list(const std::vector<Hole*>* holes) {_holes = holes;}

    /**
     * When constructing a PSLG, it is often not possible to do
     * so implicitly through the ordering of the points.  You
     * can use the segments vector to specify the segments explicitly,
     * Ex: unit square numbered counter-clockwise starting from origin
     * segments[0] = (0,1)
     * segments[1] = (1,2)
     * segments[2] = (2,3)
     * segments[3] = (3,0)
     * Note: for this case, you could use the implicit ordering!
     */
    std::vector<std::pair<unsigned int, unsigned int> > segments;

  private:
    /**
     * Reference to the mesh which is to be created by triangle.
     */
    UnstructuredMesh& _mesh;

    /**
     * A pointer to a vector of Hole*s.  If this is NULL, there
     * are no holes!
     */
    const std::vector<Hole*>* _holes;

    /**
     * The type of elements to generate.  (Defaults to
     * TRI3).
     */
    ElemType _elem_type;

    /**
     * The desired area for the elements in the resulting mesh.
     */
    Real _desired_area;

    /**
     * Minimum angle in triangles
     */
    Real _minimum_angle;

    /**
     * Additional flags to be passed to triangle
     */
    std::string _extra_flags;

    /**
     * The type of triangulation to perform: choices are:
     * convex hull
     * PSLG
     */
    TriangulationType _triangulation_type;

    /**
     * Flag which tells whether or not to insert additional nodes
     * before triangulation.  This can sometimes be used to "de-regularize"
     * the resulting triangulation.
     */
    bool _insert_extra_points;

    /**
     * Flag which tells whether we should smooth the mesh after
     * it is generated.  True by default.
     */
    bool _smooth_after_generating;

    /**
     * Triangle only operates on serial meshes.
     */
    MeshSerializer _serializer;
  };

} // namespace libMesh



#endif // LIBMESH_HAVE_TRIANGLE

#endif // ifndef LIBMESH_MESH_TRIANGLE_INTERFACE_H
