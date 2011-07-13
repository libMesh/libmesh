// $Id$
 
// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
  
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


#ifndef __mesh_triangle_support_h__
#define __mesh_triangle_support_h__

// C++ includes
#include <vector>

// Local Includes
#include "libmesh_config.h"
#include "enum_elem_type.h" // For ElemType declaration below
#include "point.h"
#include "libmesh.h"
#include "mesh_output.h" // for MeshSerializer... could this get its own header?
#include "unstructured_mesh.h"

#ifdef LIBMESH_HAVE_TRIANGLE

namespace libMesh
{

// Note: libmesh_common.h defines REAL, which is required by triangle.
// Therefore, we need to include it first.
#include "libmesh_common.h"

typedef Real REAL;

/**
 * A special namespace for wrapping the standard Triangle API,
 * as well as some helper functions for initializing/destroying
 * the structs triangle uses to communicate.
 */
namespace Triangle {
extern "C" {
#include "triangle.h"
}

  enum IO_Type {
    INPUT  = 0,
    OUTPUT = 1,
    BOTH   = 2};
  
  /**
   * Initializes the fields of t to NULL/0 as necessary.
   * This is helpful for preventing the access of uninitialized
   * memory when working with C, which has no constructors or
   * destructors.
   */
  void init(triangulateio& t);
    
  /**
   * Frees any memory which has been dynamically allocated by
   * Triangle.  Note the following facts:
   * 1) Triangle does not free any memory itself
   * 2) It is always safe to call free on a NULL pointer.
   *
   * However, triangle *does* shallow-copy (for example)
   * the holelist pointer from the input to output struct **without**
   * performing a deep copy of the holelist itself.  Therefore, double-free
   * will occur without additional care!
   */
  void destroy(triangulateio& t, IO_Type);

  /**
   * Copies triangulation data computed by triange from a triangulateio object
   * to a LibMesh mesh.   This routine is used internally by the
   * MeshTools::Generation::build_delaunay_square(...) and
   * MeshTools::Generation::build_delaunay_square_with_hole(...) routines.
   */
  void copy_tri_to_mesh(const triangulateio& triangle_data_input,
			UnstructuredMesh& mesh_output,
			const ElemType type);
    
  
}










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
  TriangleInterface(UnstructuredMesh& mesh)
    : _mesh(mesh),
      _holes(NULL),
      _elem_type(TRI3),
      _desired_area(0.1),
      _triangulation_type(GENERATE_CONVEX_HULL),
      _insert_extra_points(false),
      _smooth_after_generating(true),
      _serializer(_mesh)
  {}

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
      GENERATE_CONVEX_HULL = 0,
      /**
       * Uses the triangle library to first generate a convex hull from the set
       * of points passed in, and then triangulate this set of points.  This
       * is probably the most common type of usage.
       */

       
      PSLG = 1,
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
       
      INVALID_TRIANGULATION_TYPE
      /**
       * Does nothing, used as a default value.
       */
    };

  /**
   * The hole class and its several subclasses define the interface
   * and functionality of a "hole" which appears in a 2D mesh.  The
   * definitions of these classes appear below.
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
   * Sets and/or gets the desired triangle area.
   */
  Real& desired_area() {return _desired_area;}

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







/**
 * An abstract class for defining a 2-dimensional hole.  We assume that
 * the connectivity of the hole is implicit in the numbering of the points,
 * i.e. node 0 is connected to node 1, node 1 is connected to node 2, etc,
 * and the last node "wraps around" to connect back to node 0.
 */
class TriangleInterface::Hole
{
public:
  /**
   * Constructor
   */
  Hole() {}

  /**
   * Destructor
   */
  virtual ~Hole() {}

  /**
   * The number of geometric points which define the hole.
   */
  virtual unsigned int n_points() const = 0;

  /**
   * Return the nth point defining the hole.
   */
  virtual Point point(const unsigned int n) const = 0;

  /**
   * Return an (arbitrary) point which lies inside the hole.
   */
  virtual Point inside() const = 0;
};





    
/**
 * A concrete instantiation of the Hole class that describes polygonal
 * (triangular, square, pentagonal, ...) holes.
 */
class TriangleInterface::PolygonHole : public TriangleInterface::Hole
{
public:
  /**
   * Constructor specifying the center, radius, and number of
   * points which comprise the hole.  The points will all lie
   * on a circle of radius r.
   */
  PolygonHole(Point center, Real radius, unsigned int n_points) :
    _center(center),
    _radius(radius),
    _n_points(n_points) {}

  /**
   * Default Constructor, does not set any values
   */
  // PolygonHole() {}

  virtual unsigned int n_points() const { return _n_points; }

  virtual Point point(const unsigned int n) const
  {
    // The nth point lies at the angle theta = 2 * pi * n / _n_points
    const Real theta = static_cast<Real>(n) * 2.0 * libMesh::pi / static_cast<Real>(_n_points);
	
    return Point(_center(0) + _radius*std::cos(theta), // x=r*cos(theta)
		 _center(1) + _radius*std::sin(theta), // y=r*sin(theta)
		 0.);
  }

  /**
   * The center of the hole is definitely inside.
   */
  virtual Point inside() const { return _center;  }
      
private:
  /**
   * (x,y) location of the center of the hole
   */
  Point _center;

  /**
   * circular hole radius
   */
  Real _radius;

  /**
   * number of points used to describe the hole.  The actual
   * points can be generated knowing the center and radius.
   * For example, n_points=3 would generate a triangular hole.
   */
  unsigned int _n_points;
};





    

/**
 * Another concrete instantiation of the hole, this one should
 * be sufficiently general for most non-polygonal purposes.  The user
 * supplies, at the time of construction, a reference to a vector
 * of Points which defines the hole (in order of connectivity) and
 * an arbitrary Point which lies inside the hole. 
 */
class TriangleInterface::ArbitraryHole : public TriangleInterface::Hole
{
public:
  ArbitraryHole(const Point center,
		const std::vector<Point>& points)
    : 	_center(center),
	_points(points)
  {}

  /**
   * Required public Hole interface:
   */
      
  /**
   * The number of geometric points which define the hole.
   */
  virtual unsigned int n_points() const { return _points.size(); }

  /**
   * Return the nth point defining the hole.
   */
  virtual Point point(const unsigned int n) const
  {
    libmesh_assert (n < _points.size());
    return _points[n];
  }

  /**
   * Return an (arbitrary) point which lies inside the hole.
   */
  virtual Point inside() const { return _center;  }

private:
  /**
   * arbitrary (x,y) location inside the hole
   */
  Point _center;

  /**
   * Reference to the vector of points which makes up
   * the hole.
   */
  const std::vector<Point>& _points;
};





/**
 * Add the triangulation routines to the MeshTools::Generation namespace.
 */
namespace MeshTools
{
  namespace Generation
  {

    /**
     * Meshes a rectangular (2D) region (with or without holes) with a
     * Delaunay triangulation.  This function internally calls the
     * triangle library written by J.R. Shewchuk.
     */
    void build_delaunay_square(UnstructuredMesh& mesh,
			       const unsigned int nx, // num. of elements in x-dir
			       const unsigned int ny, // num. of elements in y-dir
			       const Real xmin, const Real xmax,
			       const Real ymin, const Real ymax,
			       const ElemType type,
			       const std::vector<TriangleInterface::Hole*>* holes=NULL);
  }
}


} // namespace libMesh



#endif // LIBMESH_HAVE_TRIANGLE

#endif // ifndef __mesh_triangle_support_h__
