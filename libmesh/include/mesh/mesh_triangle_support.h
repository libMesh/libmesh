// $Id: mesh_triangle_support.h,v 1.7 2006-10-24 17:52:58 jwpeterson Exp $
 
// The libMesh Finite Element Library.
// Copyright (C) 2002  Benjamin S. Kirk, John W. Peterson
  
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
#include "elem_type.h" // For ElemType declaration below
#include "point.h"
#include "libmesh.h"

// Forward Declarations
class Mesh;

#ifdef HAVE_TRIANGLE

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
			Mesh& mesh_output,
			const ElemType type);
    
  
}






/**
 * Add the triangulation routines to the MeshTools::Generation namespace.
 */
namespace MeshTools
{
  namespace Generation
  {
    /**
     * Meshes a square (2D) with a Delaunay triangulation.
     * This function internally calls the triangle library
     * written by J.R. Shewchuk.  
     */
    void build_delaunay_square(Mesh& mesh,
			       const unsigned int nx, // num. of elements in x-dir
			       const unsigned int ny, // num. of elements in y-dir
			       const Real xmin=0., const Real xmax=1.,
			       const Real ymin=0., const Real ymax=1.,
			       const ElemType type=INVALID_ELEM);
    
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
     * Meshes an arbitrary 2D domain with a Delaunay triangulation.  This
     * function internally calls the triangle library written by
     * J.R. Shewchuk.
     *
     * The different TriangulationType enums are described above.
     *
     * Setting the insert_points boolean to true only has an effect when
     * triangulating a PSLG.  See decription above.
     * 
     * If you just want to triangulate a rectangular
     * domain, you might want to try build_delaunay_square() instead.
     */
    void triangulate(Mesh& mesh,
		     const std::vector<Point>& points,
		     const Real desired_area,
		     const ElemType type=INVALID_ELEM,
		     const TriangulationType tt=INVALID_TRIANGULATION_TYPE,
		     const bool insert_points=false);


    


    // An abstract class for defining a 2-dimensional hole.  We assume that
    // the connectivity of the hole is implicit in the numbering of the points,
    // i.e. node 0 is connected to node 1, node 1 is connected to node 2, etc,
    // and the last node "wraps around" to connect back to node 0.
    class Hole
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
    class PolygonHole : public Hole
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
	
	return Point(_center(0) + _radius*cos(theta), // x=r*cos(theta)
		     _center(1) + _radius*sin(theta), // y=r*sin(theta)
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
    class ArbitraryHole : public Hole
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
	assert (n < _points.size());
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

    

    // Function which calls triangle to create a mesh on a square domain with
    // one or more circular holes cut out.
    void build_delaunay_square_with_hole(Mesh& mesh,
					 const std::vector<Hole*>& holes,
					 const unsigned int nx=10, // num. of nodes in x-dir (approximate)
					 const unsigned int ny=10, // num. of nodes in y-dir (approximate)
					 const Real xmin=-1., const Real xmax=1.,
					 const Real ymin=-1., const Real ymax=1.,
					 const ElemType type=INVALID_ELEM);



  }
}




#endif // HAVE_TRIANGLE

#endif // ifndef __mesh_triangle_support_h__
