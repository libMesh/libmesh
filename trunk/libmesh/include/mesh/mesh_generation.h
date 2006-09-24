// $Id: mesh_generation.h,v 1.12 2006-09-24 05:22:29 jwpeterson Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2005  Benjamin S. Kirk, John W. Peterson
  
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



#ifndef __mesh_generation_h__
#define __mesh_generation_h__



// C++ Includes   -----------------------------------
#include <vector>

// Local Includes -----------------------------------
// #include "libmesh_common.h" // needed for Real
#include "enum_elem_type.h" // needed for ElemType enum
#include "libmesh.h"        // needed for libMesh::invalid_uint
#include "point.h"

// forward declarations
class Mesh;



// ------------------------------------------------------------
// MeshTools::Generation namespace
namespace MeshTools
{
  /**
   * Tools for \p Mesh generation.
   *
   * \author Benjamin S. Kirk
   * \date 2004
   * \version $Revision: 1.12 $
   */
  namespace Generation
  {
    /**
     * Builds a \f$ nx \times ny \times nz \f$ (elements) cube.
     * Defaults to a unit cube (or line in 1D, square in 2D),
     * but the dimensions can be specified through the optional
     * arguments.
     */  
    void build_cube (Mesh& mesh,
		     const unsigned int nx=0,
		     const unsigned int ny=0,
		     const unsigned int nz=0,
		     const Real xmin=0., const Real xmax=1.,
		     const Real ymin=0., const Real ymax=1.,
		     const Real zmin=0., const Real zmax=1.,
		     const ElemType type=INVALID_ELEM,
		     const bool gauss_lobatto_grid=false);

    /**
     * A specialized \p build_cube() for 1D meshes
     */
    void build_line (Mesh& mesh,
                     const unsigned int nx,
                     const Real xmin=0., const Real xmax=1.,
                     const ElemType type=INVALID_ELEM,
                     const bool gauss_lobatto_grid=false);

    /**
     * A specialized \p build_cube() for 2D meshes.
     */
    void build_square (Mesh& mesh,
		       const unsigned int nx,
		       const unsigned int ny,
		       const Real xmin=0., const Real xmax=1.,
		       const Real ymin=0., const Real ymax=1.,
		       const ElemType type=INVALID_ELEM,
		       const bool gauss_lobatto_grid=false);

    /**
     * Meshes a spherical or mapped-spherical domain.
     */
    void build_sphere (Mesh& mesh,
		       const Real rad=1,
		       const unsigned int nr=2,
		       const ElemType type=INVALID_ELEM);

#ifdef HAVE_TRIANGLE

    /**
     * Meshes a square (2D) with a Delaunay triangulation.
     * This function internally calls the triangle library
     * written by J.R. Shewchuk.  Currently only works with
     * linear (Tri3) elements but could be extened to Tri6
     * as well.
     */
    void build_delaunay_square(Mesh& mesh,
			       const unsigned int nx, // num. of nodes in x-dir
			       const unsigned int ny, // num. of nodes in y-dir
			       const Real xmin=0., const Real xmax=1.,
			       const Real ymin=0., const Real ymax=1.,
			       const ElemType type=INVALID_ELEM);


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


#endif // HAVE_TRIANGLE
    
    namespace Private
    {
      /**
       * A useful inline function which replaces the #defines
       * used previously.  Not private since this is a namespace,
       * but would be if this were a class.  The first one returns
       * the proper node number for 2D elements while the second
       * one returns the node number for 3D elements.
       */
      inline
      unsigned int idx(const ElemType type,
		       const unsigned int nx,
		       const unsigned int i,
		       const unsigned int j)
      {
	switch(type)
	  {
	  case INVALID_ELEM:
	  case QUAD4:
	  case TRI3:
	    {
	      return i + j*(nx+1);
	      break;
	    }

	  case QUAD8:
	  case QUAD9:
	  case TRI6:
	    {
	      return i + j*(2*nx+1);
	      break;
	    }
	  
	  default:
	    {
	      std::cerr << "ERROR: Unrecognized 2D element type." << std::endl;
	      error();
	    }
	  }

	return libMesh::invalid_uint;
      }


    
      // Same as the function above, but for 3D elements
      inline
      unsigned int idx(const ElemType type,
		       const unsigned int nx,
		       const unsigned int ny,
		       const unsigned int i,
		       const unsigned int j,
		       const unsigned int k)
      {
	switch(type)
	  {
	  case INVALID_ELEM:
	  case HEX8:
	  case PRISM6:
	    {
	      return i + (nx+1)*(j + k*(ny+1));
	      break;
	    }

	  case HEX20:
	  case HEX27:
	  case TET4:  // TET4's are created from an initial HEX27 discretization
	  case TET10: // TET10's are created from an initial HEX27 discretization
	  case PRISM15:
	  case PRISM18:
	    {
	      return i + (2*nx+1)*(j + k*(2*ny+1));
	      break;
	    }
	  
	  default:
	    {
	      std::cerr << "ERROR: Unrecognized element type." << std::endl;
	      error();
	    }
	  }
      
	return libMesh::invalid_uint;
      }
    } // end namespace Meshtools::Generation::Private
  } // end namespace Meshtools::Generation
} // end namespace MeshTools


#endif // #define __mesh_generation_h__
