// $Id: mesh_triangle_support.C,v 1.13 2006-10-24 19:43:34 jwpeterson Exp $

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



#include "mesh_triangle_support.h"
#include "mesh.h"
#include "face_tri3.h"
#include "face_tri6.h"
#include "mesh_generation.h"
#include "mesh_smoother_laplace.h"
#include "boundary_info.h"


#ifdef HAVE_TRIANGLE

void MeshTools::Generation::triangulate(Mesh& mesh,
					const std::vector<Point>& points,
					const Real desired_area,
					const ElemType type,
					const MeshTools::Generation::TriangulationType tt,
					const bool insert_points)
{
  // If the initial PSLG is really simple, e.g. an L-shaped domain or
  // a square/rectangle, the resulting triangulation will be very
  // "structured" looking.  Sometimes this is a problem if your
  // intention is to work with an unstructured grid.  We can attempt
  // to work around this limitation by inserting midpoints into the
  // original PSLG.
  std::vector<Point> new_points = points;

  if ((tt==PSLG) && (insert_points))
    {
      new_points.resize  (2*points.size());
      
      // Insert a new point on each PSLG at some random location
      // np=index into new points vector
      // n =index into original points vector
      for (unsigned int np=0, n=0; np<new_points.size(); ++np)
	{
	  // the even entries are the original points
	  if (np%2==0)
	    {
	      new_points[np] = points[n];
	      n++;
	    }

	  else // the odd entries are the midpoints of the original PSLG segments
	    {
	      new_points[np] = 0.5*(points[n] + points[n-1]);
	    }
	}
    }

  
  // Triangle data structure for the mesh
  Triangle::triangulateio initial;
  Triangle::triangulateio final;

  // Pseudo-Constructor for the triangle io structs
  Triangle::init(initial);
  Triangle::init(final);
    
  initial.numberofpoints = new_points.size();
  initial.pointlist      = static_cast<REAL*>(std::malloc(initial.numberofpoints * 2 * sizeof(REAL)));

  if (tt==PSLG)
    {
      initial.numberofsegments = initial.numberofpoints; // n. of segments = n. of points
      initial.segmentlist      = static_cast<int*>(std::malloc(initial.numberofsegments * 2 * sizeof(int))); // 2 int per segment
    }
  
  // Copy all the point information into the triangle initial struct.
  for (unsigned int n=0, index=0; n<new_points.size(); ++n, index+=2)
    {
      initial.pointlist[index]   = new_points[n](0);
      initial.pointlist[index+1] = new_points[n](1);
    }

  // Generate the PSLG segments
  if (tt==PSLG)
    for (unsigned int n=0, index=0; n<new_points.size(); ++n, index+=2)
      {
	initial.segmentlist[index]   = n;
	initial.segmentlist[index+1] = (n==new_points.size()-1) ? 0 : n+1;
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
  // D ~ Conforming Delaunay: use this switch if you want all triangles in the mesh to be Delaunay, and not just constrained Delaunay
  // q ~  Quality mesh generation with no angles smaller than 20 degrees. An alternate minimum angle may be specified after the q
  // a ~ Imposes a maximum triangle area constraint.
  // -P  Suppresses the output .poly file. Saves disk space, but you lose the ability to maintain
  //     constraining segments  on later refinements of the mesh.
  // Create the flag strings, depends on element type
  std::ostringstream flags;

  // Default flags always used
  flags << "zBPQq";

  // Flags which are specific to the type of triangulation
  switch (tt)
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
      {
	error();
	break;
      }
      
    default:
      {
	error();
      }
    }


  // Flags specific to the type of element
  switch (type)
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
      {
	std::cerr << "ERROR: Unrecognized triangular element type." << std::endl;
	error();
      }
    }


  // Finally, add the area constraint
  flags << "a" << std::fixed << desired_area;

  // Refine the initial output to conform to the area constraint
  Triangle::triangulate(const_cast<char*>(flags.str().c_str()),
			&initial,
			&final,
			NULL); // voronoi ouput -- not used
  
  
  // Send the information computed by Triangle to the Mesh.
  Triangle::copy_tri_to_mesh(final,
			     mesh,
			     type);
      
  // To the naked eye, a few smoothing iterations usually looks better.
  LaplaceMeshSmoother(mesh).smooth(2);

    
  // Clean up.    
  Triangle::destroy(initial,      Triangle::INPUT);
  Triangle::destroy(final,        Triangle::OUTPUT);
}









// Definition of the function from the MeshTools::Generation namespace
void MeshTools::Generation::build_delaunay_square (Mesh& mesh,
						   const unsigned int nx, // n elem, x-direction
						   const unsigned int ny, // n elem, y-direction
						   const Real xmin, const Real xmax,
						   const Real ymin, const Real ymax,
						   const ElemType type)
{
  // Check for existing nodes and compatible dimension.
  assert (mesh.mesh_dimension() == 2);
  assert (nx >= 1); // need at least 1 element in x-direction
  assert (ny >= 1); // need at least 1 element in y-direction
  assert (xmin < xmax);
  assert (ymin < ymax);

  // Generate vector of points to be passed to MeshTools::Generation::triangulate
  std::vector<Point> points((nx+1)*(ny+1));

  // The x and y spacing between boundary points
  const Real delta_x = (xmax-xmin) / static_cast<Real>(nx);
  const Real delta_y = (ymax-ymin) / static_cast<Real>(ny);  

  // Top and Bottom
  unsigned int ctr=0;
  for (unsigned int p=0; p<=nx; ++p)
    {
      const Real x=xmin + p*delta_x;

      points[ctr++] = Point(x, ymin); // Bottom
      points[ctr++] = Point(x, ymax); // Top
    }

  // Left and Right Sides
  for (unsigned int p=1; p<ny; ++p)
    {
      const Real y = ymin + p*delta_y;

      points[ctr++] = Point(xmin, y); // Left
      points[ctr++] = Point(xmax, y); // Right
    }
 


  // Call MeshTools::Generation::triangulate
  MeshTools::Generation::triangulate(mesh,
				     points,
				     0.5 * (xmax-xmin)*(ymax-ymin) / static_cast<Real>(nx*ny), // desired area
				     type,
				     GENERATE_CONVEX_HULL, // type of triangulation to perform
				     false                 // do not insert any extra points
				     );

  // The mesh is now generated, but we still need to mark the boundaries
  // to be consistent with the other build_square routines.
  MeshBase::element_iterator       el     = mesh.elements_begin();
  const MeshBase::element_iterator end_el = mesh.elements_end();
  for ( ; el != end_el; ++el)
    {
      const Elem* elem = *el;
      
      for (unsigned int s=0; s<elem->n_sides(); s++)
	if (elem->neighbor(s) == NULL)
	  {
	    AutoPtr<Elem> side (elem->build_side(s));

	    // Check the location of the side's midpoint.  Since
	    // the square has straight sides, the midpoint is not
	    // on the corner and thus it is uniquely on one of the
	    // sides.
	    Point side_midpoint= 0.5*( (*side->get_node(0)) + (*side->get_node(1)) );

	    // The boundary ids are set following the same convention as Quad4 sides
	    // bottom = 0
	    // right  = 1
	    // top = 2
	    // left = 3
	    // airfoil = 4
	    short int bc_id=99;

	    // bottom
	    if      (fabs(side_midpoint(1) - ymin) < TOLERANCE)
	      bc_id=0;

	    // right
	    else if (fabs(side_midpoint(0) - xmax) < TOLERANCE)
	      bc_id=1;

	    // top
	    else if (fabs(side_midpoint(1) - ymax) < TOLERANCE)
	      bc_id=2;

	    // left
	    else if (fabs(side_midpoint(0) - xmin) < TOLERANCE)
	      bc_id=3;

	    assert (bc_id != 99);

	    // Finally, add this element's information to the boundary info object.
	    mesh.boundary_info->add_side(elem->id(), s, bc_id);
	  }
    }
  
}










void MeshTools::Generation::build_delaunay_square_with_hole(Mesh& mesh,
							    const std::vector<Hole*>& holes,
							    const unsigned int nx, // num. of nodes in x-dir (approximate)
							    const unsigned int ny, // num. of nodes in y-dir (approximate)
							    const Real xmin, const Real xmax,
							    const Real ymin, const Real ymax,
							    const ElemType type)
{
    // Triangle data structure for the mesh
    Triangle::triangulateio
      final_input,
      final_output;

    // Pseudo-Constructor for the triangle io structs
    Triangle::init(final_input);
    Triangle::init(final_output);
    
    const unsigned int n_holes = holes.size();
    assert (n_holes > 0);

    // Sanity checks for the holes
    for (unsigned int i=0; i<n_holes; ++i)
      {
	assert (holes[i] != NULL);
	assert (holes[i]->n_points() > 1);
      }
    
    // Pre-allocate space for the points and segments (one segment per point).
    unsigned int n_hole_points=0;
    for (unsigned int i=0; i<n_holes; ++i)
      n_hole_points += holes[i]->n_points();

    final_input.numberofpoints = n_hole_points + 4;     // 4 additional points make up the corners of the square
    final_input.pointlist      = static_cast<REAL*>(std::malloc(final_input.numberofpoints * 2 * sizeof(REAL)));

    final_input.numberofsegments = n_hole_points;
    final_input.segmentlist      = static_cast<int*>(std::malloc(final_input.numberofsegments * 2 * sizeof(int))); // 2 int per segment


    // Constant offset into the vector for each successive hole after the first.
    unsigned int offset=0;
    
    // For each hole, compute points and determine the segments, add them to the input struct
    for (unsigned int i=0; i<n_holes; ++i)
      {
	for (unsigned int ctr=0, h=0; h<holes[i]->n_points(); ctr+=2, ++h)
	  {
	    Point p = holes[i]->point(h);

	    const unsigned int index0 = 2*offset+ctr;
	    const unsigned int index1 = 2*offset+ctr+1;

	    // Save the x,y locations in the triangle struct.
	    final_input.pointlist[index0] = p(0);
	    final_input.pointlist[index1] = p(1);
	  }

	// Generate all the segments for this hole
	for (unsigned int ctr=0, h=0; h<holes[i]->n_points(); ctr+=2, ++h)
	  {

	    const unsigned int index0 = 2*offset+ctr;
	    const unsigned int index1 = 2*offset+ctr+1;

	    // The points
	    final_input.segmentlist[index0] = offset+h;
	    final_input.segmentlist[index1] = (h==holes[i]->n_points()-1) ? offset : offset+h+1; // wrap around
	  }

	// for (unsigned int h=0; h<2*final_input.numberofsegments; ++h)
	// std::cout << final_input.segmentlist[h] << std::endl;

	// Update the offset
	offset += holes[i]->n_points();
      }

    
    // Add the corner points
    unsigned int idx = 2*n_hole_points;

    final_input.pointlist[idx++] = xmin;
    final_input.pointlist[idx++] = ymin;
    
    final_input.pointlist[idx++] = xmax;
    final_input.pointlist[idx++] = ymin;
    
    final_input.pointlist[idx++] = xmax;
    final_input.pointlist[idx++] = ymax;
    
    final_input.pointlist[idx++] = xmin;
    final_input.pointlist[idx++] = ymax;

    // Tell the input structure about the hole(s) by giving it any point
    // which lies "inside" the hole (not necessarily the center).
    final_input.numberofholes = n_holes;
    final_input.holelist      = static_cast<REAL*>(std::malloc(final_input.numberofholes * 2 * sizeof(REAL)));
    for (unsigned int i=0, ctr=0; i<n_holes; ++i, ctr+=2)
      {
	Point inside_point = holes[i]->inside();
	final_input.holelist[ctr]   = inside_point(0);
	final_input.holelist[ctr+1] = inside_point(1);
      }
    
    // Perform the triangulation.
    // c ~ enclose convex hull with segments
    // z ~ use zero indexing
    // B ~ Suppresses boundary markers in the output
    // Q ~ run in "quiet" mode
    // p ~ Triangulates a Planar Straight Line Graph
    //     If the `p' switch is used, `segmentlist' must point to a list of     
    //     segments, `numberofsegments' must be properly set, and               
    //     `segmentmarkerlist' must either be set to NULL (in which case all    
    //     markers default to zero), or must point to a list of markers.
    // D ~ Conforming Delaunay: use this switch if you want all triangles in the mesh to be Delaunay, and not just constrained Delaunay
    // q ~  Quality mesh generation with no angles smaller than 20 degrees. An alternate minimum angle may be specified after the q
    // a ~ Imposes a maximum triangle area constraint.
    std::ostringstream flags;
    const Real desired_area = 0.5 * (xmax-xmin)*(ymax-ymin) / static_cast<Real>((nx-1)*(ny-1));
    //flags << "pczBQDqa" << std::fixed << 0.5 * (xmax-xmin)*(ymax-ymin) / static_cast<Real>((nx-1)*(ny-1));

    switch (type)
    {
      case TRI3:
        {
          flags << "pczBQDqa" << std::fixed << desired_area;
          break;
        }

      case TRI6:
        {
          flags << "pczBQDqo2a" << std::fixed << desired_area;
          break;
        }

      default:
        {
          std::cerr << "ERROR: Unrecognized triangular element type." << std::endl;
          error();
        }
    }

    Triangle::triangulate(const_cast<char*>(flags.str().c_str()),
			  &final_input,
			  &final_output,
			  NULL); // voronoi ouput -- not used

    
    // Send the information computed by Triangle to the Mesh.
    Triangle::copy_tri_to_mesh(final_output,
			       mesh,
			       type);
      
    // To the naked eye, a few smoothing iterations usually looks better.
    LaplaceMeshSmoother(mesh).smooth(2);

    
    // here();
    
    Triangle::destroy(final_input, Triangle::INPUT);

    // here();
    
    Triangle::destroy(final_output, Triangle::OUTPUT);

    // here();
}






// Init helper routine defined in the Triangle namespace
void Triangle::init(Triangle::triangulateio& t)
{
  t.pointlist                    = static_cast<REAL*>(NULL);
  t.pointattributelist           = static_cast<REAL*>(NULL);
  t.pointmarkerlist              = static_cast<int* >(NULL);
  t.numberofpoints               = 0 ;
  t.numberofpointattributes      = 0 ;                                   

  t.trianglelist                 = static_cast<int* >(NULL);
  t.triangleattributelist        = static_cast<REAL*>(NULL);
  t.trianglearealist             = static_cast<REAL*>(NULL);
  t.neighborlist                 = static_cast<int* >(NULL);
  t.numberoftriangles            = 0;
  t.numberofcorners              = 0;
  t.numberoftriangleattributes   = 0;
  
  t.segmentlist                  = static_cast<int* >(NULL);
  t.segmentmarkerlist            = static_cast<int* >(NULL);
  t.numberofsegments             = 0;

  t.holelist                     = static_cast<REAL*>(NULL);
  t.numberofholes                = 0;

  t.regionlist                   = static_cast<REAL*>(NULL);
  t.numberofregions              = 0;
  
  t.edgelist                     = static_cast<int* >(NULL);
  t.edgemarkerlist               = static_cast<int* >(NULL);
  t.normlist                     = static_cast<REAL*>(NULL);
  t.numberofedges                = 0;
}






// Destroy helper routine defined in the Triangle namespace
void Triangle::destroy(Triangle::triangulateio& t, Triangle::IO_Type io_type)
{
  std::free (t.pointlist                    );
  std::free (t.pointattributelist           );
  std::free (t.pointmarkerlist              );
  std::free (t.trianglelist                 );
  std::free (t.triangleattributelist        );
  std::free (t.trianglearealist             );
  std::free (t.neighborlist                 );
  std::free (t.segmentlist                  );
  std::free (t.segmentmarkerlist            );

  // Only attempt to free these when t was used as an input struct!
  if (io_type==INPUT)
    {
      std::free (t.holelist                     );
      std::free (t.regionlist                   );
    }
  
  std::free (t.edgelist                     );
  std::free (t.edgemarkerlist               );
  std::free (t.normlist                     );

  // Reset
  // Triangle::init(t);
}






void Triangle::copy_tri_to_mesh(const triangulateio& triangle_data_input,
				Mesh& mesh_output,
				const ElemType type)
{
  // Transfer the information into the LibMesh mesh.
  mesh_output.clear();
  
  // Node information
  for (int i=0, c=0; c<triangle_data_input.numberofpoints; i+=2, ++c)
    {
      mesh_output.add_point( Point(triangle_data_input.pointlist[i],
				   triangle_data_input.pointlist[i+1]) );
    }

  // Element information
  for (int i=0; i<triangle_data_input.numberoftriangles; ++i)
    {
      switch (type)
	{
        case TRI3:
	  {
	    Elem* elem = mesh_output.add_elem (new Tri3);

	    for (unsigned int n=0; n<3; ++n)
	      elem->set_node(n) = mesh_output.node_ptr(triangle_data_input.trianglelist[i*3 + n]);

	    break;
	  }

        case TRI6:
	  {
	    Elem* elem = mesh_output.add_elem (new Tri6);

	    // Triangle number TRI6 nodes in a different way to libMesh
	    elem->set_node(0) = mesh_output.node_ptr(triangle_data_input.trianglelist[i*6 + 0]);
	    elem->set_node(1) = mesh_output.node_ptr(triangle_data_input.trianglelist[i*6 + 1]);
	    elem->set_node(2) = mesh_output.node_ptr(triangle_data_input.trianglelist[i*6 + 2]);
	    elem->set_node(3) = mesh_output.node_ptr(triangle_data_input.trianglelist[i*6 + 5]);
	    elem->set_node(4) = mesh_output.node_ptr(triangle_data_input.trianglelist[i*6 + 3]);
	    elem->set_node(5) = mesh_output.node_ptr(triangle_data_input.trianglelist[i*6 + 4]);

	    break;
	  }

        default:
	  {
	    std::cerr << "ERROR: Unrecognized triangular element type." << std::endl;
	    error();
	  }
	}
    }

  // Prepare mesh for usage.
  mesh_output.prepare_for_use();
}



#endif // HAVE_TRIANGLE








