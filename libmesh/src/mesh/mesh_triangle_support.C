// $Id: mesh_triangle_support.C,v 1.11 2006-09-24 05:22:29 jwpeterson Exp $

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



#ifdef HAVE_TRIANGLE


// Definition of the function from the MeshTools::Generation namespace
void MeshTools::Generation::build_delaunay_square (Mesh& mesh,
						   const unsigned int nx,
						   const unsigned int ny,
						   const Real xmin, const Real xmax,
						   const Real ymin, const Real ymax,
						   const ElemType type)
{
  // Check for existing nodes and compatible dimension.
  assert (mesh.mesh_dimension() == 2);
  assert (nx >= 2);
  assert (ny >= 2);
  assert (xmin < xmax);
  assert (ymin < ymax);

  // Declare and initialize the Triangle interface structs
  Triangle::triangulateio input, intermediate, final;
  Triangle::init(input);
  Triangle::init(intermediate);
  Triangle::init(final);

  
  // Compute the desired final area of the triangles, based on the
  // area of the domain and the requested number of nodes.
  const Real desired_area = 0.5 * (xmax-xmin)*(ymax-ymin) / static_cast<Real>((nx-1)*(ny-1));


  // Create the flag strings, depends on element type
  std::ostringstream flags_initial, flags_final;

  switch (type)
    {
      case TRI3:
        {
          flags_initial << "czBQ";
          flags_final << "przBPQa" << std::fixed << desired_area;
          break;
        }

      case TRI6:
        {
          flags_initial << "czBQo2";
          flags_final << "przBPQo2a" << std::fixed << desired_area;
          break;
        }

      default:
        {
          std::cerr << "ERROR: Unrecognized triangular element type." << std::endl;
          error();
        }
    }


  // Allocate memory for the initial points.  Stick to malloc here
  // so that all the accompanying 'destroy's work as well.
  input.numberofpoints = 2*(nx+ny-2);
  input.pointlist      = static_cast<REAL*>(std::malloc(input.numberofpoints * 2 * sizeof(REAL)));

  // The x and y spacing between boundary points
  const Real delta_x = (xmax-xmin) / static_cast<Real>(nx-1);
  const Real delta_y = (ymax-ymin) / static_cast<Real>(ny-1);  

  // Top and Bottom Sides
  for (unsigned int i=0, n=0; n<nx; i+=4, ++n)
    {
      // Bottom
      input.pointlist[i]   = xmin + n*delta_x; // x
      input.pointlist[i+1] = ymin;             // y

      // Top
      input.pointlist[i+2] = input.pointlist[i]; // x
      input.pointlist[i+3] = ymax;  // y
    }

  // Left and Right Sides
  for (unsigned int i=4*nx, n=1; n<ny-1; i+=4, ++n)
    {
      // Left
      input.pointlist[i]   = xmin;             // x
      input.pointlist[i+1] = ymin + n*delta_y; // y

      // Right
      input.pointlist[i+2] = xmax;     // x
      input.pointlist[i+3] = input.pointlist[i+1];  // y
    }


  // Note, if instead of putting a bunch of nodes on the boundary, you
  // start with just 4 at the corners, this always results in a
  // perfectly symmetric (boring) triangulation of the square, with
  // all similar triangles.

  // Perform initial triangulation.
  Triangle::triangulate(const_cast<char*>(flags_initial.str().c_str()), // gives the desired char* 
			&input,
			&intermediate,
			NULL);

  // Final refined Triangle object.
  final.pointlist    = static_cast<Real*>(NULL);
  final.trianglelist = static_cast<int* >(NULL);

  // Perform final triangulation with area constraint
  Triangle::triangulate(const_cast<char*>(flags_final.str().c_str()),
			&intermediate,
			&final,
			static_cast<Triangle::triangulateio*>(NULL));

  // Send the information computed by Triangle to the Mesh.
  Triangle::copy_tri_to_mesh(final,
			     mesh,
			     type);

  // Run mesh through a few Laplace smoothing steps.
  LaplaceMeshSmoother s (mesh);
  s.smooth(2);

  
  // Clean up triangle data structures
  Triangle::destroy(input,        Triangle::INPUT );
  Triangle::destroy(intermediate, Triangle::BOTH  );
  Triangle::destroy(final,        Triangle::OUTPUT);
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








