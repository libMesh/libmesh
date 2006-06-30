// $Id: mesh_triangle_support.C,v 1.9 2006-06-30 15:33:13 jwpeterson Exp $

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
#include "mesh_generation.h"
#include "mesh_smoother_laplace.h"



#ifdef HAVE_TRIANGLE


// Definition of the function from the MeshTools::Generation namespace
void MeshTools::Generation::build_delaunay_square (Mesh& mesh,
						   const unsigned int nx,
						   const unsigned int ny,
						   const Real xmin, const Real xmax,
						   const Real ymin, const Real ymax)
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
  Triangle::triangulate(const_cast<char*>("czBQ"), // gives the desired char* 
			&input,
			&intermediate,
			NULL);

  // Final refined Triangle object.
  final.pointlist    = static_cast<Real*>(NULL);
  final.trianglelist = static_cast<int* >(NULL);

  // Perform final triangulation with area constraint
  std::ostringstream flags;
  flags << "przBPQa" << std::fixed << desired_area;
  Triangle::triangulate(const_cast<char*>(flags.str().c_str()),
			&intermediate,
			&final,
			static_cast<Triangle::triangulateio*>(NULL));

  // Transfer the information into the LibMesh mesh.
  mesh.clear();
  
  // Node information
  for (int i=0, c=0; c<final.numberofpoints; i+=2, ++c)
    mesh.add_point( Point(final.pointlist[i],
			  final.pointlist[i+1]) );
  
  // Element information
  for (int i=0; i<final.numberoftriangles; ++i)
    {
      Elem* elem = mesh.add_elem (new Tri3);

      for (unsigned int n=0; n<3; ++n)
	elem->set_node(n) = mesh.node_ptr(final.trianglelist[i*3 + n]);
    }

  // Clean up triangle data structures
  Triangle::destroy(input,        Triangle::INPUT );
  Triangle::destroy(intermediate, Triangle::BOTH  );
  Triangle::destroy(final,        Triangle::OUTPUT);
  
  // Prepare mesh for usage.
  mesh.prepare_for_use();
  
  // Run mesh through a few Laplace smoothing steps.
  LaplaceMeshSmoother s (mesh);
  s.smooth(2);

}




void MeshTools::Generation::build_delaunay_square_with_hole(Mesh& mesh,
							    const std::vector<Hole>& holes,
							    const unsigned int nx, // num. of nodes in x-dir (approximate)
							    const unsigned int ny, // num. of nodes in y-dir (approximate)
							    const Real xmin, const Real xmax,
							    const Real ymin, const Real ymax)
{
    // Triangle data structure for the mesh
    Triangle::triangulateio
      final_input,
      final_output;

    // Pseudo-Constructor for the triangle io structs
    Triangle::init(final_input);
    Triangle::init(final_output);

    
    const unsigned int n_holes = holes.size();

    // Pre-allocate space for the points and segments (one segment per point).
    unsigned int n_hole_points=0;
    for (unsigned int i=0; i<n_holes; ++i)
      n_hole_points += holes[i].n_points;

    final_input.numberofpoints = n_hole_points + 4;     // 4 additional points make up the corners of the square
    final_input.pointlist      = static_cast<REAL*>(std::malloc(final_input.numberofpoints * 2 * sizeof(REAL)));

    final_input.numberofsegments = n_hole_points;
    final_input.segmentlist      = static_cast<int*>(std::malloc(final_input.numberofsegments * 2 * sizeof(int))); // 2 int per segment


    // Constant offset into the vector for each successive hole after the first.
    unsigned int offset=0;
    
    // For each hole, compute points and determine the segments, add them to the input struct
    for (unsigned int i=0; i<n_holes; ++i)
      {
	// The change in angle between each point.
	const Real dtheta = 2.0*libMesh::pi / static_cast<Real>(holes[i].n_points);
	// std::cout << "dtheta=" << dtheta << std::endl;
	
	// Generate the points for the holes
	for (unsigned int ctr=0, h=0; h<holes[i].n_points; ctr+=2, ++h)
	  {
	    const Real theta = static_cast<Real>(h)*dtheta;

	    const unsigned int index0 = 2*offset+ctr;
	    const unsigned int index1 = 2*offset+ctr+1;

	    // The points
	    final_input.pointlist[index0] = holes[i].center(0) + holes[i].radius*cos(theta); // x=r*cos(theta)
	    final_input.pointlist[index1] = holes[i].center(1) + holes[i].radius*sin(theta); // y=r*sin(theta)
	  }

	// Generate all the segments for this hole
	for (unsigned int ctr=0, h=0; h<holes[i].n_points; ctr+=2, ++h)
	  {

	    const unsigned int index0 = 2*offset+ctr;
	    const unsigned int index1 = 2*offset+ctr+1;

	    // The points
	    final_input.segmentlist[index0] = offset+h;
	    final_input.segmentlist[index1] = (h==holes[i].n_points-1) ? offset : offset+h+1; // wrap around
	  }

	// for (unsigned int h=0; h<2*final_input.numberofsegments; ++h)
	// std::cout << final_input.segmentlist[h] << std::endl;

	// Update the offset
	offset += holes[i].n_points;
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

    // Tell the input structure about the hole(s)
    final_input.numberofholes = n_holes;
    final_input.holelist      = static_cast<REAL*>(std::malloc(final_input.numberofholes * 2 * sizeof(REAL)));
    for (unsigned int i=0, ctr=0; i<n_holes; ++i, ctr+=2)
      {
	final_input.holelist[ctr]   = holes[i].center(0);
	final_input.holelist[ctr+1] = holes[i].center(1);
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
    flags << "pczBQDqa" << std::fixed << 0.5 * (xmax-xmin)*(ymax-ymin) / static_cast<Real>((nx-1)*(ny-1));
    Triangle::triangulate(const_cast<char*>(flags.str().c_str()),
			  &final_input,
			  &final_output,
			  NULL); // voronoi ouput -- not used

    
    // Transfer the information into the LibMesh mesh.
    // Note: this should be its own function!!
    mesh.clear();
  
    // Node information
    for (int i=0, c=0; c<final_output.numberofpoints; i+=2, ++c)
      mesh.add_point( Point(final_output.pointlist[i],
			    final_output.pointlist[i+1]) );
  
    // Element information
    for (int i=0; i<final_output.numberoftriangles; ++i)
      {
	Elem* elem = mesh.add_elem (new Tri3);
	
	for (unsigned int n=0; n<3; ++n)
	  elem->set_node(n) = mesh.node_ptr(final_output.trianglelist[i*3 + n]);
      }


    // here();
    
    // Prepare mesh for usage.
    mesh.prepare_for_use();

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
#endif // HAVE_TRIANGLE







// // Constructor
// TriangleMeshInterface::TriangleMeshInterface (Mesh& mesh) :
//   _mesh (mesh)
// {}

// void TriangleMeshInterface::triangulate()
// {
//   // Check for existing nodes and compatible dimension.
//   assert (_mesh.n_nodes() != 0);
//   assert (_mesh.mesh_dimension() == 2);
  
//   // Construct the triangle interface struct for the input.
//   // You have to be really careful with the memory management
//   // of these things.  That might be a good thing to abstract
//   // away in the future.  AFAIK, triangle will only allocate
//   // memory (through malloc) if it finds a NULL pointer for
//   // one of its variables.
//   // Otherwise, it assumes there is enough space.
//   Triangle::triangulateio input;
  
//   // Set the number of nodes and set x,y values in Triangle's
//   // data structure.
//   input.numberofpoints = _mesh.n_nodes();
//   input.numberofpointattributes = 0;
  
//   // Allocate memory using a local vector, and
//   // let the input node list point to the beginning of the xy vector.
//   std::vector<Real> xy(2*_mesh.n_nodes());
//   for (unsigned int i=0, c=0; c<_mesh.n_nodes(); i+=2, ++c)
//     {
//       const Point& p = _mesh.point(c);
//       xy[i]   = p(0);
//       xy[i+1] = p(1);
//     }
//   input.pointlist = &xy[0];


//   // Intermediate triangulation step.
//   Triangle::triangulateio intermediate;
//   intermediate.pointlist       = static_cast<Real*>(NULL);
//   intermediate.trianglelist    = static_cast<int* >(NULL);
//   intermediate.segmentlist     = static_cast<int* >(NULL);
//   intermediate.pointmarkerlist = static_cast<int* >(NULL);
  
//   // Perform initial triangulation.
//   Triangle::triangulate("czBQ",
// 			&input,
// 			&intermediate,
// 			static_cast<Triangle::triangulateio*>(NULL));



//   // Refinement triangulation step.
//   Triangle::triangulateio final;
//   final.pointlist    = static_cast<Real*>(NULL);
//   final.trianglelist = static_cast<int* >(NULL);

//   // Perform final triangulation with area constraint
//   Triangle::triangulate("przBPQa0.002",
// 			&intermediate,
// 			&final,
// 			static_cast<Triangle::triangulateio*>(NULL));

//   // Transfer the information into the LibMesh mesh.
//   _mesh.clear();
  
//   // Node information
//   for (int i=0, c=0; c<final.numberofpoints; i+=2, ++c)
//     _mesh.add_point( Point(final.pointlist[i],
// 			   final.pointlist[i+1]) );
  
//   // Element information
//   for (int i=0; i<final.numberoftriangles; ++i)
//     {
//       Elem* elem = _mesh.add_elem (new Tri3);

//       for (unsigned int n=0; n<3; ++n)
// 	elem->set_node(n) = _mesh.node_ptr(final.trianglelist[i*3 + n]);
//     }

//   // Free vectors in the ouput struct
//   free (intermediate.pointlist             );
//   free (intermediate.trianglelist          );
//   free (intermediate.segmentlist           );
//   free (intermediate.pointmarkerlist       );
//   free (final.pointlist               );
//   free (final.trianglelist            );
  
//   // Prepare mesh for usage.
//   _mesh.prepare_for_use();
  
// }

