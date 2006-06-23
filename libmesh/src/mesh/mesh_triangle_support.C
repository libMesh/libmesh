// $Id: mesh_triangle_support.C,v 1.8 2006-06-23 16:49:48 jwpeterson Exp $

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

