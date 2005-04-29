#include "mesh_triangle_support.h"
#include "mesh.h"
#include "elem.h"
#include "face_tri3.h"
#include "mesh_generation.h"
#include "mesh_smoother_laplace.h"

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
  
  // The desired final area of the triangles.
  const Real desired_area = 0.5 * (xmax-xmin)*(ymax-ymin) / static_cast<Real>((nx-1)*(ny-1));

  // Construct a vector of x,y point locations which define the
  // points on the boundary of the square.
  const unsigned int n_bndry_pts = 2*(nx+ny-2);
  std::vector<Real> xy(2*n_bndry_pts);

  // The x and y spacing between boundary points
  const Real delta_x = (xmax-xmin) / static_cast<Real>(nx-1);
  const Real delta_y = (ymax-ymin) / static_cast<Real>(ny-1);  

  // std::cout << "desired area = " << desired_area << std::endl;
  
  // Top and Bottom Sides
  for (unsigned int i=0, n=0; n<nx; i+=4, ++n)
    {
      // Bottom
      xy[i]   = xmin + n*delta_x; // x
      xy[i+1] = ymin;             // y

      // Top
      xy[i+2] = xy[i]; // x
      xy[i+3] = ymax;  // y
    }

  // Left and Right Sides
  for (unsigned int i=4*nx, n=1; n<ny-1; i+=4, ++n)
    {
      // Left
      xy[i]   = xmin;             // x
      xy[i+1] = ymin + n*delta_y; // y

      // Right
      xy[i+2] = xmax;     // x
      xy[i+3] = xy[i+1];  // y
    }



//   // Instead of putting a bunch of nodes on the boundary, start with
//   // just 4.  This always results in a uniform triangulation of the square,
//   // with all similar triangles.
//   std::vector<Real> xy(8);
//   xy[0] = xmin; xy[1] = ymin;
//   xy[2] = xmax; xy[3] = ymin;
//   xy[4] = xmin; xy[5] = ymax;
//   xy[6] = xmax; xy[7] = ymax;

  
  //for (unsigned int i=0; i<xy.size(); i+=2)
  //  std::cout << "(" << xy[i] << "," << xy[i+1] << ")" << std::endl;

  // The Triangle input object.
  Triangle::triangulateio input;
  input.numberofpoints          = xy.size() / 2;
  input.pointlist               = &xy[0];        // points to our local vector, no need to free later
  input.numberofpointattributes = 0;

  
  // Intermediate Trinagle object.
  Triangle::triangulateio intermediate;
  intermediate.pointlist       = static_cast<Real*>(NULL);
  intermediate.trianglelist    = static_cast<int* >(NULL);
  intermediate.segmentlist     = static_cast<int* >(NULL);
  intermediate.pointmarkerlist = static_cast<int* >(NULL);
  
  // Perform initial triangulation.
  Triangle::triangulate("czBQ",
			&input,
			&intermediate,
			static_cast<Triangle::triangulateio*>(NULL));

  // Final refined Triangle object.
  Triangle::triangulateio final;
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

  // Free vectors in the ouput struct
  free (intermediate.pointlist             );
  free (intermediate.trianglelist          );
  free (intermediate.segmentlist           );
  free (intermediate.pointmarkerlist       );
  free (final.pointlist               );
  free (final.trianglelist            );
  
  // Prepare mesh for usage.
  mesh.prepare_for_use();
  
  // Run mesh through a few Laplace smoothing steps.
  LaplaceMeshSmoother s (mesh);
  s.smooth(3);
}






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

