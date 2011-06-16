/* $Id$ */

/* The Next Great Finite Element Library. */
/* Copyright (C) 2003  Benjamin S. Kirk */

/* This library is free software; you can redistribute it and/or */
/* modify it under the terms of the GNU Lesser General Public */
/* License as published by the Free Software Foundation; either */
/* version 2.1 of the License, or (at your option) any later version. */

/* This library is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU */
/* Lesser General Public License for more details. */

/* You should have received a copy of the GNU Lesser General Public */
/* License along with this library; if not, write to the Free Software */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */



 // <h1>Example 28 - Meshing with LibMesh's TetGen and Triangle Interfaces</h1>
 //
 // LibMesh provides interfaces to both Triangle and TetGen for generating 
 // Delaunay triangulations and tetrahedralizations in two and three dimensions
 // (respectively).

// Local header files
#include "mesh.h"
#include "mesh_triangle_support.h"
#include "mesh_generation.h"
#include "elem.h"
#include "mesh_tetgen_support.h"
#include "node.h"
#include "face_tri3.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

// Major functions called by main
void triangulate_domain();
void tetrahedralize_domain();

// Convenient data structure for holding/passing Mesh parameters
// for the tetrahedralize function.
struct MeshParams
{
  Real domain_x;
  Real domain_y;
  Real domain_z;
  Real hole_radius;
  Real hole_xpos;
  Real hole_ypos;
  Real hole_zpos;
};

// Helper routines for tetrahedralize_domain()
void add_sphere_surface_points_to_mesh(MeshBase& mesh, 
				       Real radius, 
				       unsigned refinements, 
				       Point sphere_center);

void add_boundary_tris(MeshBase& mesh, MeshParams& mesh_params);



// Begin the main program.
int main (int argc, char** argv)
{
  // Initialize libMesh and any dependent libaries, like in example 2.
  LibMeshInit init (argc, argv);

  // 1.) 2D triangulation of L-shaped domain with three holes of different shape
  triangulate_domain();
  
  // 2.) 3D tetrahedralization of rectangular domain with spherical hole
  tetrahedralize_domain();
  
  return 0;
}




void triangulate_domain()
{
  // Use typedefs for slightly less typing.
  typedef TriangleInterface::Hole Hole;
  typedef TriangleInterface::PolygonHole PolygonHole;
  typedef TriangleInterface::ArbitraryHole ArbitraryHole;

  // Libmesh mesh that will eventually be created.
  Mesh mesh(2);
    
  // The points which make up the L-shape:
  mesh.add_point(Point( 0. ,  0.));
  mesh.add_point(Point( 0. , -1.));
  mesh.add_point(Point(-1. , -1.));
  mesh.add_point(Point(-1. ,  1.));
  mesh.add_point(Point( 1. ,  1.));
  mesh.add_point(Point( 1. ,  0.));

  // Declare the TriangleInterface object.  This is where
  // we can set parameters of the triangulation and where the
  // actual triangulate function lives.
  TriangleInterface t(mesh);

  // Customize the variables for the triangulation
  t.desired_area()       = .01;

  // A Planar Straight Line Graph (PSLG) is essentially a list
  // of segments which have to exist in the final triangulation.
  // For an L-shaped domain, Triangle will compute the convex
  // hull of boundary points if we do not specify the PSLG.
  // The PSLG algorithm is also required for triangulating domains
  // containing holes
  t.triangulation_type() = TriangleInterface::PSLG;

    
  // Define holes...
    
  // hole_1 is a circle (discretized by 50 points)
  PolygonHole hole_1(Point(-0.5,  0.5), // center
		     0.25,              // radius
		     50);               // n. points

  // hole_2 is itself a triangle
  PolygonHole hole_2(Point(0.5, 0.5),   // center
		     0.1,               // radius
		     3);                // n. points

  // hole_3 is an ellipse of 100 points which we define here
  Point ellipse_center(-0.5,  -0.5);
  const unsigned int n_ellipse_points=100;
  std::vector<Point> ellipse_points(n_ellipse_points);
  const Real
    dtheta = 2*libMesh::pi / static_cast<Real>(n_ellipse_points),
    a = .1,
    b = .2;

  for (unsigned int i=0; i<n_ellipse_points; ++i)
    ellipse_points[i]= Point(ellipse_center(0)+a*cos(i*dtheta),
			     ellipse_center(1)+b*sin(i*dtheta));
    
  ArbitraryHole hole_3(ellipse_center, ellipse_points);
	
  // Create the vector of Hole*'s ...
  std::vector<Hole*> holes;
  holes.push_back(&hole_1);
  holes.push_back(&hole_2);
  holes.push_back(&hole_3);
	
  // ... and attach it to the triangulator object
  t.attach_hole_list(&holes);

  // Triangulate!
  t.triangulate();

  // Write the result to file
  mesh.write("delaunay_l_shaped_hole.e");
}



void tetrahedralize_domain()
{
  // The algorithm is broken up into N steps: we must generate a convex hull
  // of TRI3 surfaces defining the boundary before we can tetrahedralize its
  // interior.

  // 1.) build_sphere() generates the points on the surface of a sphere
  // 2.) Use Tetgen to generate the convex hull of the sphere's surface
  // 3.) Generate elements comprising the outer boundary of the domain.
  // 4.) Tetrahedralize the interior of the domain
  
  // The mesh we will eventually generate
  Mesh mesh(3);

  // Object to hold all the mesh parameters
  MeshParams mp;

  mp.domain_x = 1.;
  mp.domain_y = 1.;
  mp.domain_z = 1.;
  mp.hole_radius = .25;
  mp.hole_xpos = .5;
  mp.hole_ypos = .5;
  mp.hole_zpos = .5;
  
  // A Point object used to mark the interior of the spherical hole
  Point hole_center(mp.hole_xpos, 
		    mp.hole_ypos, 
		    mp.hole_zpos);
  
  // 1.) Use build_sphere() to generate surface points on the sphere
  add_sphere_surface_points_to_mesh(mesh, 
				    mp.hole_radius, 
				    /*refinements=*/3, 
				    hole_center);

  // 2.) Use Tetgen to generate the convex hull of the sphere's surface
  TetGenMeshInterface t(mesh);
  t.pointset_convexhull(); 
  
  // 3.) Generate elements comprising the outer boundary of the domain.
  add_boundary_tris(mesh, mp);

  // Also update neighbor information so that TetGen can verify there is
  // a convex hull.
  mesh.find_neighbors();

  // 4.) Tetrahedralize the interior of the domain with a hole cut out.
  Node hole_node(hole_center);
  std::vector<Node*> hole(1);
  hole[0] = &hole_node;

  // 0 means "use TetGen default value"
  Real quality_constraint = 0.;

  // The volume constraint determines the max-allowed tetrahedral
  // volume in the Mesh.  TetGen will split cells which are larger than
  // this size
  Real volume_constraint = .001;

  t.triangulate_conformingDelaunayMesh_carvehole(hole, 
						 quality_constraint, 
						 volume_constraint);
  
  // Find neighbors, etc in preparation for writing out the Mesh
  mesh.prepare_for_use();
  
  // Finally, write out the result
  mesh.write("sphere_hole.e");
}





void add_sphere_surface_points_to_mesh(MeshBase& mesh, 
				       Real radius, 
				       unsigned refinements, 
				       Point sphere_center)
{
  // n_refinements=3 leads to 386 nodes on the sphere
  // n_refinements=4 leads to 1538 nodes on the sphere
      
  // Note: to use build_sphere, you currently must disable the Laplace mesh smoothing
  // in that function, it currently dies with a weird error if we try to do it...
  Mesh sphere_mesh(3);
  MeshTools::Generation::build_sphere(sphere_mesh, radius, refinements, HEX8);
  //sphere_mesh.write("sphere_mesh.e");

  std::set<unsigned> boundary_node_ids;

  // Loop over elements, find sides on boundary, add node IDs to boundary node set
  {
    MeshBase::element_iterator it = sphere_mesh.elements_begin();
    const MeshBase::element_iterator end = sphere_mesh.elements_end();
    for ( ; it != end; ++it) 
      {
	Elem* elem = *it;
	  
	for (unsigned s=0; s<elem->n_sides(); ++s)
	  if (elem->neighbor(s) == NULL)
	    {
	      // Add the elements of this side to the set
	      AutoPtr<Elem> side = elem->side(s);
		
	      for (unsigned n=0; n<side->n_nodes(); ++n)
		boundary_node_ids.insert( side->node(n) );
	    }
      }

  }

  // For all the boundary_node_ids, add them to the real Mesh
  {
    MeshBase::node_iterator it = sphere_mesh.nodes_begin();
    const MeshBase::node_iterator end = sphere_mesh.nodes_end();
    for ( ; it != end; ++it) 
      {
	Node* node = *it;
	    
	if ( boundary_node_ids.find(node->id()) != boundary_node_ids.end() )
	  {
	    mesh.add_point (*node + sphere_center);

	    // Debugging
	    // std::cout << "Adding point: " << (*node + translation_vector);
	  }
      }
  }
}







void add_boundary_tris(MeshBase& mesh, MeshParams& mesh_params)
{
  // This function uses build_cube and tetgen to generate a convex hull of
  // the *outer* boundary, then adds those new points and elements to the
  // input mesh, thereby creating an outer boundary for the domain.
  Mesh outer_mesh(3);

  MeshTools::Generation::build_cube(outer_mesh,
				    1,1,1, // 1 element in each direction => 1 total element
				    0., mesh_params.domain_x,
				    0., mesh_params.domain_y,
				    0., mesh_params.domain_z,
				    HEX8);
  
  // The pointset_convexhull() algorithm will ignore the single Hex8
  // that's sitting in the Mesh, and just construct the triangulation
  // of the convex hull.
  TetGenMeshInterface t(outer_mesh);
  t.pointset_convexhull(); 
  
  // Now add all nodes from the boundary of the outer_mesh to the input mesh.
  // There should only be 8!  Store the IDs they are assigned by the input mesh.
  std::vector<unsigned> new_node_ids;
  new_node_ids.reserve(8);
  {
    MeshBase::node_iterator       node_it  = outer_mesh.nodes_begin();
    const MeshBase::node_iterator node_end = outer_mesh.nodes_end();
    for (; node_it != node_end; ++node_it)
      {
	// Grab a pointer to the node
	Node* old_node = *node_it;

	// Check that old_node is sequentially numbered?
	if (old_node->id() != new_node_ids.size())
	  {
	    libMesh::err << "Error, node IDs from build cube are not sequential!" << std::endl;
	    libmesh_error();
	  }

	// Add it to the input mesh, and catch the returned pointer
	Node* new_node = mesh.add_point (*old_node);

	// Store that ID so we can get a pointer to it later.
	new_node_ids.push_back(new_node->id());
      }
  }

  // Check to be sure that we didn't add too many/few nodes
  if (new_node_ids.size() != 8)
    libMesh::err << "Too many nodes found in outer_mesh!" << std::endl;
  
  // For each TRI3 element in the outer_mesh, add it to the input Mesh 
  // with proper node assignments
  {
    MeshBase::element_iterator       el     = outer_mesh.elements_begin();
    const MeshBase::element_iterator end_el = outer_mesh.elements_end();
    
    for (; el != end_el; ++el)
      {
	Elem* old_elem = *el;

	if (old_elem->type() == TRI3)
	  {
	    Elem* new_elem = mesh.add_elem(new Tri3);

	    // Assign nodes in new elements.  Since this is an example,
	    // we'll do it in several steps.
	    for (unsigned i=0; i<old_elem->n_nodes(); ++i)
	      {
		// Node ID in the temporary mesh
		unsigned old_node_id = old_elem->node(i);

		// Mapping to node ID in input mesh
		unsigned new_node_id = new_node_ids[old_node_id];

		// Node pointer assigned from input mesh
		new_elem->set_node(i) = mesh.node_ptr(new_node_id);
	      }
	  }
      }
  }
}
