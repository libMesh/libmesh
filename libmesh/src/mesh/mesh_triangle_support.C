// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2007  Benjamin S. Kirk, John W. Peterson

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
#include "unstructured_mesh.h"
#include "face_tri3.h"
#include "face_tri6.h"
#include "mesh_generation.h"
#include "mesh_smoother_laplace.h"
#include "boundary_info.h"


#ifdef HAVE_TRIANGLE

// Triangulates a 2D rectangular region with or without holes
void MeshTools::Generation::build_delaunay_square(UnstructuredMesh& mesh,
						  const unsigned int nx, // num. of elements in x-dir
						  const unsigned int ny, // num. of elements in y-dir
						  const Real xmin, const Real xmax,
						  const Real ymin, const Real ymax,
						  const ElemType type,
						  const std::vector<TriangleInterface::Hole*>* holes)
{
  // Check for existing nodes and compatible dimension.
  assert (mesh.mesh_dimension() == 2);
  assert (nx >= 1); // need at least 1 element in x-direction
  assert (ny >= 1); // need at least 1 element in y-direction
  assert (xmin < xmax);
  assert (ymin < ymax);

  // Clear out any data which may have been in the Mesh
  mesh.clear();
  
  // The x and y spacing between boundary points
  const Real delta_x = (xmax-xmin) / static_cast<Real>(nx);
  const Real delta_y = (ymax-ymin) / static_cast<Real>(ny);  

  // Bottom
  for (unsigned int p=0; p<=nx; ++p)
    mesh.add_point(Point(xmin + p*delta_x, ymin));
  
  // Right side
  for (unsigned int p=1; p<ny; ++p)
    mesh.add_point(Point(xmax, ymin + p*delta_y)); 

  // Top
  for (unsigned int p=0; p<=nx; ++p)
    mesh.add_point(Point(xmax - p*delta_x, ymax));

  // Left side
  for (unsigned int p=1; p<ny; ++p)
    mesh.add_point(Point(xmin,  ymax - p*delta_y)); 
		   
  // Be sure we added as many points as we thought we did
  assert (mesh.n_nodes() == 2*(nx+ny));

  // Construct the Triangle Interface object
  TriangleInterface t(mesh);

  // Set custom variables for the triangulation
  t.desired_area()       = 0.5 * (xmax-xmin)*(ymax-ymin) / static_cast<Real>(nx*ny);
  t.triangulation_type() = TriangleInterface::PSLG;
  t.elem_type()          = type;

  if (holes != NULL)
    t.attach_hole_list(holes);

  // Triangulate!
  t.triangulate();

  // The mesh is now generated, but we still need to mark the boundaries
  // to be consistent with the other build_square routines.  Note that only
  // the external boundaries of the square are given boundary ids, the holes
  // are not numbered.
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
	    // hole = 4
	    short int bc_id=4;

	    // bottom
	    if      (std::fabs(side_midpoint(1) - ymin) < TOLERANCE)
	      bc_id=0;

	    // right
	    else if (std::fabs(side_midpoint(0) - xmax) < TOLERANCE)
	      bc_id=1;

	    // top
	    else if (std::fabs(side_midpoint(1) - ymax) < TOLERANCE)
	      bc_id=2;

	    // left
	    else if (std::fabs(side_midpoint(0) - xmin) < TOLERANCE)
	      bc_id=3;

	    // If the point is not on any of the external boundaries, it
	    // is on one of the holes....
	    
	    // Finally, add this element's information to the boundary info object.
	    mesh.boundary_info->add_side(elem->id(), s, bc_id);
	  }
    }
  
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
				UnstructuredMesh& mesh_output,
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






// Function definitions for the TriangleInterface class
void TriangleInterface::triangulate()
{
  // Will the triangulation have holes?
  const bool have_holes = ((_holes != NULL) && (!_holes->empty()));
  
  // Currently, we don't support specifying user-defined segments
  // *and* having holes
  if (have_holes)
    assert (this->segments.empty());
  
  // If the initial PSLG is really simple, e.g. an L-shaped domain or
  // a square/rectangle, the resulting triangulation may be very
  // "structured" looking.  Sometimes this is a problem if your
  // intention is to work with an "unstructured" looking grid.  We can
  // attempt to work around this limitation by inserting midpoints
  // into the original PSLG.  Inserting additional points into a
  // set of points meant to be a convex hull usually makes less sense.

  // May or may not need to insert new points ...
  if ((_triangulation_type==PSLG) && (_insert_extra_points))
    {
      // Make a copy of the original points from the Mesh
      std::vector<Point> original_points (_mesh.n_nodes());
      
      MeshBase::node_iterator       node_it  = _mesh.nodes_begin();
      const MeshBase::node_iterator node_end = _mesh.nodes_end();
      
      for (unsigned int ctr=0; node_it != node_end; ++node_it)
	original_points[ctr++] = **node_it;
      
      // Clear out the mesh
      _mesh.clear();
      
      // Insert a new point on each PSLG at some random location
      // np=index into new points vector
      // n =index into original points vector
      for (unsigned int np=0, n=0; np<2*original_points.size(); ++np)
	{
	  // the even entries are the original points
	  if (np%2==0)
	    _mesh.add_point(original_points[n++]);

	  else // the odd entries are the midpoints of the original PSLG segments
	    _mesh.add_point (0.5*(original_points[n] + original_points[n-1]));
	}
    }
  
  // Regardless of whether we added additional points, the set of points to
  // triangulate is now sitting in the mesh.

  // If the holes vector is non-NULL (and non-empty) we need to determine
  // the number of additional points which the holes will add to the
  // triangulation.
  unsigned int n_hole_points = 0;

  if (have_holes)
    {
      for (unsigned int i=0; i<_holes->size(); ++i)
	n_hole_points += (*_holes)[i]->n_points();
    }
  
  // Triangle data structure for the mesh
  Triangle::triangulateio initial;
  Triangle::triangulateio final;

  // Pseudo-Constructor for the triangle io structs
  Triangle::init(initial);
  Triangle::init(final);
    
  initial.numberofpoints = _mesh.n_nodes() + n_hole_points;
  initial.pointlist      = static_cast<REAL*>(std::malloc(initial.numberofpoints * 2 * sizeof(REAL)));

  if (_triangulation_type==PSLG)
    {
      // Implicit segment ordering: One segment per point, including hole points
      if (this->segments.empty())
	initial.numberofsegments = initial.numberofpoints; 

      // User-defined segment ordering: One segment per entry in the segments vector
      else
	initial.numberofsegments = this->segments.size();
    }
  
  else if (_triangulation_type==GENERATE_CONVEX_HULL)
    initial.numberofsegments = n_hole_points; // One segment for each hole point

  // Debugging
  // std::cout << "Number of segments set to: " << initial.numberofsegments << std::endl;
  
  // Allocate space for the segments (2 int per segment)
  if (initial.numberofsegments > 0)
    {
      initial.segmentlist = static_cast<int*> (std::malloc(initial.numberofsegments * 2 * sizeof(int)));
    }  


  // Copy all the holes' points and segments into the triangle struct.

  // The hole_offset is a constant offset into the points vector which points
  // past the end of the last hole point added.
  unsigned int hole_offset=0;
  
  if (have_holes)
    for (unsigned int i=0; i<_holes->size(); ++i)
      {
	for (unsigned int ctr=0, h=0; h<(*_holes)[i]->n_points(); ctr+=2, ++h)
	  {
	    Point p = (*_holes)[i]->point(h);

	    const unsigned int index0 = 2*hole_offset+ctr;
	    const unsigned int index1 = 2*hole_offset+ctr+1;

	    // Save the x,y locations in the triangle struct.
	    initial.pointlist[index0] = p(0);
	    initial.pointlist[index1] = p(1);

	    // Set the points which define the segments
	    initial.segmentlist[index0] = hole_offset+h;
	    initial.segmentlist[index1] = (h==(*_holes)[i]->n_points()-1) ? hole_offset : hole_offset+h+1; // wrap around
	  }
	
	// Update the hole_offset for the next hole
	hole_offset += (*_holes)[i]->n_points();
      }

  
  // Copy all the non-hole points and segments into the triangle struct.
  for (unsigned int ctr=0, n=0; n<_mesh.n_nodes(); ctr+=2, ++n)
    {
      const unsigned int index0 = 2*hole_offset+ctr;
      const unsigned int index1 = 2*hole_offset+ctr+1;
      
      initial.pointlist[index0] = _mesh.point(n)(0);
      initial.pointlist[index1] = _mesh.point(n)(1);

      // If the user requested a PSLG, the non-hole points are also segments
      if (_triangulation_type==PSLG)
	{
	  // Use implicit ordering to define segments
	  if (this->segments.empty())
	    {
	      initial.segmentlist[index0] = hole_offset+n;
	      initial.segmentlist[index1] = (n==_mesh.n_nodes()-1) ? hole_offset : hole_offset+n+1; // wrap around
	    }
	}
    }

  
  // If the user provided it, use his ordering to define the segments
  for (unsigned int ctr=0, s=0; s<this->segments.size(); ctr+=2, ++s)
    {
      const unsigned int index0 = 2*hole_offset+ctr;
      const unsigned int index1 = 2*hole_offset+ctr+1;
      
      initial.segmentlist[index0] = hole_offset + this->segments[s].first;
      initial.segmentlist[index1] = hole_offset + this->segments[s].second;
    }



  // Tell the input struct about the holes
  if (have_holes)
    {
      initial.numberofholes = _holes->size();
      initial.holelist      = static_cast<REAL*>(std::malloc(initial.numberofholes * 2 * sizeof(REAL)));
      for (unsigned int i=0, ctr=0; i<_holes->size(); ++i, ctr+=2)
	{
	  Point inside_point = (*_holes)[i]->inside();
	  initial.holelist[ctr]   = inside_point(0);
	  initial.holelist[ctr+1] = inside_point(1);
	}
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
  // D ~ Conforming Delaunay: use this switch if you want all triangles
  //     in the mesh to be Delaunay, and not just constrained Delaunay
  // q ~  Quality mesh generation with no angles smaller than 20 degrees.
  //      An alternate minimum angle may be specified after the q
  // a ~ Imposes a maximum triangle area constraint.
  // -P  Suppresses the output .poly file. Saves disk space, but you lose the ability to maintain
  //     constraining segments on later refinements of the mesh.
  // Create the flag strings, depends on element type
std::ostringstream flags;

// Default flags always used
flags << "zBPQq";

  // Flags which are specific to the type of triangulation
  switch (_triangulation_type)
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
  switch (_elem_type)
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


  // If we do have holes and the user asked to GENERATE_CONVEX_HULL,
  // need to add the p flag so the triangulation respects those segments.
  if ((_triangulation_type==GENERATE_CONVEX_HULL) && (have_holes))
    flags << "p";
  
  // Finally, add the area constraint
  flags << "a" << std::fixed << _desired_area;

  // Refine the initial output to conform to the area constraint
  Triangle::triangulate(const_cast<char*>(flags.str().c_str()),
			&initial,
			&final,
			NULL); // voronoi ouput -- not used
  
  
  // Send the information computed by Triangle to the Mesh.
  Triangle::copy_tri_to_mesh(final,
			     _mesh,
			     _elem_type);
      
  // To the naked eye, a few smoothing iterations usually looks better,
  // so we do this by default unless the user says not to.
  if (this->_smooth_after_generating)
    LaplaceMeshSmoother(_mesh).smooth(2);

    
  // Clean up.    
  Triangle::destroy(initial,      Triangle::INPUT);
  Triangle::destroy(final,        Triangle::OUTPUT);
  
}







#endif // HAVE_TRIANGLE








