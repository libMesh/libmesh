// $Id: mesh_generation.C,v 1.32 2005-01-28 21:29:50 benkirk Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2004  Benjamin S. Kirk, John W. Peterson
  
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



// C++ includes
#include <cmath> // for std::sqrt


// Local includes
#include "mesh_generation.h"
#include "mesh.h"
#include "elem.h"
// #include "mesh_refinement.h"
#include "edge_edge2.h"
#include "edge_edge3.h"
#include "face_tri3.h"
#include "face_tri6.h"
#include "face_quad4.h"
#include "face_quad8.h"
#include "face_quad9.h"
#include "cell_hex8.h"
#include "cell_hex20.h"
#include "cell_hex27.h"
#include "cell_prism6.h"
#include "cell_prism15.h"
#include "cell_prism18.h"
#include "libmesh_logging.h"


// ------------------------------------------------------------
// MeshTools::Generation function for mesh generation
void MeshTools::Generation::build_cube(Mesh& mesh,
				       const unsigned int nx,
				       const unsigned int ny,
				       const unsigned int nz,
				       const Real xmin, const Real xmax,
				       const Real ymin, const Real ymax,
				       const Real zmin, const Real zmax,
				       const ElemType type)
{
  START_LOG("build_cube()", "MeshTools::Generation");

  // Clear the mesh and start from scratch
  mesh.clear();
  
  switch (mesh.mesh_dimension())
    {

      
      //---------------------------------------------------------------------
      // Build a 1D line
    case 1:      
      {
	assert (nx != 0);

	mesh.reserve_nodes(nx+1);
	mesh.reserve_elem (nx);

	for (unsigned int i=0; i<=nx; i++)
	  mesh.add_point (Point(static_cast<Real>(i)/static_cast<Real>(nx), 
				0, 
				0));
	
	for (unsigned int i=0; i<nx; i++)
	  {
	    Elem* elem = mesh.add_elem (new Edge2);
	    elem->set_node(0) = mesh.node_ptr(i);
	    elem->set_node(1) = mesh.node_ptr(i+1);
	  }


	// Scale the nodal positions
	for (unsigned int p=0; p<mesh.n_nodes(); p++)
	  mesh.node(p)(0) = (mesh.node(p)(0))*(xmax-xmin) + xmin;
	
	break;
      }







      
      
      
      //---------------------------------------------------------------------
      // Build a 2D quadrilateral
    case 2:
      {
	assert (nx != 0);
	assert (ny != 0);

	if ((type == INVALID_ELEM) ||
	    (type == QUAD4))
	  {
	    // Build a QUAD4

	    mesh.reserve_nodes( (nx+1)*(ny+1) );
	    mesh.reserve_elem (nx*ny);
	    
	    for (unsigned int j=0; j<=ny; j++)
	      for (unsigned int i=0; i<=nx; i++)
		mesh.add_point (Point(static_cast<Real>(i)/static_cast<Real>(nx),
				      static_cast<Real>(j)/static_cast<Real>(ny),
				      0));
	    
#define G(i,j) ( (i) + (j)*(nx+1) )
	
	    for (unsigned int j=0; j<ny; j++)
	      for (unsigned int i=0; i<nx; i++)
		{
		  Elem* elem = mesh.add_elem(new Quad4);

		  elem->set_node(0) = mesh.node_ptr(G(i,j)    );
		  elem->set_node(1) = mesh.node_ptr(G(i+1,j)  );
		  elem->set_node(2) = mesh.node_ptr(G(i+1,j+1));
		  elem->set_node(3) = mesh.node_ptr(G(i,j+1)  );

		  if (j == 0)
		    mesh.boundary_info.add_side(elem, 0, 0);
		  
		  if (j == (ny-1))
		    mesh.boundary_info.add_side(elem, 2, 2);
		  
		  if (i == 0)
		    mesh.boundary_info.add_side(elem, 3, 3);
		  
		  if (i == (nx-1))
		    mesh.boundary_info.add_side(elem, 1, 1);
		}
	    
#undef G
	  }


	
	else if (type == TRI3)
	  {
	    // Build a TRI3
	
	    mesh.reserve_nodes( (nx+1)*(ny+1) );
	    mesh.reserve_elem (2*nx*ny);
	    
	    for (unsigned int j=0; j<=ny; j++)
	      for (unsigned int i=0; i<=nx; i++)
		mesh.add_point (Point(static_cast<Real>(i)/static_cast<Real>(nx),
				      static_cast<Real>(j)/static_cast<Real>(ny),
				      0));
	    
#define G(i,j) ( (i) + (j)*(nx+1) )
	
	    for (unsigned int j=0; j<ny; j++)
	      for (unsigned int i=0; i<nx; i++)
		{
		  Elem* elem = NULL;
		  
		  elem = mesh.add_elem(new Tri3);

		  elem->set_node(0) = mesh.node_ptr(G(i,j)    );
		  elem->set_node(1) = mesh.node_ptr(G(i+1,j)  );
		  elem->set_node(2) = mesh.node_ptr(G(i+1,j+1));
		  
		  elem = mesh.add_elem(new Tri3);

		  elem->set_node(0) = mesh.node_ptr(G(i,j)    );
		  elem->set_node(1) = mesh.node_ptr(G(i+1,j+1));
		  elem->set_node(2) = mesh.node_ptr(G(i,j+1)  );
		}
	    
#undef G
	  }


	
	else if ((type == QUAD8) ||
		 (type == QUAD9))
	  {
	    // Build a Quad8 or a QUAD9
	
	    mesh.reserve_nodes( (2*nx+1)*(2*ny+1) );
	    mesh.reserve_elem (nx*ny);
	    
	    for (unsigned int j=0; j<=(2*ny); j++)
	      for (unsigned int i=0; i<=(2*nx); i++)
		mesh.add_point (Point(static_cast<Real>(i)/static_cast<Real>(2*nx),
				      static_cast<Real>(j)/static_cast<Real>(2*ny),
				      0));
	    
#define G(i,j) ( (i) + (j)*(2*nx+1) )
	
	    for (unsigned int j=0; j<(2*ny); j += 2)
	      for (unsigned int i=0; i<(2*nx); i += 2)
		{
		  Elem* elem = (type == QUAD8) ?
		    mesh.add_elem(new Quad8) :
		    mesh.add_elem(new Quad9);
		  

		  elem->set_node(0) = mesh.node_ptr(G(i,j)    );
		  elem->set_node(1) = mesh.node_ptr(G(i+2,j)  );
		  elem->set_node(2) = mesh.node_ptr(G(i+2,j+2));
		  elem->set_node(3) = mesh.node_ptr(G(i,j+2)  );
		  elem->set_node(4) = mesh.node_ptr(G(i+1,j)  );
		  elem->set_node(5) = mesh.node_ptr(G(i+2,j+1));
		  elem->set_node(6) = mesh.node_ptr(G(i+1,j+2));
		  elem->set_node(7) = mesh.node_ptr(G(i,j+1)  );
		  if (type == QUAD9)
		    elem->set_node(8) = mesh.node_ptr(G(i+1,j+1));
		  

		  if (j == 0)
		    mesh.boundary_info.add_side(elem, 0, 0);
		  
		  if (j == 2*(ny-1))
		    mesh.boundary_info.add_side(elem, 2, 2);
		  
		  if (i == 0)
		    mesh.boundary_info.add_side(elem, 3, 3);
		  
		  if (i == 2*(nx-1))
		    mesh.boundary_info.add_side(elem, 1, 1);
		}
	    
#undef G
	  }

	else if (type == TRI6)
	  {
	    // Build a TRI6
	
	    mesh.reserve_nodes( (2*nx+1)*(2*ny+1) );
	    mesh.reserve_elem(2*nx*ny);
	    
	    for (unsigned int j=0; j<=(2*ny); j++)
	      for (unsigned int i=0; i<=(2*nx); i++)
		mesh.add_point (Point(static_cast<Real>(i)/static_cast<Real>(2*nx),
				      static_cast<Real>(j)/static_cast<Real>(2*ny),
				      0));
	    
#define G(i,j) ( (i) + (j)*(2*nx+1) )
	
	    for (unsigned int j=0; j<(2*ny); j += 2)
	      for (unsigned int i=0; i<(2*nx); i += 2)
		{
		  Elem* elem = NULL;
		  
		  elem = mesh.add_elem(new Tri6);

		  elem->set_node(0) = mesh.node_ptr(G(i,j)    );
		  elem->set_node(1) = mesh.node_ptr(G(i+2,j)  );
		  elem->set_node(2) = mesh.node_ptr(G(i+2,j+2));
		  elem->set_node(3) = mesh.node_ptr(G(i+1,j)  );
		  elem->set_node(4) = mesh.node_ptr(G(i+2,j+1));
		  elem->set_node(5) = mesh.node_ptr(G(i+1,j+1));

		  
		  elem = mesh.add_elem(new Tri6);

		  elem->set_node(0) = mesh.node_ptr(G(i,j)    );
		  elem->set_node(1) = mesh.node_ptr(G(i+2,j+2));
		  elem->set_node(2) = mesh.node_ptr(G(i,j+2)  );
		  elem->set_node(3) = mesh.node_ptr(G(i+1,j+1));
		  elem->set_node(4) = mesh.node_ptr(G(i+1,j+2));
		  elem->set_node(5) = mesh.node_ptr(G(i,j+1)  );
		}
	    
#undef G
	  }

	
	else
	  {
	    error();
	  }
	
	// Scale the nodal positions
	for (unsigned int p=0; p<mesh.n_nodes(); p++)
	  {
	    mesh.node(p)(0) = (mesh.node(p)(0))*(xmax-xmin) + xmin;
	    mesh.node(p)(1) = (mesh.node(p)(1))*(ymax-ymin) + ymin;
	  }

	
	break;
      }








      


      //---------------------------------------------------------------------
      // Build a 3D hexahedral
    case 3:
      {
	assert (nx != 0);
	assert (ny != 0);
	assert (nz != 0);


	if ((type == INVALID_ELEM) ||
	    (type == HEX8) ||
	    (type == PRISM6))
	  {
	    // Build a HEX8
	    mesh.reserve_nodes( (nx+1)*(ny+1)*(nz+1) );
	    
	    if (type == PRISM6)
	      mesh.reserve_elem(2*nx*ny*nz);
	    else
	      mesh.reserve_elem(nx*ny*nz);
	    
	    for (unsigned int k=0; k<=nz; k++)
	      for (unsigned int j=0; j<=ny; j++)
		for (unsigned int i=0; i<=nx; i++)
		  mesh.add_point(Point(static_cast<Real>(i)/static_cast<Real>(nx),
				       static_cast<Real>(j)/static_cast<Real>(ny),
				       static_cast<Real>(k)/static_cast<Real>(nz)));
	    
#define G(i,j,k) ( (i) + (nx+1)*((j) + (k)*(ny+1)) )
	    
	    for (unsigned int k=0; k<nz; k++)
	      for (unsigned int j=0; j<ny; j++)
		for (unsigned int i=0; i<nx; i++)
		  {
		    if ((type == INVALID_ELEM) ||
			(type == HEX8))
		      {
			Elem* elem = mesh.add_elem(new Hex8);
			
			elem->set_node(0) = mesh.node_ptr(G(i,j,k)      );
			elem->set_node(1) = mesh.node_ptr(G(i+1,j,k)    );
			elem->set_node(2) = mesh.node_ptr(G(i+1,j+1,k)  );
			elem->set_node(3) = mesh.node_ptr(G(i,j+1,k)    );
			elem->set_node(4) = mesh.node_ptr(G(i,j,k+1)    );
			elem->set_node(5) = mesh.node_ptr(G(i+1,j,k+1)  );
			elem->set_node(6) = mesh.node_ptr(G(i+1,j+1,k+1));
			elem->set_node(7) = mesh.node_ptr(G(i,j+1,k+1)  );
			
			if (k == 0)
			  mesh.boundary_info.add_side(elem, 0, 0);
			
			if (k == (nz-1))
			  mesh.boundary_info.add_side(elem, 5, 5);
			
			if (j == 0)
			  mesh.boundary_info.add_side(elem, 1, 1);

			if (j == (ny-1))
			  mesh.boundary_info.add_side(elem, 3, 3);
 			
			if (i == 0)
			  mesh.boundary_info.add_side(elem, 4, 4);
 			
			if (i == (nx-1))
			  mesh.boundary_info.add_side(elem, 2, 2);
		      }

		    else if (type == PRISM6)
		      {
			Elem* elem = NULL;
			elem = mesh.add_elem(new Prism6);

			elem->set_node(0) = mesh.node_ptr(G(i,j,k)      );
			elem->set_node(1) = mesh.node_ptr(G(i+1,j,k)    );
			elem->set_node(2) = mesh.node_ptr(G(i,j+1,k)    );
			elem->set_node(3) = mesh.node_ptr(G(i,j,k+1)    );
			elem->set_node(4) = mesh.node_ptr(G(i+1,j,k+1)  );
			elem->set_node(5) = mesh.node_ptr(G(i,j+1,k+1)  );

			elem = mesh.add_elem(new Prism6);

			elem->set_node(0) = mesh.node_ptr(G(i+1,j,k)    );
			elem->set_node(1) = mesh.node_ptr(G(i+1,j+1,k)  );
			elem->set_node(2) = mesh.node_ptr(G(i,j+1,k)    );
			elem->set_node(3) = mesh.node_ptr(G(i+1,j,k+1)  );
			elem->set_node(4) = mesh.node_ptr(G(i+1,j+1,k+1));
			elem->set_node(5) = mesh.node_ptr(G(i,j+1,k+1)  );
		      }

		    else
		      error();
		  }
#undef G

	  }
	

	else if ((type == HEX20) ||
		 (type == HEX27) ||
		 (type == PRISM15) ||
		 (type == PRISM18))
	  {
	    // Build a HEX20 or HEX27 or PRISM15 or PRISM18
	    mesh.reserve_nodes( (2*nx+1)*(2*ny+1)*(2*nz+1) );
	    
	    if ((type == HEX20) || (type == HEX27))
	      mesh.reserve_elem(nx*ny*nz);
	    else
	      mesh.reserve_elem(2*nx*ny*nz);
	    
	    for (unsigned int k=0; k<=(2*nz); k++)
	      for (unsigned int j=0; j<=(2*ny); j++)
		for (unsigned int i=0; i<=(2*nx); i++)
		  mesh.add_point(Point(static_cast<Real>(i)/static_cast<Real>(2*nx),
				       static_cast<Real>(j)/static_cast<Real>(2*ny),
				       static_cast<Real>(k)/static_cast<Real>(2*nz)));
	    
#define G(i,j,k) ( (i) + (2*nx+1)*((j) + (k)*(2*ny+1)) )
	    
	    for (unsigned int k=0; k<(2*nz); k += 2)
	      for (unsigned int j=0; j<(2*ny); j += 2)
		for (unsigned int i=0; i<(2*nx); i += 2)
		  {
		    if ((type == HEX20) ||
			(type == HEX27))
		      {
			Elem* elem = (type == HEX20) ?
			  mesh.add_elem(new Hex20) :
			  mesh.add_elem(new Hex27);
		    
			elem->set_node(0)  = mesh.node_ptr(G(i,  j,  k)  );
			elem->set_node(1)  = mesh.node_ptr(G(i+2,j,  k)  );
			elem->set_node(2)  = mesh.node_ptr(G(i+2,j+2,k)  );
			elem->set_node(3)  = mesh.node_ptr(G(i,  j+2,k)  );
			elem->set_node(4)  = mesh.node_ptr(G(i,  j,  k+2));
			elem->set_node(5)  = mesh.node_ptr(G(i+2,j,  k+2));
			elem->set_node(6)  = mesh.node_ptr(G(i+2,j+2,k+2));
			elem->set_node(7)  = mesh.node_ptr(G(i,  j+2,k+2));
			elem->set_node(8)  = mesh.node_ptr(G(i+1,j,  k)  );
			elem->set_node(9)  = mesh.node_ptr(G(i+2,j+1,k)  );
			elem->set_node(10) = mesh.node_ptr(G(i+1,j+2,k)  );
			elem->set_node(11) = mesh.node_ptr(G(i,  j+1,k)  );
			elem->set_node(12) = mesh.node_ptr(G(i,  j,  k+1));
			elem->set_node(13) = mesh.node_ptr(G(i+2,j,  k+1));
			elem->set_node(14) = mesh.node_ptr(G(i+2,j+2,k+1));
			elem->set_node(15) = mesh.node_ptr(G(i,  j+2,k+1));
			elem->set_node(16) = mesh.node_ptr(G(i+1,j,  k+2));
			elem->set_node(17) = mesh.node_ptr(G(i+2,j+1,k+2));
			elem->set_node(18) = mesh.node_ptr(G(i+1,j+2,k+2));
			elem->set_node(19) = mesh.node_ptr(G(i,  j+1,k+2));
			if (type == HEX27)
			  {
			    elem->set_node(20) = mesh.node_ptr(G(i+1,j+1,k)  );
			    elem->set_node(21) = mesh.node_ptr(G(i+1,j,  k+1));
			    elem->set_node(22) = mesh.node_ptr(G(i+2,j+1,k+1));
			    elem->set_node(23) = mesh.node_ptr(G(i+1,j+2,k+1));
			    elem->set_node(24) = mesh.node_ptr(G(i,  j+1,k+1));
			    elem->set_node(25) = mesh.node_ptr(G(i+1,j+1,k+2));
			    elem->set_node(26) = mesh.node_ptr(G(i+1,j+1,k+1));
			  }
			
			
			if (k == 0)
			  mesh.boundary_info.add_side(elem, 0, 0);
			
			if (k == 2*(nz-1))
			  mesh.boundary_info.add_side(elem, 5, 5);
			
			if (j == 0)
			  mesh.boundary_info.add_side(elem, 1, 1);

			if (j == 2*(ny-1))
			  mesh.boundary_info.add_side(elem, 3, 3);
 			
			if (i == 0)
			  mesh.boundary_info.add_side(elem, 4, 4);
 			
			if (i == 2*(nx-1))
			  mesh.boundary_info.add_side(elem, 2, 2);
		      }

		    else if ((type == PRISM15) ||
			     (type == PRISM18))
		      {
			Elem* elem = NULL;
			elem = ((type == PRISM15) ?
				mesh.add_elem(new Prism15) :
				mesh.add_elem(new Prism18));

			elem->set_node(0)  = mesh.node_ptr(G(i,  j,  k)  );
			elem->set_node(1)  = mesh.node_ptr(G(i+2,j,  k)  );
			elem->set_node(2)  = mesh.node_ptr(G(i,  j+2,k)  );
			elem->set_node(3)  = mesh.node_ptr(G(i,  j,  k+2));
			elem->set_node(4)  = mesh.node_ptr(G(i+2,j,  k+2));
			elem->set_node(5)  = mesh.node_ptr(G(i,  j+2,k+2));
			elem->set_node(6)  = mesh.node_ptr(G(i+1,j,  k)  );
			elem->set_node(7)  = mesh.node_ptr(G(i+1,j+1,k)  );
			elem->set_node(8)  = mesh.node_ptr(G(i,  j+1,k)  );
			elem->set_node(9)  = mesh.node_ptr(G(i,  j,  k+1));
			elem->set_node(10) = mesh.node_ptr(G(i+2,j,  k+1));
			elem->set_node(11) = mesh.node_ptr(G(i,  j+2,k+1));
			elem->set_node(12) = mesh.node_ptr(G(i+1,j,  k+2));
			elem->set_node(13) = mesh.node_ptr(G(i+1,j+1,k+2));
			elem->set_node(14) = mesh.node_ptr(G(i,  j+1,k+2));
			if (type == PRISM18)
			  {
			    elem->set_node(15) = mesh.node_ptr(G(i+1,j,  k+1));
			    elem->set_node(16) = mesh.node_ptr(G(i+1,j+1,k+1));
			    elem->set_node(17) = mesh.node_ptr(G(i,  j+1,k+1));
			  }

			
			elem = ((type == PRISM15) ?
				mesh.add_elem(new Prism15) :
				mesh.add_elem(new Prism18));
			
			elem->set_node(0)  = mesh.node_ptr(G(i+2,j,k)     );
			elem->set_node(1)  = mesh.node_ptr(G(i+2,j+2,k)   );
			elem->set_node(2)  = mesh.node_ptr(G(i,j+2,k)     );
			elem->set_node(3)  = mesh.node_ptr(G(i+2,j,k+2)   );
			elem->set_node(4)  = mesh.node_ptr(G(i+2,j+2,k+2) );
			elem->set_node(5)  = mesh.node_ptr(G(i,j+2,k+2)   );
			elem->set_node(6)  = mesh.node_ptr(G(i+2,j+1,k)  );
			elem->set_node(7)  = mesh.node_ptr(G(i+1,j+2,k)  );
			elem->set_node(8)  = mesh.node_ptr(G(i+1,j+1,k)  );
			elem->set_node(9)  = mesh.node_ptr(G(i+2,j,k+1)  );
			elem->set_node(10) = mesh.node_ptr(G(i+2,j+2,k+1));
			elem->set_node(11) = mesh.node_ptr(G(i,j+2,k+1)  );
			elem->set_node(12) = mesh.node_ptr(G(i+2,j+1,k+2));
			elem->set_node(13) = mesh.node_ptr(G(i+1,j+2,k+2));
			elem->set_node(14) = mesh.node_ptr(G(i+1,j+1,k+2));
			if (type == PRISM18)
			  {
			    elem->set_node(15)  = mesh.node_ptr(G(i+2,j+1,k+1));
			    elem->set_node(16)  = mesh.node_ptr(G(i+1,j+2,k+1));
			    elem->set_node(17)  = mesh.node_ptr(G(i+1,j+1,k+1));
			  }
		      }

		    else
		      error();
		  }
#undef G

	  }

	else
	  {
	    error();
	  }
	
	// Scale the nodal positions
	for (unsigned int p=0; p<mesh.n_nodes(); p++)
	  {
	    mesh.node(p)(0) = (mesh.node(p)(0))*(xmax-xmin) + xmin;
	    mesh.node(p)(1) = (mesh.node(p)(1))*(ymax-ymin) + ymin;
	    mesh.node(p)(2) = (mesh.node(p)(2))*(zmax-zmin) + zmin;
	  }

	
	break;
      }

    default:
      {
	error();
      }
    }  

  STOP_LOG("build_cube()", "MeshTools::Generation");


  
  // Done building the mesh.  Now prepare it for use.
  mesh.prepare_for_use ();  
}



void MeshTools::Generation::build_square (Mesh& mesh,
					  const unsigned int nx,
					  const unsigned int ny,
					  const Real xmin, const Real xmax,
					  const Real ymin, const Real ymax,
					  const ElemType type)
{
  // This method only makes sense in 2D!
  assert (mesh.mesh_dimension() == 2);

  // Call the build_cube() member to actually do the work for us.
  build_cube (mesh,
	      nx, ny, 0,
	      xmin, xmax,
	      ymin, ymax,
	      0., 0.,
	      type);
}



void MeshTools::Generation::build_sphere (Mesh& /*mesh*/,
					  const Real /*rad*/,
					  const unsigned int /*nr*/,
					  const ElemType /*type*/)
{
  error();
//   assert (mesh.mesh_dimension() != 1);

//   assert (rad > 0.);

//   START_LOG("build_sphere()", "MeshTools::Generation");

//   // Clear the mesh and start from scratch
//   mesh.clear();
  
//   // Sphere is centered at origin by default
//   const Point cent;

//   const Sphere sphere (cent, rad);


  
//   switch (mesh_dimension())
//     {

//       //-----------------------------------------------------------------
//       // Build a circle in two dimensions
//     case 2:
//       {

// #ifdef ENABLE_AMR

// 	const Real sqrt_2     = sqrt(2.);
//         const Real rad_2      = .25*rad;
//         const Real rad_sqrt_2 = rad/sqrt_2;
  

// 	// Linear elements
// 	if ((type == INVALID_ELEM) ||
// 	    (type == TRI3) ||
// 	    (type == QUAD4))
// 	  {

// 	    mesh.reserve_nodes(8);
	    
// 	    // Point 0
// 	    node_ptr(0) = Node::build(-rad_2,-rad_2, 0, 0).release();
	    
// 	    // Point 1
// 	    node_ptr(1) = Node::build( rad_2,-rad_2, 0, 1).release();
	    
// 	    // Point 2
// 	    node_ptr(2) = Node::build( rad_2, rad_2, 0, 2).release();

// 	    // Point 3
// 	    node_ptr(3) = Node::build(-rad_2, rad_2, 0, 3).release();
	    
// 	    // Point 4
// 	    node_ptr(4) = Node::build(-rad_sqrt_2,-rad_sqrt_2, 0, 4).release();
	    
// 	    // Point 5
// 	    node_ptr(5) = Node::build( rad_sqrt_2,-rad_sqrt_2, 0, 5).release();
	    
// 	    // Point 6
// 	    node_ptr(6) = Node::build( rad_sqrt_2, rad_sqrt_2, 0, 6).release();
	    
// 	    // Point 7
// 	    node_ptr(7) = Node::build(-rad_sqrt_2, rad_sqrt_2, 0, 7).release();

// 	    // Build the elements
// 	    mesh.reserve_elem(5);
	    
// 	    for (unsigned int e=0; e<n_elem(); e++)
// 	      _elements[e] = new Quad4;
	    
// 	    // Element 0
// 	    elem(0)->set_node(0) = mesh.node_ptr(0);
// 	    elem(0)->set_node(1) = mesh.node_ptr(1);
// 	    elem(0)->set_node(2) = mesh.node_ptr(2);
// 	    elem(0)->set_node(3) = mesh.node_ptr(3);
	    
// 	    // Element 1
// 	    elem(1)->set_node(0) = mesh.node_ptr(4);
// 	    elem(1)->set_node(1) = mesh.node_ptr(0);
// 	    elem(1)->set_node(2) = mesh.node_ptr(3);
// 	    elem(1)->set_node(3) = mesh.node_ptr(7);
	    
// 	    // Element 2
// 	    elem(2)->set_node(0) = mesh.node_ptr(4);
// 	    elem(2)->set_node(1) = mesh.node_ptr(5);
// 	    elem(2)->set_node(2) = mesh.node_ptr(1);
// 	    elem(2)->set_node(3) = mesh.node_ptr(0);
	    
// 	    // Element 3
// 	    elem(3)->set_node(0) = mesh.node_ptr(1);
// 	    elem(3)->set_node(1) = mesh.node_ptr(5);
// 	    elem(3)->set_node(2) = mesh.node_ptr(6);
// 	    elem(3)->set_node(3) = mesh.node_ptr(2);
	    
// 	    // Element 4
// 	    elem(4)->set_node(0) = mesh.node_ptr(3);
// 	    elem(4)->set_node(1) = mesh.node_ptr(2);
// 	    elem(4)->set_node(2) = mesh.node_ptr(6);
// 	    elem(4)->set_node(3) = mesh.node_ptr(7);
// 	  }


// 	// Quadratic elements
// 	else if ((type == QUAD8) ||
// 		 (type == QUAD9) ||
// 		 (type == TRI6))
// 	  {
	    
// 	    mesh.reserve_nodes(25);

// 	    // Point 0
// 	    node_ptr(0)  = Node::build(-rad_2,-rad_2, 0, 0).release();
	    
// 	    // Point 1
// 	    node_ptr(1)  = Node::build( rad_2,-rad_2, 0, 1).release();
	    
// 	    // Point 2
// 	    node_ptr(2)  = Node::build( rad_2, rad_2, 0, 2).release();

// 	    // Point 3
// 	    node_ptr(3)  = Node::build(-rad_2, rad_2, 0, 3).release();
	    
// 	    // Point 4
// 	    node_ptr(4)  = Node::build(-rad_sqrt_2,-rad_sqrt_2, 0, 4).release();
	     
// 	    // Point 5
// 	    node_ptr(5)  = Node::build( rad_sqrt_2,-rad_sqrt_2, 0, 5).release();
	    
// 	    // Point 6
// 	    node_ptr(6)  = Node::build( rad_sqrt_2, rad_sqrt_2, 0, 6).release();
	    
// 	    // Point 7
// 	    node_ptr(7)  = Node::build(-rad_sqrt_2, rad_sqrt_2, 0, 7).release();

// 	    // Point 8
// 	    node_ptr(8)  = Node::build((point(0) + point(1))/2., 8).release();

// 	    // Point 9
// 	    node_ptr(9)  = Node::build((point(1) + point(2))/2., 9).release();
	    
// 	    // Point 10
// 	    node_ptr(10) = Node::build((point(2) + point(3))/2., 10).release();
	    
// 	    // Point 11
// 	    node_ptr(11) = Node::build((point(0) + point(3))/2., 11).release();
	    
// 	    // Point 12
// 	    node_ptr(12) = Node::build((point(0) + point(2))/2., 12).release();

// 	    // Point 13
// 	    node_ptr(13) = Node::build( 0,-rad, 0, 13).release();

// 	    // Point 14
// 	    node_ptr(14) = Node::build( rad, 0, 0, 14).release();

// 	    // Point 15
// 	    node_ptr(15) = Node::build( 0, rad, 0, 15).release();

// 	    // Point 16
// 	    node_ptr(16) = Node::build(-rad, 0, 0, 16).release();

// 	    // Point 17
// 	    node_ptr(17) = Node::build((point(8) + point(13))/2., 17).release();

// 	    // Point 18
// 	    node_ptr(18) = Node::build((point(9) + point(14))/2., 18).release();

// 	    // Point 19
// 	    node_ptr(19) = Node::build((point(10) + point(15))/2., 19).release();
	    
// 	    // Point 20
// 	    node_ptr(20) = Node::build((point(11) + point(16))/2., 20).release();
	      
// 	    // Point 21
// 	    node_ptr(21) = Node::build((point(0) + point(4))/2., 21).release();

// 	    // Point 22
// 	    node_ptr(22) = Node::build((point(1) + point(5))/2., 22).release();

// 	    // Point 23
// 	    node_ptr(23) = Node::build((point(2) + point(6))/2., 23).release();

// 	    // Point 24
// 	    node_ptr(24) = Node::build((point(3) + point(7))/2., 24).release();

	    
// 	    // Build the elements
// 	    mesh.reserve_elem(5);
	    
// 	    if ((type == QUAD9) ||
// 		(type == TRI6))
// 	      for (unsigned int e=0; e<n_elem(); e++)
// 		_elements[e] = new Quad9;
// 	    else
// 	      for (unsigned int e=0; e<n_elem(); e++)
// 		_elements[e] = new Quad8;
		
// 	    // Element 0
// 	    elem(0)->set_node(0) = mesh.node_ptr(0);
// 	    elem(0)->set_node(1) = mesh.node_ptr(1);
// 	    elem(0)->set_node(2) = mesh.node_ptr(2);
// 	    elem(0)->set_node(3) = mesh.node_ptr(3);
// 	    elem(0)->set_node(4) = mesh.node_ptr(8);
// 	    elem(0)->set_node(5) = mesh.node_ptr(9);
// 	    elem(0)->set_node(6) = mesh.node_ptr(10);
// 	    elem(0)->set_node(7) = mesh.node_ptr(11);
// 	    if (type != QUAD8) elem(0)->set_node(8) = mesh.node_ptr(12);

// 	    // Element 1
// 	    elem(1)->set_node(0) = mesh.node_ptr(4);
// 	    elem(1)->set_node(1) = mesh.node_ptr(0);
// 	    elem(1)->set_node(2) = mesh.node_ptr(3);
// 	    elem(1)->set_node(3) = mesh.node_ptr(7);
// 	    elem(1)->set_node(4) = mesh.node_ptr(21);
// 	    elem(1)->set_node(5) = mesh.node_ptr(11);
// 	    elem(1)->set_node(6) = mesh.node_ptr(24);
// 	    elem(1)->set_node(7) = mesh.node_ptr(16);
// 	    if (type != QUAD8) elem(1)->set_node(8) = mesh.node_ptr(20);

// 	    // Element 2
// 	    elem(2)->set_node(0) = mesh.node_ptr(4);
// 	    elem(2)->set_node(1) = mesh.node_ptr(5);
// 	    elem(2)->set_node(2) = mesh.node_ptr(1);
// 	    elem(2)->set_node(3) = mesh.node_ptr(0);
// 	    elem(2)->set_node(4) = mesh.node_ptr(13);
// 	    elem(2)->set_node(5) = mesh.node_ptr(22);
// 	    elem(2)->set_node(6) = mesh.node_ptr(8);
// 	    elem(2)->set_node(7) = mesh.node_ptr(21);
// 	    if (type != QUAD8) elem(2)->set_node(8) = mesh.node_ptr(17);

// 	    // Element 3
// 	    elem(3)->set_node(0) = mesh.node_ptr(1);
// 	    elem(3)->set_node(1) = mesh.node_ptr(5);
// 	    elem(3)->set_node(2) = mesh.node_ptr(6);
// 	    elem(3)->set_node(3) = mesh.node_ptr(2);
// 	    elem(3)->set_node(4) = mesh.node_ptr(22);
// 	    elem(3)->set_node(5) = mesh.node_ptr(14);
// 	    elem(3)->set_node(6) = mesh.node_ptr(23);
// 	    elem(3)->set_node(7) = mesh.node_ptr(9);
// 	    if (type != QUAD8) elem(3)->set_node(8) = mesh.node_ptr(18);

// 	    // Element 4
// 	    elem(4)->set_node(0) = mesh.node_ptr(3);
// 	    elem(4)->set_node(1) = mesh.node_ptr(2);
// 	    elem(4)->set_node(2) = mesh.node_ptr(6);
// 	    elem(4)->set_node(3) = mesh.node_ptr(7);
// 	    elem(4)->set_node(4) = mesh.node_ptr(10);
// 	    elem(4)->set_node(5) = mesh.node_ptr(23);
// 	    elem(4)->set_node(6) = mesh.node_ptr(15);
// 	    elem(4)->set_node(7) = mesh.node_ptr(24);
// 	    if (type != QUAD8) elem(4)->set_node(8) = mesh.node_ptr(19);

// 	  }
	

// 	else
// 	  error();


// 	// Now we have the beginnings of a sphere.
// 	// Add some more elements
// 	for (unsigned int r=0; r<nr; r++)
// 	  {
// 	    MeshRefinement mesh_refinement (*this);
	    
// 	    mesh_refinement.uniformly_refine(1);
	    
// 	    for (unsigned int e=0; e<n_elem(); e++)
// 	      if (elem(e)->active())
// 		for (unsigned int s=0; s<elem(e)->n_sides(); s++)
// 		  if (elem(e)->neighbor(s) == NULL)
// 		    {
// 		      AutoPtr<Elem> side(elem(e)->build_side(s));
		      
// 		      for (unsigned int n=0; n<side->n_nodes(); n++)
// 			side->point(n) =
// 			  sphere.closest_point(side->point(n));
// 		    }
// 	  }

	
// 	// Copy only the active elements to the current mesh
// 	{
//  	  std::vector<Elem*> new_elements;

// 	  new_elements.reserve (n_active_elem());
	  
// 	  for (unsigned int e=0; e<n_elem(); e++)
// 	    {
// 	      // Add only the active elements.   We can't just copy
// 	      // the pointers, because that would break the tree.
// 	      // We'll allocate new elements instead.
// 	      if (elem(e)->active())
// 		{
// 		  new_elements.push_back(Elem::build(elem(e)->type()).release());

// 		  for (unsigned int n=0; n<elem(e)->n_nodes(); n++)
// 		    new_elements.back()->set_node(n) =
// 		      elem(e)->get_node(n);
// 		}
	      
// 	      delete _elements[e);
// 	      _elements[e] = NULL;
// 	    }
	  
// 	  // Copy the new elements
// 	  _elements = new_elements;
// 	}
	
// 	// Possibly convert all the elements to triangles
// 	if ((type == TRI6) ||
// 	    (type == TRI3))
// 	  all_tri();

// #else

// 	std::cout << "Building a circle in 2D only works with AMR." << std::endl;
// 	error();

// #endif

	
// 	break;
//       }


      
//       //-----------------------------------------------------------------
//       // Build a sphere in three dimensions
//     case 3:
//       {
// 	const unsigned int nx=nr, ny=nr, nz=nr;

// 	build_cube(nx,ny,nz,
// 		   -rad/sqrt(3.), rad/sqrt(3.),
// 		   -rad/sqrt(3.), rad/sqrt(3.),
// 		   -rad/sqrt(3.), rad/sqrt(3.),
// 		   type);
	
// 	if ((type == INVALID_ELEM) ||
// 	    (type == HEX8))
// 	  {
	    
// #define G(i,j,k) ( (i) + (nx+1)*((j) + (k)*(ny+1)) )

// 	    // Pop IJ boundary to the sphere
// 	    for (unsigned int i=0; i<=nx; i++)
// 	      for (unsigned int j=0; j<=ny; j++)
// 		{
// 		  unsigned int k = 0;
		
// 		  node(G(i,j,k)) =
// 		    sphere.closest_point( point(G(i,j,k)) );
		  
// 		  k = nz;
		
// 		  node(G(i,j,k)) =
// 		    sphere.closest_point( point(G(i,j,k)) );
// 		}
	    
// 	    // Pop JK boundary to the sphere
// 	    for (unsigned int j=0; j<=ny; j++)
// 	      for (unsigned int k=0; k<=nz; k++)
// 		{
// 		  unsigned int i = 0;
		
// 		  node(G(i,j,k)) =
// 		    sphere.closest_point( point(G(i,j,k)) );
		
// 		  i = nx;
		
// 		  node(G(i,j,k)) =
// 		    sphere.closest_point( point(G(i,j,k)) );
// 		}
	    
// 	    // Pop IK boundary to the sphere
// 	    for (unsigned int i=0; i<=nx; i++)
// 	      for (unsigned int k=0; k<=nz; k++)
// 		{
// 		  unsigned int j = 0;
		
// 		  node(G(i,j,k)) =
// 		    sphere.closest_point( point(G(i,j,k)) );
		
// 		  j = ny;
		
// 		  node(G(i,j,k)) =
// 		    sphere.closest_point( point(G(i,j,k)) );
// 		}
	    
// 	    // Handle internal nodes
// 	    for (unsigned int k=1; k<nz; k++)
// 	      for (unsigned int j=1; j<ny; j++)
// 		for (unsigned int i=1; i<nx; i++)
// 		  {
// 		    const Real xmin = point(G(0 ,j, k))(0);
// 		    const Real xmax = point(G(nx,j, k))(0);
// 		    const Real ymin = point(G(i, 0, k))(1);
// 		    const Real ymax = point(G(i,ny, k))(1);
// 		    const Real zmin = point(G(i, j, 0))(2);
// 		    const Real zmax = point(G(i, j,nz))(2);

// 		    node(G(i,j,k))(0) = xmin +
// 		      (xmax - xmin)*static_cast<Real>(i)/static_cast<Real>(nx);
		  
// 		    node(G(i,j,k))(1) = ymin +
// 		      (ymax - ymin)*static_cast<Real>(j)/static_cast<Real>(ny);
		  
// 		    node(G(i,j,k))(2) = zmin +
// 		      (zmax - zmin)*static_cast<Real>(k)/static_cast<Real>(nz);
// 		  }

// 	    // Do some smoothing steps.
// 	    for (unsigned int l=0; l<10; l++)
// 	      for (unsigned int k=1; k<nz; k++)
// 		for (unsigned int j=1; j<ny; j++)
// 		  for (unsigned int i=1; i<nx; i++)
// 		    node(G(i,j,k)) = (point(G(i-1,j,  k  )) +
// 				      point(G(i+1,j,  k  )) +
// 				      point(G(i,  j-1,k  )) +
// 				      point(G(i,  j+1,k  )) +
// 				      point(G(i,  j,  k-1)) +
// 				      point(G(i,  j,  k+1)))/6.;

// #undef G
// 	  }
	
// 	else if ((type == HEX20) ||
// 		 (type == HEX27))
// 	  {
	    
// #define G(i,j,k) ( (i) + (2*nx+1)*((j) + (k)*(2*ny+1)) )

// 	    // Pop IJ boundary to the sphere
// 	    for (unsigned int i=0; i<=2*nx; i++)
// 	      for (unsigned int j=0; j<=2*ny; j++)
// 		{
// 		  unsigned int k = 0;
		
// 		  node(G(i,j,k)) =
// 		    sphere.closest_point( point(G(i,j,k)) );
		  
// 		  k = 2*nz;
		
// 		  node(G(i,j,k)) =
// 		    sphere.closest_point( point(G(i,j,k)) );
// 		}
	    
// 	    // Pop JK boundary to the sphere
// 	    for (unsigned int j=0; j<=2*ny; j++)
// 	      for (unsigned int k=0; k<=2*nz; k++)
// 		{
// 		  unsigned int i = 0;
		
// 		  node(G(i,j,k)) =
// 		    sphere.closest_point( point(G(i,j,k)) );
		
// 		  i = 2*nx;
		
// 		  node(G(i,j,k)) =
// 		    sphere.closest_point( point(G(i,j,k)) );
// 		}
	    
// 	    // Pop IK boundary to the sphere
// 	    for (unsigned int i=0; i<=2*nx; i++)
// 	      for (unsigned int k=0; k<=2*nz; k++)
// 		{
// 		  unsigned int j = 0;
		
// 		  node(G(i,j,k)) =
// 		    sphere.closest_point( point(G(i,j,k)) );
		
// 		  j = 2*ny;
		
// 		  node(G(i,j,k)) =
// 		    sphere.closest_point( point(G(i,j,k)) );
// 		}
	    
// 	    // Handle internal nodes
// 	    for (unsigned int k=1; k<2*nz; k++)
// 	      for (unsigned int j=1; j<2*ny; j++)
// 		for (unsigned int i=1; i<2*nx; i++)
// 		  {
// 		    const Real xmin = point(G(0,   j, k))(0);
// 		    const Real xmax = point(G(2*nx,j, k))(0);
// 		    const Real ymin = point(G(i,   0, k))(1);
// 		    const Real ymax = point(G(i,2*ny, k))(1);
// 		    const Real zmin = point(G(i, j,   0))(2);
// 		    const Real zmax = point(G(i, j,2*nz))(2);

// 		    node(G(i,j,k))(0) = xmin +
// 		      (xmax - xmin)*static_cast<Real>(i)/static_cast<Real>(2*nx);
		  
// 		    node(G(i,j,k))(1) = ymin +
// 		      (ymax - ymin)*static_cast<Real>(j)/static_cast<Real>(2*ny);
		  
// 		    node(G(i,j,k))(2) = zmin +
// 		      (zmax - zmin)*static_cast<Real>(k)/static_cast<Real>(2*nz);
// 		  }

// 	    // Do some smoothing steps.
// 	    for (unsigned int l=0; l<10; l++)
// 	      for (unsigned int k=1; k<2*nz; k++)
// 		for (unsigned int j=1; j<2*ny; j++)
// 		  for (unsigned int i=1; i<2*nx; i++)
// 		    node(G(i,j,k)) = (point(G(i-1,j,  k  )) +
// 				      point(G(i+1,j,  k  )) +
// 				      point(G(i,  j-1,k  )) +
// 				      point(G(i,  j+1,k  )) +
// 				      point(G(i,  j,  k-1)) +
// 				      point(G(i,  j,  k+1)))/6.;

// #undef G
// 	  }
	
// 	else
// 	  error();



// 	break;
//       }


      
//     default:
//       error();
//     }

  
//   STOP_LOG("build_sphere()", "MeshTools::Generation");


  
//   // Done building the mesh.  Now prepare it for use.
//   mesh.prepare_for_use();
}
