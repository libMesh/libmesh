// $Id: mesh.C,v 1.1.1.1 2003-01-10 16:17:48 libmesh Exp $

// The Next Great Finite Element Library.
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



// C++ includes
#include <iostream>
#include <algorithm>
#include <math.h>


// Local includes
#include "mesh.h"
#include "boundary_info.h"
#include "mesh_refinement.h"
#include "elem.h"
#include "cell.h"
#include "point.h"
#include "fe.h"
#include "sphere.h"

// Temporary 1D element includes
#include "edge_edge2.h"
#include "edge_edge3.h"

// Temporary 2D element includes
#include "face_tri3.h"
#include "face_tri6.h"
#include "face_quad4.h"
#include "face_quad8.h"
#include "face_quad9.h"

// Temporary 3D element includes
#include "cell_hex8.h"
#include "cell_hex20.h"
#include "cell_hex27.h"

#ifdef HAVE_SFCURVES
// prototype for SFC code
namespace sfc {
extern "C" {
#include "sfcurves.h"
}
}
#endif



// ------------------------------------------------------------
// Mesh class member functions
Mesh::Mesh(unsigned int d,
	   unsigned int pid) :
  MeshBase(d, pid),
#ifdef ENABLE_AMR
  boundary_info(d,*this),
  mesh_refinement(*this)
#else
  boundary_info(d,*this)
#endif
{
};



Mesh::~Mesh()
{
  boundary_info.clear();
  MeshBase::clear();
};



void Mesh::clear()
{
#ifdef ENABLE_AMR
  
  mesh_refinement.clear();
  
#endif
  
  boundary_info.clear();

  MeshBase::clear();
};



void Mesh::build_cube(const unsigned int nx,
		      const unsigned int ny,
		      const unsigned int nz,
		      const real xmin, const real xmax,
		      const real ymin, const real ymax,
		      const real zmin, const real zmax,
		      const ElemType type)
{
  switch (mesh_dimension())
    {

      
      //---------------------------------------------------------------------
      // Build a 1D line
    case 1:      
      {
	assert (nx != 0);
	
	_vertices.resize(nx+1);
	_elements.resize(nx);

	for (unsigned int i=0; i<=nx; i++)
	  {
	    _vertices[i](0) = (real) ((real) i)/((real) nx);
	  };
	
	for (unsigned int i=0; i<nx; i++)
	  {
	    _elements[i] = new Edge2;
	    elem(i)->node(0) = i;
	    elem(i)->node(1) = i+1;
	  };


	// Scale the nodal positions
	for (unsigned int p=0; p<n_nodes(); p++)
	  {
	    _vertices[p](0) = (_vertices[p](0))*(xmax-xmin) + xmin;
	  }
	
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
	
	    _vertices.resize( (nx+1)*(ny+1) );
	    _elements.resize(nx*ny);
	    
	    unsigned int p=0;
	    
	    for (unsigned int j=0; j<=ny; j++)
	      for (unsigned int i=0; i<=nx; i++)
		{
		  _vertices[p](0) = (real) ((real) i)/((real) nx);
		  _vertices[p](1) = (real) ((real) j)/((real) ny);

		  p++;
		};
	    
	    unsigned int e=0;
	    
#define G(i,j) ( (i) + (j)*(nx+1) )
	
	    for (unsigned int j=0; j<ny; j++)
	      for (unsigned int i=0; i<nx; i++)
		{
		  _elements[e] = new Quad4;

		  elem(e)->node(0) = G(i,j);
		  elem(e)->node(1) = G(i+1,j);
		  elem(e)->node(2) = G(i+1,j+1);
		  elem(e)->node(3) = G(i,j+1);

		  if (j == 0)
		    boundary_info.add_side(elem(e), 0, 0);
		  
		  else if (j == (ny-1))
		    boundary_info.add_side(elem(e), 2, 2);
		  
		  if (i == 0)
		    boundary_info.add_side(elem(e), 3, 3);
		  
		  else if (i == (nx-1))
		    boundary_info.add_side(elem(e), 1, 1);
		  
		  e++;
		};
	    
#undef G
	  }


	
	else if (type == TRI3)
	  {
	    // Build a TRI3
	
	    _vertices.resize( (nx+1)*(ny+1) );
	    _elements.resize(2*nx*ny);
	    
	    unsigned int p=0;
	    
	    for (unsigned int j=0; j<=ny; j++)
	      for (unsigned int i=0; i<=nx; i++)
		{
		  _vertices[p](0) = (real) ((real) i)/((real) nx);
		  _vertices[p](1) = (real) ((real) j)/((real) ny);

		  p++;
		};
	    
	    unsigned int e=0;
	    
#define G(i,j) ( (i) + (j)*(nx+1) )
	
	    for (unsigned int j=0; j<ny; j++)
	      for (unsigned int i=0; i<nx; i++)
		{
		  _elements[e] = new Tri3;

		  elem(e)->node(0) = G(i,j);
		  elem(e)->node(1) = G(i+1,j);
		  elem(e)->node(2) = G(i+1,j+1);
		  
		  e++;
		  
		  _elements[e] = new Tri3;

		  elem(e)->node(0) = G(i,j);
		  elem(e)->node(1) = G(i+1,j+1);
		  elem(e)->node(2) = G(i,j+1);
		  
		  e++;
		};
	    
#undef G
	  }


	
	else if ((type == QUAD8) ||
		 (type == QUAD9))
	  {
	    // Build a Quad8 or a QUAD9
	
	    _vertices.resize( (2*nx+1)*(2*ny+1) );
	    _elements.resize(nx*ny);
	    
	    unsigned int p=0;
	    
	    for (unsigned int j=0; j<=(2*ny); j++)
	      for (unsigned int i=0; i<=(2*nx); i++)
		{
		  _vertices[p](0) = (real) ((real) i)/((real) (2*nx));
		  _vertices[p](1) = (real) ((real) j)/((real) (2*ny));

		  p++;
		};
	    
	    unsigned int e=0;
	    
#define G(i,j) ( (i) + (j)*(2*nx+1) )
	
	    for (unsigned int j=0; j<(2*ny); j += 2)
	      for (unsigned int i=0; i<(2*nx); i += 2)
		{		
		  if (type == QUAD8)
		    _elements[e] = new Quad8;
		  else
		    _elements[e] = new Quad9;


		  elem(e)->node(0) = G(i,j);
		  elem(e)->node(1) = G(i+2,j);
		  elem(e)->node(2) = G(i+2,j+2);
		  elem(e)->node(3) = G(i,j+2);
		  elem(e)->node(4) = G(i+1,j);
		  elem(e)->node(5) = G(i+2,j+1);
		  elem(e)->node(6) = G(i+1,j+2);
		  elem(e)->node(7) = G(i,j+1);
		  if (type == QUAD9)
		    elem(e)->node(8) = G(i+1,j+1);
		  

		  if (j == 0)
		    boundary_info.add_side(elem(e), 0, 0);
		  
		  else if (j == (ny-1))
		    boundary_info.add_side(elem(e), 2, 2);
		  
		  if (i == 0)
		    boundary_info.add_side(elem(e), 3, 3);
		  
		  else if (i == (nx-1))
		    boundary_info.add_side(elem(e), 1, 1);
		  
		  e++;
		};
	    
#undef G
	  }

	else if (type == TRI6)
	  {
	    // Build a TRI6
	
	    _vertices.resize( (2*nx+1)*(2*ny+1) );
	    _elements.resize(2*nx*ny);
	    
	    unsigned int p=0;
	    
	    for (unsigned int j=0; j<=(2*ny); j++)
	      for (unsigned int i=0; i<=(2*nx); i++)
		{
		  _vertices[p](0) = (real) ((real) i)/((real) (2*nx));
		  _vertices[p](1) = (real) ((real) j)/((real) (2*ny));

		  p++;
		};
	    
	    unsigned int e=0;
	    
#define G(i,j) ( (i) + (j)*(2*nx+1) )
	
	    for (unsigned int j=0; j<(2*ny); j += 2)
	      for (unsigned int i=0; i<(2*nx); i += 2)
		{
		  _elements[e] = new Tri6;

		  elem(e)->node(0) = G(i,j);
		  elem(e)->node(1) = G(i+2,j);
		  elem(e)->node(2) = G(i+2,j+2);
		  elem(e)->node(3) = G(i+1,j);
		  elem(e)->node(4) = G(i+2,j+1);
		  elem(e)->node(5) = G(i+1,j+1);
		  
		  e++;
		  
		  _elements[e] = new Tri6;

		  elem(e)->node(0) = G(i,j);
		  elem(e)->node(1) = G(i+2,j+2);
		  elem(e)->node(2) = G(i,j+2);
		  elem(e)->node(3) = G(i+1,j+1);
		  elem(e)->node(4) = G(i+1,j+2);
		  elem(e)->node(5) = G(i,j+1);
		  
		  e++;
		  
		};
	    
#undef G
	  }

	
	else
	  {
	    error();
	  }
	
	// Scale the nodal positions
	for (unsigned int p=0; p<n_nodes(); p++)
	  {
	    _vertices[p](0) = (_vertices[p](0))*(xmax-xmin) + xmin;
	    _vertices[p](1) = (_vertices[p](1))*(ymax-ymin) + ymin;
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
	    (type == HEX8))
	  {
	    // Build a HEX8
	    _vertices.resize( (nx+1)*(ny+1)*(nz+1) );
	    _elements.resize(nx*ny*nz);
	    
	    unsigned int p=0;
	    
	    for (unsigned int k=0; k<=nz; k++)
	      for (unsigned int j=0; j<=ny; j++)
		for (unsigned int i=0; i<=nx; i++)
		  {
		    _vertices[p](0) = (real) ((real) i)/((real) nx);
		    _vertices[p](1) = (real) ((real) j)/((real) ny);
		    _vertices[p](2) = (real) ((real) k)/((real) nz);
		    
		    p++;
		  };
	    
	    unsigned int e=0;
	    
#define G(i,j,k) ( (i) + (nx+1)*((j) + (k)*(ny+1)) )
	    
	    for (unsigned int k=0; k<nz; k++)
	      for (unsigned int j=0; j<ny; j++)
		for (unsigned int i=0; i<nx; i++)
		  {
		    _elements[e] = new Hex8;
		    
		    elem(e)->node(0) = G(i,j,k);
		    elem(e)->node(1) = G(i+1,j,k);
		    elem(e)->node(2) = G(i+1,j+1,k);
		    elem(e)->node(3) = G(i,j+1,k);
		    elem(e)->node(4) = G(i,j,k+1);
		    elem(e)->node(5) = G(i+1,j,k+1);
		    elem(e)->node(6) = G(i+1,j+1,k+1);
		    elem(e)->node(7) = G(i,j+1,k+1);
		    
		    e++;
		  };
#undef G

	  }
	
	else if ((type == HEX20) ||
		 (type == HEX27))
	  {
	    // Build a HEX20 or HEX27
	    _vertices.resize( (2*nx+1)*(2*ny+1)*(2*nz+1) );
	    _elements.resize(nx*ny*nz);
	    
	    unsigned int p=0;
	    
	    for (unsigned int k=0; k<=(2*nz); k++)
	      for (unsigned int j=0; j<=(2*ny); j++)
		for (unsigned int i=0; i<=(2*nx); i++)
		  {
		    _vertices[p](0) = (real) ((real) i)/((real) (2*nx));
		    _vertices[p](1) = (real) ((real) j)/((real) (2*ny));
		    _vertices[p](2) = (real) ((real) k)/((real) (2*nz));
		    
		    p++;
		  };
	    
	    unsigned int e=0;
	    
#define G(i,j,k) ( (i) + (2*nx+1)*((j) + (k)*(2*ny+1)) )
	    
	    for (unsigned int k=0; k<(2*nz); k += 2)
	      for (unsigned int j=0; j<(2*ny); j += 2)
		for (unsigned int i=0; i<(2*nx); i += 2)
		  {
		    if (type == HEX20)
		      _elements[e] = new Hex20;
		    else
		      _elements[e] = new Hex27;
		    
		    elem(e)->node(0)  = G(i,  j,  k);
		    elem(e)->node(1)  = G(i+2,j,  k);
		    elem(e)->node(2)  = G(i+2,j+2,k);
		    elem(e)->node(3)  = G(i,  j+2,k);
		    elem(e)->node(4)  = G(i,  j,  k+2);
		    elem(e)->node(5)  = G(i+2,j,  k+2);
		    elem(e)->node(6)  = G(i+2,j+2,k+2);
		    elem(e)->node(7)  = G(i,  j+2,k+2);
		    elem(e)->node(8)  = G(i+1,j,  k);
		    elem(e)->node(9)  = G(i+2,j+1,k);
		    elem(e)->node(10) = G(i+1,j+2,k);
		    elem(e)->node(11) = G(i,  j+1,k);
		    elem(e)->node(12) = G(i,  j,  k+1);
		    elem(e)->node(13) = G(i+2,j,  k+1);
		    elem(e)->node(14) = G(i+2,j+2,k+1);
		    elem(e)->node(15) = G(i,  j+2,k+1);
		    elem(e)->node(16) = G(i+1,j,  k+2);
		    elem(e)->node(17) = G(i+2,j+1,k+2);
		    elem(e)->node(18) = G(i+1,j+2,k+2);
		    elem(e)->node(19) = G(i,  j+1,k+2);
		    if (type == HEX27)
		      {
			elem(e)->node(20) = G(i+1,j+1,k);
			elem(e)->node(21) = G(i+1,j,  k+1);
			elem(e)->node(22) = G(i+2,j+1,k+1);
			elem(e)->node(23) = G(i+1,j+2,k+1);
			elem(e)->node(24) = G(i,  j+1,k+1);
			elem(e)->node(25) = G(i+1,j+1,k+2);
			elem(e)->node(26) = G(i+1,j+1,k+1);
		      };
		    
		    e++;
		  };
#undef G

	  }

	else
	  {
	    error();
	  }
	
	// Scale the nodal positions
	for (unsigned int p=0; p<n_nodes(); p++)
	  {
	    _vertices[p](0) = (_vertices[p](0))*(xmax-xmin) + xmin;
	    _vertices[p](1) = (_vertices[p](1))*(ymax-ymin) + ymin;
	    _vertices[p](2) = (_vertices[p](2))*(zmax-zmin) + zmin;
	  }

	
	break;
      }








      
    default:
      {
	error();
      }
    };  
};



void Mesh::build_sphere (const real rad,
			 const unsigned int nr,
			 const ElemType type)
{
  assert (mesh_dimension() != 1);

  assert (rad > 0.);

  const Point cent;

  const Sphere sphere (cent, rad);


  
  switch (mesh_dimension())
    {

      //-----------------------------------------------------------------
      // Build a circle in two dimensions
    case 2:
      {

#ifdef ENABLE_AMR

	const real sqrt_2     = sqrt(2.);
        const real rad_2      = .25*rad;
        const real rad_sqrt_2 = rad/sqrt_2;
  

	// Linear elements
	if ((type == INVALID_ELEM) ||
	    (type == TRI3) ||
	    (type == QUAD4))
	  {

	    _vertices.resize(8);

	    // Point 0
	    _vertices[0](0) = -rad_2;
	    _vertices[0](1) = -rad_2;
	    
	    // Point 1
	    _vertices[1](0) =  rad_2;
	    _vertices[1](1) = -rad_2;
	    
	    // Point 2
	    _vertices[2](0) =  rad_2;
	    _vertices[2](1) =  rad_2;
	    
	    // Point 3
	    _vertices[3](0) = -rad_2;
	    _vertices[3](1) =  rad_2;
	    
	    // Point 4
	    _vertices[4](0) = -rad_sqrt_2;
	    _vertices[4](1) = -rad_sqrt_2;
	    
	    // Point 5
	    _vertices[5](0) =  rad_sqrt_2;
	    _vertices[5](1) = -rad_sqrt_2;
	    
	    // Point 6
	    _vertices[6](0) =  rad_sqrt_2;
	    _vertices[6](1) =  rad_sqrt_2;
	    
	    // Point 7
	    _vertices[7](0) = -rad_sqrt_2;
	    _vertices[7](1) =  rad_sqrt_2;

	    // Build the elements
	    _elements.resize(5);
	    
	    for (unsigned int e=0; e<n_elem(); e++)
	      _elements[e] = Elem::build(QUAD4);
	    
	    // Element 0
	    _elements[0]->node(0) = 0;
	    _elements[0]->node(1) = 1;
	    _elements[0]->node(2) = 2;
	    _elements[0]->node(3) = 3;
	    
	    // Element 1
	    _elements[1]->node(0) = 4;
	    _elements[1]->node(1) = 0;
	    _elements[1]->node(2) = 3;
	    _elements[1]->node(3) = 7;
	    
	    // Element 2
	    _elements[2]->node(0) = 4;
	    _elements[2]->node(1) = 5;
	    _elements[2]->node(2) = 1;
	    _elements[2]->node(3) = 0;
	    
	    // Element 3
	    _elements[3]->node(0) = 1;
	    _elements[3]->node(1) = 5;
	    _elements[3]->node(2) = 6;
	    _elements[3]->node(3) = 2;
	    
	    // Element 4
	    _elements[4]->node(0) = 3;
	    _elements[4]->node(1) = 2;
	    _elements[4]->node(2) = 6;
	    _elements[4]->node(3) = 7;
	  }


	// Quadratic elements
	else if ((type == QUAD8) ||
		 (type == QUAD9) ||
		 (type == TRI6))
	  {
	    
	    _vertices.resize(25);

	    // Point 0
	    _vertices[0](0) = -rad_2;
	    _vertices[0](1) = -rad_2;
	    
	    // Point 1
	    _vertices[1](0) =  rad_2;
	    _vertices[1](1) = -rad_2;
	    
	    // Point 2
	    _vertices[2](0) =  rad_2;
	    _vertices[2](1) =  rad_2;
	    
	    // Point 3
	    _vertices[3](0) = -rad_2;
	    _vertices[3](1) =  rad_2;
	    
	    // Point 4
	    _vertices[4](0) = -rad_sqrt_2;
	    _vertices[4](1) = -rad_sqrt_2;
	    
	    // Point 5
	    _vertices[5](0) =  rad_sqrt_2;
	    _vertices[5](1) = -rad_sqrt_2;
	    
	    // Point 6
	    _vertices[6](0) =  rad_sqrt_2;
	    _vertices[6](1) =  rad_sqrt_2;
	    
	    // Point 7
	    _vertices[7](0) = -rad_sqrt_2;
	    _vertices[7](1) =  rad_sqrt_2;

	    // Point 8
	    _vertices[8] = (_vertices[0] + _vertices[1])/2.;

	    // Point 9
	    _vertices[9] = (_vertices[1] + _vertices[2])/2.;
	    
	    // Point 10
	    _vertices[10] = (_vertices[2] + _vertices[3])/2.;
	    
	    // Point 11
	    _vertices[11] = (_vertices[0] + _vertices[3])/2.;
	    
	    // Point 12
	    _vertices[12] = (_vertices[0] + _vertices[2])/2.;

	    // Point 13
	    _vertices[13](0) =  0.;
	    _vertices[13](1) = -rad;

	    // Point 14
	    _vertices[14](0) =  rad;
	    _vertices[14](1) =  0;

	    // Point 15
	    _vertices[15](0) =  0.;
	    _vertices[15](1) =  rad;

	    // Point 16
	    _vertices[16](0) = -rad;
	    _vertices[16](1) =  0.;

	    // Point 17
	    _vertices[17] = (_vertices[8] + _vertices[13])/2.;

	    // Point 18
	    _vertices[18] = (_vertices[9] + _vertices[14])/2.;

	    // Point 19
	    _vertices[19] = (_vertices[10] + _vertices[15])/2.;

	    // Point 20
	    _vertices[20] = (_vertices[11] + _vertices[16])/2.;

	    // Point 21
	    _vertices[21] = (_vertices[0] + _vertices[4])/2.;

	    // Point 22
	    _vertices[22] = (_vertices[1] + _vertices[5])/2.;

	    // Point 23
	    _vertices[23] = (_vertices[2] + _vertices[6])/2.;

	    // Point 24
	    _vertices[24] = (_vertices[3] + _vertices[7])/2.;

	    
	    // Build the elements
	    _elements.resize(5);
	    
	    if ((type == QUAD9) ||
		(type == TRI6))
	      for (unsigned int e=0; e<n_elem(); e++)
		_elements[e] = Elem::build(QUAD9);
	    else
	      for (unsigned int e=0; e<n_elem(); e++)
		_elements[e] = Elem::build(QUAD8);
		
	    // Element 0
	    _elements[0]->node(0) = 0;
	    _elements[0]->node(1) = 1;
	    _elements[0]->node(2) = 2;
	    _elements[0]->node(3) = 3;
	    _elements[0]->node(4) = 8;
	    _elements[0]->node(5) = 9;
	    _elements[0]->node(6) = 10;
	    _elements[0]->node(7) = 11;
	    if (type != QUAD8) _elements[0]->node(8) = 12;

	    // Element 1
	    _elements[1]->node(0) = 4;
	    _elements[1]->node(1) = 0;
	    _elements[1]->node(2) = 3;
	    _elements[1]->node(3) = 7;
	    _elements[1]->node(4) = 21;
	    _elements[1]->node(5) = 11;
	    _elements[1]->node(6) = 24;
	    _elements[1]->node(7) = 16;
	    if (type != QUAD8) _elements[1]->node(8) = 20;

	    // Element 2
	    _elements[2]->node(0) = 4;
	    _elements[2]->node(1) = 5;
	    _elements[2]->node(2) = 1;
	    _elements[2]->node(3) = 0;
	    _elements[2]->node(4) = 13;
	    _elements[2]->node(5) = 22;
	    _elements[2]->node(6) = 8;
	    _elements[2]->node(7) = 21;
	    if (type != QUAD8) _elements[2]->node(8) = 17;

	    // Element 3
	    _elements[3]->node(0) = 1;
	    _elements[3]->node(1) = 5;
	    _elements[3]->node(2) = 6;
	    _elements[3]->node(3) = 2;
	    _elements[3]->node(4) = 22;
	    _elements[3]->node(5) = 14;
	    _elements[3]->node(6) = 23;
	    _elements[3]->node(7) = 9;
	    if (type != QUAD8) _elements[3]->node(8) = 18;

	    // Element 4
	    _elements[4]->node(0) = 3;
	    _elements[4]->node(1) = 2;
	    _elements[4]->node(2) = 6;
	    _elements[4]->node(3) = 7;
	    _elements[4]->node(4) = 10;
	    _elements[4]->node(5) = 23;
	    _elements[4]->node(6) = 15;
	    _elements[4]->node(7) = 24;
	    if (type != QUAD8) _elements[4]->node(8) = 19;

	  }
	

	else
	  error();


	// Now we have the beginnings of a sphere.
	// Add some more elements
	for (unsigned int r=0; r<nr; r++)
	  {
	    mesh_refinement.uniformly_refine(1);
	    
	    for (unsigned int e=0; e<n_elem(); e++)
	      if (elem(e)->active())
		for (unsigned int s=0; s<elem(e)->n_sides(); s++)
		  if (elem(e)->neighbor(s) == NULL)
		    {
#ifndef __IBMCPP__
		      const std::auto_ptr<Elem> side = elem(e)->build_side(s);
#else
		      const std::auto_ptr<Elem> side(elem(e)->build_side(s));
#endif		      
		      for (unsigned int n=0; n<side->n_nodes(); n++)
			    _vertices[side->node(n)] =
			      sphere.closest_point(vertex(side->node(n)));
		    };
	  };


	// Copy only the active elements to the current mesh
	{
	  mesh_refinement.clear();

	  std::vector<Elem*> new_elements(n_active_elem());
	  
	  unsigned int ne=0;
	  
	  for (unsigned int e=0; e<n_elem(); e++)
	    {
	      if (elem(e)->active())
		{
		  // Build the new elements
		  new_elements[ne] = Elem::build(elem(e)->type());
		  
		  // Copy the nodes
		  new_elements[ne]->_nodes = elem(e)->_nodes;
		  
		  ne++;
		};
	      
	      // Delete the old element
	      delete elem(e);
	    };
	  
	  // Copy the new elements
	  _elements = new_elements;
	};
	
	
	// Possibly convert all the elements to triangles
	if ((type == TRI6) ||
	    (type == TRI3))
	  all_tri();
	
	// Finally, find the face neighbors
	find_neighbors();



#else

	std::cout << "Building a circle in 2D only works with AMR." << std::endl;
	error();

#endif

	
	return;
      };


      
      //-----------------------------------------------------------------
      // Build a sphere in three dimensions
    case 3:
      {
	const unsigned int nx=nr, ny=nr, nz=nr;

	build_cube(nx,ny,nz,
		   -rad/sqrt(3.), rad/sqrt(3.),
		   -rad/sqrt(3.), rad/sqrt(3.),
		   -rad/sqrt(3.), rad/sqrt(3.),
		   type);
	
	if ((type == INVALID_ELEM) ||
	    (type == HEX8))
	  {
	    
#define G(i,j,k) ( (i) + (nx+1)*((j) + (k)*(ny+1)) )

	    // Pop IJ boundary to the sphere
	    for (unsigned int i=0; i<=nx; i++)
	      for (unsigned int j=0; j<=ny; j++)
		{
		  unsigned int k = 0;
		
		  _vertices[G(i,j,k)] =
		    sphere.closest_point( vertex(G(i,j,k)) );
		  
		  k = nz;
		
		_vertices[G(i,j,k)] =
		  sphere.closest_point( vertex(G(i,j,k)) );
	      };
	    
	    // Pop JK boundary to the sphere
	    for (unsigned int j=0; j<=ny; j++)
	      for (unsigned int k=0; k<=nz; k++)
		{
		  unsigned int i = 0;
		
		  _vertices[G(i,j,k)] =
		    sphere.closest_point( vertex(G(i,j,k)) );
		
		  i = nx;
		
		  _vertices[G(i,j,k)] =
		    sphere.closest_point( vertex(G(i,j,k)) );
	      };
	    
	    // Pop IK boundary to the sphere
	    for (unsigned int i=0; i<=nx; i++)
	      for (unsigned int k=0; k<=nz; k++)
		{
		  unsigned int j = 0;
		
		  _vertices[G(i,j,k)] =
		    sphere.closest_point( vertex(G(i,j,k)) );
		
		  j = ny;
		
		  _vertices[G(i,j,k)] =
		    sphere.closest_point( vertex(G(i,j,k)) );
	      };
	    
	    // Handle internal nodes
	    for (unsigned int k=1; k<nz; k++)
	      for (unsigned int j=1; j<ny; j++)
		for (unsigned int i=1; i<nx; i++)
		{
		  const real xmin = vertex(G(0 ,j, k))(0);
		  const real xmax = vertex(G(nx,j, k))(0);
		  const real ymin = vertex(G(i, 0, k))(1);
		  const real ymax = vertex(G(i,ny, k))(1);
		  const real zmin = vertex(G(i, j, 0))(2);
		  const real zmax = vertex(G(i, j,nz))(2);

		  _vertices[G(i,j,k)](0) = xmin +
		    (xmax - xmin)*((real) ((real) i)/((real) nx));
		  
		  _vertices[G(i,j,k)](1) = ymin +
		    (ymax - ymin)*((real) ((real) j)/((real) ny));
		  
		  _vertices[G(i,j,k)](2) = zmin +
		    (zmax - zmin)*((real) ((real) k)/((real) nz));
		};

	    // Do some smoothing steps.
	    for (unsigned int l=0; l<10; l++)
	      for (unsigned int k=1; k<nz; k++)
		for (unsigned int j=1; j<ny; j++)
		  for (unsigned int i=1; i<nx; i++)
		    _vertices[G(i,j,k)] = (vertex(G(i-1,j,  k  )) +
					  vertex(G(i+1,j,  k  )) +
					  vertex(G(i,  j-1,k  )) +
					  vertex(G(i,  j+1,k  )) +
					  vertex(G(i,  j,  k-1)) +
					  vertex(G(i,  j,  k+1)))/6.;

#undef G
	  }
	
	else if ((type == HEX20) ||
		 (type == HEX27))
	  {
	    
#define G(i,j,k) ( (i) + (2*nx+1)*((j) + (k)*(2*ny+1)) )

	    // Pop IJ boundary to the sphere
	    for (unsigned int i=0; i<=2*nx; i++)
	      for (unsigned int j=0; j<=2*ny; j++)
		{
		  unsigned int k = 0;
		
		  _vertices[G(i,j,k)] =
		    sphere.closest_point( vertex(G(i,j,k)) );
		  
		  k = 2*nz;
		
		_vertices[G(i,j,k)] =
		  sphere.closest_point( vertex(G(i,j,k)) );
	      };
	    
	    // Pop JK boundary to the sphere
	    for (unsigned int j=0; j<=2*ny; j++)
	      for (unsigned int k=0; k<=2*nz; k++)
		{
		  unsigned int i = 0;
		
		  _vertices[G(i,j,k)] =
		    sphere.closest_point( vertex(G(i,j,k)) );
		
		  i = 2*nx;
		
		  _vertices[G(i,j,k)] =
		    sphere.closest_point( vertex(G(i,j,k)) );
	      };
	    
	    // Pop IK boundary to the sphere
	    for (unsigned int i=0; i<=2*nx; i++)
	      for (unsigned int k=0; k<=2*nz; k++)
		{
		  unsigned int j = 0;
		
		  _vertices[G(i,j,k)] =
		    sphere.closest_point( vertex(G(i,j,k)) );
		
		  j = 2*ny;
		
		  _vertices[G(i,j,k)] =
		    sphere.closest_point( vertex(G(i,j,k)) );
	      };
	    
	    // Handle internal nodes
	    for (unsigned int k=1; k<2*nz; k++)
	      for (unsigned int j=1; j<2*ny; j++)
		for (unsigned int i=1; i<2*nx; i++)
		{
		  const real xmin = vertex(G(0,   j, k))(0);
		  const real xmax = vertex(G(2*nx,j, k))(0);
		  const real ymin = vertex(G(i,   0, k))(1);
		  const real ymax = vertex(G(i,2*ny, k))(1);
		  const real zmin = vertex(G(i, j,   0))(2);
		  const real zmax = vertex(G(i, j,2*nz))(2);

		  _vertices[G(i,j,k)](0) = xmin +
		    (xmax - xmin)*((real) ((real) i)/((real) 2*nx));
		  
		  _vertices[G(i,j,k)](1) = ymin +
		    (ymax - ymin)*((real) ((real) j)/((real) 2*ny));
		  
		  _vertices[G(i,j,k)](2) = zmin +
		    (zmax - zmin)*((real) ((real) k)/((real) 2*nz));
		};

	    // Do some smoothing steps.
	    for (unsigned int l=0; l<10; l++)
	      for (unsigned int k=1; k<2*nz; k++)
		for (unsigned int j=1; j<2*ny; j++)
		  for (unsigned int i=1; i<2*nx; i++)
		    _vertices[G(i,j,k)] = (vertex(G(i-1,j,  k  )) +
					  vertex(G(i+1,j,  k  )) +
					  vertex(G(i,  j-1,k  )) +
					  vertex(G(i,  j+1,k  )) +
					  vertex(G(i,  j,  k-1)) +
					  vertex(G(i,  j,  k+1)))/6.;

#undef G
	  }
	
	else
	  error();



	break;
      };


      
    default:
      error();
    };
};



void Mesh::read(const std::string name)
{
  _perf_log.start_event("read()");

  
  // Read the file based on extension
  {
    if (name.rfind(".mat") < name.size())
      read_matlab (name);
    
    else if (name.rfind(".ucd") < name.size())
      read_ucd (name);

    else if (name.rfind(".exd") < name.size())
      read_exd (name);

    else if (name.rfind(".xda") < name.size())
      read_xdr (name);

    else if ((name.rfind(".off")  < name.size()) ||
	     (name.rfind(".ogl")  < name.size()) ||
	     (name.rfind(".oogl") < name.size()))
      read_off(name);

    else if ((name.rfind(".xdr")  < name.size()) ||
	     (name.rfind(".0000") < name.size()))
      read_xdr_binary (name);

    else if (name.rfind(".mesh") < name.size())
      read_shanee (name);

    else if (name.rfind(".unv") < name.size())
      read_unv (name);

    else
      {
	std::cerr << " ERROR: Unrecognized file extension: " << name
		  << "\n   I understand the following:\n\n"
		  << "     *.mat  -- Matlab triangular ASCII file\n"
		  << "     *.ucd  -- AVS's ASCII UCD format\n"
		  << "     *.mesh -- Ben's \"shanee\" format\n"
		  << "     *.off  -- OOGL OFF surface format\n"
		  << "     *.exd  -- Sandia's ExodusII format\n"
		  << "     *.xda  -- Internal ASCII format\n"
		  << "     *.xdr  -- Internal binary format,\n"
		  << "               compatible with XdrMGF\n"
		  << "     *.unv  -- I-deas format\n"
		  << std::endl;
	error();

      }    
  };
  _perf_log.stop_event("read()");
};



void Mesh::write(const std::string name)
{
  _perf_log.start_event("write()");
  
  // Write the file based on extension
  {
    if (name.rfind(".dat") < name.size())
      write_tecplot (name);
    
    else if (name.rfind(".plt") < name.size())
      write_tecplot_binary (name);

    else if (name.rfind(".ucd") < name.size())
      write_ucd (name);

    else if (name.rfind(".gmv") < name.size())
      {
	if (n_subdomains() > 1)
	  write_gmv_binary(name, NULL, NULL, true);
	else
	  write_gmv_binary(name);
      }


    else if (name.rfind(".ugrid") < name.size())
      write_diva (name);
    
    else if (name.rfind(".xda") < name.size())
      write_xdr (name);
    
    else if (name.rfind(".xdr") < name.size())
      write_xdr_binary (name);

    else
      {
	std::cerr << " ERROR: Unrecognized file extension: " << name
		  << "\n   I understand the following:\n\n"
		  << "     *.dat   -- Tecplot ASCII file\n"
		  << "     *.plt   -- Tecplot binary file\n"
		  << "     *.ucd   -- AVS's ASCII UCD format\n"
		  << "     *.ugrid -- Kelly's DIVA ASCII format\n"
		  << "     *.gmv   -- LANL's GMV (General Mesh Viewer) format\n"
		  << "     *.xda   -- Internal ASCII format\n"
		  << "     *.xdr   -- Internal binary format,\n"
		  << "               compatible with XdrMGF\n"
		  << std::endl
		  << "\n Exiting without writing output\n";
      }    
  };
  
  _perf_log.stop_event("write()");
};



void Mesh::write(const std::string name,
		 std::vector<number>& v,
		 std::vector<std::string>& vn)
{
  _perf_log.start_event("write()");

  // Write the file based on extension
  {
    if (name.rfind(".dat") < name.size())
      write_tecplot (name, &v, &vn);
    
    else if (name.rfind(".plt") < name.size())
      write_tecplot_binary (name, &v, &vn);

    else if (name.rfind(".gmv") < name.size())
      {
	if (n_subdomains() > 1)
	  write_gmv_binary(name, &v, &vn, true);
	else
	  write_gmv_binary(name, &v, &vn);
      }    
    else
      {
	std::cerr << " ERROR: Unrecognized file extension: " << name
		  << "\n   I understand the following:\n\n"
		  << "     *.dat  -- Tecplot ASCII file\n"
		  << "     *.plt  -- Tecplot binary file\n"
		  << "     *.gmv  -- LANL's GMV (General Mesh Viewer) format\n"
		  << "\n Exiting without writing output\n";
      }
  };

  _perf_log.stop_event("write()");
};



#ifdef ENABLE_AMR

void Mesh::trim_unused_elements(std::set<unsigned int>& unused_elements)
{
  /**
   * Anything we clear in this routiune
   * will invalidate the unknowing boundary
   * mesh, so we need to clear it.  It must
   * be recreated before reuse.  
   */
  boundary_info.boundary_mesh.clear();
  
  
  /**
   * Trim the unused elements
   */
  {
    // We don't really need this in the
    // current implementation
    unused_elements.clear();

    // for the time being we make a copy
    // of the elements vector since the pointers
    // are relatively small.  Note that this is
    // not _necessary_, but it should be
    // less expensive than repeated calls
    // to std::vector<>::erase()    
    std::vector<Elem*> new_elements;
    
    new_elements.resize(n_elem());

    unsigned int ne=0;
    
    for (unsigned int e=0; e<n_elem(); e++)
      if (elem(e) != NULL)
	new_elements[ne++] = elem(e); 

    new_elements.resize(ne);
    
    _elements = new_elements;

#ifdef DEBUG

    for (unsigned int e=0; e<n_elem(); e++)
      assert (elem(e) != NULL);
    
#endif
    
  };
};

#endif
