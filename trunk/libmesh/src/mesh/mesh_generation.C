// $Id: mesh_generation.C,v 1.3 2003-01-21 19:24:37 benkirk Exp $

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
#include <math.h>


// Local includes
#include "mesh.h"
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



// ------------------------------------------------------------
// Mesh class member functions for mesh generation
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
	
	_nodes.resize(nx+1);
	_elements.resize(nx);

	for (unsigned int i=0; i<=nx; i++)
	  {
	    node_ptr(i) = Node::build((real) ((real) i)/((real) nx), 0, 0, i);
	  };
	
	for (unsigned int i=0; i<nx; i++)
	  {
	    _elements[i] = new Edge2;
	    elem(i)->set_node(0) = _nodes[i];
	    elem(i)->set_node(1) = _nodes[i+1];
	  };


	// Scale the nodal positions
	for (unsigned int p=0; p<n_nodes(); p++)
	  {
	    node(p)(0) = (node(p)(0))*(xmax-xmin) + xmin;
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
	
	    _nodes.resize( (nx+1)*(ny+1) );
	    _elements.resize(nx*ny);
	    
	    unsigned int p=0;
	    
	    for (unsigned int j=0; j<=ny; j++)
	      for (unsigned int i=0; i<=nx; i++)
		node_ptr(p) = Node::build((real) ((real) i)/((real) nx),
					  (real) ((real) j)/((real) ny),
					  0,
					  p++);
	    
	    unsigned int e=0;
	    
#define G(i,j) ( (i) + (j)*(nx+1) )
	
	    for (unsigned int j=0; j<ny; j++)
	      for (unsigned int i=0; i<nx; i++)
		{
		  _elements[e] = new Quad4;

		  elem(e)->set_node(0) = _nodes[G(i,j)    ];
		  elem(e)->set_node(1) = _nodes[G(i+1,j)  ];
		  elem(e)->set_node(2) = _nodes[G(i+1,j+1)];
		  elem(e)->set_node(3) = _nodes[G(i,j+1)  ];

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
	
	    _nodes.resize( (nx+1)*(ny+1) );
	    _elements.resize(2*nx*ny);
	    
	    unsigned int p=0;
	    
	    for (unsigned int j=0; j<=ny; j++)
	      for (unsigned int i=0; i<=nx; i++)
		node_ptr(p) = Node::build((real) ((real) i)/((real) nx),
					  (real) ((real) j)/((real) ny),
					  0,
					  p++);
	    
	    unsigned int e=0;
	    
#define G(i,j) ( (i) + (j)*(nx+1) )
	
	    for (unsigned int j=0; j<ny; j++)
	      for (unsigned int i=0; i<nx; i++)
		{
		  _elements[e] = new Tri3;

		  elem(e)->set_node(0) = _nodes[G(i,j)    ];
		  elem(e)->set_node(1) = _nodes[G(i+1,j)  ];
		  elem(e)->set_node(2) = _nodes[G(i+1,j+1)];
		  
		  e++;
		  
		  _elements[e] = new Tri3;

		  elem(e)->set_node(0) = _nodes[G(i,j)    ];
		  elem(e)->set_node(1) = _nodes[G(i+1,j+1)];
		  elem(e)->set_node(2) = _nodes[G(i,j+1)  ];
		  
		  e++;
		};
	    
#undef G
	  }


	
	else if ((type == QUAD8) ||
		 (type == QUAD9))
	  {
	    // Build a Quad8 or a QUAD9
	
	    _nodes.resize( (2*nx+1)*(2*ny+1) );
	    _elements.resize(nx*ny);
	    
	    unsigned int p=0;
	    
	    for (unsigned int j=0; j<=(2*ny); j++)
	      for (unsigned int i=0; i<=(2*nx); i++)
		node_ptr(p) = Node::build((real) ((real) i)/((real) (2*nx)),
					  (real) ((real) j)/((real) (2*ny)),
					  0,
					  p++);
	    
	    unsigned int e=0;
	    
#define G(i,j) ( (i) + (j)*(2*nx+1) )
	
	    for (unsigned int j=0; j<(2*ny); j += 2)
	      for (unsigned int i=0; i<(2*nx); i += 2)
		{		
		  if (type == QUAD8)
		    _elements[e] = new Quad8;
		  else
		    _elements[e] = new Quad9;


		  elem(e)->set_node(0) = _nodes[G(i,j)    ];
		  elem(e)->set_node(1) = _nodes[G(i+2,j)  ];
		  elem(e)->set_node(2) = _nodes[G(i+2,j+2)];
		  elem(e)->set_node(3) = _nodes[G(i,j+2)  ];
		  elem(e)->set_node(4) = _nodes[G(i+1,j)  ];
		  elem(e)->set_node(5) = _nodes[G(i+2,j+1)];
		  elem(e)->set_node(6) = _nodes[G(i+1,j+2)];
		  elem(e)->set_node(7) = _nodes[G(i,j+1)  ];
		  if (type == QUAD9)
		    elem(e)->set_node(8) = _nodes[G(i+1,j+1)];
		  

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
	
	    _nodes.resize( (2*nx+1)*(2*ny+1) );
	    _elements.resize(2*nx*ny);
	    
	    unsigned int p=0;
	    
	    for (unsigned int j=0; j<=(2*ny); j++)
	      for (unsigned int i=0; i<=(2*nx); i++)
		node_ptr(p) = Node::build((real) ((real) i)/((real) (2*nx)),
					  (real) ((real) j)/((real) (2*ny)),
					  0,
					  p++);
	    
	    unsigned int e=0;
	    
#define G(i,j) ( (i) + (j)*(2*nx+1) )
	
	    for (unsigned int j=0; j<(2*ny); j += 2)
	      for (unsigned int i=0; i<(2*nx); i += 2)
		{
		  _elements[e] = new Tri6;

		  elem(e)->set_node(0) = _nodes[G(i,j)    ];
		  elem(e)->set_node(1) = _nodes[G(i+2,j)  ];
		  elem(e)->set_node(2) = _nodes[G(i+2,j+2)];
		  elem(e)->set_node(3) = _nodes[G(i+1,j)  ];
		  elem(e)->set_node(4) = _nodes[G(i+2,j+1)];
		  elem(e)->set_node(5) = _nodes[G(i+1,j+1)];
		  
		  e++;
		  
		  _elements[e] = new Tri6;

		  elem(e)->set_node(0) = _nodes[G(i,j)    ];
		  elem(e)->set_node(1) = _nodes[G(i+2,j+2)];
		  elem(e)->set_node(2) = _nodes[G(i,j+2)  ];
		  elem(e)->set_node(3) = _nodes[G(i+1,j+1)];
		  elem(e)->set_node(4) = _nodes[G(i+1,j+2)];
		  elem(e)->set_node(5) = _nodes[G(i,j+1)  ];
		  
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
	    node(p)(0) = (node(p)(0))*(xmax-xmin) + xmin;
	    node(p)(1) = (node(p)(1))*(ymax-ymin) + ymin;
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
	    _nodes.resize( (nx+1)*(ny+1)*(nz+1) );
	    _elements.resize(nx*ny*nz);
	    
	    unsigned int p=0;
	    
	    for (unsigned int k=0; k<=nz; k++)
	      for (unsigned int j=0; j<=ny; j++)
		for (unsigned int i=0; i<=nx; i++)
		  node_ptr(p) = Node::build((real) ((real) i)/((real) nx),
					    (real) ((real) j)/((real) ny),
					    (real) ((real) k)/((real) nz),
					    p++);
	    
	    unsigned int e=0;
	    
#define G(i,j,k) ( (i) + (nx+1)*((j) + (k)*(ny+1)) )
	    
	    for (unsigned int k=0; k<nz; k++)
	      for (unsigned int j=0; j<ny; j++)
		for (unsigned int i=0; i<nx; i++)
		  {
		    _elements[e] = new Hex8;
		    
		    elem(e)->set_node(0) = _nodes[G(i,j,k)      ];
		    elem(e)->set_node(1) = _nodes[G(i+1,j,k)    ];
		    elem(e)->set_node(2) = _nodes[G(i+1,j+1,k)  ];
		    elem(e)->set_node(3) = _nodes[G(i,j+1,k)    ];
		    elem(e)->set_node(4) = _nodes[G(i,j,k+1)    ];
		    elem(e)->set_node(5) = _nodes[G(i+1,j,k+1)  ];
		    elem(e)->set_node(6) = _nodes[G(i+1,j+1,k+1)];
		    elem(e)->set_node(7) = _nodes[G(i,j+1,k+1)  ];
		    
		    e++;
		  };
#undef G

	  }
	
	else if ((type == HEX20) ||
		 (type == HEX27))
	  {
	    // Build a HEX20 or HEX27
	    _nodes.resize( (2*nx+1)*(2*ny+1)*(2*nz+1) );
	    _elements.resize(nx*ny*nz);
	    
	    unsigned int p=0;
	    
	    for (unsigned int k=0; k<=(2*nz); k++)
	      for (unsigned int j=0; j<=(2*ny); j++)
		for (unsigned int i=0; i<=(2*nx); i++)
		  node_ptr(p) = Node::build((real) ((real) i)/((real) (2*nx)),
					    (real) ((real) j)/((real) (2*ny)),
					    (real) ((real) k)/((real) (2*nz)),
					    p++);
	    
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
		    
		    elem(e)->set_node(0)  = _nodes[G(i,  j,  k)  ];
		    elem(e)->set_node(1)  = _nodes[G(i+2,j,  k)  ];
		    elem(e)->set_node(2)  = _nodes[G(i+2,j+2,k)  ];
		    elem(e)->set_node(3)  = _nodes[G(i,  j+2,k)  ];
		    elem(e)->set_node(4)  = _nodes[G(i,  j,  k+2)];
		    elem(e)->set_node(5)  = _nodes[G(i+2,j,  k+2)];
		    elem(e)->set_node(6)  = _nodes[G(i+2,j+2,k+2)];
		    elem(e)->set_node(7)  = _nodes[G(i,  j+2,k+2)];
		    elem(e)->set_node(8)  = _nodes[G(i+1,j,  k)  ];
		    elem(e)->set_node(9)  = _nodes[G(i+2,j+1,k)  ];
		    elem(e)->set_node(10) = _nodes[G(i+1,j+2,k)  ];
		    elem(e)->set_node(11) = _nodes[G(i,  j+1,k)  ];
		    elem(e)->set_node(12) = _nodes[G(i,  j,  k+1)];
		    elem(e)->set_node(13) = _nodes[G(i+2,j,  k+1)];
		    elem(e)->set_node(14) = _nodes[G(i+2,j+2,k+1)];
		    elem(e)->set_node(15) = _nodes[G(i,  j+2,k+1)];
		    elem(e)->set_node(16) = _nodes[G(i+1,j,  k+2)];
		    elem(e)->set_node(17) = _nodes[G(i+2,j+1,k+2)];
		    elem(e)->set_node(18) = _nodes[G(i+1,j+2,k+2)];
		    elem(e)->set_node(19) = _nodes[G(i,  j+1,k+2)];
		    if (type == HEX27)
		      {
			elem(e)->set_node(20) = _nodes[G(i+1,j+1,k)  ];
			elem(e)->set_node(21) = _nodes[G(i+1,j,  k+1)];
			elem(e)->set_node(22) = _nodes[G(i+2,j+1,k+1)];
			elem(e)->set_node(23) = _nodes[G(i+1,j+2,k+1)];
			elem(e)->set_node(24) = _nodes[G(i,  j+1,k+1)];
			elem(e)->set_node(25) = _nodes[G(i+1,j+1,k+2)];
			elem(e)->set_node(26) = _nodes[G(i+1,j+1,k+1)];
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
	    node(p)(0) = (node(p)(0))*(xmax-xmin) + xmin;
	    node(p)(1) = (node(p)(1))*(ymax-ymin) + ymin;
	    node(p)(2) = (node(p)(2))*(zmax-zmin) + zmin;
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

	    _nodes.resize(8);
	    
	    // Point 0
	    node_ptr(0) = Node::build(-rad_2,-rad_2, 0, 0);
	    
	    // Point 1
	    node_ptr(1) = Node::build( rad_2,-rad_2, 0, 1);
	    
	    // Point 2
	    node_ptr(2) = Node::build( rad_2, rad_2, 0, 2);

	    // Point 3
	    node_ptr(3) = Node::build(-rad_2, rad_2, 0, 3);
	    
	    // Point 4
	    node_ptr(4) = Node::build(-rad_sqrt_2,-rad_sqrt_2, 0, 4);
	    
	    // Point 5
	    node_ptr(5) = Node::build( rad_sqrt_2,-rad_sqrt_2, 0, 5);
	    
	    // Point 6
	    node_ptr(6) = Node::build( rad_sqrt_2, rad_sqrt_2, 0, 6);
	    
	    // Point 7
	    node_ptr(7) = Node::build(-rad_sqrt_2, rad_sqrt_2, 0, 7);

	    // Build the elements
	    _elements.resize(5);
	    
	    for (unsigned int e=0; e<n_elem(); e++)
	      _elements[e] = Elem::build(QUAD4);
	    
	    // Element 0
	    elem(0)->set_node(0) = _nodes[0];
	    elem(0)->set_node(1) = _nodes[1];
	    elem(0)->set_node(2) = _nodes[2];
	    elem(0)->set_node(3) = _nodes[3];
	    
	    // Element 1
	    elem(1)->set_node(0) = _nodes[4];
	    elem(1)->set_node(1) = _nodes[0];
	    elem(1)->set_node(2) = _nodes[3];
	    elem(1)->set_node(3) = _nodes[7];
	    
	    // Element 2
	    elem(2)->set_node(0) = _nodes[4];
	    elem(2)->set_node(1) = _nodes[5];
	    elem(2)->set_node(2) = _nodes[1];
	    elem(2)->set_node(3) = _nodes[0];
	    
	    // Element 3
	    elem(3)->set_node(0) = _nodes[1];
	    elem(3)->set_node(1) = _nodes[5];
	    elem(3)->set_node(2) = _nodes[6];
	    elem(3)->set_node(3) = _nodes[2];
	    
	    // Element 4
	    elem(4)->set_node(0) = _nodes[3];
	    elem(4)->set_node(1) = _nodes[2];
	    elem(4)->set_node(2) = _nodes[6];
	    elem(4)->set_node(3) = _nodes[7];
	  }


	// Quadratic elements
	else if ((type == QUAD8) ||
		 (type == QUAD9) ||
		 (type == TRI6))
	  {
	    
	    _nodes.resize(25);

	    // Point 0
	    node_ptr(0)  = Node::build(-rad_2,-rad_2, 0, 0);
	    
	    // Point 1
	    node_ptr(1)  = Node::build( rad_2,-rad_2, 0, 1);
	    
	    // Point 2
	    node_ptr(2)  = Node::build( rad_2, rad_2, 0, 2);

	    // Point 3
	    node_ptr(3)  = Node::build(-rad_2, rad_2, 0, 3);
	    
	    // Point 4
	    node_ptr(4)  = Node::build(-rad_sqrt_2,-rad_sqrt_2, 0, 4);
	     
	    // Point 5
	    node_ptr(5)  = Node::build( rad_sqrt_2,-rad_sqrt_2, 0, 5);
	    
	    // Point 6
	    node_ptr(6)  = Node::build( rad_sqrt_2, rad_sqrt_2, 0, 6);
	    
	    // Point 7
	    node_ptr(7)  = Node::build(-rad_sqrt_2, rad_sqrt_2, 0, 7);

	    // Point 8
	    node_ptr(8)  = Node::build((point(0) + point(1))/2., 8);

	    // Point 9
	    node_ptr(9)  = Node::build((point(1) + point(2))/2., 9);
	    
	    // Point 10
	    node_ptr(10) = Node::build((point(2) + point(3))/2., 10);
	    
	    // Point 11
	    node_ptr(11) = Node::build((point(0) + point(3))/2., 11);
	    
	    // Point 12
	    node_ptr(12) = Node::build((point(0) + point(2))/2., 12);

	    // Point 13
	    node_ptr(13) = Node::build( 0,-rad, 0, 13);

	    // Point 14
	    node_ptr(14) = Node::build( rad, 0, 0, 14);

	    // Point 15
	    node_ptr(15) = Node::build( 0, rad, 0, 15);

	    // Point 16
	    node_ptr(16) = Node::build(-rad, 0, 0, 16);

	    // Point 17
	    node_ptr(17) = Node::build((point(8) + point(13))/2., 17);

	    // Point 18
	    node_ptr(18) = Node::build((point(9) + point(14))/2., 18);

	    // Point 19
	    node_ptr(19) = Node::build((point(10) + point(15))/2., 19);
	    
	    // Point 20
	    node_ptr(20) = Node::build((point(11) + point(16))/2., 20);
	      
	    // Point 21
	    node_ptr(21) = Node::build((point(0) + point(4))/2., 21);

	    // Point 22
	    node_ptr(22) = Node::build((point(1) + point(5))/2., 22);

	    // Point 23
	    node_ptr(23) = Node::build((point(2) + point(6))/2., 23);

	    // Point 24
	    node_ptr(24) = Node::build((point(3) + point(7))/2., 24);

	    
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
	    elem(0)->set_node(0) = _nodes[0];
	    elem(0)->set_node(1) = _nodes[1];
	    elem(0)->set_node(2) = _nodes[2];
	    elem(0)->set_node(3) = _nodes[3];
	    elem(0)->set_node(4) = _nodes[8];
	    elem(0)->set_node(5) = _nodes[9];
	    elem(0)->set_node(6) = _nodes[10];
	    elem(0)->set_node(7) = _nodes[11];
	    if (type != QUAD8) elem(0)->set_node(8) = _nodes[12];

	    // Element 1
	    elem(1)->set_node(0) = _nodes[4];
	    elem(1)->set_node(1) = _nodes[0];
	    elem(1)->set_node(2) = _nodes[3];
	    elem(1)->set_node(3) = _nodes[7];
	    elem(1)->set_node(4) = _nodes[21];
	    elem(1)->set_node(5) = _nodes[11];
	    elem(1)->set_node(6) = _nodes[24];
	    elem(1)->set_node(7) = _nodes[16];
	    if (type != QUAD8) elem(1)->set_node(8) = _nodes[20];

	    // Element 2
	    elem(2)->set_node(0) = _nodes[4];
	    elem(2)->set_node(1) = _nodes[5];
	    elem(2)->set_node(2) = _nodes[1];
	    elem(2)->set_node(3) = _nodes[0];
	    elem(2)->set_node(4) = _nodes[13];
	    elem(2)->set_node(5) = _nodes[22];
	    elem(2)->set_node(6) = _nodes[8];
	    elem(2)->set_node(7) = _nodes[21];
	    if (type != QUAD8) elem(2)->set_node(8) = _nodes[17];

	    // Element 3
	    elem(3)->set_node(0) = _nodes[1];
	    elem(3)->set_node(1) = _nodes[5];
	    elem(3)->set_node(2) = _nodes[6];
	    elem(3)->set_node(3) = _nodes[2];
	    elem(3)->set_node(4) = _nodes[22];
	    elem(3)->set_node(5) = _nodes[14];
	    elem(3)->set_node(6) = _nodes[23];
	    elem(3)->set_node(7) = _nodes[9];
	    if (type != QUAD8) elem(3)->set_node(8) = _nodes[18];

	    // Element 4
	    elem(4)->set_node(0) = _nodes[3];
	    elem(4)->set_node(1) = _nodes[2];
	    elem(4)->set_node(2) = _nodes[6];
	    elem(4)->set_node(3) = _nodes[7];
	    elem(4)->set_node(4) = _nodes[10];
	    elem(4)->set_node(5) = _nodes[23];
	    elem(4)->set_node(6) = _nodes[15];
	    elem(4)->set_node(7) = _nodes[24];
	    if (type != QUAD8) elem(4)->set_node(8) = _nodes[19];

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
		      AutoPtr<Elem> side(elem(e)->build_side(s));
		      
		      for (unsigned int n=0; n<side->n_nodes(); n++)
			side->point(n) =
			  sphere.closest_point(side->point(n));
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
		  // Build the new elements.  We can't just copy
		  // the pointers since we must delete the parents to
		  // avoid a memory leak.  Deleting the parents will
		  // also delete the children, so we need to build a
		  // whole new element for this to work properly.
		  // Clear?
		  new_elements[ne] = Elem::build(elem(e)->type());
		  
		  for (unsigned int n=0; n<elem(e)->n_nodes(); n++)
		    new_elements[ne]->set_node(n) =
		      elem(e)->get_node(n);
		  
		  ne++;
		};

	      // Delete the element
	      delete elem(e);
	      _elements[e] = NULL;
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
		
		  node(G(i,j,k)) =
		    sphere.closest_point( point(G(i,j,k)) );
		  
		  k = nz;
		
		  node(G(i,j,k)) =
		    sphere.closest_point( point(G(i,j,k)) );
		};
	    
	    // Pop JK boundary to the sphere
	    for (unsigned int j=0; j<=ny; j++)
	      for (unsigned int k=0; k<=nz; k++)
		{
		  unsigned int i = 0;
		
		  node(G(i,j,k)) =
		    sphere.closest_point( point(G(i,j,k)) );
		
		  i = nx;
		
		  node(G(i,j,k)) =
		    sphere.closest_point( point(G(i,j,k)) );
		};
	    
	    // Pop IK boundary to the sphere
	    for (unsigned int i=0; i<=nx; i++)
	      for (unsigned int k=0; k<=nz; k++)
		{
		  unsigned int j = 0;
		
		  node(G(i,j,k)) =
		    sphere.closest_point( point(G(i,j,k)) );
		
		  j = ny;
		
		  node(G(i,j,k)) =
		    sphere.closest_point( point(G(i,j,k)) );
		};
	    
	    // Handle internal nodes
	    for (unsigned int k=1; k<nz; k++)
	      for (unsigned int j=1; j<ny; j++)
		for (unsigned int i=1; i<nx; i++)
		  {
		    const real xmin = point(G(0 ,j, k))(0);
		    const real xmax = point(G(nx,j, k))(0);
		    const real ymin = point(G(i, 0, k))(1);
		    const real ymax = point(G(i,ny, k))(1);
		    const real zmin = point(G(i, j, 0))(2);
		    const real zmax = point(G(i, j,nz))(2);

		    node(G(i,j,k))(0) = xmin +
		      (xmax - xmin)*((real) ((real) i)/((real) nx));
		  
		    node(G(i,j,k))(1) = ymin +
		      (ymax - ymin)*((real) ((real) j)/((real) ny));
		  
		    node(G(i,j,k))(2) = zmin +
		      (zmax - zmin)*((real) ((real) k)/((real) nz));
		  };

	    // Do some smoothing steps.
	    for (unsigned int l=0; l<10; l++)
	      for (unsigned int k=1; k<nz; k++)
		for (unsigned int j=1; j<ny; j++)
		  for (unsigned int i=1; i<nx; i++)
		    node(G(i,j,k)) = (point(G(i-1,j,  k  )) +
				      point(G(i+1,j,  k  )) +
				      point(G(i,  j-1,k  )) +
				      point(G(i,  j+1,k  )) +
				      point(G(i,  j,  k-1)) +
				      point(G(i,  j,  k+1)))/6.;

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
		
		  node(G(i,j,k)) =
		    sphere.closest_point( point(G(i,j,k)) );
		  
		  k = 2*nz;
		
		  node(G(i,j,k)) =
		    sphere.closest_point( point(G(i,j,k)) );
		};
	    
	    // Pop JK boundary to the sphere
	    for (unsigned int j=0; j<=2*ny; j++)
	      for (unsigned int k=0; k<=2*nz; k++)
		{
		  unsigned int i = 0;
		
		  node(G(i,j,k)) =
		    sphere.closest_point( point(G(i,j,k)) );
		
		  i = 2*nx;
		
		  node(G(i,j,k)) =
		    sphere.closest_point( point(G(i,j,k)) );
		};
	    
	    // Pop IK boundary to the sphere
	    for (unsigned int i=0; i<=2*nx; i++)
	      for (unsigned int k=0; k<=2*nz; k++)
		{
		  unsigned int j = 0;
		
		  node(G(i,j,k)) =
		    sphere.closest_point( point(G(i,j,k)) );
		
		  j = 2*ny;
		
		  node(G(i,j,k)) =
		    sphere.closest_point( point(G(i,j,k)) );
		};
	    
	    // Handle internal nodes
	    for (unsigned int k=1; k<2*nz; k++)
	      for (unsigned int j=1; j<2*ny; j++)
		for (unsigned int i=1; i<2*nx; i++)
		  {
		    const real xmin = point(G(0,   j, k))(0);
		    const real xmax = point(G(2*nx,j, k))(0);
		    const real ymin = point(G(i,   0, k))(1);
		    const real ymax = point(G(i,2*ny, k))(1);
		    const real zmin = point(G(i, j,   0))(2);
		    const real zmax = point(G(i, j,2*nz))(2);

		    node(G(i,j,k))(0) = xmin +
		      (xmax - xmin)*((real) ((real) i)/((real) 2*nx));
		  
		    node(G(i,j,k))(1) = ymin +
		      (ymax - ymin)*((real) ((real) j)/((real) 2*ny));
		  
		    node(G(i,j,k))(2) = zmin +
		      (zmax - zmin)*((real) ((real) k)/((real) 2*nz));
		  };

	    // Do some smoothing steps.
	    for (unsigned int l=0; l<10; l++)
	      for (unsigned int k=1; k<2*nz; k++)
		for (unsigned int j=1; j<2*ny; j++)
		  for (unsigned int i=1; i<2*nx; i++)
		    node(G(i,j,k)) = (point(G(i-1,j,  k  )) +
				      point(G(i+1,j,  k  )) +
				      point(G(i,  j-1,k  )) +
				      point(G(i,  j+1,k  )) +
				      point(G(i,  j,  k-1)) +
				      point(G(i,  j,  k+1)))/6.;

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
