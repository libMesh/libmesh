// $Id: cell_pyramid.C,v 1.1.1.1 2003-01-10 16:17:48 libmesh Exp $

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

// Local includes
#include "mesh.h"

// Temporary includes
#include "cell_pyramid.h"
#include "face_tri3.h"
#include "face_quad4.h"




// ------------------------------------------------------------
// Pyramid class member functions
Elem Pyramid::side (const unsigned int i) const
{
  assert (i < n_sides());
  assert (_nodes.size() == n_nodes());


  
  switch (i)
    {
    case 0:  // triangular face 1
      {
	Elem face(3);

	face.node(0) = node(0);
	face.node(1) = node(1);
	face.node(2) = node(4);
	
	return face;
      }
    case 1:  // triangular face 2
      {
	Elem face(3);

	face.node(0) = node(1);
	face.node(1) = node(2);
	face.node(2) = node(4);
	
	return face;
      }
    case 2:  // triangular face 3
      {
	Elem face(3);

	face.node(0) = node(2);
	face.node(1) = node(3);
	face.node(2) = node(4);
	
	return face;
      }
    case 3:  // triangular face 4
      {
	Elem face(3);

	face.node(0) = node(3);
	face.node(1) = node(0);
	face.node(2) = node(4);
	
	return face;
      }
    case 4:  // the quad face at z=0
      {
	Elem face(4);
	
	face.node(0) = node(0);
	face.node(1) = node(3);
	face.node(2) = node(2);
	face.node(3) = node(1);

	return face;
      }
    default:
      {
	error();
	Elem face(0);
	return face;
      }
    };

  // We'll never get here.
  error();
  Elem face(0);

  return face;
};
