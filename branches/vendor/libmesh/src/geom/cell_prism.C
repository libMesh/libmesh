// $Id: cell_prism.C,v 1.1.1.1 2003-01-10 16:17:48 libmesh Exp $

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

#include "cell_prism.h"
#include "face_tri3.h"
#include "face_quad4.h"




// ------------------------------------------------------------
// Prism class member functions
Elem Prism::side (const unsigned int i) const
{
  assert (i < n_sides());
  assert (_nodes.size() == n_nodes());


  
  Elem faceq(4);
  Elem facet(3);

  switch (i)
    {
    case 0:  // the triangular face at z=0
      {
	facet.node(0) = node(0);
	facet.node(1) = node(2);
	facet.node(2) = node(1);

	return facet;
      }
    case 1:  // the quad face at y=0
      {
	faceq.node(0) = node(0);
	faceq.node(1) = node(1);
	faceq.node(2) = node(4);
	faceq.node(3) = node(3);
	
	return faceq;
      }
    case 2:  // the other quad face
      {
	faceq.node(0) = node(1);
	faceq.node(1) = node(2);
	faceq.node(2) = node(5);
	faceq.node(3) = node(4);

	return faceq;
      }
    case 3: // the quad face at x=0
      {
	faceq.node(0) = node(2);
	faceq.node(1) = node(0);
	faceq.node(2) = node(3);
	faceq.node(3) = node(5);
	
	return faceq;
      }
    case 4: // the triangular face at z=1
      {
	facet.node(0) = node(3);
	facet.node(1) = node(4);
	facet.node(2) = node(5);

	return facet;
      }
    default:
      {
	error();
	return facet;
      }
    };

  // We'll never get here.
  error();

  return facet;
};
