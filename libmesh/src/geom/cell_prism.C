// $Id: cell_prism.C,v 1.5 2003-01-24 17:24:43 jwpeterson Exp $

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
#include "cell_prism.h"




// ------------------------------------------------------------
// Prism class member functions
AutoPtr<Elem> Prism::side (const unsigned int i) const
{
  assert (i < n_sides());


  
  AutoPtr<Elem> faceq(Elem::build(QUAD4));
  AutoPtr<Elem> facet(Elem::build(TRI3));

  switch (i)
    {
    case 0:  // the triangular face at z=0
      {
	facet->set_node(0) = get_node(0);
	facet->set_node(1) = get_node(2);
	facet->set_node(2) = get_node(1);

	return facet;
      }
    case 1:  // the quad face at y=0
      {
	faceq->set_node(0) = get_node(0);
	faceq->set_node(1) = get_node(1);
	faceq->set_node(2) = get_node(4);
	faceq->set_node(3) = get_node(3);
	
	return faceq;
      }
    case 2:  // the other quad face
      {
	faceq->set_node(0) = get_node(1);
	faceq->set_node(1) = get_node(2);
	faceq->set_node(2) = get_node(5);
	faceq->set_node(3) = get_node(4);

	return faceq;
      }
    case 3: // the quad face at x=0
      {
	faceq->set_node(0) = get_node(2);
	faceq->set_node(1) = get_node(0);
	faceq->set_node(2) = get_node(3);
	faceq->set_node(3) = get_node(5);
	
	return faceq;
      }
    case 4: // the triangular face at z=1
      {
	facet->set_node(0) = get_node(3);
	facet->set_node(1) = get_node(4);
	facet->set_node(2) = get_node(5);

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
