// $Id: cell_pyramid.C,v 1.11 2003-09-02 18:02:42 benkirk Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002-2003  Benjamin S. Kirk, John W. Peterson
  
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
#include "cell_pyramid.h"
#include "face_tri3.h"
#include "face_quad4.h"


// ------------------------------------------------------------
// Pyramid class member functions
unsigned int Pyramid::key (const unsigned int s) const
{  
  assert (s < this->n_sides());

  
  switch (s)
    {
    case 0:  // triangular face 1

      return
	this->compute_key (this->node(0),
			   this->node(1),
			   this->node(4));
      	
    case 1:  // triangular face 2

      return
	this->compute_key (this->node(1),
			   this->node(2),
			   this->node(4));
	
    case 2:  // triangular face 3

      return
	this->compute_key (this->node(2),
			   this->node(3),
			   this->node(4));
      	
    case 3:  // triangular face 4

      return
	this->compute_key (this->node(3),
			   this->node(0),
			   this->node(4));
      
    case 4:  // the quad face at z=0

      return
	this->compute_key (this->node(0),
			   this->node(3),
			   this->node(2),
			   this->node(1));
    }

  // We'll never get here.
  error();
  return 0;
}



AutoPtr<Elem> Pyramid::side (const unsigned int i) const
{
  assert (i < this->n_sides());


  
  switch (i)
    {
    case 0:  // triangular face 1
      {
	AutoPtr<Elem> face(new Tri3); 

	face->set_node(0) = this->get_node(0);
	face->set_node(1) = this->get_node(1);
	face->set_node(2) = this->get_node(4);
	
	return face;
      }
    case 1:  // triangular face 2
      {
	AutoPtr<Elem> face(new Tri3);

	face->set_node(0) = this->get_node(1);
	face->set_node(1) = this->get_node(2);
	face->set_node(2) = this->get_node(4);
	
	return face;
      }
    case 2:  // triangular face 3
      {
	AutoPtr<Elem> face(new Tri3);

	face->set_node(0) = this->get_node(2);
	face->set_node(1) = this->get_node(3);
	face->set_node(2) = this->get_node(4);
	
	return face;
      }
    case 3:  // triangular face 4
      {
	AutoPtr<Elem> face(new Tri3);

	face->set_node(0) = this->get_node(3);
	face->set_node(1) = this->get_node(0);
	face->set_node(2) = this->get_node(4);
	
	return face;
      }
    case 4:  // the quad face at z=0
      {
	AutoPtr<Elem> face(new Quad4);
	
	face->set_node(0) = this->get_node(0);
	face->set_node(1) = this->get_node(3);
	face->set_node(2) = this->get_node(2);
	face->set_node(3) = this->get_node(1);

	return face;
      }
    default:
      {
	error();
	AutoPtr<Elem> face(NULL);
	return face;
      }
    }

  // We'll never get here.
  error();
  AutoPtr<Elem> face(NULL);

  return face;
}
