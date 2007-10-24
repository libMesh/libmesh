// $Id: cell_tet.C,v 1.14 2004-10-25 21:49:25 benkirk Exp $

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

// Local includes
#include "cell_tet.h"
#include "face_tri3.h"



// ------------------------------------------------------------
// Tet class member functions
unsigned int Tet::key (const unsigned int s) const
{
  assert (s < this->n_sides());

  switch (s)
    {
    case 0:
      return
	this->compute_key (this->node(0),
			   this->node(2),
			   this->node(1));
      
    case 1:
      return
	this->compute_key (this->node(0),
			   this->node(1),
			   this->node(3));

    case 2:
      return
	this->compute_key (this->node(1),
			   this->node(2),
			   this->node(3));

    case 3:
      return
	this->compute_key (this->node(2),
			   this->node(0),
			   this->node(3));	
    }

  // We'll never get here.
  error();
  return 0;
}



AutoPtr<Elem> Tet::side (const unsigned int i) const
{
  assert (i < this->n_sides());


  
  AutoPtr<Elem> face(new Tri3);
  
  switch (i)
    {
    case 0:
      {
	face->set_node(0) = this->get_node(0);
	face->set_node(1) = this->get_node(2);
	face->set_node(2) = this->get_node(1);

	return face;
      }
    case 1:
      {
	face->set_node(0) = this->get_node(0);
	face->set_node(1) = this->get_node(1);
	face->set_node(2) = this->get_node(3);

	return face;
      }
    case 2:
      {
	face->set_node(0) = this->get_node(1);
	face->set_node(1) = this->get_node(2);
	face->set_node(2) = this->get_node(3);

	return face;
      }
    case 3:
      {
	face->set_node(0) = this->get_node(2);
	face->set_node(1) = this->get_node(0);
	face->set_node(2) = this->get_node(3);
	
	return face;
      }
    default:
      {
	error();
      }
    }

  // We'll never get here.
  error();
  return face;
}



Real Tet::quality(const ElemQuality) const
{
  return -1.0; // Not implemented
}




std::pair<Real, Real> Tet::qual_bounds (const ElemQuality q) const
{
  std::pair<Real, Real> bounds;
  
  switch (q)
    {

    case ASPECT_RATIO_BETA:
    case ASPECT_RATIO_GAMMA:
      bounds.first  = 1.;
      bounds.second = 3.;
      break;
      
    case SIZE:
    case SHAPE:
      bounds.first  = 0.2;
      bounds.second = 1.;
      break;

    case CONDITION:
      bounds.first  = 1.;
      bounds.second = 3.;
      break;
      
    case DISTORTION:
      bounds.first  = 0.6;
      bounds.second = 1.;
      break;  

    case JACOBIAN:
      bounds.first  = 0.5;
      bounds.second = 1.414;
      break;
      
    default:
      std::cout << "Warning: Invalid quality measure chosen." << std::endl;
      bounds.first  = -1;
      bounds.second = -1;
    }

  return bounds;
}
