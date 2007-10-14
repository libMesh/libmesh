// $Id: cell_tet.C,v 1.17 2006-11-02 13:50:36 jwpeterson Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2005  Benjamin S. Kirk, John W. Peterson
  
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



AutoPtr<DofObject> Tet::side (const unsigned int i) const
{
  assert (i < this->n_sides());


  
  Elem* face = new Tri3;
  
  switch (i)
    {
    case 0:
      {
	face->set_node(0) = this->get_node(0);
	face->set_node(1) = this->get_node(2);
	face->set_node(2) = this->get_node(1);

        AutoPtr<DofObject> ap_face(face);
	return ap_face;
      }
    case 1:
      {
	face->set_node(0) = this->get_node(0);
	face->set_node(1) = this->get_node(1);
	face->set_node(2) = this->get_node(3);

        AutoPtr<DofObject> ap_face(face);
	return ap_face;
      }
    case 2:
      {
	face->set_node(0) = this->get_node(1);
	face->set_node(1) = this->get_node(2);
	face->set_node(2) = this->get_node(3);

        AutoPtr<DofObject> ap_face(face);
	return ap_face;
      }
    case 3:
      {
	face->set_node(0) = this->get_node(2);
	face->set_node(1) = this->get_node(0);
	face->set_node(2) = this->get_node(3);
	
        AutoPtr<DofObject> ap_face(face);
	return ap_face;
      }
    default:
      {
	error();
      }
    }

  // We'll never get here.
  error();
  AutoPtr<DofObject> ap_face(face);
  return ap_face;
}



Real Tet::quality(const ElemQuality q) const
{
  return Elem::quality(q); // Not implemented
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
