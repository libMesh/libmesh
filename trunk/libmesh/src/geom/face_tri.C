// $Id: face_tri.C,v 1.16 2005-02-22 22:17:39 jwpeterson Exp $

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
#include "face_tri.h"
#include "edge_edge2.h"





// ------------------------------------------------------------
// Tri class member functions
unsigned int Tri::key (const unsigned int s) const
{
  assert (s < this->n_sides());

  switch (s)
    {
    case 0:
      return
	this->compute_key (this->node(0),
			   this->node(1));
	
    case 1:
      return
	this->compute_key (this->node(1),
			   this->node(2));
    case 2:
      return
	this->compute_key (this->node(2),
			   this->node(0));
    }

  
  // We will never get here...  Look at the code above.
  error();
  return 0;
}



AutoPtr<Elem> Tri::side (const unsigned int i) const
{
  assert (i < this->n_sides());

  AutoPtr<Elem> edge(new Edge2);

  switch (i)
    {
    case 0:
      {
	edge->set_node(0) = this->get_node(0);
	edge->set_node(1) = this->get_node(1);
	
	return edge;
      }
    case 1:
      {
	edge->set_node(0) = this->get_node(1);
	edge->set_node(1) = this->get_node(2);
	
	return edge;
      }
    case 2:
      {
	edge->set_node(0) = this->get_node(2);
	edge->set_node(1) = this->get_node(0);
	
	return edge;
      }
    default:
      {
	error();
      }
    }

  
  // We will never get here...  Look at the code above.
  error();
  return edge;
}






Real Tri::quality (const ElemQuality) const
{
  return 0.; // not implemented
}






std::pair<Real, Real> Tri::qual_bounds (const ElemQuality q) const
{
  std::pair<Real, Real> bounds;
  
  switch (q)
    {

    case MAX_ANGLE:
      bounds.first  = 60.;
      bounds.second = 90.;
      break;
      
    case MIN_ANGLE:
      bounds.first  = 30.;
      bounds.second = 60.;
      break;

    case CONDITION:
      bounds.first  = 1.;
      bounds.second = 1.3;
      break;

    case JACOBIAN:
      bounds.first  = 0.5;
      bounds.second = 1.155;
      break;
     
    case SIZE:
    case SHAPE:
      bounds.first  = 0.25;
      bounds.second = 1.;
      break;   

    case DISTORTION:
      bounds.first  = 0.6;
      bounds.second = 1.;
      break;
      
    default:
      std::cout << "Warning: Invalid quality measure chosen." << std::endl;
      bounds.first  = -1;
      bounds.second = -1;
    }

  return bounds;
}
