// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2007  Benjamin S. Kirk, John W. Peterson
  
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


// Local includes
#include "libmesh_config.h"
#ifdef ENABLE_INFINITE_ELEMENTS

// C++ includes

// Local includes cont'd
#include "face_inf_quad.h"
#include "edge_edge2.h"
#include "edge_inf_edge2.h"


// ------------------------------------------------------------
// InfQuad class member functions
unsigned int InfQuad::key (const unsigned int s) const
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
			   this->node(3));	

    case 2:

      return
	this->compute_key (this->node(0),
			   this->node(2));
    }

  // We will never get here...  Look at the code above.
  error();
  return 0;
}



AutoPtr<DofObject> InfQuad::side (const unsigned int i) const
{
  assert (i < this->n_sides());


  switch (i)
    {
    case 0:
      {
	// base face
	Edge2* edge = new Edge2;

	edge->set_node(0) = this->get_node(0);
	edge->set_node(1) = this->get_node(1);
	
	AutoPtr<DofObject> ap(edge);  return ap;
      }

    case 1:
      {
	// adjacent to another infinite element
        InfEdge2* edge = new InfEdge2;

	edge->set_node(0) = this->get_node(1);
	edge->set_node(1) = this->get_node(3);	

	AutoPtr<DofObject> ap(edge);  return ap;
      }

    case 2:
      {
	// adjacent to another infinite element	
	InfEdge2* edge = new InfEdge2;

	edge->set_node(0) = this->get_node(0); // be aware of swapped nodes,
	edge->set_node(1) = this->get_node(2); // compared to conventional side numbering

	AutoPtr<DofObject> ap(edge);  return ap;
      }

    default:
      {
	error();
	AutoPtr<DofObject> ap(NULL);  return ap;
      }
    }

  // We will never get here...  Look at the code above.
  error();
  AutoPtr<DofObject> ap(NULL);  return ap;
}





Real InfQuad::quality (const ElemQuality) const
{
  return 0.; // Not implemented
}




std::pair<Real, Real> InfQuad::qual_bounds (const ElemQuality q) const
{
  std::pair<Real, Real> bounds;
  
  switch (q)
    {

    case ASPECT_RATIO:
      bounds.first  = 1.;
      bounds.second = 4.;
      break;
      
    case SKEW:
      bounds.first  = 0.;
      bounds.second = 0.5;
      break;

    case TAPER:
      bounds.first  = 0.;
      bounds.second = 0.7;
      break;

    case WARP:
      bounds.first  = 0.9;
      bounds.second = 1.;
      break;
     
    case STRETCH:
      bounds.first  = 0.25;
      bounds.second = 1.;
      break;

    case MIN_ANGLE:
      bounds.first  = 45.;
      bounds.second = 90.;
      break;
      
    case MAX_ANGLE:
      bounds.first  = 90.;
      bounds.second = 135.;
      break;
      
    case CONDITION:
      bounds.first  = 1.;
      bounds.second = 4.;
      break;

    case JACOBIAN:
      bounds.first  = 0.5;
      bounds.second = 1.;
      break;
      
    case SHEAR:
    case SHAPE:
    case SIZE:
      bounds.first  = 0.3;
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



#endif // ifdef ENABLE_INFINITE_ELEMENTS
