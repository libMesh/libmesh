// $Id: face_inf_quad.C,v 1.1 2003-03-11 00:47:45 ddreyer Exp $

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


// Local includes
#include "mesh_config.h"
#ifdef ENABLE_INFINITE_ELEMENTS

// C++ includes

// Local includes cont'd
#include "face_inf_quad.h"
#include "edge_edge2.h"
#include "edge_inf_edge2.h"


// ------------------------------------------------------------
// InfQuad class member functions
AutoPtr<Elem> InfQuad::side (const unsigned int i) const
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
	
	AutoPtr<Elem> ap(edge);  return ap;
      }

    case 1:
      {
	// adjacent to another infinite element
        InfEdge2* edge = new InfEdge2;

	edge->set_node(0) = this->get_node(1);
	edge->set_node(1) = this->get_node(2);	

	AutoPtr<Elem> ap(edge);  return ap;
      }

    case 2:
      {
	// adjacent to another infinite element	
	InfEdge2* edge = new InfEdge2;

	edge->set_node(0) = this->get_node(0); // be aware of swapped nodes,
	edge->set_node(1) = this->get_node(3); // compared to conventional side numbering

	AutoPtr<Elem> ap(edge);  return ap;
      }

    default:
      {
	error();
	AutoPtr<Elem> ap(NULL);  return ap;
      }
    }

  // We will never get here...  Look at the code above.
  error();
  AutoPtr<Elem> ap(NULL);  return ap;
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



#ifdef ENABLE_AMR

const unsigned int InfQuad::_side_children_matrix[3][2] =
{
  // note different storage scheme
  {0,  1}, // 2 side-0 children
  {1, 42}, // 1 side-1 children
  {0, 42}  // 1 side-2 children
};

#endif



#endif // ifdef ENABLE_INFINITE_ELEMENTS
