// $Id: face_quad.C,v 1.12 2003-08-18 14:44:52 ddreyer Exp $

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
#include <iterator>
#include <algorithm>

// Local includes
#include "face_quad.h"
#include "edge_edge2.h"


// ------------------------------------------------------------
// Quad class member functions
unsigned int Quad::key (const unsigned int s) const
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
			   this->node(3));
	
    case 3:
      return
	this->compute_key (this->node(3),
			   this->node(0));	
    }


  // We will never get here...  Look at the code above.
  error();  
  return 0;
}



AutoPtr<Elem> Quad::side (const unsigned int i) const
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
	edge->set_node(1) = this->get_node(3);
	
	return edge;
      }
    case 3:
      {
	edge->set_node(0) = this->get_node(3);
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





Real Quad::quality (const ElemQuality) const
{
  return 0.0; // Not implemented
}






std::pair<Real, Real> Quad::qual_bounds (const ElemQuality q) const
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




const unsigned short int Quad::_second_order_adjacent_vertices[4][2] = 
{
  {0, 1}, // vertices adjacent to node 4 
  {1, 2}, // vertices adjacent to node 5 
  {2, 3}, // vertices adjacent to node 6 
  {0, 3}  // vertices adjacent to node 7 
};





#ifdef ENABLE_AMR

const unsigned int Quad::_side_children_matrix[4][2] =
{
  {0, 1}, // side-0 children
  {1, 3}, // side-1 children
  {2, 3}, // side-2 children
  {0, 2}  // side-3 children
};

#endif
