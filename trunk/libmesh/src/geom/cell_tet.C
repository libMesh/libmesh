// $Id: cell_tet.C,v 1.5 2003-01-24 17:24:43 jwpeterson Exp $

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
#include "cell_tet.h"




// ------------------------------------------------------------
// Tet class member functions
AutoPtr<Elem> Tet::side (const unsigned int i) const
{
  assert (i < n_sides());


  
  AutoPtr<Elem> face(Elem::build(TRI3));
  
  switch (i)
    {
    case 0:
      {
	face->set_node(0) = get_node(0);
	face->set_node(1) = get_node(2);
	face->set_node(2) = get_node(1);

	return face;
      }
    case 1:
      {
	face->set_node(0) = get_node(0);
	face->set_node(1) = get_node(1);
	face->set_node(2) = get_node(3);

	return face;
      }
    case 2:
      {
	face->set_node(0) = get_node(1);
	face->set_node(1) = get_node(2);
	face->set_node(2) = get_node(3);

	return face;
      }
    case 3:
      {
	face->set_node(0) = get_node(2);
	face->set_node(1) = get_node(0);
	face->set_node(2) = get_node(3);
	
	return face;
      }
    default:
      {
	error();
      }
    };

  // We'll never get here.
  error();
  return face;
};



real Tet::quality(const ElemQuality) const
{
  return -1.0; // Not implemented
}




std::pair<real, real> Tet::qual_bounds (const ElemQuality q) const
{
  std::pair<real, real> bounds;
  
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
};
