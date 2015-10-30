// $Id: cell_tet.C,v 1.1.1.1 2003-01-10 16:17:48 libmesh Exp $

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
#include "cell_tet.h"
#include "face_tri3.h"




// ------------------------------------------------------------
// Tet class member functions
Elem Tet::side (const unsigned int i) const
{
  assert (i < n_sides());
  assert (_nodes.size() == n_nodes());

  Elem face(3);
  
  switch (i)
    {
    case 0:
      {
	face.node(0) = node(0);
	face.node(1) = node(2);
	face.node(2) = node(1);

	return face;
      }
    case 1:
      {
	face.node(0) = node(0);
	face.node(1) = node(1);
	face.node(2) = node(3);

	return face;
      }
    case 2:
      {
	face.node(0) = node(1);
	face.node(1) = node(2);
	face.node(2) = node(3);

	return face;
      }
    case 3:
      {
	face.node(0) = node(2);
	face.node(1) = node(0);
	face.node(2) = node(3);
	
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



real Tet::quality(const MeshBase& mesh, const ElemQuality q) const
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
