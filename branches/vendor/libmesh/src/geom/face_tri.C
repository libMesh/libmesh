// $Id: face_tri.C,v 1.1.1.1 2003-01-10 16:17:48 libmesh Exp $

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
#include "face_tri.h"
#include "mesh_base.h"
#include "mesh_refinement.h"
#include "fe.h"






// ------------------------------------------------------------
// Tri class member functions
Elem Tri::side (const unsigned int i) const
{
  assert (i < n_sides());

  Elem edge(2);

  switch (i)
    {
    case 0:
      {
	edge.node(0) = node(0);
	edge.node(1) = node(1);
	
	return edge;
      }
    case 1:
      {
	edge.node(0) = node(1);
	edge.node(1) = node(2);
	
	return edge;
      }
    case 2:
      {
	edge.node(0) = node(2);
	edge.node(1) = node(0);
	
	return edge;
      }
    default:
      {
	error();
      }
    };

  
  // We will never get here...  Look at the code above.
  return edge;
};






real Tri::quality (const MeshBase& mesh, const ElemQuality q) const
{
  return 0.; // not implemented
}






std::pair<real, real> Tri::qual_bounds (const ElemQuality q) const
{
  std::pair<real, real> bounds;
  
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


