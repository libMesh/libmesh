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

// C++ includes

// Local includes
#include "face_tri.h"
#include "edge_edge2.h"





// ------------------------------------------------------------
// Tri class member functions
unsigned int Tri::key (const unsigned int s) const
{
  libmesh_assert (s < this->n_sides());

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
  libmesh_error();
  return 0;
}



AutoPtr<DofObject> Tri::side (const unsigned int i) const
{
  libmesh_assert (i < this->n_sides());

  Elem* edge = new Edge2;

  switch (i)
    {
    case 0:
      {
	edge->set_node(0) = this->get_node(0);
	edge->set_node(1) = this->get_node(1);
	
        AutoPtr<DofObject> ap_edge(edge);
	return ap_edge;
      }
    case 1:
      {
	edge->set_node(0) = this->get_node(1);
	edge->set_node(1) = this->get_node(2);
	
        AutoPtr<DofObject> ap_edge(edge);
	return ap_edge;
      }
    case 2:
      {
	edge->set_node(0) = this->get_node(2);
	edge->set_node(1) = this->get_node(0);
	
        AutoPtr<DofObject> ap_edge(edge);
	return ap_edge;
      }
    default:
      {
	libmesh_error();
      }
    }

  
  // We will never get here...  Look at the code above.
  libmesh_error();
  AutoPtr<DofObject> ap_edge(edge);
  return ap_edge;
}



bool Tri::is_child_on_side(const unsigned int c,
                           const unsigned int s) const
{
  libmesh_assert (c < this->n_children());
  libmesh_assert (s < this->n_sides());

  return (c == s || c == (s+1)%3);
}



Real Tri::quality (const ElemQuality q) const
{
  switch (q)
    {
      
      /**
       * Source: Netgen, meshtool.cpp, TriangleQualityInst
       */
    case DISTORTION:
    case STRETCH:
      {
        const Node* p1 = this->get_node(0);
        const Node* p2 = this->get_node(1);
        const Node* p3 = this->get_node(2);

        Point v1 = (*p2) - (*p1);
        Point v2 = (*p3) - (*p1);
        Point v3 = (*p3) - (*p2);
        const Real l1 = v1.size();
        const Real l2 = v2.size();
        const Real l3 = v3.size();

        // if one length is 0, quality is quite bad!
        if ((l1 <=0.) || (l2 <= 0.) || (l3 <= 0.))
          return 0.;

        const Real s1 = std::sin(std::acos(v1*v2/l1/l2)/2.);
        v1 *= -1;
        const Real s2 = std::sin(std::acos(v1*v3/l1/l3)/2.);
        const Real s3 = std::sin(std::acos(v2*v3/l2/l3)/2.);
        
        return 8. * s1 * s2 * s3;
        
      }
    default:
      return Elem::quality(q);
    }
  
  /**
   * I don't know what to do for this metric. 
   * Maybe the base class knows.  We won't get
   * here because of the defualt case above.
   */
  return Elem::quality(q);
  
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
