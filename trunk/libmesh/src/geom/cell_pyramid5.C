// $Id: cell_pyramid5.C,v 1.21 2006-06-20 15:20:26 jwpeterson Exp $

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
#include "side.h"
#include "cell_pyramid5.h"
#include "edge_edge2.h"
#include "face_tri3.h"
#include "face_quad4.h"




// ------------------------------------------------------------
// Pyramid5 class static member initializations
const unsigned int Pyramid5::side_nodes_map[5][4] =
{
  {0, 1, 4, 99}, // Side 0
  {1, 2, 4, 99}, // Side 1
  {2, 3, 4, 99}, // Side 2
  {3, 0, 4, 99}, // Side 3
  {0, 3, 2,  1}  // Side 4
};

const unsigned int Pyramid5::edge_nodes_map[8][2] =
{
  {0, 1}, // Side 0
  {1, 2}, // Side 1
  {2, 3}, // Side 2
  {0, 3}, // Side 3
  {0, 4}, // Side 4
  {1, 4}, // Side 5
  {2, 4}, // Side 6
  {3, 4}  // Side 7
};



// ------------------------------------------------------------
// Pyramid5 class member functions

bool Pyramid5::is_vertex(const unsigned int) const
{
  return true;
}

bool Pyramid5::is_edge(const unsigned int) const
{
  return false;
}

bool Pyramid5::is_face(const unsigned int) const
{
  return false;
}

bool Pyramid5::is_node_on_side(const unsigned int n,
			       const unsigned int s) const
{
  assert(s < n_sides());
  for (unsigned int i = 0; i != 4; ++i)
    if (side_nodes_map[s][i] == n)
      return true;
  return false;
}

bool Pyramid5::is_node_on_edge(const unsigned int n,
			       const unsigned int e) const
{
  assert(e < n_edges());
  for (unsigned int i = 0; i != 2; ++i)
    if (edge_nodes_map[e][i] == n)
      return true;
  return false;
}



bool Pyramid5::has_affine_map() const
{
  return (this->point(3) - this->point(0) ==
      this->point(2) - this->point(1));
}



AutoPtr<Elem> Pyramid5::build_side (const unsigned int i) const
{
  assert (i < this->n_sides());


  
  switch (i)
    {
    case 0:  // triangular face 1
      {
	AutoPtr<Elem> face(new Side<Tri3,Pyramid5>(this,i));

// 	face->set_node(0) = this->get_node(0);
// 	face->set_node(1) = this->get_node(1);
// 	face->set_node(2) = this->get_node(4);

	return face;
      }
    case 1:  // triangular face 2
      {
	AutoPtr<Elem> face(new Side<Tri3,Pyramid5>(this,i));

// 	face->set_node(0) = this->get_node(1);
// 	face->set_node(1) = this->get_node(2);
// 	face->set_node(2) = this->get_node(4);

	return face;
      }
    case 2:  // triangular face 3
      {
	AutoPtr<Elem> face(new Side<Tri3,Pyramid5>(this,i));

// 	face->set_node(0) = this->get_node(2);
// 	face->set_node(1) = this->get_node(3);
// 	face->set_node(2) = this->get_node(4);

	return face;
      }
    case 3:  // triangular face 4
      {
	AutoPtr<Elem> face(new Side<Tri3,Pyramid5>(this,i));

// 	face->set_node(0) = this->get_node(3);
// 	face->set_node(1) = this->get_node(0);
// 	face->set_node(2) = this->get_node(4);

	return face;
      }
    case 4:  // the quad face at z=0
      {
	AutoPtr<Elem> face(new Side<Quad4,Pyramid5>(this,i));

// 	face->set_node(0) = this->get_node(0);
// 	face->set_node(1) = this->get_node(3);
// 	face->set_node(2) = this->get_node(2);
// 	face->set_node(3) = this->get_node(1);

	return face;
      }
    default:
      {
	error();
      }
    }

  // We'll never get here.
  error();

  AutoPtr<Elem> ap(NULL);  return ap;
}



AutoPtr<Elem> Pyramid5::build_edge (const unsigned int i) const
{
  assert (i < this->n_edges());

  return AutoPtr<Elem>(new SideEdge<Edge2,Pyramid5>(this,i));
}



void Pyramid5::connectivity(const unsigned int sc,
			    const IOPackage iop,
			    std::vector<unsigned int>& conn) const
{
  assert (_nodes != NULL);
  assert (sc < this->n_sub_elem());
  assert (iop != INVALID_IO_PACKAGE);

    switch (iop)
    {
    case TECPLOT:
      {
	conn.resize(8);
	conn[0] = this->node(0)+1;
	conn[1] = this->node(1)+1;
	conn[2] = this->node(2)+1;
	conn[3] = this->node(3)+1;
	conn[4] = this->node(4)+1;
	conn[5] = this->node(4)+1;
	conn[6] = this->node(4)+1;
	conn[7] = this->node(4)+1;
	return;
      }

    case VTK:
      {
	conn.resize(5);
	conn[0] = this->node(0);
	conn[1] = this->node(1);
	conn[2] = this->node(2);
	conn[3] = this->node(3);
	conn[4] = this->node(4);
	return;
      }
      
    default:
      error();
    }

    error();
}


Real Pyramid5::volume () const
{
  // The formula for computing the volume of a pyramid is simple: (1/3)*A_b*h
  // where A_b is the area of the base, and h is the vertical height.
  // The formula does not depend on the shape of the base or the position of
  // the apex relative to the base.  For example, the same formula works for
  // a cone.  The trick is to find out what "h" is, for that we need to know
  // the perpendicular distance from the point which forms the apex to the
  // plane which forms the base of the pyramid.

  // The equation of a plane is Ax + By + Cz + D = 0.  Given 3 points on the
  // plane, we can find A,B,C,D:
  Node* node0 = this->get_node(0);
  Node* node1 = this->get_node(1); 
  Node* node2 = this->get_node(2);
  //Node* node3 = this->get_node(3);
  Node* apex  = this->get_node(4);
			      
  Point x  ((*node0)(0), (*node1)(0), (*node2)(0));
  Point y  ((*node0)(1), (*node1)(1), (*node2)(1));
  Point z  ((*node0)(2), (*node1)(2), (*node2)(2));
  Point one(1.,1.,1.);

  const Real A = one * (  y.cross(z));
  const Real B =  x  * (one.cross(z));
  const Real C =  x  * (y.cross(one));
  const Real D = -x  * (y.cross(  z));

  // To find the distance from the base to the apex node, the formula is
  // simple once you have A,B,C,D.  Note that h<0 since the apex is on the
  // side of the plane opposite the normal vector.
  const Real h = (A*(*apex)(0) + B*(*apex)(1) + C*(*apex)(2) + D) / sqrt(A*A + B*B + C*C);
  assert (h<0.);
  
  // Area of the base: We can just use area formula for the Quad4...
  const Real base_area = this->build_side(4)->volume();
  assert (base_area>0.);
  
  // Finally, ready to return the volume!
  return (-1./3.)*base_area*h;
    
}
