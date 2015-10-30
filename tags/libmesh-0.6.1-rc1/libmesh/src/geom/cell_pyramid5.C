// $Id: cell_pyramid5.C,v 1.26 2007-10-03 22:09:24 roystgnr Exp $

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
//  Point v = this->point(3) - this->point(0);
//  return (v.relative_fuzzy_equals(this->point(2) - this->point(1)));
  return false;
}



AutoPtr<Elem> Pyramid5::build_side (const unsigned int i,
				    bool proxy) const
{
  assert (i < this->n_sides());

  if (proxy)
    {
      switch (i)
	{
	case 0:  
	case 1:
	case 2:
	case 3:
	  {
	    AutoPtr<Elem> face(new Side<Tri3,Pyramid5>(this,i));
	    return face;
	  }

	case 4:
	  {
	    AutoPtr<Elem> face(new Side<Quad4,Pyramid5>(this,i));
	    return face;
	  }

	default:
	  {
	    error();
	  }
	}
    }

  else
    {
      switch (i)
	{
	case 0:  // triangular face 1
	  {
	    AutoPtr<Elem> face(new Tri3);

	    face->set_node(0) = this->get_node(0);
	    face->set_node(1) = this->get_node(1);
	    face->set_node(2) = this->get_node(4);

	    return face;
	  }
	case 1:  // triangular face 2
	  {
	    AutoPtr<Elem> face(new Tri3);

	    face->set_node(0) = this->get_node(1);
	    face->set_node(1) = this->get_node(2);
	    face->set_node(2) = this->get_node(4);

	    return face;
	  }
	case 2:  // triangular face 3
	  {
	    AutoPtr<Elem> face(new Tri3);

	    face->set_node(0) = this->get_node(2);
	    face->set_node(1) = this->get_node(3);
	    face->set_node(2) = this->get_node(4);

	    return face;
	  }
	case 3:  // triangular face 4
	  {
	    AutoPtr<Elem> face(new Tri3);

	    face->set_node(0) = this->get_node(3);
	    face->set_node(1) = this->get_node(0);
	    face->set_node(2) = this->get_node(4);

	    return face;
	  }
	case 4:  // the quad face at z=0
	  {
	    AutoPtr<Elem> face(new Quad4);

	    face->set_node(0) = this->get_node(0);
	    face->set_node(1) = this->get_node(3);
	    face->set_node(2) = this->get_node(2);
	    face->set_node(3) = this->get_node(1);

	    return face;
	  }
	default:
	  {
	    error();
	  }
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
	conn[0] = this->node(3);
	conn[1] = this->node(2);
	conn[2] = this->node(1);
	conn[3] = this->node(0);
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
  // The pyramid with a bilinear base has volume given by the
  // formula in: "Calculation of the Volume of a General Hexahedron
  // for Flow Predictions", AIAA Journal v.23, no.6, 1984, p.954-
  //
  // Note: the volume returned by summing the Jacobian times the
  // quadrature weights over all the quadrature points for the
  // pyramid element does *not* give the correct volume (the formula
  // in this function gives the correct volume).  This implies there
  // is something wrong with the first-order Lagrange shape functions
  // or the quadrature rules for the Pyramid5.
  Node* node0 = this->get_node(0);
  Node* node1 = this->get_node(1); 
  Node* node2 = this->get_node(2);
  Node* node3 = this->get_node(3);
  Node* node4 = this->get_node(4);

  // Construct Various edge and diagonal vectors
  Point v40 ( *node0 - *node4 );
  Point v13 ( *node3 - *node1 );
  Point v02 ( *node2 - *node0 );
  Point v03 ( *node3 - *node0 );
  Point v01 ( *node1 - *node0 );
  
  // Finally, ready to return the volume!
  return (1./6.)*(v40*(v13.cross(v02))) + (1./12.)*(v02*(v01.cross(v03)));
}
