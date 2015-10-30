// $Id: cell_prism15.C,v 1.15 2005-05-11 18:31:16 roystgnr Exp $

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
#include "cell_prism15.h"
#include "edge_edge3.h"
#include "face_quad8.h"
#include "face_tri6.h"



// ------------------------------------------------------------
// Prism15 class static member initializations
const unsigned int Prism15::side_nodes_map[5][8] =
{
  {0, 2, 1,  8,  7,  6, 99, 99}, // Side 0
  {0, 1, 4,  3,  6, 10, 12,  9}, // Side 1
  {1, 2, 5,  4,  7, 11, 13, 10}, // Side 2
  {2, 0, 3,  5,  8,  9, 14, 11}, // Side 3
  {3, 4, 5, 12, 13, 14, 99, 99}  // Side 4
};

const unsigned int Prism15::edge_nodes_map[9][3] =
{
  {0, 1, 6},  // Side 0
  {1, 2, 7},  // Side 1
  {0, 2, 8},  // Side 2
  {0, 3, 9},  // Side 3
  {1, 4, 10}, // Side 4
  {2, 5, 11}, // Side 5
  {3, 4, 12}, // Side 6
  {4, 5, 13}, // Side 7
  {3, 5, 14}  // Side 8
};


// ------------------------------------------------------------
// Prism15 class member functions

bool Prism15::is_vertex(const unsigned int i) const
{
  if (i < 6)
    return true;
  return false;
}

bool Prism15::is_edge(const unsigned int i) const
{
  if (i < 6)
    return false;
  return true;
}

bool Prism15::is_face(const unsigned int) const
{
  return false;
}

bool Prism15::is_node_on_side(const unsigned int n,
			      const unsigned int s) const
{
  assert(s < n_sides());
  for (unsigned int i = 0; i != 8; ++i)
    if (side_nodes_map[s][i] == n)
      return true;
  return false;
}

bool Prism15::is_node_on_edge(const unsigned int n,
			      const unsigned int e) const
{
  assert(e < n_edges());
  for (unsigned int i = 0; i != 3; ++i)
    if (edge_nodes_map[e][i] == n)
      return true;
  return false;
}

AutoPtr<Elem> Prism15::build_side (const unsigned int i) const
{
  assert (i < this->n_sides());


  
  switch (i)
    {
    case 0:  // the triangular face at z=-1
      {
	AutoPtr<Elem> face(new Side<Tri6,Prism15>(this,i));

// 	face->set_node(0) = this->get_node(0);
// 	face->set_node(1) = this->get_node(2);
// 	face->set_node(2) = this->get_node(1);
// 	face->set_node(3) = this->get_node(8);
// 	face->set_node(4) = this->get_node(7);
// 	face->set_node(5) = this->get_node(6);

	return face;
      }
    case 1:  // the quad face at y=0
      {
	AutoPtr<Elem> face(new Side<Quad8,Prism15>(this,i));
	
// 	face->set_node(0) = this->get_node(0);
// 	face->set_node(1) = this->get_node(1);
// 	face->set_node(2) = this->get_node(4);
// 	face->set_node(3) = this->get_node(3);
// 	face->set_node(4) = this->get_node(6);
// 	face->set_node(5) = this->get_node(10);
// 	face->set_node(6) = this->get_node(12);
// 	face->set_node(7) = this->get_node(9);
	
	return face;
      }
    case 2:  // the other quad face
      {
	AutoPtr<Elem> face(new Side<Quad8,Prism15>(this,i));

// 	face->set_node(0) = this->get_node(1);
// 	face->set_node(1) = this->get_node(2);
// 	face->set_node(2) = this->get_node(5);
// 	face->set_node(3) = this->get_node(4);
// 	face->set_node(4) = this->get_node(7);
// 	face->set_node(5) = this->get_node(11);
// 	face->set_node(6) = this->get_node(13);
// 	face->set_node(7) = this->get_node(10);

	return face;
      }
    case 3: // the quad face at x=0
      {
	AutoPtr<Elem> face(new Side<Quad8,Prism15>(this,i));

// 	face->set_node(0) = this->get_node(2);
// 	face->set_node(1) = this->get_node(0);
// 	face->set_node(2) = this->get_node(3);
// 	face->set_node(3) = this->get_node(5);
// 	face->set_node(4) = this->get_node(8);
// 	face->set_node(5) = this->get_node(9);
// 	face->set_node(6) = this->get_node(14);
// 	face->set_node(7) = this->get_node(11);
	
	return face;
      }
    case 4: // the triangular face at z=1
      {
	AutoPtr<Elem> face(new Side<Tri6,Prism15>(this,i));

// 	face->set_node(0) = this->get_node(3);
// 	face->set_node(1) = this->get_node(4);
// 	face->set_node(2) = this->get_node(5);
// 	face->set_node(3) = this->get_node(12);
// 	face->set_node(4) = this->get_node(13);
// 	face->set_node(5) = this->get_node(14);

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


AutoPtr<Elem> Prism15::build_edge (const unsigned int i) const
{
  assert (i < this->n_edges());

  return AutoPtr<Elem>(new SideEdge<Edge3,Prism15>(this,i));
}


void Prism15::connectivity(const unsigned int sc,
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
	conn[3] = this->node(2)+1;
	conn[4] = this->node(3)+1;
	conn[5] = this->node(4)+1;
	conn[6] = this->node(5)+1;
	conn[7] = this->node(5)+1;
	return;
      }

    case VTK:
      {
	conn.resize(6);
	conn[0] = this->node(0);
	conn[1] = this->node(2);
	conn[2] = this->node(1);
	conn[3] = this->node(3);
	conn[4] = this->node(5);
	conn[5] = this->node(4);
	return;
      }

    default:
      error();
    }

  error();

}




unsigned short int Prism15::second_order_adjacent_vertex (const unsigned int n,
							  const unsigned int v) const
{ 
  assert (n >= this->n_vertices());
  assert (n <  this->n_nodes());
  assert (v < 2);
  return _second_order_adjacent_vertices[n-this->n_vertices()][v]; 
}


