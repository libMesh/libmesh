// $Id: cell_prism15.C,v 1.4 2003-08-07 19:25:31 ddreyer Exp $

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
#include "cell_prism15.h"
#include "face_quad8.h"
#include "face_tri6.h"


// ------------------------------------------------------------
// Prism15 class member functions
AutoPtr<Elem> Prism15::build_side (const unsigned int i) const
{
  assert (i < this->n_sides());


  
  switch (i)
    {
    case 0:  // the triangular face at z=-1
      {
	AutoPtr<Elem> face(new Tri6);

	face->set_node(0) = this->get_node(0);
	face->set_node(1) = this->get_node(2);
	face->set_node(2) = this->get_node(1);
	face->set_node(3) = this->get_node(8);
	face->set_node(4) = this->get_node(7);
	face->set_node(5) = this->get_node(6);

	return face;
      }
    case 1:  // the quad face at y=0
      {
	AutoPtr<Elem> face(new Quad8);
	
	face->set_node(0) = this->get_node(0);
	face->set_node(1) = this->get_node(1);
	face->set_node(2) = this->get_node(4);
	face->set_node(3) = this->get_node(3);
	face->set_node(4) = this->get_node(6);
	face->set_node(5) = this->get_node(10);
	face->set_node(6) = this->get_node(12);
	face->set_node(7) = this->get_node(9);
	
	return face;
      }
    case 2:  // the other quad face
      {
	AutoPtr<Elem> face(new Quad8);

	face->set_node(0) = this->get_node(1);
	face->set_node(1) = this->get_node(2);
	face->set_node(2) = this->get_node(5);
	face->set_node(3) = this->get_node(4);
	face->set_node(4) = this->get_node(7);
	face->set_node(5) = this->get_node(11);
	face->set_node(6) = this->get_node(13);
	face->set_node(7) = this->get_node(10);

	return face;
      }
    case 3: // the quad face at x=0
      {
	AutoPtr<Elem> face(new Quad8);

	face->set_node(0) = this->get_node(2);
	face->set_node(1) = this->get_node(0);
	face->set_node(2) = this->get_node(3);
	face->set_node(3) = this->get_node(5);
	face->set_node(4) = this->get_node(8);
	face->set_node(5) = this->get_node(9);
	face->set_node(6) = this->get_node(14);
	face->set_node(7) = this->get_node(11);
	
	return face;
      }
    case 4: // the triangular face at z=1
      {
	AutoPtr<Elem> face(new Tri6);

	face->set_node(0) = this->get_node(3);
	face->set_node(1) = this->get_node(4);
	face->set_node(2) = this->get_node(5);
	face->set_node(3) = this->get_node(12);
	face->set_node(4) = this->get_node(13);
	face->set_node(5) = this->get_node(14);

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



const std::vector<unsigned int> Prism15::tecplot_connectivity(const unsigned int sc) const
{

  assert (_nodes != NULL);
  assert (sc < this->n_sub_elem());

  std::vector<unsigned int> conn(8);
  
  conn[0] = this->node(0)+1;
  conn[1] = this->node(1)+1;
  conn[2] = this->node(2)+1;
  conn[3] = this->node(2)+1;
  conn[4] = this->node(3)+1;
  conn[5] = this->node(4)+1;
  conn[6] = this->node(5)+1;
  conn[7] = this->node(5)+1;

  return conn;
}



void Prism15::vtk_connectivity(const unsigned int sc,
			       std::vector<unsigned int> *conn) const
{
  assert (_nodes != NULL);
  assert (sc < this->n_sub_elem());
  
  if (conn == NULL)
    conn = new std::vector<unsigned int>;

  conn->resize(6);

  (*conn)[0] = this->node(0);
  (*conn)[1] = this->node(2);
  (*conn)[2] = this->node(1);
  (*conn)[3] = this->node(3);
  (*conn)[4] = this->node(5);
  (*conn)[5] = this->node(4);

  return;
}




unsigned int Prism15::second_order_adjacent_vertex (const unsigned int n,
						    const unsigned int v) const
{ 
  assert (n >= this->n_vertices());
  assert (n <  this->n_nodes());
  return _second_order_adjacent_vertices[n-this->n_vertices()][v]; 
}



const unsigned int Prism15::_second_order_adjacent_vertices[9][2] = 
{
  { 0,  1}, // vertices adjacent to node 6 
  { 1,  2}, // vertices adjacent to node 7 
  { 0,  2}, // vertices adjacent to node 8 

  { 0,  3}, // vertices adjacent to node 9 
  { 1,  4}, // vertices adjacent to node 10 
  { 2,  5}, // vertices adjacent to node 11

  { 3,  4}, // vertices adjacent to node 12
  { 4,  5}, // vertices adjacent to node 13
  { 3,  5}  // vertices adjacent to node 14
};
