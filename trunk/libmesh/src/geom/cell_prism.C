// $Id: cell_prism.C,v 1.11 2003-08-18 14:44:52 ddreyer Exp $

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
#include "cell_prism.h"
#include "face_quad4.h"
#include "face_tri3.h"


// ------------------------------------------------------------
// Prism class member functions
unsigned int Prism::key (const unsigned int s) const
{
  assert (s < this->n_sides());

  switch (s)
    {
    case 0:  // the triangular face at z=0

      return	
	this->compute_key (this->node(0),
			   this->node(2),
			   this->node(1));
      
    case 1:  // the quad face at y=0

      return
	this->compute_key (this->node(0),
			   this->node(1),
			   this->node(4),
			   this->node(3));
	
    case 2:  // the other quad face

      return
	this->compute_key (this->node(1),
			   this->node(2),
			   this->node(5),
			   this->node(4));

    case 3: // the quad face at x=0

      return
	this->compute_key (this->node(2),
			   this->node(0),
			   this->node(3),
			   this->node(5));	
    case 4: // the triangular face at z=1

      return
	this->compute_key (this->node(3),
			   this->node(4),
			   this->node(5));
    }
  
  // We'll never get here.
  error();
  return 0;
}



AutoPtr<Elem> Prism::side (const unsigned int i) const
{
  assert (i < this->n_sides());


  
  AutoPtr<Elem> faceq(new Quad4);
  AutoPtr<Elem> facet(new Tri3);

  switch (i)
    {
    case 0:  // the triangular face at z=0
      {
	facet->set_node(0) = this->get_node(0);
	facet->set_node(1) = this->get_node(2);
	facet->set_node(2) = this->get_node(1);

	return facet;
      }
    case 1:  // the quad face at y=0
      {
	faceq->set_node(0) = this->get_node(0);
	faceq->set_node(1) = this->get_node(1);
	faceq->set_node(2) = this->get_node(4);
	faceq->set_node(3) = this->get_node(3);
	
	return faceq;
      }
    case 2:  // the other quad face
      {
	faceq->set_node(0) = this->get_node(1);
	faceq->set_node(1) = this->get_node(2);
	faceq->set_node(2) = this->get_node(5);
	faceq->set_node(3) = this->get_node(4);

	return faceq;
      }
    case 3: // the quad face at x=0
      {
	faceq->set_node(0) = this->get_node(2);
	faceq->set_node(1) = this->get_node(0);
	faceq->set_node(2) = this->get_node(3);
	faceq->set_node(3) = this->get_node(5);
	
	return faceq;
      }
    case 4: // the triangular face at z=1
      {
	facet->set_node(0) = this->get_node(3);
	facet->set_node(1) = this->get_node(4);
	facet->set_node(2) = this->get_node(5);

	return facet;
      }
    default:
      {
	error();
	return facet;
      }
    }

  // We'll never get here.
  error();
  return facet;
}




const unsigned short int Prism::_second_order_adjacent_vertices[9][2] = 
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



#ifdef ENABLE_AMR

const unsigned int Prism::_side_children_matrix[5][4] =
{
  {0, 1, 2, 3}, // side-0 children
  {0, 1, 4, 5}, // side-1 children
  {1, 2, 5, 6}, // side-2 children
  {0, 2, 4, 6}, // side-3 children
  {4, 5, 6, 7}  // side-4 children
};

#endif
