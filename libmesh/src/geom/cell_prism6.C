// $Id: cell_prism6.C,v 1.22 2005-05-06 17:06:58 roystgnr Exp $

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
#include "cell_prism6.h"
#include "face_quad4.h"
#include "face_tri3.h"



// ------------------------------------------------------------
// Prism6 class static member initializations
const unsigned int Prism6::side_nodes_map[5][4] =
{
  {0, 2, 1, 99}, // Side 0
  {0, 1, 4,  3}, // Side 1
  {1, 2, 5,  4}, // Side 2
  {2, 0, 3,  5}, // Side 3
  {3, 4, 5, 99}  // Side 4
};

const unsigned int Prism6::edge_nodes_map[9][2] =
{
  {0, 1}, // Side 0
  {1, 2}, // Side 1
  {0, 2}, // Side 2
  {0, 3}, // Side 3
  {1, 4}, // Side 4
  {2, 5}, // Side 5
  {3, 4}, // Side 6
  {4, 5}, // Side 7
  {3, 5}  // Side 8
};


// ------------------------------------------------------------
// Prism6 class member functions

bool Prism6::is_vertex(const unsigned int) const
{
  return true;
}

bool Prism6::is_edge(const unsigned int) const
{
  return false;
}

bool Prism6::is_face(const unsigned int) const
{
  return false;
}

bool Prism6::is_node_on_side(const unsigned int n,
			     const unsigned int s) const
{
  assert(s < n_sides());
  for (unsigned int i = 0; i != 4; ++i)
    if (side_nodes_map[s][i] == n)
      return true;
  return false;
}

bool Prism6::is_node_on_edge(const unsigned int n,
			     const unsigned int e) const
{
  assert(e < n_edges());
  for (unsigned int i = 0; i != 2; ++i)
    if (edge_nodes_map[e][i] == n)
      return true;
  return false;
}

AutoPtr<Elem> Prism6::build_side (const unsigned int i) const
{
  assert (i < this->n_sides());

  switch (i)
    {
    case 0:  // the triangular face at z=-1
      {
	AutoPtr<Elem> face(new Side<Tri3,Prism6>(this,i));

// 	face->set_node(0) = this->get_node(0);
// 	face->set_node(1) = this->get_node(2);
// 	face->set_node(2) = this->get_node(1);

	return face;
      }
    case 1:  // the quad face at y=0
      {
	AutoPtr<Elem> face(new Side<Quad4,Prism6>(this,i));
	
// 	face->set_node(0) = this->get_node(0);
// 	face->set_node(1) = this->get_node(1);
// 	face->set_node(2) = this->get_node(4);
// 	face->set_node(3) = this->get_node(3);
	
	return face;
      }
    case 2:  // the other quad face
      {
	AutoPtr<Elem> face(new Side<Quad4,Prism6>(this,i));

// 	face->set_node(0) = this->get_node(1);
// 	face->set_node(1) = this->get_node(2);
// 	face->set_node(2) = this->get_node(5);
// 	face->set_node(3) = this->get_node(4);

	return face;
      }
    case 3: // the quad face at x=0
      {
	AutoPtr<Elem> face(new Side<Quad4,Prism6>(this,i));

// 	face->set_node(0) = this->get_node(2);
// 	face->set_node(1) = this->get_node(0);
// 	face->set_node(2) = this->get_node(3);
// 	face->set_node(3) = this->get_node(5);
	
	return face;
      }
    case 4: // the triangular face at z=1
      {
	AutoPtr<Elem> face(new Side<Tri3,Prism6>(this,i));

// 	face->set_node(0) = this->get_node(3);
// 	face->set_node(1) = this->get_node(4);
// 	face->set_node(2) = this->get_node(5);

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




void Prism6::connectivity(const unsigned int sc,
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



#ifdef ENABLE_AMR

const float Prism6::_embedding_matrix[8][6][6] =
{
  // embedding matrix for child 0
  {
    //  0     1     2     3     4     5
    { 1.0,  0.0,  0.0,  0.0,  0.0,  0.0}, // 0
    { 0.5,  0.5,  0.0,  0.0,  0.0,  0.0}, // 1
    { 0.5,  0.0,  0.5,  0.0,  0.0,  0.0}, // 2
    { 0.5,  0.0,  0.0,  0.5,  0.0,  0.0}, // 3
    { .25,  .25,  0.0,  .25,  .25,  0.0}, // 4
    { .25,  0.0,  .25,  .25,  0.0,  .25}  // 5
  },

  // embedding matrix for child 1
  {
    //  0     1     2     3     4     5
    { 0.5,  0.5,  0.0,  0.0,  0.0,  0.0}, // 0
    { 0.0,  1.0,  0.0,  0.0,  0.0,  0.0}, // 1
    { 0.0,  0.5,  0.5,  0.0,  0.0,  0.0}, // 2
    { .25,  .25,  0.0,  .25,  .25,  0.0}, // 3
    { 0.0,  0.5,  0.0,  0.0,  0.5,  0.0}, // 4
    { 0.0,  .25,  .25,  0.0,  .25,  .25}  // 5
  },

  // embedding matrix for child 2
  {
    //  0     1     2     3     4     5
    { 0.5,  0.0,  0.5,  0.0,  0.0,  0.0}, // 0
    { 0.0,  0.5,  0.5,  0.0,  0.0,  0.0}, // 1
    { 0.0,  0.0,  1.0,  0.0,  0.0,  0.0}, // 2
    { .25,  0.0,  .25,  .25,  0.0,  .25}, // 3
    { 0.0,  .25,  .25,  0.0,  .25,  .25}, // 4
    { 0.0,  0.0,  0.5,  0.0,  0.0,  0.5}  // 5
  },

  // embedding matrix for child 3
  {
    //  0     1     2     3     4     5
    { 0.5,  0.5,  0.0,  0.0,  0.0,  0.0}, // 0
    { 0.0,  0.5,  0.5,  0.0,  0.0,  0.0}, // 1
    { 0.5,  0.0,  0.5,  0.0,  0.0,  0.0}, // 2
    { .25,  .25,  0.0,  .25,  .25,  0.0}, // 3
    { 0.0,  .25,  .25,  0.0,  .25,  .25}, // 4
    { .25,  0.0,  .25,  .25,  0.0,  .25}  // 5
  },

  // embedding matrix for child 4
  {
    //  0     1     2     3     4     5
    { 0.5,  0.0,  0.0,  0.5,  0.0,  0.0}, // 0
    { .25,  .25,  0.0,  .25,  .25,  0.0}, // 1
    { .25,  0.0,  .25,  .25,  0.0,  .25}, // 2
    { 0.0,  0.0,  0.0,  1.0,  0.0,  0.0}, // 3
    { 0.0,  0.0,  0.0,  0.5,  0.5,  0.0}, // 4
    { 0.0,  0.0,  0.0,  0.5,  0.0,  0.5}  // 5
  },

  // embedding matrix for child 5
  {
    //  0     1     2     3     4     5
    { .25,  .25,  0.0,  .25,  .25,  0.0}, // 0
    { 0.0,  0.5,  0.0,  0.0,  0.5,  0.0}, // 1
    { 0.0,  .25,  .25,  0.0,  .25,  .25}, // 2
    { 0.0,  0.0,  0.0,  0.5,  0.5,  0.0}, // 3
    { 0.0,  0.0,  0.0,  0.0,  1.0,  0.0}, // 4
    { 0.0,  0.0,  0.0,  0.0,  0.5,  0.5}  // 5
  },

  // embedding matrix for child 6
  {
    //  0     1     2     3     4     5
    { .25,  0.0,  .25,  .25,  0.0,  .25}, // 0
    { 0.0,  .25,  .25,  0.0,  .25,  .25}, // 1
    { 0.0,  0.0,  0.5,  0.0,  0.0,  0.5}, // 2
    { 0.0,  0.0,  0.0,  0.5,  0.0,  0.5}, // 3
    { 0.0,  0.0,  0.0,  0.0,  0.5,  0.5}, // 4
    { 0.0,  0.0,  0.0,  0.0,  0.0,  1.0}  // 5
  },

  // embedding matrix for child 7
  {
    //  0     1     2     3     4     5
    { .25,  .25,  0.0,  .25,  .25,  0.0}, // 0
    { 0.0,  .25,  .25,  0.0,  .25,  .25}, // 1
    { .25,  0.0,  .25,  .25,  0.0,  .25}, // 2
    { 0.0,  0.0,  0.0,  0.5,  0.5,  0.0}, // 3
    { 0.0,  0.0,  0.0,  0.0,  0.5,  0.5}, // 4
    { 0.0,  0.0,  0.0,  0.5,  0.0,  0.5}  // 5
  }
};

#endif
