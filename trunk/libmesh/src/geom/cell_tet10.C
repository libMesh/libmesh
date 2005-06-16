// $Id: cell_tet10.C,v 1.23 2005-06-16 23:03:52 roystgnr Exp $

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
#include "cell_tet10.h"
#include "edge_edge3.h"
#include "face_tri6.h"



// ------------------------------------------------------------
// Tet10 class static member initializations
const unsigned int Tet10::side_nodes_map[4][6] =
{
  {0, 2, 1, 6, 5, 4}, // Side 0
  {0, 1, 3, 4, 8, 7}, // Side 1
  {1, 2, 3, 5, 9, 8}, // Side 2
  {2, 0, 3, 6, 7, 9}  // Side 3
};

const unsigned int Tet10::edge_nodes_map[6][3] =
{
  {0, 1, 4}, // Side 0
  {1, 2, 5}, // Side 1
  {0, 2, 6}, // Side 2
  {0, 3, 7}, // Side 3
  {1, 3, 8}, // Side 4
  {2, 3, 9}  // Side 5
};



// ------------------------------------------------------------
// Tet10 class member functions

bool Tet10::is_vertex(const unsigned int i) const
{
  if (i < 4)
    return true;
  return false;
}

bool Tet10::is_edge(const unsigned int i) const
{
  if (i < 4)
    return false;
  return true;
}

bool Tet10::is_face(const unsigned int) const
{
  return false;
}

bool Tet10::is_node_on_side(const unsigned int n,
			    const unsigned int s) const
{
  assert(s < n_sides());
  for (unsigned int i = 0; i != 6; ++i)
    if (side_nodes_map[s][i] == n)
      return true;
  return false;
}

bool Tet10::is_node_on_edge(const unsigned int n,
			    const unsigned int e) const
{
  assert(e < n_edges());
  for (unsigned int i = 0; i != 3; ++i)
    if (edge_nodes_map[e][i] == n)
      return true;
  return false;
}



bool Tet10::has_affine_map() const
{
  // Make sure edges are straight
  if ((this->point(0) + this->point(1))/2 != this->point(4))
    return false;
  if ((this->point(1) + this->point(2))/2 != this->point(5))
    return false;
  if ((this->point(2) + this->point(0))/2 != this->point(6))
    return false;
  if ((this->point(3) + this->point(0))/2 != this->point(7))
    return false;
  if ((this->point(3) + this->point(1))/2 != this->point(8))
    return false;
  if ((this->point(3) + this->point(2))/2 != this->point(9))
    return false;
  return true;
}



AutoPtr<Elem> Tet10::build_side (const unsigned int i) const
{
  assert (i < this->n_sides());

  AutoPtr<Elem> ap(new Side<Tri6,Tet10>(this,i));
  return ap;
  
//   AutoPtr<Elem> face(new Tri6);
  
//   switch (i)
//     {
//     case 0:
//       {
// 	face->set_node(0) = this->get_node(0);
// 	face->set_node(1) = this->get_node(2);
// 	face->set_node(2) = this->get_node(1);
// 	face->set_node(3) = this->get_node(6);
// 	face->set_node(4) = this->get_node(5);
// 	face->set_node(5) = this->get_node(4);

// 	return face;
//       }
//     case 1:
//       {
// 	face->set_node(0) = this->get_node(0);
// 	face->set_node(1) = this->get_node(1);
// 	face->set_node(2) = this->get_node(3);
// 	face->set_node(3) = this->get_node(4);
// 	face->set_node(4) = this->get_node(8);
// 	face->set_node(5) = this->get_node(7);

// 	return face;
//       }
//     case 2:
//       {
// 	face->set_node(0) = this->get_node(1);
// 	face->set_node(1) = this->get_node(2);
// 	face->set_node(2) = this->get_node(3);
// 	face->set_node(3) = this->get_node(5);
// 	face->set_node(4) = this->get_node(9);
// 	face->set_node(5) = this->get_node(8);

// 	return face;
//       }
//     case 3:
//       {
// 	face->set_node(0) = this->get_node(2);
// 	face->set_node(1) = this->get_node(0);
// 	face->set_node(2) = this->get_node(3);
// 	face->set_node(3) = this->get_node(6);
// 	face->set_node(4) = this->get_node(7);
// 	face->set_node(5) = this->get_node(9);

// 	return face;
//       }
//     default:
//       {
// 	error();
//       }
//     }

//   // We'll never get here.
//   error();
//   return face;
}



AutoPtr<Elem> Tet10::build_edge (const unsigned int i) const
{
  assert (i < this->n_edges());

  return AutoPtr<Elem>(new SideEdge<Edge3,Tet10>(this,i));
}



void Tet10::connectivity(const unsigned int sc,
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
	switch (sc)
	  {
      
      
	    // Linear sub-tet 0
	  case 0:
      
	    conn[0] = this->node(0)+1;
	    conn[1] = this->node(4)+1;
	    conn[2] = this->node(6)+1;
	    conn[3] = this->node(6)+1;
	    conn[4] = this->node(7)+1;
	    conn[5] = this->node(7)+1;
	    conn[6] = this->node(7)+1;
	    conn[7] = this->node(7)+1;

	    return;

	    // Linear sub-tet 1
	  case 1:
      
	    conn[0] = this->node(4)+1;
	    conn[1] = this->node(1)+1;
	    conn[2] = this->node(5)+1;
	    conn[3] = this->node(5)+1;
	    conn[4] = this->node(8)+1;
	    conn[5] = this->node(8)+1;
	    conn[6] = this->node(8)+1;
	    conn[7] = this->node(8)+1;

	    return;

	    // Linear sub-tet 2
	  case 2:
      
	    conn[0] = this->node(5)+1;
	    conn[1] = this->node(2)+1;
	    conn[2] = this->node(6)+1;
	    conn[3] = this->node(6)+1;
	    conn[4] = this->node(9)+1;
	    conn[5] = this->node(9)+1;
	    conn[6] = this->node(9)+1;
	    conn[7] = this->node(9)+1;

	    return;

	    // Linear sub-tet 3
	  case 3:
      
	    conn[0] = this->node(7)+1;
	    conn[1] = this->node(8)+1;
	    conn[2] = this->node(9)+1;
	    conn[3] = this->node(9)+1;
	    conn[4] = this->node(3)+1;
	    conn[5] = this->node(3)+1;
	    conn[6] = this->node(3)+1;
	    conn[7] = this->node(3)+1;

	    return;

	    // Linear sub-tet 4
	  case 4:
      
	    conn[0] = this->node(4)+1;
	    conn[1] = this->node(8)+1;
	    conn[2] = this->node(6)+1;
	    conn[3] = this->node(6)+1;
	    conn[4] = this->node(7)+1;
	    conn[5] = this->node(7)+1;
	    conn[6] = this->node(7)+1;
	    conn[7] = this->node(7)+1;

	    return;

	    // Linear sub-tet 5
	  case 5:
      
	    conn[0] = this->node(4)+1;
	    conn[1] = this->node(5)+1;
	    conn[2] = this->node(6)+1;
	    conn[3] = this->node(6)+1;
	    conn[4] = this->node(8)+1;
	    conn[5] = this->node(8)+1;
	    conn[6] = this->node(8)+1;
	    conn[7] = this->node(8)+1;

	    return;

	    // Linear sub-tet 6
	  case 6:
      
	    conn[0] = this->node(5)+1;
	    conn[1] = this->node(9)+1;
	    conn[2] = this->node(6)+1;
	    conn[3] = this->node(6)+1;
	    conn[4] = this->node(8)+1;
	    conn[5] = this->node(8)+1;
	    conn[6] = this->node(8)+1;
	    conn[7] = this->node(8)+1;

	    return;

	    // Linear sub-tet 7
	  case 7:
      
	    conn[0] = this->node(7)+1;
	    conn[1] = this->node(6)+1;
	    conn[2] = this->node(9)+1;
	    conn[3] = this->node(9)+1;
	    conn[4] = this->node(8)+1;
	    conn[5] = this->node(8)+1;
	    conn[6] = this->node(8)+1;
	    conn[7] = this->node(8)+1;

	    return;


	  default:

	    error();
	  }
      }

    case VTK:
      {
	conn.resize(4);
	switch (sc)
	  {            
	    // Linear sub-tet 0
	  case 0:
      
	    conn[0] = this->node(0);
	    conn[1] = this->node(4);
	    conn[2] = this->node(6);
	    conn[3] = this->node(7);

	    return;

	    // Linear sub-tet 1
	  case 1:
      
	    conn[0] = this->node(4);
	    conn[1] = this->node(1);
	    conn[2] = this->node(5);
	    conn[3] = this->node(8);

	    return;

	    // Linear sub-tet 2
	  case 2:
      
	    conn[0] = this->node(5);
	    conn[1] = this->node(2);
	    conn[2] = this->node(6);
	    conn[3] = this->node(9);

	    return;

	    // Linear sub-tet 3
	  case 3:
      
	    conn[0] = this->node(7);
	    conn[1] = this->node(8);
	    conn[2] = this->node(9);
	    conn[3] = this->node(3);

	    return;

	    // Linear sub-tet 4
	  case 4:
      
	    conn[0] = this->node(4);
	    conn[1] = this->node(8);
	    conn[2] = this->node(6);
	    conn[3] = this->node(7);

	    return;

	    // Linear sub-tet 5
	  case 5:
      
	    conn[0] = this->node(4);
	    conn[1] = this->node(5);
	    conn[2] = this->node(6);
	    conn[3] = this->node(8);

	    return;

	    // Linear sub-tet 6
	  case 6:
      
	    conn[0] = this->node(5);
	    conn[1] = this->node(9);
	    conn[2] = this->node(6);
	    conn[3] = this->node(8);

	    return;

	    // Linear sub-tet 7
	  case 7:
      
	    conn[0] = this->node(7);
	    conn[1] = this->node(6);
	    conn[2] = this->node(9);
	    conn[3] = this->node(8);

	    return;


	  default:

	    error();
	  }
      }

    default:
      error();
    }

  error();
}




unsigned short int Tet10::second_order_adjacent_vertex (const unsigned int n,
							const unsigned int v) const
{ 
  assert (n >= this->n_vertices());
  assert (n <  this->n_nodes());
  assert (v < 2);
  return _second_order_adjacent_vertices[n-this->n_vertices()][v]; 
}



const unsigned short int Tet10::_second_order_adjacent_vertices[6][2] = 
{
  {0, 1}, // vertices adjacent to node 4 
  {1, 2}, // vertices adjacent to node 5 
  {0, 2}, // vertices adjacent to node 6 
  {0, 3}, // vertices adjacent to node 7 
  {1, 3}, // vertices adjacent to node 8 
  {2, 3}  // vertices adjacent to node 9 
};





#ifdef ENABLE_AMR

const float Tet10::_embedding_matrix[8][10][10] =
{
  // embedding matrix for child 0
  {
    //    0      1      2      3      4      5      6      7      8      9
    {    1.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.}, // 0
    {    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.,    0.,    0.}, // 1
    {    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.}, // 2
    {    0.,    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.,    0.}, // 3
    { 0.375,-0.125,    0.,    0.,  0.75,    0.,    0.,    0.,    0.,    0.}, // 4
    {    0.,-0.125,-0.125,    0.,   0.5,  0.25,   0.5,    0.,    0.,    0.}, // 5
    { 0.375,    0.,-0.125,    0.,    0.,    0.,  0.75,    0.,    0.,    0.}, // 6
    { 0.375,    0.,    0.,-0.125,    0.,    0.,    0.,  0.75,    0.,    0.}, // 7
    {    0.,-0.125,    0.,-0.125,   0.5,    0.,    0.,   0.5,  0.25,    0.}, // 8
    {    0.,    0.,-0.125,-0.125,    0.,    0.,   0.5,   0.5,    0.,  0.25}  // 9
  },

  // embedding matrix for child 1
  {
    //    0      1      2      3      4      5      6      7      8      9
    {    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.,    0.,    0.}, // 0
    {    0.,    1.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.}, // 1
    {    0.,    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.,    0.}, // 2
    {    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.}, // 3
    {-0.125, 0.375,    0.,    0.,  0.75,    0.,    0.,    0.,    0.,    0.}, // 4
    {    0., 0.375,-0.125,    0.,    0.,  0.75,    0.,    0.,    0.,    0.}, // 5
    {-0.125,    0.,-0.125,    0.,   0.5,   0.5,  0.25,    0.,    0.,    0.}, // 6
    {-0.125,    0.,    0.,-0.125,   0.5,    0.,    0.,  0.25,   0.5,    0.}, // 7
    {    0., 0.375,    0.,-0.125,    0.,    0.,    0.,    0.,  0.75,    0.}, // 8
    {    0.,    0.,-0.125,-0.125,    0.,   0.5,    0.,    0.,   0.5,  0.25}  // 9
  },

  // embedding matrix for child 2
  {
    //    0      1      2      3      4      5      6      7      8      9
    {    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.}, // 0
    {    0.,    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.,    0.}, // 1
    {    0.,    0.,    1.,    0.,    0.,    0.,    0.,    0.,    0.,    0.}, // 2
    {    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    1.}, // 3
    {-0.125,-0.125,    0.,    0.,  0.25,   0.5,   0.5,    0.,    0.,    0.}, // 4
    {    0.,-0.125, 0.375,    0.,    0.,  0.75,    0.,    0.,    0.,    0.}, // 5
    {-0.125,    0., 0.375,    0.,    0.,    0.,  0.75,    0.,    0.,    0.}, // 6
    {-0.125,    0.,    0.,-0.125,    0.,    0.,   0.5,  0.25,    0.,   0.5}, // 7
    {    0.,-0.125,    0.,-0.125,    0.,   0.5,    0.,    0.,  0.25,   0.5}, // 8
    {    0.,    0., 0.375,-0.125,    0.,    0.,    0.,    0.,    0.,  0.75}  // 9
  },

  // embedding matrix for child 3
  {
    //    0      1      2      3      4      5      6      7      8      9
    {    0.,    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.,    0.}, // 0
    {    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.}, // 1
    {    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    1.}, // 2
    {    0.,    0.,    0.,    1.,    0.,    0.,    0.,    0.,    0.,    0.}, // 3
    {-0.125,-0.125,    0.,    0.,  0.25,    0.,    0.,   0.5,   0.5,    0.}, // 4
    {    0.,-0.125,-0.125,    0.,    0.,  0.25,    0.,    0.,   0.5,   0.5}, // 5
    {-0.125,    0.,-0.125,    0.,    0.,    0.,  0.25,   0.5,    0.,   0.5}, // 6
    {-0.125,    0.,    0., 0.375,    0.,    0.,    0.,  0.75,    0.,    0.}, // 7
    {    0.,-0.125,    0., 0.375,    0.,    0.,    0.,    0.,  0.75,    0.}, // 8
    {    0.,    0.,-0.125, 0.375,    0.,    0.,    0.,    0.,    0.,  0.75}  // 9
  },

  // embedding matrix for child 4
  {
    //    0      1      2      3      4      5      6      7      8      9
    {    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.,    0.,    0.}, // 0
    {    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.}, // 1
    {    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.}, // 2
    {    0.,    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.,    0.}, // 3
    {-0.125,    0.,    0.,-0.125,   0.5,    0.,    0.,  0.25,   0.5,    0.}, // 4
    {-0.125,-0.125,-0.125,-0.125,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25}, // 5
    {    0.,-0.125,-0.125,    0.,   0.5,  0.25,   0.5,    0.,    0.,    0.}, // 6
    {    0.,-0.125,    0.,-0.125,   0.5,    0.,    0.,   0.5,  0.25,    0.}, // 7
    {-0.125,-0.125,    0.,    0.,  0.25,    0.,    0.,   0.5,   0.5,    0.}, // 8
    {    0.,    0.,-0.125,-0.125,    0.,    0.,   0.5,   0.5,    0.,  0.25}  // 9
  },

  // embedding matrix for child 5
  {
    //    0      1      2      3      4      5      6      7      8      9
    {    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.,    0.,    0.}, // 0
    {    0.,    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.,    0.}, // 1
    {    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.}, // 2
    {    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.}, // 3
    {-0.125,    0.,-0.125,    0.,   0.5,   0.5,  0.25,    0.,    0.,    0.}, // 4
    {-0.125,-0.125,    0.,    0.,  0.25,   0.5,   0.5,    0.,    0.,    0.}, // 5
    {    0.,-0.125,-0.125,    0.,   0.5,  0.25,   0.5,    0.,    0.,    0.}, // 6
    {-0.125,    0.,    0.,-0.125,   0.5,    0.,    0.,  0.25,   0.5,    0.}, // 7
    {    0.,    0.,-0.125,-0.125,    0.,   0.5,    0.,    0.,   0.5,  0.25}, // 8
    {-0.125,-0.125,-0.125,-0.125,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25}  // 9
  },

  // embedding matrix for child 6
  {
    //    0      1      2      3      4      5      6      7      8      9
    {    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.}, // 0
    {    0.,    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.,    0.}, // 1
    {    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    1.}, // 2
    {    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.}, // 3
    {-0.125,-0.125,    0.,    0.,  0.25,   0.5,   0.5,    0.,    0.,    0.}, // 4
    {    0.,-0.125,    0.,-0.125,    0.,   0.5,    0.,    0.,  0.25,   0.5}, // 5
    {-0.125,    0.,    0.,-0.125,    0.,    0.,   0.5,  0.25,    0.,   0.5}, // 6
    {-0.125,-0.125,-0.125,-0.125,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25}, // 7
    {    0.,    0.,-0.125,-0.125,    0.,   0.5,    0.,    0.,   0.5,  0.25}, // 8
    {    0.,-0.125,-0.125,    0.,    0.,  0.25,    0.,    0.,   0.5,   0.5}  // 9
  },

  // embedding matrix for child 7
  {
    //    0      1      2      3      4      5      6      7      8      9
    {    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.}, // 0
    {    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.}, // 1
    {    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    1.}, // 2
    {    0.,    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.,    0.}, // 3
    {-0.125,-0.125,-0.125,-0.125,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25}, // 4
    {    0.,-0.125,-0.125,    0.,    0.,  0.25,    0.,    0.,   0.5,   0.5}, // 5
    {-0.125,    0.,    0.,-0.125,    0.,    0.,   0.5,  0.25,    0.,   0.5}, // 6
    {    0.,    0.,-0.125,-0.125,    0.,    0.,   0.5,   0.5,    0.,  0.25}, // 7
    {-0.125,-0.125,    0.,    0.,  0.25,    0.,    0.,   0.5,   0.5,    0.}, // 8
    {-0.125,    0.,-0.125,    0.,    0.,    0.,  0.25,   0.5,    0.,   0.5}  // 9
  }
};

#endif // #ifdef ENABLE_AMR
