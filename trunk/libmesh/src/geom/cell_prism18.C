// $Id: cell_prism18.C,v 1.20 2006-10-13 03:05:32 roystgnr Exp $

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
#include "cell_prism18.h"
#include "edge_edge3.h"
#include "face_quad9.h"
#include "face_tri6.h"



// ------------------------------------------------------------
// Prism18 class static member initializations
const unsigned int Prism18::side_nodes_map[5][9] =
{
  {0, 2, 1,  8,  7,  6, 99, 99, 99}, // Side 0
  {0, 1, 4,  3,  6, 10, 12,  9, 15}, // Side 1
  {1, 2, 5,  4,  7, 11, 13, 10, 16}, // Side 2
  {2, 0, 3,  5,  8,  9, 14, 11, 17}, // Side 3
  {3, 4, 5, 12, 13, 14, 99, 99, 99}  // Side 4
};

const unsigned int Prism18::edge_nodes_map[9][3] =
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
// Prism18 class member functions

bool Prism18::is_vertex(const unsigned int i) const
{
  if (i < 6)
    return true;
  return false;
}

bool Prism18::is_edge(const unsigned int i) const
{
  if (i < 6)
    return false;
  if (i > 14)
    return false;
  return true;
}

bool Prism18::is_face(const unsigned int i) const
{
  if (i > 14)
    return true;
  return false;
}

bool Prism18::is_node_on_side(const unsigned int n,
			      const unsigned int s) const
{
  assert(s < n_sides());
  for (unsigned int i = 0; i != 9; ++i)
    if (side_nodes_map[s][i] == n)
      return true;
  return false;
}

bool Prism18::is_node_on_edge(const unsigned int n,
			      const unsigned int e) const
{
  assert(e < n_edges());
  for (unsigned int i = 0; i != 3; ++i)
    if (edge_nodes_map[e][i] == n)
      return true;
  return false;
}



bool Prism18::has_affine_map() const
{
  // Make sure z edges are affine
  Point v = this->point(3) - this->point(0);
  if (!v.relative_fuzzy_equals(this->point(4) - this->point(1)) ||
      !v.relative_fuzzy_equals(this->point(5) - this->point(2)))
    return false;
  // Make sure edges are straight
  v /= 2;
  if (!v.relative_fuzzy_equals(this->point(9) - this->point(0)) ||
      !v.relative_fuzzy_equals(this->point(10) - this->point(1)) ||
      !v.relative_fuzzy_equals(this->point(11) - this->point(2)) ||
      !v.relative_fuzzy_equals(this->point(15) - this->point(6)) ||
      !v.relative_fuzzy_equals(this->point(16) - this->point(7)) ||
      !v.relative_fuzzy_equals(this->point(17) - this->point(8)))
    return false;
  v = (this->point(1) - this->point(0))/2;
  if (!v.relative_fuzzy_equals(this->point(6) - this->point(0)) ||
      !v.relative_fuzzy_equals(this->point(12) - this->point(3)))
    return false;
  v = (this->point(2) - this->point(0))/2;
  if (!v.relative_fuzzy_equals(this->point(8) - this->point(0)) ||
      !v.relative_fuzzy_equals(this->point(14) - this->point(3)))
    return false;
  v = (this->point(2) - this->point(1))/2;
  if (!v.relative_fuzzy_equals(this->point(7) - this->point(1)) ||
      !v.relative_fuzzy_equals(this->point(13) - this->point(4)))
    return false;
  return true;
}



AutoPtr<Elem> Prism18::build_side (const unsigned int i) const
{
  assert (i < this->n_sides());


  
  switch (i)
    {
    case 0:  // the triangular face at z=-1
      {
	AutoPtr<Elem> face(new Side<Tri6,Prism18>(this,i));

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
	AutoPtr<Elem> face(new Side<Quad9,Prism18>(this,i));
	
// 	face->set_node(0) = this->get_node(0);
// 	face->set_node(1) = this->get_node(1);
// 	face->set_node(2) = this->get_node(4);
// 	face->set_node(3) = this->get_node(3);
// 	face->set_node(4) = this->get_node(6);
// 	face->set_node(5) = this->get_node(10);
// 	face->set_node(6) = this->get_node(12);
// 	face->set_node(7) = this->get_node(9);
// 	face->set_node(8) = this->get_node(15);
	
	return face;
      }
    case 2:  // the other quad face
      {
	AutoPtr<Elem> face(new Side<Quad9,Prism18>(this,i));

// 	face->set_node(0) = this->get_node(1);
// 	face->set_node(1) = this->get_node(2);
// 	face->set_node(2) = this->get_node(5);
// 	face->set_node(3) = this->get_node(4);
// 	face->set_node(4) = this->get_node(7);
// 	face->set_node(5) = this->get_node(11);
// 	face->set_node(6) = this->get_node(13);
// 	face->set_node(7) = this->get_node(10);
// 	face->set_node(8) = this->get_node(16);

	return face;
      }
    case 3: // the quad face at x=0
      {
	AutoPtr<Elem> face(new Side<Quad9,Prism18>(this,i));

// 	face->set_node(0) = this->get_node(2);
// 	face->set_node(1) = this->get_node(0);
// 	face->set_node(2) = this->get_node(3);
// 	face->set_node(3) = this->get_node(5);
// 	face->set_node(4) = this->get_node(8);
// 	face->set_node(5) = this->get_node(9);
// 	face->set_node(6) = this->get_node(14);
// 	face->set_node(7) = this->get_node(11);
// 	face->set_node(8) = this->get_node(17);
	
	return face;
      }
    case 4: // the triangular face at z=1
      {
	AutoPtr<Elem> face(new Side<Tri6,Prism18>(this,i));

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



AutoPtr<Elem> Prism18::build_edge (const unsigned int i) const
{
  assert (i < this->n_edges());

  return AutoPtr<Elem>(new SideEdge<Edge3,Prism18>(this,i));
}



void Prism18::connectivity(const unsigned int sc,
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
      
	  case 0:
	    {
	      conn[0] = this->node(0)+1;
	      conn[1] = this->node(6)+1;
	      conn[2] = this->node(8)+1;
	      conn[3] = this->node(8)+1;
	      conn[4] = this->node(9)+1;
	      conn[5] = this->node(15)+1;
	      conn[6] = this->node(17)+1;
	      conn[7] = this->node(17)+1;
	
	      return;
	    }

	  case 1:
	    {
	      conn[0] = this->node(6)+1;
	      conn[1] = this->node(1)+1;
	      conn[2] = this->node(7)+1;
	      conn[3] = this->node(7)+1;
	      conn[4] = this->node(15)+1;
	      conn[5] = this->node(10)+1;
	      conn[6] = this->node(16)+1;
	      conn[7] = this->node(16)+1;
	
	      return;
	    }
      
	  case 2:
	    {
	      conn[0] = this->node(8)+1;
	      conn[1] = this->node(7)+1;
	      conn[2] = this->node(2)+1;
	      conn[3] = this->node(2)+1;
	      conn[4] = this->node(17)+1;
	      conn[5] = this->node(16)+1;
	      conn[6] = this->node(11)+1;
	      conn[7] = this->node(11)+1;
	
	      return;
	    }
      
	  case 3:
	    {
	      conn[0] = this->node(6)+1;
	      conn[1] = this->node(7)+1;
	      conn[2] = this->node(8)+1;
	      conn[3] = this->node(8)+1;
	      conn[4] = this->node(15)+1;
	      conn[5] = this->node(16)+1;
	      conn[6] = this->node(17)+1;
	      conn[7] = this->node(17)+1;
	
	      return;
	    }

	  case 4:
	    {
	      conn[0] = this->node(9)+1;
	      conn[1] = this->node(15)+1;
	      conn[2] = this->node(17)+1;
	      conn[3] = this->node(17)+1;
	      conn[4] = this->node(3)+1;
	      conn[5] = this->node(12)+1;
	      conn[6] = this->node(14)+1;
	      conn[7] = this->node(14)+1;
	
	      return;
	    }

	  case 5:
	    {
	      conn[0] = this->node(15)+1;
	      conn[1] = this->node(10)+1;
	      conn[2] = this->node(16)+1;
	      conn[3] = this->node(16)+1;
	      conn[4] = this->node(12)+1;
	      conn[5] = this->node(4)+1;
	      conn[6] = this->node(13)+1;
	      conn[7] = this->node(13)+1;
	
	      return;
	    }

	  case 6:
	    {
	      conn[0] = this->node(17)+1;
	      conn[1] = this->node(16)+1;
	      conn[2] = this->node(11)+1;
	      conn[3] = this->node(11)+1;
	      conn[4] = this->node(14)+1;
	      conn[5] = this->node(13)+1;
	      conn[6] = this->node(5)+1;
	      conn[7] = this->node(5)+1;
	
	      return;
	    }

	  case 7:
	    {
	      conn[0] = this->node(15)+1;
	      conn[1] = this->node(16)+1;
	      conn[2] = this->node(17)+1;
	      conn[3] = this->node(17)+1;
	      conn[4] = this->node(12)+1;
	      conn[5] = this->node(13)+1;
	      conn[6] = this->node(14)+1;
	      conn[7] = this->node(14)+1;
	
	      return;
	    }

	  default:
	    error();
	  }

      }

    case VTK:
      {
	conn.resize(6);
	switch (sc)
	  {
      
	  case 0:
	    {
	      conn[0] = this->node(0);
	      conn[1] = this->node(6);
	      conn[2] = this->node(8);
	      conn[3] = this->node(9);
	      conn[4] = this->node(15);
	      conn[5] = this->node(17);
	
	      return;
	    }

	  case 1:
	    {
	      conn[0] = this->node(6);
	      conn[1] = this->node(1);
	      conn[2] = this->node(7);
	      conn[3] = this->node(15);
	      conn[4] = this->node(10);
	      conn[5] = this->node(16);
	
	      return;
	    }
      
	  case 2:
	    {
	      conn[0] = this->node(8);
	      conn[1] = this->node(7);
	      conn[2] = this->node(2);
	      conn[3] = this->node(17);
	      conn[4] = this->node(16);
	      conn[5] = this->node(11);
	
	      return;
	    }
      
	  case 3:
	    {
	      conn[0] = this->node(6);
	      conn[1] = this->node(7);
	      conn[2] = this->node(8);
	      conn[3] = this->node(15);
	      conn[4] = this->node(16);
	      conn[5] = this->node(17);
	
	      return;
	    }

	  case 4:
	    {
	      conn[0] = this->node(9);
	      conn[1] = this->node(15);
	      conn[2] = this->node(17);
	      conn[3] = this->node(3);
	      conn[4] = this->node(12);
	      conn[5] = this->node(14);
	
	      return;
	    }

	  case 5:
	    {
	      conn[0] = this->node(15);
	      conn[1] = this->node(10);
	      conn[2] = this->node(16);
	      conn[3] = this->node(12);
	      conn[4] = this->node(4);
	      conn[5] = this->node(13);
	
	      return;
	    }

	  case 6:
	    {
	      conn[0] = this->node(17);
	      conn[1] = this->node(16);
	      conn[2] = this->node(11);
	      conn[3] = this->node(14);
	      conn[4] = this->node(13);
	      conn[5] = this->node(5);
	
	      return;
	    }

	  case 7:
	    {
	      conn[0] = this->node(15);
	      conn[1] = this->node(16);
	      conn[2] = this->node(17);
	      conn[3] = this->node(12);
	      conn[4] = this->node(13);
	      conn[5] = this->node(14);
	
	      return;
	    }

	  default:
	    error();
	  }

      }

    default:
      error();
    }

  error();

}




unsigned int Prism18::n_second_order_adjacent_vertices (const unsigned int n) const
{
  switch (n)
    {
      case 6:
      case 7:
      case 8:
      case 9:
      case 10:
      case 11:
      case 12:
      case 13:
      case 14:
	return 2;

      case 15:
      case 16:
      case 17:
	return 4;

      default:
	error();
    }
  error();
  return static_cast<unsigned int>(-1);
}





unsigned short int Prism18::second_order_adjacent_vertex (const unsigned int n,
							  const unsigned int v) const
{ 
  assert (n >= this->n_vertices());
  assert (n <  this->n_nodes());

  switch (n)
    {
      /*
       * These nodes are unique to \p Prism18,
       * let our _remaining_... matrix handle
       * this.
       */
      case 15:
      case 16:
      case 17:
      {
	assert (v < 4);
	return _remaining_second_order_adjacent_vertices[n-15][v]; 
      }

      /*
       * All other second-order nodes (6,...,14) are
       * identical with Prism15 and are therefore
       * delegated to the _second_order matrix of
       * \p Prism
       */
      default:
      {
	assert (v < 2);
	return _second_order_adjacent_vertices[n-this->n_vertices()][v]; 
      }

    }

  error();
  return libMesh::invalid_uint;
}



const unsigned short int Prism18::_remaining_second_order_adjacent_vertices[3][4] = 
{
  { 0,  1,  3,  4}, // vertices adjacent to node 15
  { 1,  2,  4,  5}, // vertices adjacent to node 16
  { 0,  2,  3,  5}  // vertices adjacent to node 17
};







#ifdef ENABLE_AMR

const float Prism18::_embedding_matrix[8][18][18] =
{
  // embedding matrix for child 0
  {
    //       0         1         2         3         4         5         6         7         8         9        10        11        12        13        14        15        16        17
    {       1.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 0
    {       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 1
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 2
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 3
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.}, // 4
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.}, // 5
    {    0.375,   -0.125,       0.,       0.,       0.,       0.,     0.75,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 6
    {       0.,   -0.125,   -0.125,       0.,       0.,       0.,      0.5,     0.25,      0.5,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 7
    {    0.375,       0.,   -0.125,       0.,       0.,       0.,       0.,       0.,     0.75,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 8
    {    0.375,       0.,       0.,   -0.125,       0.,       0.,       0.,       0.,       0.,     0.75,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 9
    {       0.,       0.,       0.,       0.,       0.,       0.,    0.375,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,       0.,     0.75,       0.,       0.}, // 10
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,    0.375,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,       0.,     0.75}, // 11
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,    0.375,   -0.125,       0.,       0.,       0.,       0.,     0.75,       0.,       0.}, // 12
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,   -0.125,       0.,       0.,       0.,      0.5,     0.25,      0.5}, // 13
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,    0.375,       0.,   -0.125,       0.,       0.,       0.,       0.,       0.,     0.75}, // 14
    { 0.140625,-0.046875,       0.,-0.046875, 0.015625,       0.,  0.28125,       0.,       0.,  0.28125, -0.09375,       0., -0.09375,       0.,       0.,   0.5625,       0.,       0.}, // 15
    {       0.,-0.046875,-0.046875,       0., 0.015625, 0.015625,   0.1875,  0.09375,   0.1875,       0., -0.09375, -0.09375,  -0.0625, -0.03125,  -0.0625,    0.375,   0.1875,    0.375}, // 16
    { 0.140625,       0.,-0.046875,-0.046875,       0., 0.015625,       0.,       0.,  0.28125,  0.28125,       0., -0.09375,       0.,       0., -0.09375,       0.,       0.,   0.5625}  // 17
  },

  // embedding matrix for child 1
  {
    //       0         1         2         3         4         5         6         7         8         9        10        11        12        13        14        15        16        17
    {       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 0
    {       0.,       1.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 1
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 2
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.}, // 3
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 4
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.}, // 5
    {   -0.125,    0.375,       0.,       0.,       0.,       0.,     0.75,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 6
    {       0.,    0.375,   -0.125,       0.,       0.,       0.,       0.,     0.75,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 7
    {   -0.125,       0.,   -0.125,       0.,       0.,       0.,      0.5,      0.5,     0.25,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 8
    {       0.,       0.,       0.,       0.,       0.,       0.,    0.375,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,       0.,     0.75,       0.,       0.}, // 9
    {       0.,    0.375,       0.,       0.,   -0.125,       0.,       0.,       0.,       0.,       0.,     0.75,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 10
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,    0.375,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,       0.,     0.75,       0.}, // 11
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,    0.375,       0.,       0.,       0.,       0.,     0.75,       0.,       0.}, // 12
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,    0.375,   -0.125,       0.,       0.,       0.,       0.,     0.75,       0.}, // 13
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,   -0.125,       0.,       0.,       0.,      0.5,      0.5,     0.25}, // 14
    {-0.046875, 0.140625,       0., 0.015625,-0.046875,       0.,  0.28125,       0.,       0., -0.09375,  0.28125,       0., -0.09375,       0.,       0.,   0.5625,       0.,       0.}, // 15
    {       0., 0.140625,-0.046875,       0.,-0.046875, 0.015625,       0.,  0.28125,       0.,       0.,  0.28125, -0.09375,       0., -0.09375,       0.,       0.,   0.5625,       0.}, // 16
    {-0.046875,       0.,-0.046875, 0.015625,       0., 0.015625,   0.1875,   0.1875,  0.09375, -0.09375,       0., -0.09375,  -0.0625,  -0.0625, -0.03125,    0.375,    0.375,   0.1875}  // 17
  },

  // embedding matrix for child 2
  {
    //       0         1         2         3         4         5         6         7         8         9        10        11        12        13        14        15        16        17
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 0
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 1
    {       0.,       0.,       1.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 2
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.}, // 3
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.}, // 4
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.,       0.,       0.}, // 5
    {   -0.125,   -0.125,       0.,       0.,       0.,       0.,     0.25,      0.5,      0.5,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 6
    {       0.,   -0.125,    0.375,       0.,       0.,       0.,       0.,     0.75,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 7
    {   -0.125,       0.,    0.375,       0.,       0.,       0.,       0.,       0.,     0.75,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 8
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,    0.375,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,       0.,     0.75}, // 9
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,    0.375,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,       0.,     0.75,       0.}, // 10
    {       0.,       0.,    0.375,       0.,       0.,   -0.125,       0.,       0.,       0.,       0.,       0.,     0.75,       0.,       0.,       0.,       0.,       0.,       0.}, // 11
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,   -0.125,       0.,       0.,       0.,       0.,     0.25,      0.5,      0.5}, // 12
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,    0.375,       0.,       0.,       0.,       0.,     0.75,       0.}, // 13
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,    0.375,       0.,       0.,       0.,       0.,       0.,     0.75}, // 14
    {-0.046875,-0.046875,       0., 0.015625, 0.015625,       0.,  0.09375,   0.1875,   0.1875, -0.09375, -0.09375,       0., -0.03125,  -0.0625,  -0.0625,   0.1875,    0.375,    0.375}, // 15
    {       0.,-0.046875, 0.140625,       0., 0.015625,-0.046875,       0.,  0.28125,       0.,       0., -0.09375,  0.28125,       0., -0.09375,       0.,       0.,   0.5625,       0.}, // 16
    {-0.046875,       0., 0.140625, 0.015625,       0.,-0.046875,       0.,       0.,  0.28125, -0.09375,       0.,  0.28125,       0.,       0., -0.09375,       0.,       0.,   0.5625}  // 17
  },

  // embedding matrix for child 3
  {
    //       0         1         2         3         4         5         6         7         8         9        10        11        12        13        14        15        16        17
    {       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 0
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 1
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 2
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.}, // 3
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.}, // 4
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.}, // 5
    {   -0.125,       0.,   -0.125,       0.,       0.,       0.,      0.5,      0.5,     0.25,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 6
    {   -0.125,   -0.125,       0.,       0.,       0.,       0.,     0.25,      0.5,      0.5,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 7
    {       0.,   -0.125,   -0.125,       0.,       0.,       0.,      0.5,     0.25,      0.5,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 8
    {       0.,       0.,       0.,       0.,       0.,       0.,    0.375,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,       0.,     0.75,       0.,       0.}, // 9
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,    0.375,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,       0.,     0.75,       0.}, // 10
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,    0.375,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,       0.,     0.75}, // 11
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,   -0.125,       0.,       0.,       0.,      0.5,      0.5,     0.25}, // 12
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,   -0.125,       0.,       0.,       0.,       0.,     0.25,      0.5,      0.5}, // 13
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,   -0.125,       0.,       0.,       0.,      0.5,     0.25,      0.5}, // 14
    {-0.046875,       0.,-0.046875, 0.015625,       0., 0.015625,   0.1875,   0.1875,  0.09375, -0.09375,       0., -0.09375,  -0.0625,  -0.0625, -0.03125,    0.375,    0.375,   0.1875}, // 15
    {-0.046875,-0.046875,       0., 0.015625, 0.015625,       0.,  0.09375,   0.1875,   0.1875, -0.09375, -0.09375,       0., -0.03125,  -0.0625,  -0.0625,   0.1875,    0.375,    0.375}, // 16
    {       0.,-0.046875,-0.046875,       0., 0.015625, 0.015625,   0.1875,  0.09375,   0.1875,       0., -0.09375, -0.09375,  -0.0625, -0.03125,  -0.0625,    0.375,   0.1875,    0.375}  // 17
  },

  // embedding matrix for child 4
  {
    //       0         1         2         3         4         5         6         7         8         9        10        11        12        13        14        15        16        17
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 0
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.}, // 1
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.}, // 2
    {       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 3
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.,       0.}, // 4
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.}, // 5
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,    0.375,   -0.125,       0.,       0.,       0.,       0.,     0.75,       0.,       0.}, // 6
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,   -0.125,       0.,       0.,       0.,      0.5,     0.25,      0.5}, // 7
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,    0.375,       0.,   -0.125,       0.,       0.,       0.,       0.,       0.,     0.75}, // 8
    {   -0.125,       0.,       0.,    0.375,       0.,       0.,       0.,       0.,       0.,     0.75,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 9
    {       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,       0.,       0.,       0.,       0.,    0.375,       0.,       0.,     0.75,       0.,       0.}, // 10
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,       0.,       0.,       0.,       0.,    0.375,       0.,       0.,     0.75}, // 11
    {       0.,       0.,       0.,    0.375,   -0.125,       0.,       0.,       0.,       0.,       0.,       0.,       0.,     0.75,       0.,       0.,       0.,       0.,       0.}, // 12
    {       0.,       0.,       0.,       0.,   -0.125,   -0.125,       0.,       0.,       0.,       0.,       0.,       0.,      0.5,     0.25,      0.5,       0.,       0.,       0.}, // 13
    {       0.,       0.,       0.,    0.375,       0.,   -0.125,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,     0.75,       0.,       0.,       0.}, // 14
    {-0.046875, 0.015625,       0., 0.140625,-0.046875,       0., -0.09375,       0.,       0.,  0.28125, -0.09375,       0.,  0.28125,       0.,       0.,   0.5625,       0.,       0.}, // 15
    {       0., 0.015625, 0.015625,       0.,-0.046875,-0.046875,  -0.0625, -0.03125,  -0.0625,       0., -0.09375, -0.09375,   0.1875,  0.09375,   0.1875,    0.375,   0.1875,    0.375}, // 16
    {-0.046875,       0., 0.015625, 0.140625,       0.,-0.046875,       0.,       0., -0.09375,  0.28125,       0., -0.09375,       0.,       0.,  0.28125,       0.,       0.,   0.5625}  // 17
  },

  // embedding matrix for child 5
  {
    //       0         1         2         3         4         5         6         7         8         9        10        11        12        13        14        15        16        17
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.}, // 0
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 1
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.}, // 2
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.,       0.}, // 3
    {       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 4
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.}, // 5
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,    0.375,       0.,       0.,       0.,       0.,     0.75,       0.,       0.}, // 6
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,    0.375,   -0.125,       0.,       0.,       0.,       0.,     0.75,       0.}, // 7
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,   -0.125,       0.,       0.,       0.,      0.5,      0.5,     0.25}, // 8
    {       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,       0.,       0.,       0.,       0.,    0.375,       0.,       0.,     0.75,       0.,       0.}, // 9
    {       0.,   -0.125,       0.,       0.,    0.375,       0.,       0.,       0.,       0.,       0.,     0.75,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 10
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,       0.,       0.,       0.,       0.,    0.375,       0.,       0.,     0.75,       0.}, // 11
    {       0.,       0.,       0.,   -0.125,    0.375,       0.,       0.,       0.,       0.,       0.,       0.,       0.,     0.75,       0.,       0.,       0.,       0.,       0.}, // 12
    {       0.,       0.,       0.,       0.,    0.375,   -0.125,       0.,       0.,       0.,       0.,       0.,       0.,       0.,     0.75,       0.,       0.,       0.,       0.}, // 13
    {       0.,       0.,       0.,   -0.125,       0.,   -0.125,       0.,       0.,       0.,       0.,       0.,       0.,      0.5,      0.5,     0.25,       0.,       0.,       0.}, // 14
    { 0.015625,-0.046875,       0.,-0.046875, 0.140625,       0., -0.09375,       0.,       0., -0.09375,  0.28125,       0.,  0.28125,       0.,       0.,   0.5625,       0.,       0.}, // 15
    {       0.,-0.046875, 0.015625,       0., 0.140625,-0.046875,       0., -0.09375,       0.,       0.,  0.28125, -0.09375,       0.,  0.28125,       0.,       0.,   0.5625,       0.}, // 16
    { 0.015625,       0., 0.015625,-0.046875,       0.,-0.046875,  -0.0625,  -0.0625, -0.03125, -0.09375,       0., -0.09375,   0.1875,   0.1875,  0.09375,    0.375,    0.375,   0.1875}  // 17
  },

  // embedding matrix for child 6
  {
    //       0         1         2         3         4         5         6         7         8         9        10        11        12        13        14        15        16        17
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.}, // 0
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.}, // 1
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.,       0.,       0.}, // 2
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.}, // 3
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.}, // 4
    {       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.}, // 5
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,   -0.125,       0.,       0.,       0.,       0.,     0.25,      0.5,      0.5}, // 6
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,    0.375,       0.,       0.,       0.,       0.,     0.75,       0.}, // 7
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,    0.375,       0.,       0.,       0.,       0.,       0.,     0.75}, // 8
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,       0.,       0.,       0.,       0.,    0.375,       0.,       0.,     0.75}, // 9
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,       0.,       0.,       0.,       0.,    0.375,       0.,       0.,     0.75,       0.}, // 10
    {       0.,       0.,   -0.125,       0.,       0.,    0.375,       0.,       0.,       0.,       0.,       0.,     0.75,       0.,       0.,       0.,       0.,       0.,       0.}, // 11
    {       0.,       0.,       0.,   -0.125,   -0.125,       0.,       0.,       0.,       0.,       0.,       0.,       0.,     0.25,      0.5,      0.5,       0.,       0.,       0.}, // 12
    {       0.,       0.,       0.,       0.,   -0.125,    0.375,       0.,       0.,       0.,       0.,       0.,       0.,       0.,     0.75,       0.,       0.,       0.,       0.}, // 13
    {       0.,       0.,       0.,   -0.125,       0.,    0.375,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,     0.75,       0.,       0.,       0.}, // 14
    { 0.015625, 0.015625,       0.,-0.046875,-0.046875,       0., -0.03125,  -0.0625,  -0.0625, -0.09375, -0.09375,       0.,  0.09375,   0.1875,   0.1875,   0.1875,    0.375,    0.375}, // 15
    {       0., 0.015625,-0.046875,       0.,-0.046875, 0.140625,       0., -0.09375,       0.,       0., -0.09375,  0.28125,       0.,  0.28125,       0.,       0.,   0.5625,       0.}, // 16
    { 0.015625,       0.,-0.046875,-0.046875,       0., 0.140625,       0.,       0., -0.09375, -0.09375,       0.,  0.28125,       0.,       0.,  0.28125,       0.,       0.,   0.5625}  // 17
  },

  // embedding matrix for child 7
  {
    //       0         1         2         3         4         5         6         7         8         9        10        11        12        13        14        15        16        17
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.}, // 0
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.}, // 1
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.}, // 2
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.,       0.}, // 3
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.,       0.}, // 4
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       1.,       0.,       0.,       0.}, // 5
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,   -0.125,       0.,       0.,       0.,      0.5,      0.5,     0.25}, // 6
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,   -0.125,       0.,       0.,       0.,       0.,     0.25,      0.5,      0.5}, // 7
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,   -0.125,       0.,       0.,       0.,      0.5,     0.25,      0.5}, // 8
    {       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,       0.,       0.,       0.,       0.,    0.375,       0.,       0.,     0.75,       0.,       0.}, // 9
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,       0.,       0.,       0.,       0.,    0.375,       0.,       0.,     0.75,       0.}, // 10
    {       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -0.125,       0.,       0.,       0.,       0.,       0.,    0.375,       0.,       0.,     0.75}, // 11
    {       0.,       0.,       0.,   -0.125,       0.,   -0.125,       0.,       0.,       0.,       0.,       0.,       0.,      0.5,      0.5,     0.25,       0.,       0.,       0.}, // 12
    {       0.,       0.,       0.,   -0.125,   -0.125,       0.,       0.,       0.,       0.,       0.,       0.,       0.,     0.25,      0.5,      0.5,       0.,       0.,       0.}, // 13
    {       0.,       0.,       0.,       0.,   -0.125,   -0.125,       0.,       0.,       0.,       0.,       0.,       0.,      0.5,     0.25,      0.5,       0.,       0.,       0.}, // 14
    { 0.015625,       0., 0.015625,-0.046875,       0.,-0.046875,  -0.0625,  -0.0625, -0.03125, -0.09375,       0., -0.09375,   0.1875,   0.1875,  0.09375,    0.375,    0.375,   0.1875}, // 15
    { 0.015625, 0.015625,       0.,-0.046875,-0.046875,       0., -0.03125,  -0.0625,  -0.0625, -0.09375, -0.09375,       0.,  0.09375,   0.1875,   0.1875,   0.1875,    0.375,    0.375}, // 16
    {       0., 0.015625, 0.015625,       0.,-0.046875,-0.046875,  -0.0625, -0.03125,  -0.0625,       0., -0.09375, -0.09375,   0.1875,  0.09375,   0.1875,    0.375,   0.1875,    0.375}  // 17
  }
};

#endif
