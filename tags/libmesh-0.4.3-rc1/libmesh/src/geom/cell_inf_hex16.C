// $Id: cell_inf_hex16.C,v 1.22 2004-07-14 19:23:18 jwpeterson Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2004  Benjamin S. Kirk, John W. Peterson
  
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

// Local includes
#include "libmesh_config.h"

#ifdef ENABLE_INFINITE_ELEMENTS

// C++ includes

// Local includes cont'd
#include "cell_inf_hex16.h"
#include "face_quad8.h"
#include "face_inf_quad6.h"


// ------------------------------------------------------------
// InfHex16 class member functions
AutoPtr<Elem> InfHex16::build_side (const unsigned int i) const
{
  assert (i < this->n_sides());

  // Think of a unit cube: (-1,1) x (-1,1)x (1,1)
  switch (i)
    {
    case 0: // the base face
      {
	AutoPtr<Elem> face(new Quad8);

	// Only here, the face element's normal points inward
	face->set_node(0) = this->get_node(0);
	face->set_node(1) = this->get_node(1);
	face->set_node(2) = this->get_node(2);
	face->set_node(3) = this->get_node(3);
	face->set_node(4) = this->get_node(8);
	face->set_node(5) = this->get_node(9);
	face->set_node(6) = this->get_node(10);
	face->set_node(7) = this->get_node(11);

	return face;
      }

    case 1:  // connecting to another infinite element
      {
	AutoPtr<Elem> face(new InfQuad6);

	face->set_node(0) = this->get_node(0);
	face->set_node(1) = this->get_node(1);
	face->set_node(2) = this->get_node(4);
	face->set_node(3) = this->get_node(5);
	face->set_node(4) = this->get_node(8);
	face->set_node(5) = this->get_node(12);

	return face;
      }

    case 2:  // connecting to another infinite element
      {
	AutoPtr<Elem> face(new InfQuad6);

	face->set_node(0) = this->get_node(1);
	face->set_node(1) = this->get_node(2);
	face->set_node(2) = this->get_node(5);
	face->set_node(3) = this->get_node(6);
	face->set_node(4) = this->get_node(9);
	face->set_node(5) = this->get_node(13);

	return face;
      }

    case 3:  // connecting to another infinite element
      {
	AutoPtr<Elem> face(new InfQuad6);
	
	face->set_node(0) = this->get_node(2);
	face->set_node(1) = this->get_node(3);
	face->set_node(2) = this->get_node(6);
	face->set_node(3) = this->get_node(7);
	face->set_node(4) = this->get_node(10);
	face->set_node(5) = this->get_node(14);

	return face;
      }

    case 4:  // connecting to another infinite element
      {
	AutoPtr<Elem> face(new InfQuad6);

	face->set_node(0) = this->get_node(3);
	face->set_node(1) = this->get_node(0);
	face->set_node(2) = this->get_node(7);
	face->set_node(3) = this->get_node(4);
	face->set_node(4) = this->get_node(11);
	face->set_node(5) = this->get_node(15);

	return face;
      }

    default:
      {
	error();
	AutoPtr<Elem> ap(NULL);  return ap;
      }
    }

  // We'll never get here.
  error();
  AutoPtr<Elem> ap(NULL);  return ap;
}



void InfHex16::connectivity(const unsigned int sc,
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
	switch (sc)
	  {
	  case 0:

	    conn[0] = this->node(0)+1;
	    conn[1] = this->node(1)+1;
	    conn[2] = this->node(2)+1;
	    conn[3] = this->node(3)+1;
	    conn[4] = this->node(4)+1;
	    conn[5] = this->node(5)+1;
	    conn[6] = this->node(6)+1;
	    conn[7] = this->node(7)+1;
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



// void InfHex16::tecplot_connectivity(const unsigned int sc,
// 				    std::vector<unsigned int>& conn) const
// {
//   assert (_nodes != NULL);
//   assert (sc < this->n_sub_elem());

//   // std::vector<unsigned int> conn(8);
//   conn.resize(8);


//   switch (sc)
//     {
//     case 0:

//       conn[0] = this->node(0)+1;
//       conn[1] = this->node(1)+1;
//       conn[2] = this->node(2)+1;
//       conn[3] = this->node(3)+1;
//       conn[4] = this->node(4)+1;
//       conn[5] = this->node(5)+1;
//       conn[6] = this->node(6)+1;
//       conn[7] = this->node(7)+1;

//       return;

      
//     default:
//       error();
      
//     }

//   error();
// }




unsigned short int InfHex16::second_order_adjacent_vertex (const unsigned int n,
							   const unsigned int v) const
{ 
  assert (n >= this->n_vertices());
  assert (n <  this->n_nodes());
  assert (v <  2);
  // note that the _second_order_adjacent_vertices matrix is
  // stored in \p InfHex
  return _second_order_adjacent_vertices[n-this->n_vertices()][v]; 
}





#ifdef ENABLE_AMR

const float InfHex16::_embedding_matrix[4][16][16] =
{
  // embedding matrix for child 0
  {
    //          0           1           2           3           4           5           6           7           8           9          10          11          12          13          14          15 th parent Node
    {         1.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 0th child N.
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 1
    {       -0.25,      -0.25,      -0.25,      -0.25,        0.0,        0.0,        0.0,        0.0,        0.5,        0.5,        0.5,        0.5,        0.0,        0.0,        0.0,        0.0}, // 2
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0}, // 3
    {         0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 4
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0}, // 5
    {         0.0,        0.0,        0.0,        0.0,      -0.25,      -0.25,      -0.25,      -0.25,        0.0,        0.0,        0.0,        0.0,        0.5,        0.5,        0.5,        0.5}, // 6
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0}, // 7
    {       0.375,     -0.125,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 8
    {     -0.1875,    -0.1875,    -0.1875,    -0.1875,        0.0,        0.0,        0.0,        0.0,       0.75,      0.375,       0.25,      0.375,        0.0,        0.0,        0.0,        0.0}, // 9
    {     -0.1875,    -0.1875,    -0.1875,    -0.1875,        0.0,        0.0,        0.0,        0.0,      0.375,       0.25,      0.375,       0.75,        0.0,        0.0,        0.0,        0.0}, // 10
    {       0.375,        0.0,        0.0,     -0.125,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0,        0.0,        0.0}, // 11
    {         0.0,        0.0,        0.0,        0.0,      0.375,     -0.125,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0,        0.0}, // 12
    {         0.0,        0.0,        0.0,        0.0,    -0.1875,    -0.1875,    -0.1875,    -0.1875,        0.0,        0.0,        0.0,        0.0,       0.75,      0.375,       0.25,      0.375}, // 13
    {         0.0,        0.0,        0.0,        0.0,    -0.1875,    -0.1875,    -0.1875,    -0.1875,        0.0,        0.0,        0.0,        0.0,      0.375,       0.25,      0.375,       0.75}, // 14
    {         0.0,        0.0,        0.0,        0.0,      0.375,        0.0,        0.0,     -0.125,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75}  // 15
  },

  // embedding matrix for child 1
  {
    //          0           1           2           3           4           5           6           7           8           9          10          11          12          13          14          15 th parent Node
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 0th child N.
    {         0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 1
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 2
    {       -0.25,      -0.25,      -0.25,      -0.25,        0.0,        0.0,        0.0,        0.0,        0.5,        0.5,        0.5,        0.5,        0.0,        0.0,        0.0,        0.0}, // 3
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0}, // 4
    {         0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 5
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0}, // 6
    {         0.0,        0.0,        0.0,        0.0,      -0.25,      -0.25,      -0.25,      -0.25,        0.0,        0.0,        0.0,        0.0,        0.5,        0.5,        0.5,        0.5}, // 7
    {      -0.125,      0.375,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 8
    {         0.0,      0.375,     -0.125,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 9
    {     -0.1875,    -0.1875,    -0.1875,    -0.1875,        0.0,        0.0,        0.0,        0.0,      0.375,       0.75,      0.375,       0.25,        0.0,        0.0,        0.0,        0.0}, // 10
    {     -0.1875,    -0.1875,    -0.1875,    -0.1875,        0.0,        0.0,        0.0,        0.0,       0.75,      0.375,       0.25,      0.375,        0.0,        0.0,        0.0,        0.0}, // 11
    {         0.0,        0.0,        0.0,        0.0,     -0.125,      0.375,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0,        0.0}, // 12
    {         0.0,        0.0,        0.0,        0.0,        0.0,      0.375,     -0.125,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0}, // 13
    {         0.0,        0.0,        0.0,        0.0,    -0.1875,    -0.1875,    -0.1875,    -0.1875,        0.0,        0.0,        0.0,        0.0,      0.375,       0.75,      0.375,       0.25}, // 14
    {         0.0,        0.0,        0.0,        0.0,    -0.1875,    -0.1875,    -0.1875,    -0.1875,        0.0,        0.0,        0.0,        0.0,       0.75,      0.375,       0.25,      0.375}  // 15
  },

  // embedding matrix for child 2
  {
    //          0           1           2           3           4           5           6           7           8           9          10          11          12          13          14          15 th parent Node
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0}, // 0th child N.
    {       -0.25,      -0.25,      -0.25,      -0.25,        0.0,        0.0,        0.0,        0.0,        0.5,        0.5,        0.5,        0.5,        0.0,        0.0,        0.0,        0.0}, // 1
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 2
    {         0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 3
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0}, // 4
    {         0.0,        0.0,        0.0,        0.0,      -0.25,      -0.25,      -0.25,      -0.25,        0.0,        0.0,        0.0,        0.0,        0.5,        0.5,        0.5,        0.5}, // 5
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0}, // 6
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 7
    {     -0.1875,    -0.1875,    -0.1875,    -0.1875,        0.0,        0.0,        0.0,        0.0,      0.375,       0.25,      0.375,       0.75,        0.0,        0.0,        0.0,        0.0}, // 8
    {     -0.1875,    -0.1875,    -0.1875,    -0.1875,        0.0,        0.0,        0.0,        0.0,       0.25,      0.375,       0.75,      0.375,        0.0,        0.0,        0.0,        0.0}, // 9
    {         0.0,        0.0,     -0.125,      0.375,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0,        0.0,        0.0,        0.0}, // 10
    {      -0.125,        0.0,        0.0,      0.375,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0,        0.0,        0.0}, // 11
    {         0.0,        0.0,        0.0,        0.0,    -0.1875,    -0.1875,    -0.1875,    -0.1875,        0.0,        0.0,        0.0,        0.0,      0.375,       0.25,      0.375,       0.75}, // 12
    {         0.0,        0.0,        0.0,        0.0,    -0.1875,    -0.1875,    -0.1875,    -0.1875,        0.0,        0.0,        0.0,        0.0,       0.25,      0.375,       0.75,      0.375}, // 13
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,     -0.125,      0.375,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0}, // 14
    {         0.0,        0.0,        0.0,        0.0,     -0.125,        0.0,        0.0,      0.375,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75}  // 15
  },

  // embedding matrix for child 3
  {
    //          0           1           2           3           4           5           6           7           8           9          10          11          12          13          14          15 th parent Node
    {       -0.25,      -0.25,      -0.25,      -0.25,        0.0,        0.0,        0.0,        0.0,        0.5,        0.5,        0.5,        0.5,        0.0,        0.0,        0.0,        0.0}, // 0th child N.
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 1
    {         0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 2
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 3
    {         0.0,        0.0,        0.0,        0.0,      -0.25,      -0.25,      -0.25,      -0.25,        0.0,        0.0,        0.0,        0.0,        0.5,        0.5,        0.5,        0.5}, // 4
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0}, // 5
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 6
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0}, // 7
    {     -0.1875,    -0.1875,    -0.1875,    -0.1875,        0.0,        0.0,        0.0,        0.0,      0.375,       0.75,      0.375,       0.25,        0.0,        0.0,        0.0,        0.0}, // 8
    {         0.0,     -0.125,      0.375,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 9
    {         0.0,        0.0,      0.375,     -0.125,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0,        0.0,        0.0,        0.0}, // 10
    {     -0.1875,    -0.1875,    -0.1875,    -0.1875,        0.0,        0.0,        0.0,        0.0,       0.25,      0.375,       0.75,      0.375,        0.0,        0.0,        0.0,        0.0}, // 11
    {         0.0,        0.0,        0.0,        0.0,    -0.1875,    -0.1875,    -0.1875,    -0.1875,        0.0,        0.0,        0.0,        0.0,      0.375,       0.75,      0.375,       0.25}, // 12
    {         0.0,        0.0,        0.0,        0.0,        0.0,     -0.125,      0.375,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0}, // 13
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,      0.375,     -0.125,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0}, // 14
    {         0.0,        0.0,        0.0,        0.0,    -0.1875,    -0.1875,    -0.1875,    -0.1875,        0.0,        0.0,        0.0,        0.0,       0.25,      0.375,       0.75,      0.375}  // 15
  }
};



#endif

#endif  // ifdef ENABLE_INFINITE_ELEMENTS
