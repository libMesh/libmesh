// $Id: cell_inf_hex18.C,v 1.17 2003-07-12 16:33:18 ddreyer Exp $

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

// Local includes
#include "mesh_config.h"

#ifdef ENABLE_INFINITE_ELEMENTS

// C++ includes

// Local includes cont'd
#include "cell_inf_hex18.h"
#include "face_quad9.h"
#include "face_inf_quad6.h"


// ------------------------------------------------------------
// InfHex18 class member functions
unsigned int InfHex18::key (const unsigned int s) const
{
  assert (s < this->n_sides());

  // Think of a unit cube: (-1,1) x (-1,1) x (1,1)
  switch (s)
    {
    case 0: // the base face

      return
	this->compute_key (this->node(16));

          
    case 1:  // the face at y = -1

      return
	this->compute_key (this->node(0),
			   this->node(1),
			   this->node(5),
			   this->node(4));
      
    case 2:  // the face at x = 1

      return
	this->compute_key (this->node(1),
			   this->node(2),
			   this->node(6),
			   this->node(5));

    case 3: // the face at y = 1

      return
	this->compute_key (this->node(2),
			   this->node(3),
			   this->node(7),
			   this->node(6));
      
    case 4: // the face at x = -1
      
      return
	this->compute_key (this->node(3),
			   this->node(0),
			   this->node(4),
			   this->node(7));
    }

  // We'll never get here.
  error();
  return 0;
}



AutoPtr<Elem> InfHex18::build_side (const unsigned int i) const
{
  assert (i < this->n_sides());

  // Think of a unit cube: (-1,1) x (-1,1) x (1,1)
  switch (i)
    {
    case 0: // the base face
      {
        AutoPtr<Elem> face(new Quad9);

	// This is the exception: all other face elements' normals
	// point outwards; but the base element's normal points inward
	face->set_node(0) = this->get_node(0);
	face->set_node(1) = this->get_node(1);
	face->set_node(2) = this->get_node(2);
	face->set_node(3) = this->get_node(3);
	face->set_node(4) = this->get_node(8);
	face->set_node(5) = this->get_node(9);
	face->set_node(6) = this->get_node(10);
	face->set_node(7) = this->get_node(11);
	face->set_node(8) = this->get_node(16);

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



const std::vector<unsigned int> InfHex18::tecplot_connectivity(const unsigned int sc) const
{
  assert (_nodes != NULL);
  assert (sc < this->n_sub_elem());

  std::vector<unsigned int> conn(8);


  switch (sc)
    {
    case 0:

      conn[0] = this->node(0)+1;
      conn[1] = this->node(8)+1;
      conn[2] = this->node(16)+1;
      conn[3] = this->node(11)+1;
      conn[4] = this->node(4)+1;
      conn[5] = this->node(12)+1;
      conn[6] = this->node(17)+1;
      conn[7] = this->node(15)+1;

      return conn;
      
    case 1:

      conn[0] = this->node(8)+1;
      conn[1] = this->node(1)+1;
      conn[2] = this->node(9)+1;
      conn[3] = this->node(16)+1;
      conn[4] = this->node(12)+1;
      conn[5] = this->node(5)+1;
      conn[6] = this->node(13)+1;
      conn[7] = this->node(17)+1;

      return conn;
      
    case 2:

      conn[0] = this->node(11)+1;
      conn[1] = this->node(16)+1;
      conn[2] = this->node(10)+1;
      conn[3] = this->node(3)+1; 
      conn[4] = this->node(15)+1;
      conn[5] = this->node(17)+1;
      conn[6] = this->node(14)+1;
      conn[7] = this->node(7)+1;

      return conn;
      
    case 3:

      conn[0] = this->node(16)+1;
      conn[1] = this->node(9)+1;
      conn[2] = this->node(2)+1;
      conn[3] = this->node(10)+1;
      conn[4] = this->node(17)+1;
      conn[5] = this->node(13)+1;
      conn[6] = this->node(6)+1;
      conn[7] = this->node(14)+1;

      return conn;
      
    default:
      error();
      
    }
  
  return conn;
}





#ifdef ENABLE_AMR

const float InfHex18::_embedding_matrix[4][18][18] =
{
  // embedding matrix for child 0
  {
    //          0           1           2           3           4           5           6           7           8           9          10          11          12          13          14          15          16          17 th parent Node
    {         1.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 0th child N.
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 1
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0}, // 2
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 3
    {         0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 4
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 5
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0}, // 6
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0}, // 7
    {       0.375,     -0.125,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 8
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,      0.375,        0.0,     -0.125,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0}, // 9
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,     -0.125,        0.0,      0.375,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0}, // 10
    {       0.375,        0.0,        0.0,     -0.125,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 11
    {         0.0,        0.0,        0.0,        0.0,      0.375,     -0.125,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0,        0.0,        0.0,        0.0}, // 12
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,      0.375,        0.0,     -0.125,        0.0,        0.0,       0.75}, // 13
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,     -0.125,        0.0,      0.375,        0.0,       0.75}, // 14
    {         0.0,        0.0,        0.0,        0.0,      0.375,        0.0,        0.0,     -0.125,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0}, // 15
    {    0.140625,  -0.046875,   0.015625,  -0.046875,        0.0,        0.0,        0.0,        0.0,    0.28125,   -0.09375,   -0.09375,    0.28125,        0.0,        0.0,        0.0,        0.0,     0.5625,        0.0}, // 16
    {         0.0,        0.0,        0.0,        0.0,   0.140625,  -0.046875,   0.015625,  -0.046875,        0.0,        0.0,        0.0,        0.0,    0.28125,   -0.09375,   -0.09375,    0.28125,        0.0,     0.5625}  // 17
  },

  // embedding matrix for child 1
  {
    //          0           1           2           3           4           5           6           7           8           9          10          11          12          13          14          15          16          17 th parent Node
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 0th child N.
    {         0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 1
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 2
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0}, // 3
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 4
    {         0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 5
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0}, // 6
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0}, // 7
    {      -0.125,      0.375,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 8
    {         0.0,      0.375,     -0.125,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 9
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,      0.375,        0.0,     -0.125,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0}, // 10
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,      0.375,        0.0,     -0.125,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0}, // 11
    {         0.0,        0.0,        0.0,        0.0,     -0.125,      0.375,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0,        0.0,        0.0,        0.0}, // 12
    {         0.0,        0.0,        0.0,        0.0,        0.0,      0.375,     -0.125,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0,        0.0,        0.0}, // 13
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,      0.375,        0.0,     -0.125,        0.0,       0.75}, // 14
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,      0.375,        0.0,     -0.125,        0.0,        0.0,       0.75}, // 15
    {   -0.046875,   0.140625,  -0.046875,   0.015625,        0.0,        0.0,        0.0,        0.0,    0.28125,    0.28125,   -0.09375,   -0.09375,        0.0,        0.0,        0.0,        0.0,     0.5625,        0.0}, // 16
    {         0.0,        0.0,        0.0,        0.0,  -0.046875,   0.140625,  -0.046875,   0.015625,        0.0,        0.0,        0.0,        0.0,    0.28125,    0.28125,   -0.09375,   -0.09375,        0.0,     0.5625}  // 17
  },

  // embedding matrix for child 2
  {
    //          0           1           2           3           4           5           6           7           8           9          10          11          12          13          14          15          16          17 th parent Node
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 0th child N.
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0}, // 1
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 2
    {         0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 3
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0}, // 4
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0}, // 5
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0}, // 6
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 7
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,     -0.125,        0.0,      0.375,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0}, // 8
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,     -0.125,        0.0,      0.375,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0}, // 9
    {         0.0,        0.0,     -0.125,      0.375,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 10
    {      -0.125,        0.0,        0.0,      0.375,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 11
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,     -0.125,        0.0,      0.375,        0.0,       0.75}, // 12
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,     -0.125,        0.0,      0.375,        0.0,        0.0,       0.75}, // 13
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,     -0.125,      0.375,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0,        0.0}, // 14
    {         0.0,        0.0,        0.0,        0.0,     -0.125,        0.0,        0.0,      0.375,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0}, // 15
    {   -0.046875,   0.015625,  -0.046875,   0.140625,        0.0,        0.0,        0.0,        0.0,   -0.09375,   -0.09375,    0.28125,    0.28125,        0.0,        0.0,        0.0,        0.0,     0.5625,        0.0}, // 16
    {         0.0,        0.0,        0.0,        0.0,  -0.046875,   0.015625,  -0.046875,   0.140625,        0.0,        0.0,        0.0,        0.0,   -0.09375,   -0.09375,    0.28125,    0.28125,        0.0,     0.5625}  // 17
  },

  // embedding matrix for child 3
  {
    //          0           1           2           3           4           5           6           7           8           9          10          11          12          13          14          15          16          17 th parent Node
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0}, // 0th child N.
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 1
    {         0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 2
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 3
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0}, // 4
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0}, // 5
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 6
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0}, // 7
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,      0.375,        0.0,     -0.125,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0}, // 8
    {         0.0,     -0.125,      0.375,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 9
    {         0.0,        0.0,      0.375,     -0.125,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 10
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,     -0.125,        0.0,      0.375,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0}, // 11
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,      0.375,        0.0,     -0.125,        0.0,       0.75}, // 12
    {         0.0,        0.0,        0.0,        0.0,        0.0,     -0.125,      0.375,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0,        0.0,        0.0}, // 13
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,      0.375,     -0.125,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0,        0.0}, // 14
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,     -0.125,        0.0,      0.375,        0.0,        0.0,       0.75}, // 15
    {    0.015625,  -0.046875,   0.140625,  -0.046875,        0.0,        0.0,        0.0,        0.0,   -0.09375,    0.28125,    0.28125,   -0.09375,        0.0,        0.0,        0.0,        0.0,     0.5625,        0.0}, // 16
    {         0.0,        0.0,        0.0,        0.0,   0.015625,  -0.046875,   0.140625,  -0.046875,        0.0,        0.0,        0.0,        0.0,   -0.09375,    0.28125,    0.28125,   -0.09375,        0.0,     0.5625}  // 17
  }
};




#endif

#endif // ifdef ENABLE_INFINITE_ELEMENTS

