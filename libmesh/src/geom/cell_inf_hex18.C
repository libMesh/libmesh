// $Id: cell_inf_hex18.C,v 1.12 2003-02-26 04:43:14 jwpeterson Exp $

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
#include "mesh.h"
#include "cell_inf_hex18.h"
#include "face_quad9.h"
#include "face_inf_quad6.h"


// ------------------------------------------------------------
// InfHex18 class member functions
AutoPtr<Elem> InfHex18::build_side (const unsigned int i) const
{
  assert (i < this->n_sides());


  
  // Think of a unit cube: (-1,1) x (-1,1) x (1,1)
  switch (i)
    {
    case 0:  // the face at z=0, base face
      {
        AutoPtr<Elem> face(new Quad9);

	face->set_node(0) = this->get_node(0);
	face->set_node(1) = this->get_node(1);
	face->set_node(2) = this->get_node(2);
	face->set_node(3) = this->get_node(3);
	face->set_node(4) = this->get_node(8);
	face->set_node(5) = this->get_node(9);
	face->set_node(6) = this->get_node(10);
	face->set_node(7) = this->get_node(11);
	face->set_node(8) = this->get_node(16);
/* old code
	face->set_node(0) = this->get_node(0);
	face->set_node(1) = this->get_node(3);
	face->set_node(2) = this->get_node(2);
	face->set_node(3) = this->get_node(1);
	face->set_node(4) = this->get_node(11);
	face->set_node(5) = this->get_node(10);
	face->set_node(6) = this->get_node(9);
	face->set_node(7) = this->get_node(8);
	face->set_node(8) = this->get_node(16);
*/

	return face;
      }

    case 1:  // the face at y = 0
      {
        AutoPtr<Elem> face(new InfQuad6);
	
	face->set_node(0) = this->get_node(0);
	face->set_node(1) = this->get_node(1);
	face->set_node(2) = this->get_node(5);
	face->set_node(3) = this->get_node(4);
	face->set_node(4) = this->get_node(8);
	face->set_node(5) = this->get_node(12);

	return face;
      }

    case 2:  // the face at x=1
      {
        AutoPtr<Elem> face(new InfQuad6);

	face->set_node(0) = this->get_node(1);
	face->set_node(1) = this->get_node(2);
	face->set_node(2) = this->get_node(6);
	face->set_node(3) = this->get_node(5);
	face->set_node(4) = this->get_node(9);
	face->set_node(5) = this->get_node(13);

	return face;
      }

    case 3: // the face at y=1
      {
        AutoPtr<Elem> face(new InfQuad6);
	
	face->set_node(0) = this->get_node(2);
	face->set_node(1) = this->get_node(3);
	face->set_node(2) = this->get_node(7);
	face->set_node(3) = this->get_node(6);
	face->set_node(4) = this->get_node(10);
	face->set_node(5) = this->get_node(14);

	return face;
      }

    case 4: // the face at x=0
      {
        AutoPtr<Elem> face(new InfQuad6);

	face->set_node(0) = this->get_node(3);
	face->set_node(1) = this->get_node(0);
	face->set_node(2) = this->get_node(4);
	face->set_node(3) = this->get_node(7);
	face->set_node(4) = this->get_node(11);
	face->set_node(5) = this->get_node(15);

	return face;
      }

    case 5: // the face at z=1
      // disable this face, since this is supposed to lie at infinity
      {
        std::cerr << "No face 5 in case of infinite elements!" << std::endl;
        error();
	AutoPtr<Elem> ap(NULL);  return ap;

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



void InfHex18::write_tecplot_connectivity(std::ostream &out) const
{
  assert (out);
  assert (_nodes != NULL);

  for (unsigned int sc=0; sc <this->n_sub_elem(); sc++)
    {
      std::vector<unsigned int> conn = tecplot_connectivity(sc);

      for (unsigned int i=0; i<conn.size(); i++)
	out << conn[i] << " ";

      out << std::endl;
    }
}



#ifdef ENABLE_AMR

const float  InfHex18::embedding_matrix[4][18][18] =
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



const unsigned int InfHex18::side_children_matrix[6][5] =
{
  // note different storage scheme
  {4,   0, 1, 2, 3}, // 4 side-0 children
  {2,   0, 1,42,42}, // 2 side-1 children
  {2,   1, 3,42,42}, // 2 side-2 children
  {2,   2, 3,42,42}, // 2 side-3 children
  {2,   0, 2,42,42}, // 2 side-4 children
  {4,   0, 1, 2, 3}  // 4 side-5 children
};



void InfHex18::refine(Mesh& mesh)
{
  assert (this->refinement_flag() == Elem::REFINE);
  assert (this->active());
  assert (_children == NULL);

  // Create my children
  {
    _children = new Elem*[this->n_children()];

    for (unsigned int c=0; c<this->n_children(); c++)
      {
	_children[c] = new InfHex18(this);
	_children[c]->set_refinement_flag() = Elem::JUST_REFINED;
      }
  }


  // Compute new nodal locations
  // and asssign nodes to children
  {
    std::vector<std::vector<Point> >  p(this->n_children());
    
    for (unsigned int c=0; c<this->n_children(); c++)
      p[c].resize(this->child(c)->n_nodes());
    

    // compute new nodal locations
    for (unsigned int c=0; c<this->n_children(); c++)
      for (unsigned int nc=0; nc<this->child(c)->n_nodes(); nc++)
	for (unsigned int n=0; n<this->n_nodes(); n++)
	  if (embedding_matrix[c][nc][n] != 0.)
	    p[c][nc].add_scaled (this->point(n), static_cast<Real>(embedding_matrix[c][nc][n]));
    
    
    // assign nodes to children & add them to the mesh
    for (unsigned int c=0; c<this->n_children(); c++)
      {
	for (unsigned int nc=0; nc<this->child(c)->n_nodes(); nc++)
	  _children[c]->set_node(nc) = mesh.mesh_refinement.add_point(p[c][nc]);

	mesh.add_elem(this->child(c), mesh.mesh_refinement.new_element_number());
      }
  }


  
  // Possibly add boundary information
  {
    for (unsigned int s=0; s<this->n_sides(); s++)
      if (this->neighbor(s) == NULL)
	{
	  const short int id = mesh.boundary_info.boundary_id(this, s);
	
	  if (id != mesh.boundary_info.invalid_id)
	    // the upper limit for sc is stored in the 0th column
	    for (unsigned int sc=1; sc <=side_children_matrix[s][0]; sc++)
	      mesh.boundary_info.add_side(this->child(side_children_matrix[s][sc]), s, id);
	}
  }


  // Un-set my refinement flag now
  this->set_refinement_flag() = Elem::DO_NOTHING;
}


#endif

#endif
