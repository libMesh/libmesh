// $Id: cell_hex8.C,v 1.19 2005-02-22 22:17:38 jwpeterson Exp $

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
#include "cell_hex8.h"
#include "face_quad4.h"




// ------------------------------------------------------------
// Hex8 class static member initializations
const unsigned int Hex8::side_nodes_map[6][4] =
{
  {0, 3, 2, 1}, // Side 0
  {0, 1, 5, 4}, // Side 1
  {1, 2, 6, 5}, // Side 2
  {2, 3, 7, 6}, // Side 3
  {3, 0, 4, 7}, // Side 4
  {4, 5, 6, 7}  // Side 5
};


// ------------------------------------------------------------
// Hex8 class member functions

bool Hex8::is_vertex(const unsigned int) const
{
  return true;
}

bool Hex8::is_edge(const unsigned int) const
{
  return false;
}

bool Hex8::is_face(const unsigned int) const
{
  return false;
}

AutoPtr<Elem> Hex8::build_side (const unsigned int i) const
{
  assert (i < this->n_sides());

  AutoPtr<Elem> ap(new Side<Quad4,Hex8>(this,i));
  return ap;
  
//   AutoPtr<Elem> face(new Quad4);

//   // Think of a unit cube: (-1,1) x (-1,1)x (-1,1)
//   switch (i)
//     {
//     case 0:  // the face at z = -1
//       {
// 	face->set_node(0) = this->get_node(0);
// 	face->set_node(1) = this->get_node(3);
// 	face->set_node(2) = this->get_node(2);
// 	face->set_node(3) = this->get_node(1);

// 	return face;
//       }
//     case 1:  // the face at y = -1
//       {
// 	face->set_node(0) = this->get_node(0);
// 	face->set_node(1) = this->get_node(1);
// 	face->set_node(2) = this->get_node(5);
// 	face->set_node(3) = this->get_node(4);
	
// 	return face;
//       }
//     case 2:  // the face at x = 1
//       {
// 	face->set_node(0) = this->get_node(1);
// 	face->set_node(1) = this->get_node(2);
// 	face->set_node(2) = this->get_node(6);
// 	face->set_node(3) = this->get_node(5);

// 	return face;
//       }
//     case 3: // the face at y = 1
//       {
// 	face->set_node(0) = this->get_node(2);
// 	face->set_node(1) = this->get_node(3);
// 	face->set_node(2) = this->get_node(7);
// 	face->set_node(3) = this->get_node(6);
	
// 	return face;
//       }
//     case 4: // the face at x = -1
//       {
// 	face->set_node(0) = this->get_node(3);
// 	face->set_node(1) = this->get_node(0);
// 	face->set_node(2) = this->get_node(4);
// 	face->set_node(3) = this->get_node(7);

// 	return face;
//       }
//     case 5: // the face at z = 1
//       {
// 	face->set_node(0) = this->get_node(4);
// 	face->set_node(1) = this->get_node(5);
// 	face->set_node(2) = this->get_node(6);
// 	face->set_node(3) = this->get_node(7);
	
// 	return face;
//       }
//     default:
//       {
// 	error();
// 	return face;
//       }
//     }

//   // We'll never get here.
//   error();
//   return face;
}




void Hex8::connectivity(const unsigned int sc,
			const IOPackage iop,
			std::vector<unsigned int>& conn) const
{
  assert (_nodes != NULL);
  assert (sc < this->n_sub_elem());
  assert (iop != INVALID_IO_PACKAGE);

  conn.resize(8);

  switch (iop)
    {
    case TECPLOT:
      {
	conn[0] = this->node(0)+1;
	conn[1] = this->node(1)+1;
	conn[2] = this->node(2)+1;
	conn[3] = this->node(3)+1;
	conn[4] = this->node(4)+1;
	conn[5] = this->node(5)+1;
	conn[6] = this->node(6)+1;
	conn[7] = this->node(7)+1;
	return;
      }

    case VTK:
      {
	conn[0] = this->node(0);
	conn[1] = this->node(1);
	conn[2] = this->node(2);
	conn[3] = this->node(3);
	conn[4] = this->node(4);
	conn[5] = this->node(5);
	conn[6] = this->node(6);
	conn[7] = this->node(7);
      }

    default:
      error();
    }

  error();
}

// void Hex8::tecplot_connectivity(const unsigned int sc,
// 				std::vector<unsigned int>& conn) const
// {
//   assert (_nodes != NULL);
//   assert (sc < this->n_sub_elem());

//   // std::vector<unsigned int> conn(8);
//   conn.resize(8);
  
//   conn[0] = this->node(0)+1;
//   conn[1] = this->node(1)+1;
//   conn[2] = this->node(2)+1;
//   conn[3] = this->node(3)+1;
//   conn[4] = this->node(4)+1;
//   conn[5] = this->node(5)+1;
//   conn[6] = this->node(6)+1;
//   conn[7] = this->node(7)+1;
// }






// void Hex8::vtk_connectivity(const unsigned int sc,
// 			    std::vector<unsigned int> *conn) const
// {
//   assert (_nodes != NULL);
//   assert (sc < this->n_sub_elem());
  
//   if (conn == NULL)
//     conn = new std::vector<unsigned int>;

//   conn->resize(8);

//   (*conn)[0] = this->node(0);
//   (*conn)[1] = this->node(1);
//   (*conn)[2] = this->node(2);
//   (*conn)[3] = this->node(3);
//   (*conn)[4] = this->node(4);
//   (*conn)[5] = this->node(5);
//   (*conn)[6] = this->node(6);
//   (*conn)[7] = this->node(7);

//   return;
// }



#ifdef ENABLE_AMR

const float Hex8::_embedding_matrix[8][8][8] =
{
  // embedding matrix for child 0
  {
    //  0     1     2     3     4     5     6     7
    { 1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0}, // 0
    { 0.5,  0.5,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0}, // 1
    { .25,  .25,  .25,  .25,  0.0,  0.0,  0.0,  0.0}, // 2
    { 0.5,  0.0,  0.0,  0.5,  0.0,  0.0,  0.0,  0.0}, // 3
    { 0.5,  0.0,  0.0,  0.0,  0.5,  0.0,  0.0,  0.0}, // 4
    { .25,  .25,  0.0,  0.0,  .25,  .25,  0.0,  0.0}, // 5
    {.125, .125, .125, .125, .125, .125, .125, .125}, // 6
    { .25,  0.0,  0.0,  .25,  .25,  0.0,  0.0,  .25}  // 7
  },

  // embedding matrix for child 1
  {
    //  0     1     2     3     4     5     6     7
    { 0.5,  0.5,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0}, // 0
    { 0.0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0}, // 1
    { 0.0,  0.5,  0.5,  0.0,  0.0,  0.0,  0.0,  0.0}, // 2
    { .25,  .25,  .25,  .25,  0.0,  0.0,  0.0,  0.0}, // 3
    { .25,  .25,  0.0,  0.0,  .25,  .25,  0.0,  0.0}, // 4
    { 0.0,  0.5,  0.0,  0.0,  0.0,  0.5,  0.0,  0.0}, // 5
    { 0.0,  .25,  .25,  0.0,  0.0,  .25,  .25,  0.0}, // 6
    {.125, .125, .125, .125, .125, .125, .125, .125}  // 7
  },

  // embedding matrix for child 2
  {
    //  0      1    2     3     4     5     6     7
    { 0.5,  0.0,  0.0,  0.5,  0.0,  0.0,  0.0,  0.0}, // 0
    { .25,  .25,  .25,  .25,  0.0,  0.0,  0.0,  0.0}, // 1
    { 0.0,  0.0,  0.5,  0.5,  0.0,  0.0,  0.0,  0.0}, // 2
    { 0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  0.0}, // 3
    { .25,  0.0,  0.0,  .25,  .25,  0.0,  0.0,  .25}, // 4
    {.125, .125, .125, .125, .125, .125, .125, .125}, // 5
    { 0.0,  0.0,  .25,  .25,  0.0,  0.0,  .25,  .25}, // 6
    { 0.0,  0.0,  0.0,  0.5,  0.0,  0.0,  0.0,  0.5}  // 7
  },

  // embedding matrix for child 3
  {
    //  0      1    2     3     4     5     6     7
    { .25,  .25,  .25,  .25,  0.0,  0.0,  0.0,  0.0}, // 0
    { 0.0,  0.5,  0.5,  0.0,  0.0,  0.0,  0.0,  0.0}, // 1
    { 0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0}, // 2
    { 0.0,  0.0,  0.5,  0.5,  0.0,  0.0,  0.0,  0.0}, // 3
    {.125, .125, .125, .125, .125, .125, .125, .125}, // 4
    { 0.0,  .25,  .25,  0.0,  0.0,  .25,  .25,  0.0}, // 5
    { 0.0,  0.0,  0.5,  0.0,  0.0,  0.0,  0.5,  0.0}, // 6
    { 0.0,  0.0,  .25,  .25,  0.0,  0.0,  .25,  .25}  // 7
  },

  // embedding matrix for child 4
  {
    //  0      1    2     3     4     5     6     7
    { 0.5,  0.0,  0.0,  0.0,  0.5,  0.0,  0.0,  0.0}, // 0
    { .25,  .25,  0.0,  0.0,  .25,  .25,  0.0,  0.0}, // 1
    {.125, .125, .125, .125, .125, .125, .125, .125}, // 2
    { .25,  0.0,  0.0,  .25,  .25,  0.0,  0.0,  .25}, // 3
    { 0.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0}, // 4
    { 0.0,  0.0,  0.0,  0.0,  0.5,  0.5,  0.0,  0.0}, // 5
    { 0.0,  0.0,  0.0,  0.0,  .25,  .25,  .25,  .25}, // 6
    { 0.0,  0.0,  0.0,  0.0,  0.5,  0.0,  0.0,  0.5}  // 7
  },

  // embedding matrix for child 5
  {
    //  0      1    2     3     4     5     6     7
    { .25,  .25,  0.0,  0.0,  .25,  .25,  0.0,  0.0}, // 0
    { 0.0,  0.5,  0.0,  0.0,  0.0,  0.5,  0.0,  0.0}, // 1
    { 0.0,  .25,  .25,  0.0,  0.0,  .25,  .25,  0.0}, // 2
    {.125, .125, .125, .125, .125, .125, .125, .125}, // 3
    { 0.0,  0.0,  0.0,  0.0,  0.5,  0.5,  0.0,  0.0}, // 4
    { 0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0}, // 5
    { 0.0,  0.0,  0.0,  0.0,  0.0,  0.5,  0.5,  0.0}, // 6
    { 0.0,  0.0,  0.0,  0.0,  .25,  .25,  .25,  .25}  // 7
  },

  // embedding matrix for child 6
  {
    //  0      1    2     3     4     5     6     7
    { .25,  0.0,  0.0,  .25,  .25,  0.0,  0.0,  .25}, // 0
    {.125, .125, .125, .125, .125, .125, .125, .125}, // 1
    { 0.0,  0.0,  .25,  .25,  0.0,  0.0,  .25,  .25}, // 2
    { 0.0,  0.0,  0.0,  0.5,  0.0,  0.0,  0.0,  0.5}, // 3
    { 0.0,  0.0,  0.0,  0.0,  0.5,  0.0,  0.0,  0.5}, // 4
    { 0.0,  0.0,  0.0,  0.0,  .25,  .25,  .25,  .25}, // 5
    { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.5,  0.5}, // 6
    { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0}  // 7
  },

  // embedding matrix for child 7
  {
    //  0      1    2     3     4     5     6     7
    {.125, .125, .125, .125, .125, .125, .125, .125}, // 0
    { 0.0,  .25,  .25,  0.0,  0.0,  .25,  .25,  0.0}, // 1
    { 0.0,  0.0,  0.5,  0.0,  0.0,  0.0,  0.5,  0.0}, // 2
    { 0.0,  0.0,  .25,  .25,  0.0,  0.0,  .25,  .25}, // 3
    { 0.0,  0.0,  0.0,  0.0,  .25,  .25,  .25,  .25}, // 4
    { 0.0,  0.0,  0.0,  0.0,  0.0,  0.5,  0.5,  0.0}, // 5
    { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  0.0}, // 6
    { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.5,  0.5}  // 7
  }
};

#endif
