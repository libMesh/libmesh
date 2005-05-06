// $Id: cell_tet4.C,v 1.22 2005-05-06 17:06:58 roystgnr Exp $

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
#include "cell_tet4.h"
#include "face_tri3.h"



// ------------------------------------------------------------
// Tet4 class static member initializations
const unsigned int Tet4::side_nodes_map[4][3] =
{
  {0, 2, 1}, // Side 0
  {0, 1, 3}, // Side 1
  {1, 2, 3}, // Side 2
  {2, 0, 3}  // Side 3
};

const unsigned int Tet4::edge_nodes_map[6][2] =
{
  {0, 1}, // Side 0
  {1, 2}, // Side 1
  {0, 2}, // Side 2
  {0, 3}, // Side 3
  {1, 3}, // Side 4
  {2, 3}  // Side 5
};


// ------------------------------------------------------------
// Tet4 class member functions

bool Tet4::is_vertex(const unsigned int) const
{
  return true;
}

bool Tet4::is_edge(const unsigned int) const
{
  return false;
}

bool Tet4::is_face(const unsigned int) const
{
  return false;
}

bool Tet4::is_node_on_edge(const unsigned int n,
			   const unsigned int e) const
{
  assert(e < n_edges());
  for (unsigned int i = 0; i != 2; ++i)
    if (edge_nodes_map[e][i] == n)
      return true;
  return false;
}

bool Tet4::is_node_on_side(const unsigned int n,
			   const unsigned int s) const
{
  assert(s < n_sides());
  for (unsigned int i = 0; i != 3; ++i)
    if (side_nodes_map[s][i] == n)
      return true;
  return false;
}

AutoPtr<Elem> Tet4::build_side (const unsigned int i) const
{
  assert (i < this->n_sides());

  AutoPtr<Elem> ap(new Side<Tri3,Tet4>(this,i));
  return ap;
  
//   AutoPtr<Elem> face(new Tri3);

//   switch (i)
//     {
//     case 0:
//       {
// 	face->set_node(0) = this->get_node(0);
// 	face->set_node(1) = this->get_node(2);
// 	face->set_node(2) = this->get_node(1);

// 	return face;
//       }
//     case 1:
//       {
// 	face->set_node(0) = this->get_node(0);
// 	face->set_node(1) = this->get_node(1);
// 	face->set_node(2) = this->get_node(3);

// 	return face;
//       }
//     case 2:
//       {
// 	face->set_node(0) = this->get_node(1);
// 	face->set_node(1) = this->get_node(2);
// 	face->set_node(2) = this->get_node(3);

// 	return face;
//       }
//     case 3:
//       {
// 	face->set_node(0) = this->get_node(2);
// 	face->set_node(1) = this->get_node(0);
// 	face->set_node(2) = this->get_node(3);
	
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


void Tet4::connectivity(const unsigned int sc,
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
	conn[5] = this->node(3)+1;
	conn[6] = this->node(3)+1;
	conn[7] = this->node(3)+1;
	return;
      }

    case VTK:
      {
	conn.resize(4);
	conn[0] = this->node(0);
	conn[1] = this->node(1);
	conn[2] = this->node(2);
	conn[3] = this->node(3);
	return;
      }

    default:
      error();
    }

  error();
}



#ifdef ENABLE_AMR

const float Tet4::_embedding_matrix[8][4][4] =
{
  // embedding matrix for child 0
  {
    // 0    1    2    3  
    {1.0, 0.0, 0.0, 0.0}, // 0
    {0.5, 0.5, 0.0, 0.0}, // 1
    {0.5, 0.0, 0.5, 0.0}, // 2
    {0.5, 0.0, 0.0, 0.5}  // 3
  },
  
  // embedding matrix for child 1
  {
    // 0    1    2    3  
    {0.5, 0.5, 0.0, 0.0}, // 0
    {0.0, 1.0, 0.0, 0.0}, // 1
    {0.0, 0.5, 0.5, 0.0}, // 2
    {0.0, 0.5, 0.0, 0.5}  // 3
  },
  
  // embedding matrix for child 2
  {
    // 0    1    2    3  
    {0.5, 0.0, 0.5, 0.0}, // 0
    {0.0, 0.5, 0.5, 0.0}, // 1
    {0.0, 0.0, 1.0, 0.0}, // 2
    {0.0, 0.0, 0.5, 0.5}  // 3
  },
  
  // embedding matrix for child 3
  {
    // 0    1    2    3  
    {0.5, 0.0, 0.0, 0.5}, // 0
    {0.0, 0.5, 0.0, 0.5}, // 1
    {0.0, 0.0, 0.5, 0.5}, // 2
    {0.0, 0.0, 0.0, 1.0}  // 3
  },
  
  // embedding matrix for child 4
  {
    // 0    1    2    3  
    {0.5, 0.5, 0.0, 0.0}, // 0
    {0.0, 0.5, 0.0, 0.5}, // 1
    {0.5, 0.0, 0.5, 0.0}, // 2
    {0.5, 0.0, 0.0, 0.5}  // 3
  },
  
  // embedding matrix for child 5
  {
    // 0    1    2    3  
    {0.5, 0.5, 0.0, 0.0}, // 0
    {0.0, 0.5, 0.5, 0.0}, // 1
    {0.5, 0.0, 0.5, 0.0}, // 2
    {0.0, 0.5, 0.0, 0.5}  // 3
  },
  
  // embedding matrix for child 6
  {
    // 0    1    2    3  
    {0.5, 0.0, 0.5, 0.0}, // 0
    {0.0, 0.5, 0.5, 0.0}, // 1
    {0.0, 0.0, 0.5, 0.5}, // 2
    {0.0, 0.5, 0.0, 0.5}  // 3
  },
  
  // embedding matrix for child 7
  {
    // 0    1    2    3  
    {0.5, 0.0, 0.5, 0.0}, // 0
    {0.0, 0.5, 0.0, 0.5}, // 1
    {0.0, 0.0, 0.5, 0.5}, // 2
    {0.5, 0.0, 0.0, 0.5}  // 3
  }  
};

#endif // #ifdef ENABLE_AMR
