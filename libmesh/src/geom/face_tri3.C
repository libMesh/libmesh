// $Id: face_tri3.C,v 1.17 2005-02-19 19:07:01 roystgnr Exp $

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

// C++ includes

// Local includes
#include "side.h"
#include "edge_edge2.h"
#include "face_tri3.h"




// ------------------------------------------------------------
// Tri3 class static member initializations
const unsigned int Tri3::side_nodes_map[3][2] =
{
  {0, 1}, // Side 0
  {1, 2}, // Side 1
  {2, 0}  // Side 2
};


#ifdef ENABLE_AMR

const float Tri3::_embedding_matrix[4][3][3] =
{
  // embedding matrix for child 0
  {
    // 0    1    2  
    {1.0, 0.0, 0.0}, // 0
    {0.5, 0.5, 0.0}, // 1
    {0.5, 0.0, 0.5}  // 2
  },

  // embedding matrix for child 1
  {
    // 0    1    2  
    {0.5, 0.5, 0.0}, // 0
    {0.0, 1.0, 0.0}, // 1
    {0.0, 0.5, 0.5}  // 2
  },

  // embedding matrix for child 2
  {
    // 0    1    2  
    {0.5, 0.0, 0.5}, // 0
    {0.0, 0.5, 0.5}, // 1
    {0.0, 0.0, 1.0}  // 2
  },

  // embedding matrix for child 3
  {
    // 0    1    2  
    {0.5, 0.5, 0.0}, // 0
    {0.0, 0.5, 0.5}, // 1
    {0.5, 0.0, 0.5}  // 2
  }
};

#endif



// ------------------------------------------------------------
// Tri3 class member functions

bool Tri3::is_vertex(const unsigned int) const
{
  return true;
}

bool Tri3::is_edge(const unsigned int) const
{
  return false;
}

bool Tri3::is_face(const unsigned int) const
{
  return false;
}

AutoPtr<Elem> Tri3::build_side (const unsigned int i) const
{
  assert (i < this->n_sides());

  AutoPtr<Elem> ap(new Side<Edge2,Tri3>(this,i));
  return ap;
  
//   Edge2* edge = new Edge2;

//   switch (i)
//     {
//     case 0:
//       {
// 	edge->set_node(0) = this->get_node(0);
// 	edge->set_node(1) = this->get_node(1);
	
// 	AutoPtr<Elem> ap(edge);  return ap;
//       }
//     case 1:
//       {
// 	edge->set_node(0) = this->get_node(1);
// 	edge->set_node(1) = this->get_node(2);
	
// 	AutoPtr<Elem> ap(edge);  return ap;
//       }
//     case 2:
//       {
// 	edge->set_node(0) = this->get_node(2);
// 	edge->set_node(1) = this->get_node(0);
	
// 	AutoPtr<Elem> ap(edge);  return ap;
//       }
//     default:
//       {
// 	error();
//       }
//     }

  
//   // We will never get here...  Look at the code above.
//   error();
//   AutoPtr<Elem> ap(NULL);  return ap;
}


void Tri3::connectivity(const unsigned int sf,
			const IOPackage iop,
			std::vector<unsigned int>& conn) const
{
  assert (sf <this->n_sub_elem());
  assert (_nodes != NULL);
  assert (iop != INVALID_IO_PACKAGE);
  
  switch (iop)
    {
    case TECPLOT:
      {
	conn.resize(4);
	conn[0] = this->node(0)+1;
	conn[1] = this->node(1)+1;
	conn[2] = this->node(2)+1;
	conn[3] = this->node(2)+1;
	return;
      }

    case VTK:
      {
	conn.resize(3);
	conn[0] = this->node(0);
	conn[1] = this->node(1);
	conn[2] = this->node(2);
	return;
      }

    default:
      error();
    }

  error();
}



// void Tri3::tecplot_connectivity(const unsigned int sf,
// 				std::vector<unsigned int>& conn) const
// {
//   assert (sf <this->n_sub_elem());
  
//   // std::vector<unsigned int> conn(4);
//   conn.resize(4);
  
//   conn[0] = this->node(0)+1;
//   conn[1] = this->node(1)+1;
//   conn[2] = this->node(2)+1;
//   conn[3] = this->node(2)+1;
// }



// void Tri3::vtk_connectivity(const unsigned int sf,
// 			    std::vector<unsigned int> *conn) const
// {
//   assert (_nodes != NULL);
//   assert (sf < this->n_sub_elem());
  
//   if (conn == NULL)
//     conn = new std::vector<unsigned int>;

//   conn->resize(3);

//   (*conn)[0] = this->node(0);
//   (*conn)[1] = this->node(1);
//   (*conn)[2] = this->node(2);

//   return;
// }
