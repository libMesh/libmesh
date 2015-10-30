// $Id: face_quad4.C,v 1.15 2004-07-14 19:23:18 jwpeterson Exp $

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
#include "edge_edge2.h"
#include "face_quad4.h"




// ------------------------------------------------------------
// Quad class static member initialization
#ifdef ENABLE_AMR

const float Quad4::_embedding_matrix[4][4][4] =
{
  // embedding matrix for child 0
  {
    // 0    1    2    3
    {1.0, 0.0, 0.0, 0.0}, // 0
    {0.5, 0.5, 0.0, 0.0}, // 1
    {.25, .25, .25, .25}, // 2
    {0.5, 0.0, 0.0, 0.5}  // 3
  },

  // embedding matrix for child 1
  {
    // 0    1    2    3
    {0.5, 0.5, 0.0, 0.0}, // 0
    {0.0, 1.0, 0.0, 0.0}, // 1
    {0.0, 0.5, 0.5, 0.0}, // 2
    {.25, .25, .25, .25}  // 3
  },

  // embedding matrix for child 2
  {
    // 0    1    2    3
    {0.5, 0.0, 0.0, 0.5}, // 0
    {.25, .25, .25, .25}, // 1
    {0.0, 0.0, 0.5, 0.5}, // 2
    {0.0, 0.0, 0.0, 1.0}  // 3
  },

  // embedding matrix for child 3
  {
    // 0    1    2    3
    {.25, .25, .25, .25}, // 0
    {0.0, 0.5, 0.5, 0.0}, // 1
    {0.0, 0.0, 1.0, 0.0}, // 2
    {0.0, 0.0, 0.5, 0.5}  // 3
  }
};

#endif





// ------------------------------------------------------------
// Quad4 class member functions
AutoPtr<Elem> Quad4::build_side (const unsigned int i) const
{
  assert (i < this->n_sides());

  switch (i)
    {
    case 0:
      {
	Edge2* edge = new Edge2;

	edge->set_node(0) = this->get_node(0);
	edge->set_node(1) = this->get_node(1);
	
	AutoPtr<Elem> ap(edge);  return ap;
      }
    case 1:
      {
	Edge2* edge = new Edge2;

	edge->set_node(0) = this->get_node(1);
	edge->set_node(1) = this->get_node(2);
	
	AutoPtr<Elem> ap(edge);  return ap;
      }
    case 2:
      {
	Edge2* edge = new Edge2;

	edge->set_node(0) = this->get_node(2);
	edge->set_node(1) = this->get_node(3);
	
	AutoPtr<Elem> ap(edge);  return ap;
      }
    case 3:
      {
	Edge2* edge = new Edge2;

	edge->set_node(0) = this->get_node(3);
	edge->set_node(1) = this->get_node(0);
	
	AutoPtr<Elem> ap(edge);  return ap;
      }
    default:
      {
	error();
      }
    }


  // We will never get here...  Look at the code above.
  AutoPtr<Elem> ap(NULL);  return ap;
}





void Quad4::connectivity(const unsigned int sf,
			 const IOPackage iop,
			 std::vector<unsigned int>& conn) const
{
  assert (sf < this->n_sub_elem());
  assert (_nodes != NULL);
  assert (iop != INVALID_IO_PACKAGE);

  // Create storage.
  conn.resize(4);

  switch (iop)
    {
    case TECPLOT:
      {
	conn[0] = this->node(0)+1;
	conn[1] = this->node(1)+1;
	conn[2] = this->node(2)+1;
	conn[3] = this->node(3)+1;
	return;
      }

    case VTK:
      {
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

// void Quad4::tecplot_connectivity(const unsigned int sf,
// 				 std::vector<unsigned int>& conn) const
// {
//   assert (sf < this->n_sub_elem());

//   //std::vector<unsigned int> conn(4);
//   conn.resize(4);

  
//   conn[0] = this->node(0)+1;
//   conn[1] = this->node(1)+1;
//   conn[2] = this->node(2)+1;
//   conn[3] = this->node(3)+1;

// }



// void Quad4::vtk_connectivity(const unsigned int sf,
// 			     std::vector<unsigned int> *conn) const
// {
//   assert (_nodes != NULL);
//   assert (sf < this->n_sub_elem());
  
//   if (conn == NULL)
//     conn = new std::vector<unsigned int>;

//   conn->resize(4);

//   (*conn)[0] = this->node(0);
//   (*conn)[1] = this->node(1);
//   (*conn)[2] = this->node(2);
//   (*conn)[3] = this->node(3);

//   return;
// }
