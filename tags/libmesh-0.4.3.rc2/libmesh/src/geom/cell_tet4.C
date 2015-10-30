// $Id: cell_tet4.C,v 1.16 2004-07-14 19:23:18 jwpeterson Exp $

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
#include "cell_tet4.h"
#include "face_tri3.h"



// ------------------------------------------------------------
// Tet4 class member functions
AutoPtr<Elem> Tet4::build_side (const unsigned int i) const
{
  assert (i < this->n_sides());

  
  
  AutoPtr<Elem> face(new Tri3);

  switch (i)
    {
    case 0:
      {
	face->set_node(0) = this->get_node(0);
	face->set_node(1) = this->get_node(2);
	face->set_node(2) = this->get_node(1);

	return face;
      }
    case 1:
      {
	face->set_node(0) = this->get_node(0);
	face->set_node(1) = this->get_node(1);
	face->set_node(2) = this->get_node(3);

	return face;
      }
    case 2:
      {
	face->set_node(0) = this->get_node(1);
	face->set_node(1) = this->get_node(2);
	face->set_node(2) = this->get_node(3);

	return face;
      }
    case 3:
      {
	face->set_node(0) = this->get_node(2);
	face->set_node(1) = this->get_node(0);
	face->set_node(2) = this->get_node(3);
	
	return face;
      }
    default:
      {
	error();
      }
    }

  // We'll never get here.
  error();  
  return face;
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





// void Tet4::tecplot_connectivity(const unsigned int sc,
// 				std::vector<unsigned int>& conn) const
// {
//   assert (_nodes != NULL);
//   assert (sc < this->n_sub_elem());

//   // std::vector<unsigned int> conn(8);
//   conn.resize(8);
  
//   conn[0] = this->node(0)+1;
//   conn[1] = this->node(1)+1;
//   conn[2] = this->node(2)+1;
//   conn[3] = this->node(2)+1;
//   conn[4] = this->node(3)+1;
//   conn[5] = this->node(3)+1;
//   conn[6] = this->node(3)+1;
//   conn[7] = this->node(3)+1;
// }



// void Tet4::vtk_connectivity(const unsigned int sc,
// 			    std::vector<unsigned int> *conn) const
// {
//   assert (_nodes != NULL);
//   assert (sc < this->n_sub_elem());
  
//   if (conn == NULL)
//     conn = new std::vector<unsigned int>;

//   conn->resize(4);

//   (*conn)[0] = this->node(0);
//   (*conn)[1] = this->node(1);
//   (*conn)[2] = this->node(2);
//   (*conn)[3] = this->node(3);

//   return;
// }



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
