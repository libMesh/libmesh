// $Id: edge_edge2.C,v 1.13 2005-02-25 19:16:37 roystgnr Exp $

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



// Local includes
#include "edge_edge2.h"


bool Edge2::is_vertex(const unsigned int) const
{
  return true;
}

bool Edge2::is_edge(const unsigned int) const
{
  return false;
}

bool Edge2::is_face(const unsigned int) const
{
  return false;
}

bool Edge2::is_node_on_side(const unsigned int n,
			    const unsigned int s) const
{
  assert(s < 2);
  return (s == n);
}

void Edge2::connectivity(const unsigned int sc,
			 const IOPackage iop,
			 std::vector<unsigned int>& conn) const
{
  assert (sc == 0);
  assert (_nodes != NULL);
  assert (sc < this->n_sub_elem());
  assert (iop != INVALID_IO_PACKAGE);
  
  // Create storage
  conn.resize(2);

  switch (iop)
    {
    case TECPLOT:
      {
	conn[0] = this->node(0)+1;
	conn[1] = this->node(1)+1;
	return;
      }

    case VTK:
      {
	conn[0] = this->node(0);
	conn[1] = this->node(1);
	return;
      }

    default:
      {
	error();
      }
    }

  error();
}




// void Edge2::tecplot_connectivity(const unsigned int se,
// 				 std::vector<unsigned int>& conn) const
// {
//   assert (se == 0);

//   //std::vector<unsigned int> conn(2);
//   conn.resize(2);

//   conn[0] = this->node(0)+1;
//   conn[1] = this->node(1)+1;
// }



// void Edge2::vtk_connectivity(const unsigned int se,
// 			     std::vector<unsigned int> *conn) const
// {
//   assert (_nodes != NULL);
//   assert (se < this->n_sub_elem());
  
//   if (conn == NULL)
//     conn = new std::vector<unsigned int>;

//   conn->resize(2);

//   (*conn)[0] = this->node(0);
//   (*conn)[1] = this->node(1);

//   return;
// }


