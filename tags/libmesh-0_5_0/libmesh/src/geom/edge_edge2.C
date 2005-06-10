// $Id: edge_edge2.C,v 1.16 2005-06-06 16:24:13 knezed01 Exp $

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


#ifdef ENABLE_AMR

const float Edge2::_embedding_matrix[2][2][2] =
{
  // embedding matrix for child 0
  {
    // 0    1    2  
    {1.0, 0.0}, // 0
    {0.5, 0.5}  // 1
  },

  // embedding matrix for child 1
  {
    // 0    1    2  
    {0.5, 0.5}, // 0
    {0.0, 1.0}  // 1
  }
};

#endif

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

bool Edge2::is_node_on_edge(const unsigned int,
			    const unsigned int e) const
{
  assert(e == 0);
  return true;
}

void Edge2::connectivity(const unsigned int sc,
			 const IOPackage iop,
			 std::vector<unsigned int>& conn) const
{
  assert (sc == 0);
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
