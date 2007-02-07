// $Id: edge_edge3.C,v 1.17 2005-06-16 23:03:52 roystgnr Exp $

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
#include "edge_edge3.h"

#ifdef ENABLE_AMR

const float Edge3::_embedding_matrix[2][3][3] =
{
  // embedding matrix for child 0
  {
    // 0    1    2  
    {1.0, 0.0, 0.0}, // left
    {0.0, 0.0, 1.0}, // right
    {3./8.,-1./8.,0.75} // middle
  },

  // embedding matrix for child 1
  {
    // 0    1    2  
    {0.0, 0.0, 1.0}, // left
    {0.0, 1.0, 0.0},  // right
    {-1./8.,3./8.,0.75} // middle
  }
};

#endif

bool Edge3::is_vertex(const unsigned int i) const
{
  if (i < 2)
    return true;
  return false;
}

bool Edge3::is_edge(const unsigned int i) const
{
  if (i < 2)
   return false;
  return true;
}

bool Edge3::is_face(const unsigned int ) const
{
  return false;
}

bool Edge3::is_node_on_side(const unsigned int n,
			    const unsigned int s) const
{
  assert(s < 2);
  assert(n < 3);
  return (s == n);
}

bool Edge3::is_node_on_edge(const unsigned int,
			    const unsigned int e) const
{
  assert(e == 0);
  return true;
}



bool Edge3::has_affine_map() const
{
  return ((this->point(0) + this->point(1))/2 == this->point(2));
}



void Edge3::connectivity(const unsigned int sc,
			 const IOPackage iop,
			 std::vector<unsigned int>& conn) const
{
  assert (sc <= 1);
  assert (sc < this->n_sub_elem());
  assert (iop != INVALID_IO_PACKAGE);

  // Create storage
  conn.resize(2);

  switch (iop)
    {
    case TECPLOT:
      {
	switch (sc)
	  {
	  case 0: 
	    conn[0] = this->node(0)+1;
	    conn[1] = this->node(2)+1;
	    return;
      
	  case 1: 
	    conn[0] = this->node(2)+1;
	    conn[1] = this->node(1)+1;
	    return;

	  default:
	    error();
	  }
      }

      
    case VTK:
      {
	switch (sc)
	  {
	  case 0: 
	    conn[0] = this->node(0);
	    conn[1] = this->node(2);
      
	    return;
      
	  case 1: 
	    conn[0] = this->node(2);
	    conn[1] = this->node(1);
      
	    return;

	  default:
	    error();
	  }
      }

    default:
      {
	error();
      }
    }
}
