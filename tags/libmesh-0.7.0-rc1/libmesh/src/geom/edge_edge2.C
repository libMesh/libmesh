// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
  
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

namespace libMesh
{


#ifdef LIBMESH_ENABLE_AMR

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
  libmesh_assert(s < 2);
  return (s == n);
}

bool Edge2::is_node_on_edge(const unsigned int,
			    const unsigned int e) const
{
  libmesh_assert(e == 0);
  return true;
}

void Edge2::connectivity(const unsigned int sc,
			 const IOPackage iop,
			 std::vector<unsigned int>& conn) const
{
  libmesh_assert (sc == 0);
  libmesh_assert (sc < this->n_sub_elem());
  libmesh_assert (iop != INVALID_IO_PACKAGE);
  
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
	libmesh_error();
      }
    }

  libmesh_error();
}


Real Edge2::volume () const
{
  // OK, so this is probably overkill, since it is equivalent to
  // Elem::hmax() for the Edge2, but here it is nonetheless...
  return (this->point(1) - this->point(0)).size();
}

} // namespace libMesh
