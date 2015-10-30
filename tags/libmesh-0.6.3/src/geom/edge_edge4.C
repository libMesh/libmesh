// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2007  Benjamin S. Kirk, John W. Peterson
  
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
#include "edge_edge4.h"

#ifdef ENABLE_AMR

const float Edge4::_embedding_matrix[2][4][4] =
{
  // embedding matrix for child 0

  {
    // 0       1       2        3        // Shape function index
    {1.0,      0.0,    0.0,     0.0},    // left,  xi = -1
    {-0.0625, -0.0625, 0.5625,  0.5625}, // right, xi = 0
    {0.3125,   0.0625, 0.9375, -0.3125}, // middle left,  xi = -2/3
    {0.0,      0.0,    1.0,     0.0}     // middle right, xi = -1/3
  },

  // embedding matrix for child 1
  {
    // 0       1        2        3        // Shape function index
    {-0.0625, -0.0625,  0.5625,  0.5625}, // left,  xi = 0
    {0.0,      1.0,     0.0,     0.0},    // right, xi = 1
    {0.0,      0.0,     0.0,     1.0},    // middle left,  xi = 1/3
    {0.0625,   0.3125, -0.3125,  0.9375}  // middle right, xi = 2/3
  }
};

#endif

bool Edge4::is_vertex(const unsigned int i) const
{
  return (i==0) || (i==1);
}

bool Edge4::is_edge(const unsigned int i) const
{
  return (i==2) || (i==3);
}

bool Edge4::is_face(const unsigned int ) const
{
  return false;
}

bool Edge4::is_node_on_side(const unsigned int n,
                            const unsigned int s) const
{
  libmesh_assert(s < 2);
  libmesh_assert(n < 4);
  return (s == n);
}

bool Edge4::is_node_on_edge(const unsigned int,
                            const unsigned int e) const
{
  libmesh_assert(e == 0);
  return true;
}



bool Edge4::has_affine_map() const
{
  if (!this->point(2).relative_fuzzy_equals
      ((this->point(0)*2. + this->point(1))/3.))
    return false;
  if (!this->point(3).relative_fuzzy_equals
      ((this->point(0) + this->point(1)*2.)/3.))
    return false;
  return true;
}



void Edge4::connectivity(const unsigned int sc,
			 const IOPackage iop,
			 std::vector<unsigned int>& conn) const
{
  libmesh_assert (sc <= 2);
  libmesh_assert (sc < this->n_sub_elem());
  libmesh_assert (iop != INVALID_IO_PACKAGE);

  // Create storage
  conn.resize(2);

  switch(iop)
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
          conn[1] = this->node(3)+1;
          return;

        case 2:
          conn[0] = this->node(3)+1;
          conn[1] = this->node(1)+1;
          return;

        default:
          libmesh_error();
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
          conn[1] = this->node(3);
          return;

        case 2:
          conn[0] = this->node(3);
          conn[1] = this->node(1);
          return;

        default:
          libmesh_error();
      }
   }
  
  default:
    libmesh_error();

  }

}
