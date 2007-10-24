// $Id: face_inf_quad4.C,v 1.14 2003-03-11 00:47:45 ddreyer Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002  Benjamin S. Kirk, John W. Peterson
  
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
#include "mesh_config.h"
#ifdef ENABLE_INFINITE_ELEMENTS


// Local includes cont'd
#include "face_inf_quad4.h"




// ------------------------------------------------------------
// InfQuad4 class member functions
#ifdef ENABLE_AMR

const float InfQuad4::_embedding_matrix[2][4][4] =
{
  // embedding matrix for child 0
  {
    // 0    1    2    3rd parent node
    {1.0, 0.0, 0.0, 0.0}, // 0th child node
    {0.5, 0.5, 0.0, 0.0}, // 1
    {0.0, 0.0, 0.5, 0.5}, // 2
    {0.0, 0.0, 0.0, 1.0}  // 3
  },

  // embedding matrix for child 1
  {
    // 0    1    2    3
    {0.5, 0.5, 0.0, 0.0}, // 0
    {0.0, 1.0, 0.0, 0.0}, // 1
    {0.0, 0.0, 1.0, 0.0}, // 2
    {0.0, 0.0, 0.5, 0.5}  // 3
  }
};

#endif



const std::vector<unsigned int> InfQuad4::tecplot_connectivity(const unsigned int sf) const
{
  assert (sf < this->n_sub_elem());
  
  std::vector<unsigned int> conn(4);

  conn[0] = this->node(0)+1;
  conn[1] = this->node(1)+1;
  conn[2] = this->node(2)+1;
  conn[3] = this->node(3)+1;

  return conn;
}



void InfQuad4::vtk_connectivity(const unsigned int,
				std::vector<unsigned int> *) const
{
  error();  // Not yet implemented
}


#endif // ifdef ENABLE_INFINITE_ELEMENTS
