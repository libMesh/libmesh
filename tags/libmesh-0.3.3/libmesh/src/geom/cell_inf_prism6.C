// $Id: cell_inf_prism6.C,v 1.14 2003-03-11 00:47:44 ddreyer Exp $

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

// C++ includes

// Local includes cont'd
#include "cell_inf_prism6.h"



// ------------------------------------------------------------
// InfPrism6 class member functions
const std::vector<unsigned int> InfPrism6::tecplot_connectivity(const unsigned int sc) const
{
  assert (_nodes != NULL);
  assert (sc < this->n_sub_elem());

  std::vector<unsigned int> conn(8);

  // guess this is a collapsed hex8
  conn[0] = this->node(0)+1;
  conn[1] = this->node(1)+1;
  conn[2] = this->node(2)+1;
  conn[3] = this->node(2)+1;
  conn[4] = this->node(3)+1;
  conn[5] = this->node(4)+1;
  conn[6] = this->node(5)+1;
  conn[7] = this->node(5)+1;

  return conn;
}





#ifdef ENABLE_AMR

const float InfPrism6::_embedding_matrix[4][6][6] =
{
  // embedding matrix for child 0
  {
    //          0           1           2           3           4           5 th parent Node
    {         1.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 0th child N.
    {         0.5,        0.5,        0.0,        0.0,        0.0,        0.0}, // 1
    {         0.5,        0.0,        0.5,        0.0,        0.0,        0.0}, // 2
    {         0.0,        0.0,        0.0,        1.0,        0.0,        0.0}, // 3
    {         0.0,        0.0,        0.0,        0.5,        0.5,        0.0}, // 4
    {         0.0,        0.0,        0.0,        0.5,        0.0,        0.5}  // 5
  },

  // embedding matrix for child 1
  {
    //          0           1           2           3           4           5 th parent Node
    {         0.5,        0.5,        0.0,        0.0,        0.0,        0.0}, // 0th child N.
    {         0.0,        1.0,        0.0,        0.0,        0.0,        0.0}, // 1
    {         0.0,        0.5,        0.5,        0.0,        0.0,        0.0}, // 2
    {         0.0,        0.0,        0.0,        0.5,        0.5,        0.0}, // 3
    {         0.0,        0.0,        0.0,        0.0,        1.0,        0.0}, // 4
    {         0.0,        0.0,        0.0,        0.0,        0.5,        0.5}  // 5
  },

  // embedding matrix for child 2
  {
    //          0           1           2           3           4           5 th parent Node
    {         0.5,        0.0,        0.5,        0.0,        0.0,        0.0}, // 0th child N.
    {         0.0,        0.5,        0.5,        0.0,        0.0,        0.0}, // 1
    {         0.0,        0.0,        1.0,        0.0,        0.0,        0.0}, // 2
    {         0.0,        0.0,        0.0,        0.5,        0.0,        0.5}, // 3
    {         0.0,        0.0,        0.0,        0.0,        0.5,        0.5}, // 4
    {         0.0,        0.0,        0.0,        0.0,        0.0,        1.0}  // 5
  },

  // embedding matrix for child 3
  {
    //          0           1           2           3           4           5 th parent Node
    {         0.5,        0.5,        0.0,        0.0,        0.0,        0.0}, // 0th child N.
    {         0.0,        0.5,        0.5,        0.0,        0.0,        0.0}, // 1
    {         0.5,        0.0,        0.5,        0.0,        0.0,        0.0}, // 2
    {         0.0,        0.0,        0.0,        0.5,        0.5,        0.0}, // 3
    {         0.0,        0.0,        0.0,        0.0,        0.5,        0.5}, // 4
    {         0.0,        0.0,        0.0,        0.5,        0.0,        0.5}  // 5
  }
};



#endif

#endif  // ifdef ENABLE_INFINITE_ELEMENTS
