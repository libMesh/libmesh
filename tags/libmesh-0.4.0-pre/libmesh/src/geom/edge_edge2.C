// $Id: edge_edge2.C,v 1.7 2003-02-26 04:43:14 jwpeterson Exp $

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
#include "edge_edge2.h"



const std::vector<unsigned int> Edge2::tecplot_connectivity(const unsigned int se) const
{
  assert (se == 0);

  std::vector<unsigned int> conn(2);

  conn[0] = this->node(0)+1;
  conn[1] = this->node(1)+1;

  return conn;
}



void Edge2::vtk_connectivity(const unsigned int se,
			     std::vector<unsigned int> *conn) const
{
  assert (_nodes != NULL);
  assert (se < this->n_sub_elem());
  
  if (conn == NULL)
    conn = new std::vector<unsigned int>;

  conn->resize(2);

  (*conn)[0] = this->node(0);
  (*conn)[1] = this->node(1);

  return;
}


