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
#include "libmesh_config.h"
#ifdef ENABLE_INFINITE_ELEMENTS


// Local includes cont'd
#include "face_inf_quad6.h"
#include "edge_edge3.h"
#include "side.h"
#include "edge_inf_edge2.h"




// ------------------------------------------------------------
// InfQuad6 class static member initializations
const unsigned int InfQuad6::side_nodes_map[3][3] =
{
  {0, 1, 4}, // Side 0
  {1, 3},    // Side 1
  {0, 2}     // Side 2
};


// ------------------------------------------------------------
// InfQuad6 class member functions

bool InfQuad6::is_vertex(const unsigned int i) const
{
  if (i < 2)
    return true;
  return false;
}

bool InfQuad6::is_edge(const unsigned int i) const
{
  if (i < 2)
    return false;
  return true;
}

bool InfQuad6::is_face(const unsigned int) const
{
  return false;
}

bool InfQuad6::is_node_on_side(const unsigned int n,
			       const unsigned int s) const
{
  assert(s < n_sides());
  for (unsigned int i = 0; i != 3; ++i)
    if (side_nodes_map[s][i] == n)
      return true;
  return false;
}

#ifdef ENABLE_AMR

const float InfQuad6::_embedding_matrix[2][6][6] =
{
  // embedding matrix for child 0
  {
    //     0       1       2       3       4       5th parent node
    {    1.0,    0.0,    0.0,    0.0,    0.0,    0.0 }, // 0th child node
    {    0.0,    0.0,    0.0,    0.0,    1.0,    0.0 }, // 1
    {    0.0,    0.0,    1.0,    0.0,    0.0,    0.0 }, // 2
    {    0.0,    0.0,    0.0,    0.0,    0.0,    1.0 }, // 3
    {  0.375, -0.125,    0.0,    0.0,   0.75,    0.0 }, // 4
    {    0.0,    0.0,  0.375, -0.125,    0.0,   0.75 }  // 5
  },

  // embedding matrix for child 1
  {
    //     0       1       2       3       4       5th parent node
    {    0.0,    0.0,    0.0,    0.0,    1.0,    0.0 }, // 0th child node
    {    0.0,    1.0,    0.0,    0.0,    0.0,    0.0 }, // 1
    {    0.0,    0.0,    0.0,    0.0,    0.0,    1.0 }, // 2
    {    0.0,    0.0,    0.0,    1.0,    0.0,    0.0 }, // 3
    { -0.125,  0.375,    0.0,    0.0,   0.75,    0.0 }, // 4
    {    0.0,    0.0, -0.125,  0.375,    0.0,   0.75 }  // 5
  }
};

#endif




AutoPtr<Elem> InfQuad6::build_side (const unsigned int i,
				    bool proxy) const
{
  // assert (i < this->n_sides());

  if (proxy)
    {
      switch (i)
	{
	case 0:
	  {
	    AutoPtr<Elem> ap(new Side<Edge3,InfQuad6>(this,i));
	    return ap;
	  }
	case 1:
	case 2:
	  {
	    AutoPtr<Elem> ap(new Side<InfEdge2,InfQuad6>(this,i));
	    return ap;
	  }
	default:
	  error();
	}
    }

  else
    {
      switch (i)
	{
	case 0:
	  {
	    Edge3* edge = new Edge3;

	    edge->set_node(0) = this->get_node(0);
	    edge->set_node(1) = this->get_node(1);
	    edge->set_node(2) = this->get_node(4);
	
	    AutoPtr<Elem> ap(edge);  return ap;
	  }

	case 1:
	  {
	    // adjacent to another infinite element	
	    InfEdge2* edge = new InfEdge2;

	    edge->set_node(0) = this->get_node(1);
	    edge->set_node(1) = this->get_node(3);

	    AutoPtr<Elem> ap(edge);  return ap;
	  }

	case 2:
	  {
	    // adjacent to another infinite element	
	    InfEdge2* edge = new InfEdge2;

	    edge->set_node(0) = this->get_node(0); // be aware of swapped nodes,
	    edge->set_node(1) = this->get_node(2); // compared to conventional side numbering

	    AutoPtr<Elem> ap(edge);  return ap;
	  }
	default:
	  {
	    error();
	    AutoPtr<Elem> ap(NULL);  return ap;
	  }
	}
    }

  // We will never get here...  
  error();
  AutoPtr<Elem> ap(NULL);  return ap;
}




void InfQuad6::connectivity(const unsigned int sf,
			    const IOPackage iop,
			    std::vector<unsigned int>& conn) const
{
  assert (sf < this->n_sub_elem());
  assert (iop != INVALID_IO_PACKAGE);

  conn.resize(4);

  switch (iop)
    {
    case TECPLOT:
      {
	switch(sf)
	  {
	  case 0:
	    // linear sub-quad 0
	    conn[0] = this->node(0)+1;
	    conn[1] = this->node(4)+1;
	    conn[2] = this->node(5)+1;
	    conn[3] = this->node(2)+1;

	    return;

	  case 1:
	    // linear sub-quad 1
	    conn[0] = this->node(4)+1;
	    conn[1] = this->node(1)+1;
	    conn[2] = this->node(3)+1;
	    conn[3] = this->node(5)+1;

	    return;
	    
	  default:
	    error();
	  }
      }
      
    default:
      error();
    }
  
  error();
}




unsigned short int InfQuad6::second_order_adjacent_vertex (const unsigned int n,
							   const unsigned int v) const
{ 
  assert (n >= this->n_vertices());
  assert (n <  this->n_nodes());
  assert (v < 2);
  return _second_order_adjacent_vertices[n-this->n_vertices()][v]; 
}



const unsigned short int InfQuad6::_second_order_adjacent_vertices[2][2] = 
{
  {0, 1}, // vertices adjacent to node 4
  {2, 3}  // vertices adjacent to node 5  
};



std::pair<unsigned short int, unsigned short int>
InfQuad6::second_order_child_vertex (const unsigned int n) const
{
  assert (n >= this->n_vertices());
  assert (n < this->n_nodes());

  return std::pair<unsigned short int, unsigned short int>
    (0, 2*n-7);
}




#endif // ifdef ENABLE_INFINITE_ELEMENTS
