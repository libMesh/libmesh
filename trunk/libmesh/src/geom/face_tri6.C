// $Id: face_tri6.C,v 1.22 2005-02-25 19:16:39 roystgnr Exp $

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

// C++ includes

// Local includes
#include "side.h"
#include "edge_edge3.h"
#include "face_tri6.h"




// ------------------------------------------------------------
// Tri6 class static member initializations
const unsigned int Tri6::side_nodes_map[3][3] =
{
  {0, 1, 3}, // Side 0
  {1, 2, 4}, // Side 1
  {2, 0, 5}  // Side 2
};


#ifdef ENABLE_AMR

const float Tri6::_embedding_matrix[4][6][6] =
{
  // embedding matrix for child 0
  {
    //  0      1      2    3    4    5
    { 1.0,   0.0,   0.0, 0.0, 0.0, 0.0}, // 0
    { 0.0,   0.0,   0.0, 1.0, 0.0, 0.0}, // 1
    { 0.0,   0.0,   0.0, 0.0, 0.0, 1.0}, // 2
    {.375, -.125,   0.0, .75, 0.0, 0.0}, // 3
    { 0.0, -.125, -.125, 0.5, .25, 0.5}, // 4
    {.375,   0.0, -.125, 0.0, 0.0, .75}  // 5
  },

  // embedding matrix for child 1
  {
    //  0      1      2    3    4    5
    {  0.0,  0.0,   0.0, 1.0, 0.0, 0.0}, // 0
    {  0.0,  1.0,   0.0, 0.0, 0.0, 0.0}, // 1
    {  0.0,  0.0,   0.0, 0.0, 1.0, 0.0}, // 2
    {-.125, .375,   0.0, .75, 0.0, 0.0}, // 3
    {  0.0, .375, -.125, 0.0, .75, 0.0}, // 4
    {-.125,  0.0, -.125, 0.5, 0.5, .25}  // 5
  },

  // embedding matrix for child 2
  {
    //  0       1     2    3    4    5
    {  0.0,   0.0,  0.0, 0.0, 0.0, 1.0}, // 0
    {  0.0,   0.0,  0.0, 0.0, 1.0, 0.0}, // 1
    {  0.0,   0.0,  1.0, 0.0, 0.0, 0.0}, // 2
    {-.125, -.125,  0.0, .25, 0.5, 0.5}, // 3
    {  0.0, -.125, .375, 0.0, .75, 0.0}, // 4
    {-.125,   0.0, .375, 0.0, 0.0, .75}  // 5
  },

  // embedding matrix for child 3
  {
    //  0       1      2    3    4    5
    {  0.0,   0.0,   0.0, 1.0, 0.0, 0.0}, // 0
    {  0.0,   0.0,   0.0, 0.0, 1.0, 0.0}, // 1
    {  0.0,   0.0,   0.0, 0.0, 0.0, 1.0}, // 2
    {-.125,   0.0, -.125, 0.5, 0.5, .25}, // 3
    {-.125, -.125,   0.0, .25, 0.5, 0.5}, // 4
    {  0.0, -.125, -.125, 0.5, .25, 0.5}  // 5
  }
};

#endif



// ------------------------------------------------------------
// Tri6 class member functions

bool Tri6::is_vertex(const unsigned int i) const
{
  if (i < 3)
    return true;
  return false;
}

bool Tri6::is_edge(const unsigned int i) const
{
  if (i < 3)
    return false;
  return true;
}

bool Tri6::is_face(const unsigned int) const
{
  return false;
}

bool Tri6::is_node_on_side(const unsigned int n,
			   const unsigned int s) const
{
  assert(s < n_sides());
  for (unsigned int i = 0; i != 3; ++i)
    if (side_nodes_map[s][i] == n)
      return true;
  return false;
}

unsigned int Tri6::key (const unsigned int s) const
{
  assert (s < this->n_sides());

  switch (s)
    {
    case 0:

      return
	this->compute_key (this->node(3));
	
    case 1:

      return
	this->compute_key (this->node(4));
	
    case 2:

      return
	this->compute_key (this->node(5));
    }

  
  // We will never get here...  Look at the code above.
  error();
  return 0;
}



AutoPtr<Elem> Tri6::build_side (const unsigned int i) const
{
  assert (i < this->n_sides());

  AutoPtr<Elem> ap(new Side<Edge3,Tri6>(this,i));
  return ap;
  
//   Edge3* edge = new Edge3;

//   switch (i)
//     {
//     case 0:
//       {
// 	edge->set_node(0) = this->get_node(0);
// 	edge->set_node(1) = this->get_node(1);
// 	edge->set_node(2) = this->get_node(3);
	
// 	AutoPtr<Elem> ap(edge);  return ap;
//       }
//     case 1:
//       {
// 	edge->set_node(0) = this->get_node(1);
// 	edge->set_node(1) = this->get_node(2);
// 	edge->set_node(2) = this->get_node(4);
	
// 	AutoPtr<Elem> ap(edge);  return ap;
//       }
//     case 2:
//       {
// 	edge->set_node(0) = this->get_node(2);
// 	edge->set_node(1) = this->get_node(0);
// 	edge->set_node(2) = this->get_node(5);
	
// 	AutoPtr<Elem> ap(edge);  return ap;
//       }
//     default:
//       {
// 	error();
//       }
//     }

  
//   // We will never get here...  Look at the code above.
//   error();
//   AutoPtr<Elem> ap(NULL);  return ap;
}


void Tri6::connectivity(const unsigned int sf,
			const IOPackage iop,
			std::vector<unsigned int>& conn) const
{
  assert (_nodes != NULL);
  assert (sf < this->n_sub_elem());
  assert (iop != INVALID_IO_PACKAGE);

  switch (iop)
    {
    case TECPLOT:
      {
	conn.resize(4);
	switch(sf)
	  {
	  case 0:
	    // linear sub-triangle 0
	    conn[0] = this->node(0)+1;
	    conn[1] = this->node(3)+1;
	    conn[2] = this->node(5)+1;
	    conn[3] = this->node(5)+1;

	    return;

	  case 1:
	    // linear sub-triangle 1
	    conn[0] = this->node(3)+1;
	    conn[1] = this->node(1)+1;
	    conn[2] = this->node(4)+1;
	    conn[3] = this->node(4)+1;

	    return;

	  case 2:
	    // linear sub-triangle 2
	    conn[0] = this->node(5)+1;
	    conn[1] = this->node(4)+1;
	    conn[2] = this->node(2)+1;
	    conn[3] = this->node(2)+1;

	    return;

	  case 3:
	    // linear sub-triangle 3
	    conn[0] = this->node(3)+1;
	    conn[1] = this->node(4)+1;
	    conn[2] = this->node(5)+1;
	    conn[3] = this->node(5)+1;

	    return;

	  default:
	    error();
	  }
      }

    case VTK:
      {
	conn.resize(3);
	switch(sf)
	  {
	  case 0:
	    // linear sub-triangle 0
	    conn[0] = this->node(0);
	    conn[1] = this->node(3);
	    conn[2] = this->node(5);

	    return;

	  case 1:
	    // linear sub-triangle 1
	    conn[0] = this->node(3);
	    conn[1] = this->node(1);
	    conn[2] = this->node(4);

	    return;

	  case 2:
	    // linear sub-triangle 2
	    conn[0] = this->node(5);
	    conn[1] = this->node(4);
	    conn[2] = this->node(2);

	    return;

	  case 3:
	    // linear sub-triangle 3
	    conn[0] = this->node(3);
	    conn[1] = this->node(4);
	    conn[2] = this->node(5);

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




// void Tri6::tecplot_connectivity(const unsigned int sf,
// 				std::vector<unsigned int>& conn) const
// {
//   assert (_nodes != NULL);
//   assert (sf < this->n_sub_elem());

//   // std::vector<unsigned int> conn(4);
//   conn.resize(4);

//   switch(sf)
//     {
//     case 0:
//       // linear sub-triangle 0
//       conn[0] = this->node(0)+1;
//       conn[1] = this->node(3)+1;
//       conn[2] = this->node(5)+1;
//       conn[3] = this->node(5)+1;

//       return;

//     case 1:
//       // linear sub-triangle 1
//       conn[0] = this->node(3)+1;
//       conn[1] = this->node(1)+1;
//       conn[2] = this->node(4)+1;
//       conn[3] = this->node(4)+1;

//       return;

//     case 2:
//       // linear sub-triangle 2
//       conn[0] = this->node(5)+1;
//       conn[1] = this->node(4)+1;
//       conn[2] = this->node(2)+1;
//       conn[3] = this->node(2)+1;

//       return;

//     case 3:
//       // linear sub-triangle 3
//       conn[0] = this->node(3)+1;
//       conn[1] = this->node(4)+1;
//       conn[2] = this->node(5)+1;
//       conn[3] = this->node(5)+1;

//       return;

//     default:
//       error();
//     }

//   error();
// }



// void Tri6::vtk_connectivity(const unsigned int sf,
// 			    std::vector<unsigned int> *conn) const
// {
//   assert (_nodes != NULL);
//   assert (sf < this->n_sub_elem());
  
//   if (conn == NULL)
//     conn = new std::vector<unsigned int>;

//   conn->resize(3);
  
//   switch(sf)
//     {
//     case 0:
//       // linear sub-triangle 0
//       (*conn)[0] = this->node(0);
//       (*conn)[1] = this->node(3);
//       (*conn)[2] = this->node(5);

//       return;

//     case 1:
//       // linear sub-triangle 1
//       (*conn)[0] = this->node(3);
//       (*conn)[1] = this->node(1);
//       (*conn)[2] = this->node(4);

//       return;

//     case 2:
//       // linear sub-triangle 2
//       (*conn)[0] = this->node(5);
//       (*conn)[1] = this->node(4);
//       (*conn)[2] = this->node(2);

//       return;

//     case 3:
//       // linear sub-triangle 3
//       (*conn)[0] = this->node(3);
//       (*conn)[1] = this->node(4);
//       (*conn)[2] = this->node(5);

//       return;

//     default:
//       error();
//     }

//   error();
  
//   return;
// }





unsigned short int Tri6::second_order_adjacent_vertex (const unsigned int n,
						       const unsigned int v) const
{ 
  assert (n >= this->n_vertices());
  assert (n <  this->n_nodes());
  assert (v < 2);
  return _second_order_adjacent_vertices[n-this->n_vertices()][v]; 
}



const unsigned short int Tri6::_second_order_adjacent_vertices[3][2] = 
{
  {0, 1}, // vertices adjacent to node 3 
  {1, 2}, // vertices adjacent to node 4 
  {0, 2}  // vertices adjacent to node 5  
};




