// $Id: face_quad8.C,v 1.20 2005-01-28 19:14:18 benkirk Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2004  Benjamin S. Kirk, John W. Peterson
  
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
#include "face_quad8.h"




// ------------------------------------------------------------
// Quad8 class static member initializations
const unsigned int Quad8::side_nodes_map[4][3] =
{
  {0, 1, 4}, // Side 0
  {1, 2, 5}, // Side 1
  {2, 3, 6}, // Side 2
  {3, 0, 7}  // Side 3
};


#ifdef ENABLE_AMR

const float Quad8::_embedding_matrix[4][8][8] =
{
  // embedding matrix for child 0
  {
    //         0           1           2           3           4           5           6           7
    {    1.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 0
    {    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000 }, // 1
    {  -0.250000,  -0.250000,  -0.250000,  -0.250000,   0.500000,   0.500000,   0.500000,   0.500000 }, // 2
    {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000 }, // 3
    {   0.375000,  -0.125000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000 }, // 4
    {  -0.187500,  -0.187500,  -0.187500,  -0.187500,   0.750000,   0.375000,   0.250000,   0.375000 }, // 5
    {  -0.187500,  -0.187500,  -0.187500,  -0.187500,   0.375000,   0.250000,   0.375000,   0.750000 }, // 6
    {   0.375000,    0.00000,    0.00000,  -0.125000,    0.00000,    0.00000,    0.00000,   0.750000 }  // 7
  },

  // embedding matrix for child 1
  {
    //         0           1           2           3           4           5           6           7
    {    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000 }, // 0
    {    0.00000,    1.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 1
    {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000 }, // 2
    {  -0.250000,  -0.250000,  -0.250000,  -0.250000,   0.500000,   0.500000,   0.500000,   0.500000 }, // 3
    {  -0.125000,   0.375000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000 }, // 4
    {    0.00000,   0.375000,  -0.125000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000 }, // 5
    {  -0.187500,  -0.187500,  -0.187500,  -0.187500,   0.375000,   0.750000,   0.375000,   0.250000 }, // 6
    {  -0.187500,  -0.187500,  -0.187500,  -0.187500,   0.750000,   0.375000,   0.250000,   0.375000 }  // 7
  },

  // embedding matrix for child 2
  {
    //         0           1           2           3           4           5           6           7
    {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000 }, // 0
    {  -0.250000,  -0.250000,  -0.250000,  -0.250000,   0.500000,   0.500000,   0.500000,   0.500000 }, // 1
    {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000 }, // 2
    {    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 3
    {  -0.187500,  -0.187500,  -0.187500,  -0.187500,   0.375000,   0.250000,   0.375000,   0.750000 }, // 4
    {  -0.187500,  -0.187500,  -0.187500,  -0.187500,   0.250000,   0.375000,   0.750000,   0.375000 }, // 5
    {    0.00000,    0.00000,  -0.125000,   0.375000,    0.00000,    0.00000,   0.750000,    0.00000 }, // 6
    {  -0.125000,    0.00000,    0.00000,   0.375000,    0.00000,    0.00000,    0.00000,   0.750000 }  // 7
  },

  // embedding matrix for child 3
  {
    //         0           1           2           3           4           5           6           7
    {  -0.250000,  -0.250000,  -0.250000,  -0.250000,   0.500000,   0.500000,   0.500000,   0.500000 }, // 0
    {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000 }, // 1
    {    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 2
    {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000 }, // 3
    {  -0.187500,  -0.187500,  -0.187500,  -0.187500,   0.375000,   0.750000,   0.375000,   0.250000 }, // 4
    {    0.00000,  -0.125000,   0.375000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000 }, // 5
    {    0.00000,    0.00000,   0.375000,  -0.125000,    0.00000,    0.00000,   0.750000,    0.00000 }, // 6
    {  -0.187500,  -0.187500,  -0.187500,  -0.187500,   0.250000,   0.375000,   0.750000,   0.375000 }  // 7
  }
};


#endif


// ------------------------------------------------------------
// Quad8 class member functions
unsigned int Quad8::key (const unsigned int s) const
{
  assert (s < this->n_sides());
  
  switch (s)
    {
    case 0:

      return
	this->compute_key (this->node(4));
	      
    case 1:

      return
	this->compute_key (this->node(5));
	
    case 2:

      return
	this->compute_key (this->node(6));
	
    case 3:

      return
	this->compute_key (this->node(7));
    }


  // We will never get here...  Look at the code above.
  error();
  return 0;
}



AutoPtr<Elem> Quad8::build_side (const unsigned int i) const
{
  assert (i < this->n_sides());

  AutoPtr<Elem> ap(new Side<Edge3,Quad8>(this,i));
  return ap;
  
//   Edge3* edge = new Edge3;

//   switch (i)
//     {
//     case 0:
//       {
// 	edge->set_node(0) = this->get_node(0);
// 	edge->set_node(1) = this->get_node(1);
// 	edge->set_node(2) = this->get_node(4);
	
// 	AutoPtr<Elem> ap(edge);  return ap;
//       }
//     case 1:
//       {
// 	edge->set_node(0) = this->get_node(1);
// 	edge->set_node(1) = this->get_node(2);
// 	edge->set_node(2) = this->get_node(5);
	
// 	AutoPtr<Elem> ap(edge);  return ap;
//       }
//     case 2:
//       {
// 	edge->set_node(0) = this->get_node(2);
// 	edge->set_node(1) = this->get_node(3);
// 	edge->set_node(2) = this->get_node(6);
	
// 	AutoPtr<Elem> ap(edge);  return ap;
//       }
//     case 3:
//       {
// 	edge->set_node(0) = this->get_node(3);
// 	edge->set_node(1) = this->get_node(0);
// 	edge->set_node(2) = this->get_node(7);
	
// 	AutoPtr<Elem> ap(edge);  return ap;
//       }
//     default:
//       {
// 	error();
//       }
//     }


//   // We will never get here...  Look at the code above.
//   AutoPtr<Elem> ap(NULL);  return ap;
}






void Quad8::connectivity(const unsigned int sf,
			 const IOPackage iop,
			 std::vector<unsigned int>& conn) const
{
  assert (_nodes != NULL);
  assert (sf < this->n_sub_elem());

  switch (iop)
    {
      // Note: TECPLOT connectivity is output as four triangles with
      // a central quadrilateral.  Therefore, the first four connectivity
      // arrays are degenerate quads (triangles in Tecplot).
    case TECPLOT:
      {
	// Create storage
	conn.resize(4);

	switch(sf)
	  {
	  case 0:
	    // linear sub-tri 0
	    conn[0] = this->node(0)+1;
	    conn[1] = this->node(4)+1;
	    conn[2] = this->node(7)+1;
	    conn[3] = this->node(7)+1;

	    return;

	  case 1:
	    // linear sub-tri 1
	    conn[0] = this->node(4)+1;
	    conn[1] = this->node(1)+1;
	    conn[2] = this->node(5)+1;
	    conn[3] = this->node(5)+1;

	    return;

	  case 2:
	    // linear sub-tri 2
	    conn[0] = this->node(5)+1;
	    conn[1] = this->node(2)+1;
	    conn[2] = this->node(6)+1;
	    conn[3] = this->node(6)+1;

	    return;

	  case 3:
	    // linear sub-tri 3
	    conn[0] = this->node(7)+1;
	    conn[1] = this->node(6)+1;
	    conn[2] = this->node(3)+1;
	    conn[3] = this->node(3)+1;

	    return;

	  case 4:
	    // linear sub-quad
	    conn[0] = this->node(4)+1;
	    conn[1] = this->node(5)+1;
	    conn[2] = this->node(6)+1;
	    conn[3] = this->node(7)+1;

	    return;

	  default:
	    error();
	  }
      }

      
      // Note: VTK connectivity is output as four triangles with
      // a central quadrilateral.  Therefore most of the connectivity
      // arrays have length three.
    case VTK:
      {
	// Create storage
	conn.resize(3);

	switch (sf)
	  {
	  case 0:
	    // linear sub-tri 0
	    conn[0] = this->node(0);
	    conn[1] = this->node(4);
	    conn[2] = this->node(7);

	    return;

	  case 1:
	    // linear sub-tri 1
	    conn[0] = this->node(4);
	    conn[1] = this->node(1);
	    conn[2] = this->node(5);

	    return;

	  case 2:
	    // linear sub-tri 2
	    conn[0] = this->node(5);
	    conn[1] = this->node(2);
	    conn[2] = this->node(6);

	    return;

	  case 3:
	    // linear sub-tri 3
	    conn[0] = this->node(7);
	    conn[1] = this->node(6);
	    conn[2] = this->node(3);

	    return;

	  case 4:
	    conn.resize(4);
      
	    // linear sub-quad
	    conn[0] = this->node(4);
	    conn[1] = this->node(5);
	    conn[2] = this->node(6);
	    conn[3] = this->node(7);

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





// void Quad8::tecplot_connectivity(const unsigned int sf,
// 				 std::vector<unsigned int>& conn) const
// {
//   assert (_nodes != NULL);
//   assert (sf < this->n_sub_elem());

//   // std::vector<unsigned int> conn(4);
//   conn.resize(4);

//   switch(sf)
//     {
//     case 0:
//       // linear sub-tri 0
//       conn[0] = this->node(0)+1;
//       conn[1] = this->node(4)+1;
//       conn[2] = this->node(7)+1;
//       conn[3] = this->node(7)+1;

//       return;

//     case 1:
//       // linear sub-tri 1
//       conn[0] = this->node(4)+1;
//       conn[1] = this->node(1)+1;
//       conn[2] = this->node(5)+1;
//       conn[3] = this->node(5)+1;

//       return;

//     case 2:
//       // linear sub-tri 2
//       conn[0] = this->node(5)+1;
//       conn[1] = this->node(2)+1;
//       conn[2] = this->node(6)+1;
//       conn[3] = this->node(6)+1;

//       return;

//     case 3:
//       // linear sub-tri 3
//       conn[0] = this->node(7)+1;
//       conn[1] = this->node(6)+1;
//       conn[2] = this->node(3)+1;
//       conn[3] = this->node(3)+1;

//       return;

//     case 4:
//       // linear sub-quad
//       conn[0] = this->node(4)+1;
//       conn[1] = this->node(5)+1;
//       conn[2] = this->node(6)+1;
//       conn[3] = this->node(7)+1;

//       return;

//     default:
//       error();
//     }

//   error();
// }



// void Quad8::vtk_connectivity(const unsigned int sf,
// 			     std::vector<unsigned int> *conn) const
// {
//   assert (_nodes != NULL);
//   assert (sf < this->n_sub_elem());
  
//   if (conn == NULL)
//     conn = new std::vector<unsigned int>;

//   conn->resize(3);

//   switch (sf)
//     {
//     case 0:
//       // linear sub-tri 0
//       (*conn)[0] = this->node(0);
//       (*conn)[1] = this->node(4);
//       (*conn)[2] = this->node(7);

//       return;

//     case 1:
//       // linear sub-tri 1
//       (*conn)[0] = this->node(4);
//       (*conn)[1] = this->node(1);
//       (*conn)[2] = this->node(5);

//       return;

//     case 2:
//       // linear sub-tri 2
//       (*conn)[0] = this->node(5);
//       (*conn)[1] = this->node(2);
//       (*conn)[2] = this->node(6);

//       return;

//     case 3:
//       // linear sub-tri 3
//       (*conn)[0] = this->node(7);
//       (*conn)[1] = this->node(6);
//       (*conn)[2] = this->node(3);

//       return;

//     case 4:
//       conn->resize(4);
      
//       // linear sub-quad
//       (*conn)[0] = this->node(4);
//       (*conn)[1] = this->node(5);
//       (*conn)[2] = this->node(6);
//       (*conn)[3] = this->node(7);

//       return;

//     default:
//       error();
//     }

//   error();

//   return;
// }



// unsigned int Quad8::vtk_element_type (const unsigned int sf) const
// {
//   if (sf == 4)
//     return 9;

//   return 5;
// }




unsigned short int Quad8::second_order_adjacent_vertex (const unsigned int n,
							const unsigned int v) const
{ 
  assert (n >= this->n_vertices());
  assert (n <  this->n_nodes());
  assert (v < 2);
  // use the matrix from \p face_quad.C
  return _second_order_adjacent_vertices[n-this->n_vertices()][v]; 
}


