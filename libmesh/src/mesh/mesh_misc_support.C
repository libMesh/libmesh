// $Id: mesh_misc_support.C,v 1.15 2004-03-23 04:47:29 jwpeterson Exp $

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
#include <iomanip>
#include <fstream>


// Local includes
#include "mesh_base.h"
#include "elem.h"
#include "face_quad4.h"
#include "face_tri3.h"




void MeshBase::read_shanee(const std::string& name)
{
  std::ifstream in (name.c_str());

  read_shanee (in);

  return;
}




void MeshBase::read_shanee (std::istream &in)
{
  /**
   * Clear any existing mesh data
   */
  clear();
  
  // Check input buffer
  assert (in.good());

  // Shanee only works in 2D
  assert (_dim == 2);


  unsigned int nNodes=0, nElem=0;
  std::string dmy;

  in >> nNodes   // Read the number of nodes from the stream
     >> nElem;    // Read the number of elements from the stream
 
  // Read the nodal coordinates. Note that UCD format always
  // stores (x,y,z), and in 2D z=0. We don't need to store this,
  // however.  So, we read in x,y,z for each node and make a point
  // in the proper way based on what dimension we're in
  {
    Real x=0., y=0., z=0.;

    // Reserve space in the nNodes vector to avoid unnecessary
    // allocations
    _nodes.resize (nNodes);
    
    
    for (unsigned int i=0; i<nNodes; i++)
      {
	in >> x   // x-coordinate value
	   >> y;  // y-coordinate value
                  // there is no z-coordinate
	
	node_ptr(i) = Node::build(x,y,z,i);
      }
  }


  
  {
    unsigned int node=0;
    _elements.resize (nElem);

    
    for (unsigned int i=0; i<nElem; i++)
      {
       _elements[i] = new Quad4;
       _elements[i]->set_id (i);

	for (unsigned int n=0; n<elem(i)->n_nodes(); n++)
	  {
	    in >> node; // read the current node
		  
	    elem(i)->set_node(n) = node_ptr(node); // assign the node		                     
	  }
      }   
  }
}





// void MeshBase::read_off(const std::string& name)
// {
//   std::ifstream in (name.c_str());

//   read_off(in);

//   return;
// }




// void MeshBase::read_off(std::istream& in)
// {
//   /**
//    * Clear any existing mesh data
//    */
//   clear();
  
//   // STL only works in 2D
//   assert (_dim == 2);

//   // Check the input buffer
//   assert (in.good());

//   unsigned int nn, ne, nf;

//   std::string label;

//   // Read the first string.  It should say "OFF"
//   in >> label;

//   assert (label == "OFF");

//   // read the number of nodes, faces, and edges
//   in >> nn >> nf >> ne;

//   // resize local vectors.
//   _nodes.resize(nn);
//   _elements.resize(nf);

  
//   Real x=0., y=0., z=0.;
  
//   // Read the nodes
//   for (unsigned int n=0; n<nn; n++)
//     {
//       assert (in.good());

//       in >> x
// 	 >> y
// 	 >> z;
      
//       node_ptr(n) = Node::build(x,y,z,n);
//     }

//   unsigned int dummy, n0, n1, n2;
  
//   // Read the triangles
//   for (unsigned int e=0; e<nf; e++)
//     {
//       assert (in.good());
      
//       _elements[e] = new Tri3;
//       _elements[e]->set_id (e);

//       // The number of nodes in the object
//       in >> dummy;

//       assert (dummy == 3);

//       in >> n0
// 	 >> n1
// 	 >> n2;
      
//       _elements[e]->set_node(0) = node_ptr(n0);
//       _elements[e]->set_node(1) = node_ptr(n1);
//       _elements[e]->set_node(2) = node_ptr(n2);
//     }  
// }
