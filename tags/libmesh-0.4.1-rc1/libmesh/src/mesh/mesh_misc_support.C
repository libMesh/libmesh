// $Id: mesh_misc_support.C,v 1.11 2003-09-02 18:02:43 benkirk Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002-2003  Benjamin S. Kirk, John W. Peterson
  
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

	for (unsigned int n=0; n<elem(i)->n_nodes(); n++)
	  {
	    in >> node; // read the current node
		  
	    elem(i)->set_node(n) = node_ptr(node); // assign the node		                     
	  }
      }   
  }
}




void MeshBase::read_matlab(const std::string& name)
{
  std::ifstream in (name.c_str());

  read_matlab(in);

  return;
}




void MeshBase::read_matlab(std::istream& in)
{

  // A VALID INPUT FILE for this type of mesh should be
  // generated in Matlab with the following steps:
  // 1.) Draw the domain and triangulate it in the GUI
  // 2.) Export the mesh to matlab using Mesh->Export Mesh
  // 3.) Create a file with this script:
  //     fid = fopen('filename', 'w');
  //     fprintf(fid, '%d %d \n', length(p), length(t));
  //     fprintf(fid, '%f %f \n', p);
  //     fprintf(fid, '%d %d %d %d \n', t);
  //     fclose(fid);

  // What's going on here? 
  // There is no standard for exporting PDE toolkit meshes
  // to files in Matlab.  When you choose "export mesh" in the GUI,
  // it returns three matrices that it likes to call
  // p, e, and t.  All meshes (as far as I can tell) that
  // come from the PDE toolkit are 2D triangle meshes.
  
  // p is the point matrix...
  // Row 1: x coordinate
  // Row 2: y coordinate
  
  // e is the edge matrix ... 
  // Row 1: starting point number          (dummy)
  // Row 2: ending point number            (dummy)
  // Row 3: starting parameter value (?)   (dummy)
  // Row 4: ending parameter value (?)     (dummy) 
  // Row 5: boundary segment number (?)    (dummy)
  // Row 6: left-hand subdomain number     (dummy)
  // Row 7: right-hand subdomain number    (dummy)

  // t is the triangle matrix ...
  // Row 1: Node number 1
  // Row 2: Node number 2
  // Row 3: Node number 3
  // Row 4: subdomain number               (dummy)

  // There are some important things to notice here:
  // o The implied ordering of the p matrix is 1..N
  // o The e matrix is entirely irrelevant in this code
  // o All of the matrices are row based

  /**
   * Clear any existing mesh data
   */
  clear();
  
  // PDE toolkit only works in 2D
  assert(_dim == 2);

  // Check the input buffer
  assert (in.good());

  unsigned int nNodes=0, nElem=0;

  in >> nNodes   // Read the number of nodes
     >> nElem;   // Read the number of elements

  // Sort of check that it worked
  assert(nNodes > 0);
  assert(nElem > 0);

  // Read the nodal coordinates
  {
    Real x=0., y=0., z=0.;

    // Resize the nodes vector
    _nodes.resize(nNodes);

    for (unsigned int i=0; i<nNodes; i++)
      {
	in >> x   // x-coordinate value
	   >> y;  // y-coordinate value

	node_ptr(i) = Node::build(x,y,z,i);
      }
  }

  // Read the elements (elements)
  {
    unsigned int node=0, dummy=0;
    
    // Resize the elements vector
    _elements.resize(nElem);

    for (unsigned int i=0; i<nElem; i++)
      {
	_elements[i] = new Tri3; // Always build a triangle
	for (unsigned int n=0; n<3; n++)  // Always read three 3 nodes
	  {
	    in >> node;
	    elem(i)->set_node(n) = node_ptr(node-1);  // Assign the node number
	  }
	
	// There is an additional subdomain number here,
	// so we read it and get rid of it!
	in >> dummy;
      }
  }
  
}




void MeshBase::read_off(const std::string& name)
{
  std::ifstream in (name.c_str());

  read_off(in);

  return;
}




void MeshBase::read_off(std::istream& in)
{
  /**
   * Clear any existing mesh data
   */
  clear();
  
  // STL only works in 2D
  assert (_dim == 2);

  // Check the input buffer
  assert (in.good());

  unsigned int nn, ne, nf;

  std::string label;

  // Read the first string.  It should say "OFF"
  in >> label;

  assert (label == "OFF");

  // read the number of nodes, faces, and edges
  in >> nn >> nf >> ne;

  // resize local vectors.
  _nodes.resize(nn);
  _elements.resize(nf);

  
  Real x=0., y=0., z=0.;
  
  // Read the nodes
  for (unsigned int n=0; n<nn; n++)
    {
      assert (in.good());

      in >> x
	 >> y
	 >> z;
      
      node_ptr(n) = Node::build(x,y,z,n);
    }

  unsigned int dummy, n0, n1, n2;
  
  // Read the triangles
  for (unsigned int e=0; e<nf; e++)
    {
      assert (in.good());
      
      _elements[e] = new Tri3;

      // The number of nodes in the object
      in >> dummy;

      assert (dummy == 3);

      in >> n0
	 >> n1
	 >> n2;
      
      _elements[e]->set_node(0) = node_ptr(n0);
      _elements[e]->set_node(1) = node_ptr(n1);
      _elements[e]->set_node(2) = node_ptr(n2);
    }  
}
