// $Id: off_io.C,v 1.1 2004-03-23 04:47:29 jwpeterson Exp $

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
#include "off_io.h"
#include "face_tri3.h"



// ------------------------------------------------------------
// OFFIO class members

void OFFIO::read(const std::string& name)
{
  std::ifstream in (name.c_str());

  read_stream(in);
}



void OFFIO::read_stream(std::istream& in)
{
  // Get a reference to the mesh
  Mesh& mesh = this->mesh();
  
  // Clear any existing mesh data
  mesh.clear();
  
  // STL only works in 2D
  assert (mesh.mesh_dimension() == 2);

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
  //  _nodes.resize(nn);
  //_elements.resize(nf);

  
  Real x=0., y=0., z=0.;
  
  // Read the nodes
  for (unsigned int n=0; n<nn; n++)
    {
      assert (in.good());

      in >> x
	 >> y
	 >> z;
      
      // node_ptr(n) = Node::build(x,y,z,n);
      mesh.add_point ( Point(x,y,z) );
    }

  unsigned int dummy, n0, n1, n2;
  
  // Read the triangles
  for (unsigned int e=0; e<nf; e++)
    {
      assert (in.good());
      
      // _elements[e] = new Tri3;
      // _elements[e]->set_id (e);
      Elem* elem = mesh.add_elem (new Tri3);

      // The number of nodes in the object
      in >> dummy;

      assert (dummy == 3);

      in >> n0
	 >> n1
	 >> n2;
      
      elem->set_node(0) = mesh.node_ptr(n0);
      elem->set_node(1) = mesh.node_ptr(n1);
      elem->set_node(2) = mesh.node_ptr(n2);
    }  
}
