// $Id: matlab_io.C,v 1.1 2004-03-23 04:17:26 jwpeterson Exp $

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
#include "matlab_io.h"
#include "face_tri3.h"

// ------------------------------------------------------------
// MatlabIO class members

void MatlabIO::read(const std::string& name)
{
  std::ifstream in (name.c_str());

  this->read_stream(in);
}


void MatlabIO::read_stream(std::istream& in)
{
  // Get a reference to the mesh
  Mesh& mesh = this->mesh();

  // Clear any existing mesh data
  mesh.clear();
  
  // PDE toolkit only works in 2D
  assert(mesh.mesh_dimension() == 2);

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

    for (unsigned int i=0; i<nNodes; i++)
      {
	in >> x   // x-coordinate value
	   >> y;  // y-coordinate value

	mesh.add_point ( Point(x,y,z) );
      }
  }

  // Read the elements (elements)
  {
    unsigned int node=0, dummy=0;
    
    for (unsigned int i=0; i<nElem; i++)
      {
	Elem* elem = mesh.add_elem (new Tri3); // Always build a triangle
	
	for (unsigned int n=0; n<3; n++)  // Always read three 3 nodes
	  {
	    in >> node;
	    elem->set_node(n) = mesh.node_ptr(node-1);  // Assign the node number
	  }
	
	// There is an additional subdomain number here,
	// so we read it and get rid of it!
	in >> dummy;
      }
  }
  
}
