// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include <fstream>

// Local includes
#include "libmesh/matlab_io.h"
#include "libmesh/mesh_base.h"
#include "libmesh/face_tri3.h"

namespace libMesh
{

// ------------------------------------------------------------
// MatlabIO class members

void MatlabIO::read(const std::string& name)
{
  std::ifstream in (name.c_str());

  this->read_stream(in);
}


void MatlabIO::read_stream(std::istream& in)
{
  // This is a serial-only process for now;
  // the Mesh should be read on processor 0 and
  // broadcast later
  libmesh_assert_equal_to (this->mesh().processor_id(), 0);

  // Get a reference to the mesh
  MeshBase& the_mesh = MeshInput<MeshBase>::mesh();

  // Clear any existing mesh data
  the_mesh.clear();

  // PDE toolkit only works in 2D
  the_mesh.set_mesh_dimension(2);

#if LIBMESH_DIM < 2
  libMesh::err << "Cannot open dimension 2 mesh file when configured without 2D support." <<
                  std::endl;
  libmesh_error();
#endif

  // Check the input buffer
  libmesh_assert (in.good());

  unsigned int nNodes=0, nElem=0;

  in >> nNodes   // Read the number of nodes
     >> nElem;   // Read the number of elements

  // Sort of check that it worked
  libmesh_assert_greater (nNodes, 0);
  libmesh_assert_greater (nElem, 0);

  // Read the nodal coordinates
  {
    Real x=0., y=0., z=0.;

    for (unsigned int i=0; i<nNodes; i++)
      {
	in >> x   // x-coordinate value
	   >> y;  // y-coordinate value

	the_mesh.add_point ( Point(x,y,z), i);
      }
  }

  // Read the elements (elements)
  {
    unsigned int node=0, dummy=0;

    for (unsigned int i=0; i<nElem; i++)
      {
	Elem* elem = new Tri3; // Always build a triangle
        elem->set_id(i);
	the_mesh.add_elem (elem);

	for (unsigned int n=0; n<3; n++)  // Always read three 3 nodes
	  {
	    in >> node;
	    elem->set_node(n) = the_mesh.node_ptr(node-1);  // Assign the node number
	  }

	// There is an additional subdomain number here,
	// so we read it and get rid of it!
	in >> dummy;
      }
  }

}

} // namespace libMesh
