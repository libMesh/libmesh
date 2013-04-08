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
#include "libmesh/off_io.h"
#include "libmesh/mesh_base.h"
#include "libmesh/edge_edge2.h"
#include "libmesh/face_tri3.h"

namespace libMesh
{



// ------------------------------------------------------------
// OFFIO class members

void OFFIO::read(const std::string& name)
{
  std::ifstream in (name.c_str());

  read_stream(in);
}



void OFFIO::read_stream(std::istream& in)
{
  // This is a serial-only process for now;
  // the Mesh should be read on processor 0 and
  // broadcast later
  libmesh_assert_equal_to (this->mesh().processor_id(), 0);

  // Get a reference to the mesh
  MeshBase& the_mesh = MeshInput<MeshBase>::mesh();

  // Clear any existing mesh data
  the_mesh.clear();

  // Check the input buffer
  libmesh_assert (in.good());

  unsigned int nn, ne, nf;

  std::string label;

  // Read the first string.  It should say "OFF"
  in >> label;

  libmesh_assert_equal_to (label, "OFF");

  // read the number of nodes, faces, and edges
  in >> nn >> nf >> ne;


  Real x=0., y=0., z=0.;

  // Read the nodes
  for (unsigned int n=0; n<nn; n++)
    {
      libmesh_assert (in.good());

      in >> x
	 >> y
	 >> z;

      the_mesh.add_point ( Point(x,y,z), n );
    }

  unsigned int nv, nid;

  // Read the elements
  for (unsigned int e=0; e<nf; e++)
    {
      libmesh_assert (in.good());

      // The number of vertices in the element
      in >> nv;

      libmesh_assert(nv == 2 || nv == 3);
      if (e == 0)
      {
        the_mesh.set_mesh_dimension(nv-1);
        if (nv == 3)
        {
#if LIBMESH_DIM < 2
          libMesh::err << "Cannot open dimension 2 mesh file when configured without 2D support." <<
                          std::endl;
          libmesh_error();
#endif
        }
      }

      Elem* elem;
      switch (nv)
      {
        case 2: elem = new Edge2; break;
        case 3: elem = new Tri3 ; break;
        default: libmesh_error();
      }

      elem->set_id(e);
      the_mesh.add_elem (elem);

      for (unsigned int i=0; i<nv; i++)
      {
        in >> nid;
        elem->set_node(i) = the_mesh.node_ptr(nid);
      }
    }
}

} // namespace libMesh
