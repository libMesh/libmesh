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


// C++ includes
#include <fstream>


// Local includes
#include "libmesh_config.h"
#include "ucd_io.h"
#include "mesh_base.h"
#include "face_quad4.h"
#include "face_tri3.h"
#include "cell_tet4.h"
#include "cell_hex8.h"
#include "cell_prism6.h"

#ifdef HAVE_GZSTREAM
# include "gzstream.h" // For reading/writing compressed streams
#endif





// ------------------------------------------------------------
// UCDIO class members
void UCDIO::read (const std::string& file_name)
{
  if (file_name.rfind(".gz") < file_name.size())
    {
#ifdef HAVE_GZSTREAM
      
      igzstream in_stream (file_name.c_str());
      this->read_implementation (in_stream);
      
#else
      
      std::cerr << "ERROR:  You must have the zlib.h header "
		<< "files and libraries to read and write "
		<< "compressed streams."
		<< std::endl;
      libmesh_error();
      
#endif
      return;      
    }
  
  else
    {
      std::ifstream in_stream (file_name.c_str());
      this->read_implementation (in_stream);
      return;
    }
}



void UCDIO::write (const std::string& file_name)
{
  if (file_name.rfind(".gz") < file_name.size())
    {
#ifdef HAVE_GZSTREAM
      
      ogzstream out_stream (file_name.c_str());
      this->write_implementation (out_stream);
      
#else
      
      std::cerr << "ERROR:  You must have the zlib.h header "
		<< "files and libraries to read and write "
		<< "compressed streams."
		<< std::endl;
      libmesh_error();
      
#endif
      return;      
    }
  
  else
    {
      std::ofstream out_stream (file_name.c_str());
      this->write_implementation (out_stream);
      return;
    }
}



void UCDIO::read_implementation (std::istream& in)
{
  // This is a serial-only process for now;
  // the Mesh should be read on processor 0 and
  // broadcast later
  assert(libMesh::processor_id() == 0);

  // Check input buffer
  assert (in.good());

  MeshBase& mesh = MeshInput<MeshBase>::mesh();
  
  // UCD doesn't work in 1D
  assert (mesh.mesh_dimension() != 1);

  this->skip_comment_lines (in, '#');
  
  unsigned int nNodes=0, nElem=0, dummy=0;

  in >> nNodes   // Read the number of nodes from the stream
     >> nElem    // Read the number of elements from the stream
     >> dummy
     >> dummy
     >> dummy;


  // Read the nodal coordinates. Note that UCD format always
  // stores (x,y,z), and in 2D z=0. We don't need to store this,
  // however.  So, we read in x,y,z for each node and make a point
  // in the proper way based on what dimension we're in
  {
    Point xyz;    
    
    for (unsigned int i=0; i<nNodes; i++)
      {
	assert (in.good());
	
	in >> dummy   // Point number
	   >> xyz(0)  // x-coordinate value
	   >> xyz(1)  // y-coordinate value
	   >> xyz(2); // z-coordinate value

	// Build the node
	mesh.add_point (xyz, i);
      }
  }


  
  // Read the elements from the stream. Notice that the UCD node-numbering
  // scheme is 1-based, and we just created a 0-based scheme above
  // (which is of course what we want). So, when we read in the nodal
  // connectivity for each element we need to take 1 off the value of
  // each node so that we get the right thing.
  {
    unsigned int material_id=0, node=0;
    std::string type;
    
    for (unsigned int i=0; i<nElem; i++)
      {
	Elem* elem = NULL;

	assert (in.good());
	
	in >> dummy        // Cell number, means nothing to us
	   >> material_id  // doesn't mean anything at present, might later
	   >> type;        // string describing cell type:
	                   // either tri, quad, tet, hex, or prism for the
	                   // obvious cases


	                   // Now read the connectivity.
	if (type == "quad")
	  elem = new Quad4;
	else if (type == "tri")
	  elem = new Tri3;
	else if (type == "hex")
	  elem = new Hex8;
	else if (type == "tet")
	  elem = new Tet4;
	else if (type == "prism")
	  elem = new Prism6;
	else
	  libmesh_error();
	
	for (unsigned int n=0; n<elem->n_nodes(); n++)
	  {
	    assert (in.good());
	    
	    in >> node; // read the current node
	    node -= 1;  // UCD is 1-based, so subtract

	    assert (node < mesh.n_nodes());
	    
	    elem->set_node(n) =
	      mesh.node_ptr(node); // assign the node
	  }

	// Add the element to the mesh
	elem->set_id(i);
	mesh.add_elem (elem);
      }
  }  
}



void UCDIO::write_implementation (std::ostream& out)
{
  assert (out.good());

  const MeshBase& mesh = MeshOutput<MeshBase>::mesh();
  
  // UCD doesn't work in 1D
  assert (mesh.mesh_dimension() != 1);
  
  // Write header to stream
  out << "# This file was generated by:\n"
      << "#\n"
      << "# $Id$\n"
      << "#\n"
      << "# For a description of the UCD format see the AVS Developer's guide.\n"
      << "#\n";

  
  // Write the mesh info
  out << mesh.n_nodes() << " "
      << mesh.n_elem()  << " "
      << " 0 0 0\n";

  // Write the coordinates
  {
//     const_node_iterator       it  (mesh.nodes_begin());
//     const const_node_iterator end (mesh.nodes_end());

    MeshBase::const_node_iterator       it  = mesh.nodes_begin();
    const MeshBase::const_node_iterator end = mesh.nodes_end();
	
    unsigned int n=1; // 1-based node number for UCD
    
    for (; it != end; ++it)
      {
	assert (out.good());
	
	out << n++ << "\t";
	(*it)->write_unformatted(out);
      }
  }

  // Write the elements
  {
    std::string type[] =
      { "edge",  "edge",  "edge",	
	"tri",   "tri",	 
	"quad",  "quad",  "quad",
	"tet",   "tet",	 
	"hex",   "hex",   "hex",
	"prism", "prism", "prism",	
	"pyramid" };
    
//     const_elem_iterator       it  (mesh.elements_begin());
//     const const_elem_iterator end (mesh.elements_end());

    MeshBase::const_element_iterator it  = mesh.elements_begin();
    const MeshBase::const_element_iterator end = mesh.elements_end();

    unsigned int e=1; // 1-based element number for UCD
    
    for (; it != end; ++it)
      {
	assert (out.good());
	
	out << e++ << " 0 " << type[(*it)->type()] << "\t"; 
	// (*it)->write_ucd_connectivity(out);
	(*it)->write_connectivity(out, UCD);
      }
  }
}
