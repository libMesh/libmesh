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
#include "libmesh/libmesh_config.h"
#include "libmesh/ucd_io.h"
#include "libmesh/mesh_base.h"
#include "libmesh/face_quad4.h"
#include "libmesh/face_tri3.h"
#include "libmesh/cell_tet4.h"
#include "libmesh/cell_hex8.h"
#include "libmesh/cell_prism6.h"

#ifdef LIBMESH_HAVE_GZSTREAM
# include "gzstream.h" // For reading/writing compressed streams
#endif


namespace libMesh
{




// ------------------------------------------------------------
// UCDIO class members
void UCDIO::read (const std::string& file_name)
{
  if (file_name.rfind(".gz") < file_name.size())
    {
#ifdef LIBMESH_HAVE_GZSTREAM

      igzstream in_stream (file_name.c_str());
      this->read_implementation (in_stream);

#else

      libMesh::err << "ERROR:  You must have the zlib.h header "
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
#ifdef LIBMESH_HAVE_GZSTREAM

      ogzstream out_stream (file_name.c_str());
      this->write_implementation (out_stream);

#else

      libMesh::err << "ERROR:  You must have the zlib.h header "
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
  libmesh_assert_equal_to (MeshOutput<MeshBase>::mesh().processor_id(), 0);

  // Check input buffer
  libmesh_assert (in.good());

  MeshBase& mesh = MeshInput<MeshBase>::mesh();

  // Keep track of what kinds of elements this file contains
  elems_of_dimension.clear();
  elems_of_dimension.resize(4, false);

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
	libmesh_assert (in.good());

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

	libmesh_assert (in.good());

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
	    libmesh_assert (in.good());

	    in >> node; // read the current node
	    node -= 1;  // UCD is 1-based, so subtract

	    libmesh_assert_less (node, mesh.n_nodes());

	    elem->set_node(n) =
	      mesh.node_ptr(node); // assign the node
	  }

        elems_of_dimension[elem->dim()] = true;

	// Add the element to the mesh
	elem->set_id(i);
	mesh.add_elem (elem);
      }

    // Set the mesh dimension to the largest encountered for an element
    for (unsigned int i=0; i!=4; ++i)
      if (elems_of_dimension[i])
        mesh.set_mesh_dimension(i);

#if LIBMESH_DIM < 3
    if (mesh.mesh_dimension() > LIBMESH_DIM)
      {
        libMesh::err << "Cannot open dimension " <<
		        mesh.mesh_dimension() <<
		        " mesh file when configured without " <<
                        mesh.mesh_dimension() << "D support." <<
                        std::endl;
        libmesh_error();
      }
#endif
  }
}



void UCDIO::write_implementation (std::ostream& out_stream)
{
  libmesh_assert (out_stream.good());

  const MeshBase& mesh = MeshOutput<MeshBase>::mesh();

  // UCD doesn't work in 1D
  libmesh_assert_not_equal_to (mesh.mesh_dimension(), 1);
  if(mesh.mesh_dimension() != 3)
    {
      libMesh::err << "Error: Can't write boundary elements for meshes of dimension less than 3"
		   << "Mesh dimension = " << mesh.mesh_dimension()
		   << std::endl;
      libmesh_error();
    }

  // Write header
  this->write_header(out_stream,mesh,mesh.n_elem(),0);

  // Write the node coordinates
  this->write_nodes(out_stream,mesh);

  // Write the elements
  this->write_interior_elems(out_stream,mesh);

  return;
}

void UCDIO::write_header(std::ostream& out_stream, const MeshBase& mesh,
			 dof_id_type n_elems, unsigned int n_vars )
{
  libmesh_assert (out_stream.good());
  // TODO: We used to print out the SVN revision here when we did keyword expansions...
  out_stream << "# For a description of the UCD format see the AVS Developer's guide.\n"
      << "#\n";

  // Write the mesh info
  out_stream << mesh.n_nodes() << " "
      << n_elems  << " "
      << n_vars << " "
      << " 0 0\n";
  return;
}

void UCDIO::write_nodes(std::ostream& out_stream, const MeshBase& mesh)
{
  MeshBase::const_node_iterator       it  = mesh.nodes_begin();
  const MeshBase::const_node_iterator end = mesh.nodes_end();

  unsigned int n=1; // 1-based node number for UCD

  // Write the node coordinates
  for (; it != end; ++it)
    {
      libmesh_assert (out_stream.good());

      out_stream << n++ << "\t";
      (*it)->write_unformatted(out_stream);
    }

  return;
}

void UCDIO::write_interior_elems(std::ostream& out_stream, const MeshBase& mesh)
{
  std::string type[] =
    { "edge",  "edge",  "edge",
      "tri",   "tri",
      "quad",  "quad",  "quad",
      "tet",   "tet",
      "hex",   "hex",   "hex",
      "prism", "prism", "prism",
      "pyramid" };

  MeshBase::const_element_iterator it  = mesh.elements_begin();
  const MeshBase::const_element_iterator end = mesh.elements_end();

  unsigned int e=1; // 1-based element number for UCD

  // Write element information
  for (; it != end; ++it)
    {
      libmesh_assert (out_stream.good());

      // PB: I believe these are the only supported ElemTypes.
      const ElemType etype = (*it)->type();
      if( (etype != TRI3) && (etype != QUAD4) &&
	  (etype != TET4) && (etype != HEX8) &&
          (etype != PRISM6) && (etype != PYRAMID5) )
	{
	  libMesh::err << "Error: Unsupported ElemType for UCDIO."
		       << std::endl;
	  libmesh_error();
	}

      out_stream << e++ << " 0 " << type[etype] << "\t";
      // (*it)->write_ucd_connectivity(out_stream);
      (*it)->write_connectivity(out_stream, UCD);
    }

  return;
}

void UCDIO::write_nodal_data(const std::string& fname,
			     const std::vector<Number>&soln,
			     const std::vector<std::string>& names)
{
  std::ofstream out_stream(fname.c_str());

  const MeshBase& mesh = MeshOutput<MeshBase>::mesh();

  // UCD doesn't work in 1D
  libmesh_assert (mesh.mesh_dimension() != 1);

  // Write header
  this->write_header(out_stream,mesh,mesh.n_elem(),names.size());

  // Write the node coordinates
  this->write_nodes(out_stream,mesh);

  // Write the elements
  this->write_interior_elems(out_stream,mesh);

  // Write the solution
  this->write_soln(out_stream,mesh,names,soln);

  return;
}

void UCDIO::write_soln(std::ostream& out_stream, const MeshBase& mesh,
		       const std::vector<std::string>& names,
		       const std::vector<Number>&soln)
{
  libmesh_assert (out_stream.good());

  // First write out how many variables and how many components per variable
  out_stream << names.size();
  for( unsigned int i = 0; i < names.size(); i++ )
    {
      libmesh_assert (out_stream.good());
      // Each named variable has only 1 component
      out_stream << " 1";
    }
  out_stream << std::endl;

  // Now write out variable names and units. Since we don't store units
  // We just write out dummy.
  for( std::vector<std::string>::const_iterator var = names.begin();
       var != names.end();
       var++)
    {
      libmesh_assert (out_stream.good());
      out_stream << (*var) << ", dummy" << std::endl;
    }

  // Now, for each node, write out the solution variables
  std::size_t nv = names.size();
  for( std::size_t n = 1; // 1-based node number for UCD
       n <= mesh.n_nodes(); n++)
    {
      libmesh_assert (out_stream.good());
      out_stream << n;

      for( std::size_t var = 0; var != nv; var++ )
	{
	  std::size_t idx = nv*(n-1) + var;

	  out_stream << " " << soln[idx];
	}
      out_stream << std::endl;
    }

  return;
}

} // namespace libMesh
